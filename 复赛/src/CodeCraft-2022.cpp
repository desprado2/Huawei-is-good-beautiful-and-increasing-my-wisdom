#include <bits/stdc++.h>
#include <sys/mman.h>
using namespace std;


#define cur_time() chrono::system_clock::now().time_since_epoch().count()

auto clock_start = cur_time();

#ifdef LOCAL
const int FIRST_STAGE_TIME_LIMIT = 10;
const int HARD_TIME_LIMIT = 40;
#else
const int FIRST_STAGE_TIME_LIMIT = 80;
const int HARD_TIME_LIMIT = 296;
/*
#define cerr 0 && cerr
#undef assert
#define assert(expr) (expr)*/
#endif

inline bool time_limit_ok(int time_limit = HARD_TIME_LIMIT) {
    return cur_time() - clock_start < time_limit * ((long long)1e9);
}

#ifdef LOCAL
#warning "LOCAL"
#define PREFIX "../"
#else
#define PREFIX "/"
#endif

const char* FILE_DEMAND = PREFIX "data/demand.csv";
const char* FILE_BANDWIDTH = PREFIX "data/site_bandwidth.csv";
const char* FILE_QOS = PREFIX "data/qos.csv";
const char* FILE_CONFIG = PREFIX "data/config.ini";
const char* FILE_OUTPUT = PREFIX "output/solution.txt";

const int M = 35; // number of customers
const int N = 135; // number of servers
const int TS = 8928; // number of time slots
const int S = 100; // number of stream kinds
const int INF = 1 << 29;
// const int THRESHOLD = 0;

int n;      // input number of customers
int m;      // input number of servers
int t;      // input number of time slots
int q;      // qos limit
int V;      // base cost
int cnt95;  // ceil (0.95 * t)
int cnt5;  // t - cnt95


string cid[M], sid[N];
vector <string> stream_id[TS];
int demand[TS][M][S], total_demand[TS];
int bandwidth[N];
bool qos[N][M]; // ping[i][j] < q

inline vector <string> getline_split (FILE* f) {
	const int BUFFER_SIZE = 1 << 16;
	static char tmp[BUFFER_SIZE];
	if (fgets(tmp, BUFFER_SIZE, f) == NULL) return {};
	size_t len = strlen(tmp);
	while (len && tmp[len - 1] <= ' ') len--;
	vector <string> ret;
	for (size_t i = 0; i < len; ) {
		size_t j = i;
		while (j < len && tmp[j] != ',') j++;
		tmp[j] = '\0';
		ret.push_back(tmp + i);
		i = j + 1;
	}
	return ret;
}

void input_data () {
	cerr << "input_data" << endl;
	n = m = t = 0;
	{
		FILE* fdemand = fopen (FILE_DEMAND, "r");
		auto ids = getline_split(fdemand);
		for (size_t i = 2; i < ids.size(); i++) {
			m++;
			cid[i - 2] = std::move(ids[i]);
		}
        string last_time = "";
		for ( ; ; ) {
			auto ret = getline_split(fdemand);
			if (!ret.size()) break;
			assert (ret.size() == m + 2u);
            string cur_time = ret[0];
            if (cur_time != last_time) last_time = cur_time, total_demand[t++] = 0;
			int tick = t - 1;
            int stream_idx = int(stream_id[tick].size());
            stream_id[tick].push_back(ret[1]);
			for (size_t i = 2; i < ret.size(); i++) {
                int cur_demand = stoi(ret[i]);
				demand[tick][i - 2][stream_idx] = cur_demand;
                total_demand[tick] += cur_demand;
			}
		}
		fclose (fdemand);
	}
	{
		FILE* fbandwidth = fopen (FILE_BANDWIDTH, "r");
		auto header = getline_split(fbandwidth);
		for ( ; ; ) {
			auto ret = getline_split(fbandwidth);
			if (!ret.size()) break;
			assert (ret.size() == 2);
			int s = n ++;
			sid[s] = std::move(ret[0]);
			bandwidth[s] = stoi(ret[1]);
		}
		fclose (fbandwidth);
	}
	{
		FILE* fconfig = fopen (FILE_CONFIG, "r");
		assert (fscanf(fconfig, "%*s qos_constraint=%d base_cost=%d", &q, &V) == 1);
		fclose (fconfig);
	}
	{
		FILE* fqos = fopen (FILE_QOS, "r");
		auto ids = getline_split(fqos);
		vector <int> c_ids (ids.size());
		for (size_t i = 1; i < ids.size(); i++) {
			int c = -1;
			for (int j = 0; j < m; j++) {
				if (ids[i] == cid[j]) {
					c = j;
					break;
				}
			}
			assert (c != -1);
			c_ids [i] = c;
		}
		for ( ; ; ) {
			auto ret = getline_split(fqos);
			if (!ret.size()) break;
			assert (ret.size() == ids.size());
			int s = -1;
			for (int j = 0; j < n; j++) {
				if (ret[0] == sid[j]) {
					s = j;
					break;
				}
			}
			assert (s != -1);
			for (size_t i = 1; i < ret.size(); i++) {
				qos[s][c_ids[i]] = stoi(ret[i]) < q;
			}
		}
		fclose (fqos);
	}
	cerr << "n = " << n << " m = " << m << " t = " << t << endl;
}

struct server_record {
    int customer, stream, flow;
};
vector<int> init_server_order;
struct FlowGraph{
	int maxFlow[N];
	vector<server_record> rec[N];

	FlowGraph(){
		reset();
		for (int s = 0; s < n; s++)
			init_server_order.push_back(s);
	}
	void reset(){
		memset(maxFlow, 0, sizeof(maxFlow));
		for (int s = 0; s < n; s++)
			rec[s].clear();
	}

	void set_max_flow(int server, int flow){
		maxFlow[server] = flow;
	}

	bool binPacking_bestFit(int tick, const vector<int>& server_order = init_server_order){
		int remain_flow[N];
		for (int s = 0; s < n; s++)
			remain_flow[s] = maxFlow[s];
		
		vector<pair<int, int>> streams; // first: customer id; second: stream id
		for (int c = 0; c < m; c++){
			int stream_num = (int)stream_id[tick].size();
			for (int stream = 0; stream < stream_num; stream++) if (demand[tick][c][stream]){
				streams.push_back({c, stream});
			}
		}

		sort(streams.begin(), streams.end(), [&](auto x, auto y){
			int D = demand[tick][x.first][x.second] - demand[tick][y.first][y.second];
			if (D) return D > 0;
			return x > y; // ***this order can be changed***
		});

		for (auto [c, stream] : streams){
			int streamFlow = demand[tick][c][stream];
			int bestFit = -1;
			for (int s = 0; s < n; s++){
				if (qos[s][c] && remain_flow[s] >= streamFlow){
					if (bestFit == -1 || remain_flow[bestFit] > remain_flow[s])
						bestFit = s;
				}
			}
			if (bestFit == -1) return false;
			remain_flow[bestFit] -= streamFlow;

			rec[bestFit].push_back({c, stream, streamFlow});
		}

		return true;
	}
};

struct Solution {
    
    struct customer_record {
        int server, stream, flow;
    };
    vector <vector<server_record>> allocation[N]; //allocation[s][tick] : {<customer, stream, flow>...}
    int flow95[N], value;
    Solution () {value = INT_MAX;}

    void output() {
        cerr << value << endl;
        freopen (FILE_OUTPUT, "w", stdout);
        for (int tick = 0; tick < t; tick ++) {
            vector < customer_record > alloc[m];
            for (int s = 0; s < n; s++) {
                for (auto r : allocation[s][tick]) {
                    alloc[r.customer].push_back({s, r.stream, r.flow});
                }
            }
            for (int c = 0; c < m; c++) {
                printf("%s:", cid[c].c_str());
                bool first = true;
                for (auto &r : alloc[c]) {
                    if (first) first = false;
                    else putchar (',');
                    printf("<%s,%s,%d>", sid[r.server].c_str(), stream_id[tick][r.stream].c_str(), r.flow);
                }
                putchar ('\n');
            }
        }
        fclose(stdout);
    }
};

void init() {
	cnt95 = ceil(0.95 * t);
	cnt5 = t - cnt95;
	for (int s = 0; s < n; s++) init_server_order.push_back(s);
}

Solution global_ans;
pthread_mutex_t global_mutex = PTHREAD_MUTEX_INITIALIZER;

struct Solver {
    mt19937 rng; // absolute hero

	FlowGraph G[TS];

    Solution ans;

	Solution getSolution(){
		Solution sol;
		double sum = 0;
		for (int s = 0; s < n; s++){
			sol.allocation[s].resize(t);
			vector<int> flows(t);
			for (int tick = 0; tick < t; tick++){
				sol.allocation[s][tick] = G[tick].rec[s];
				flows[tick] = 0;
				for (auto r : sol.allocation[s][tick]){
					flows[tick] += r.flow;
				}
			}

			sort(flows.begin(), flows.end());
			sol.flow95[s] = flows[cnt95 - 1];

			if (flows.back() > 0){
				if (sol.flow95[s] <= V) sum += V;
				else sum += (double)(sol.flow95[s] - V) * (sol.flow95[s] - V) / bandwidth[s] + sol.flow95[s];
			}
		}

		sol.value = int(sum + 0.5);
		
		return sol;
	}

    void main(int pass) {
		for (int tick = 0; tick < t; tick++){
			for (int s = 0; s < n; s++){
				G[tick].set_max_flow(s, bandwidth[s]);
			}
			G[tick].binPacking_bestFit(tick);
		}

		ans = getSolution();
    }
};

// multi-thread

void set_affinity (int core) {
    cpu_set_t mask;
	CPU_ZERO(&mask);
	CPU_SET((int)core, &mask);
	assert(!pthread_setaffinity_np(pthread_self(), sizeof(mask), &mask));
}

const int NThread = 4;
Solver worker[NThread];

void *call_worker_main(void *arg) {
    int pass = (int)(unsigned long long)arg;
	set_affinity (pass);
	worker[pass].main(pass);
    return NULL;
}

void dispatch_threads (void* (*func)(void*)) {
    pthread_t thrd[NThread]; void *thrd_ret;
	for(int i = 0; i < NThread; i++)
		pthread_create(&thrd[i], NULL, func, (void*)(1ull*i));
    for(int i = 0; i < NThread; i++) {
		pthread_join(thrd[i], &thrd_ret);
		if (global_ans.value > worker[i].ans.value)
			global_ans = worker[i].ans;
    }
}

int main() {
	input_data();
	cerr << "hello";

    init();

    dispatch_threads (call_worker_main);

    global_ans.output();
}

