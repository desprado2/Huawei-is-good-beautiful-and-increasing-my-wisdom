#include <bits/stdc++.h>
#include <sys/mman.h>
using namespace std;


#define cur_time() chrono::system_clock::now().time_since_epoch().count()

auto clock_start = cur_time();

#ifdef LOCAL
const int FIRST_STAGE_TIME_LIMIT = 10;
const int HARD_TIME_LIMIT = 40;
#else
const int FIRST_STAGE_TIME_LIMIT = 200;
const int HARD_TIME_LIMIT = 296;
#define cerr 0 && cerr
#undef assert
#define assert(expr) (expr)
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
	string fuck(tmp);
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
	cerr << "read demand.csv finished\n";
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
	cerr << "read bandwidth.csv finished\n";
	{
		FILE* fconfig = fopen (FILE_CONFIG, "r");
		assert (fscanf(fconfig, "%*s qos_constraint=%d base_cost=%d", &q, &V) == 2);
		fclose (fconfig);
	}
	cerr << "read config.ini finished\n";
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
	int remain_bw[N];
	vector<server_record> rec[N];

	struct StreamRecord{
		int customer;
		int stream;
	};
	vector<StreamRecord> remain_streams;

	FlowGraph(){
		reset();
	}
	void reset(){
		fill(remain_bw, remain_bw + n, 0);
		for (int s = 0; s < n; s++)
			rec[s].clear();
		remain_streams.clear();
	}

	/*
	O(NS log(NS))
	*/
	void init_streams(int tick){
		for (int c = 0; c < m; c++){
			int stream_num = (int)stream_id[tick].size();
			for (int stream = 0; stream < stream_num; stream++) if (demand[tick][c][stream]){
				remain_streams.push_back({c, stream});
			}
		}

		sort(remain_streams.begin(), remain_streams.end(), [&](auto x, auto y){
			int D = demand[tick][x.customer][x.stream] - demand[tick][y.customer][y.stream];
			if (D) return D < 0;
			return make_pair(x.customer, x.stream) < make_pair(y.customer, y.stream); // this order can be changed
		});
	}

	void add_bw(int server, int flow){
		remain_bw[server] += flow;
	}

	/*
	O(NMS)
	*/
	bool binPacking_bestFit(int tick, const vector<int>& server_order = init_server_order){
		/*
		目前是发现一个无法塞进箱子就立即退出
		实际上可以整个轮询一遍再退出，但是这样的话要remain_streams要改为list
		*/
		while (!remain_streams.empty()){
			auto [c, stream] = remain_streams.back();
			int streamFlow = demand[tick][c][stream];
			int bestFit = -1;
			for (int sID = 0; sID < n; sID++){
				int s = server_order[sID];
				if (qos[s][c] && remain_bw[s] >= streamFlow){
					if (bestFit == -1 || remain_bw[bestFit] > remain_bw[s])
						bestFit = s;
				}
			}
			if (bestFit == -1) return false;

			remain_bw[bestFit] -= streamFlow;
			rec[bestFit].push_back({c, stream, streamFlow});

			remain_streams.pop_back();
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
			map<int, vector<int>> alloc[m]; // alloc[customer]: {server:{stream1, stream2, ...}, server:{...}}
			for (int s = 0; s < n; s++){
				for (auto r : allocation[s][tick]){
					alloc[r.customer][s].push_back(r.stream);
				}
			}
			for (int c = 0; c < m; c++){
				printf("%s:", cid[c].c_str());
				bool first = true;
				for (auto & [server, streamVec] : alloc[c]){
					if (first) first = false;
					else putchar(',');
					printf("<%s", sid[server].c_str());
					for (auto stream : streamVec)
						printf(",%s", stream_id[tick][stream].c_str());
					putchar('>');
				}
				putchar('\n');
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

	vector<int> server_order;

	void Solver_init(){
		server_order = init_server_order;
	}
	/*
	TODO: set server order
	*/
	void set_server_order(){
	}

	vector<int> high_server;
	vector<int> burst_pool[TS / 20 + 5];
	vector<int> sorted_moments;
	void init_classification(int highServer_number){
		high_server.clear();
		for (int burst_time = 0; burst_time <= t / 20; burst_time++)
			burst_pool[burst_time].clear();

		for (int sID = 0; sID < highServer_number; sID++){
			high_server.push_back(server_order[sID]);
			burst_pool[t / 20].push_back(server_order[sID]);
		}
		for (int burst_time = t / 20; burst_time; burst_time--)
			shuffle(burst_pool[burst_time].begin(), burst_pool[burst_time].end(), rng);

		sorted_moments.clear();
		for (int tick = 0; tick < t; tick++)
			sorted_moments.push_back(tick);
		sort(sorted_moments.begin(), sorted_moments.end(), [](int x, int y){
			if (total_demand[x] != total_demand[y]) return total_demand[x] > total_demand[y];
			return stream_id[x].size() > stream_id[y].size();
		});
	}

	double server_flow95_distribution[N];
	/*
	It can be set in other ways, such as simulated annealing.
	*/
	void set_distribution(){
		for (int s = 0; s < n; s++)
			server_flow95_distribution[s] = 0;
		
		for (auto s : high_server){
			server_flow95_distribution[s] = 1.0;
		}
	}

	/*
	O((T + TN/20) * NMS)
	*/
	bool check(int midFlow){
        vector<int> pool[t / 20 + 1];
        for (int burst_time = 1; burst_time <= t / 20; burst_time++)
            pool[burst_time] = burst_pool[burst_time];

		int remain_burst[n];
		memset(remain_burst, 0, sizeof(remain_burst));
		for (auto s : high_server)
			remain_burst[s] = t / 20;
		
		for (auto tick : sorted_moments){
			bool bursted[n];
            int current_flow[n];
            memset (bursted, 0, sizeof(bursted));
			memset (current_flow, 0, sizeof(current_flow));

			G[tick].reset();
			G[tick].init_streams(tick);
			for (int s = 0; s < n; s++){
				current_flow[s] = min(bandwidth[s], int(server_flow95_distribution[s] * midFlow));
				G[tick].add_bw(s, current_flow[s]);
			}

			bool flag = G[tick].binPacking_bestFit(tick, server_order);
			while (!flag){
				/*
				The strategy of finding the new burstID can be changed
				For example, think of which servers can connect G[tick].remain_streams.back().customer
				*/
				int burstID = -1;
				for (int burst_time = t / 20; burst_time; burst_time--){
					for (int i = 0; i < (int)pool[burst_time].size(); i++){
						int s = pool[burst_time][i];
						if (!bursted[s]){
							burstID = s;
							bursted[s] = true;
							remain_burst[s]--;

							pool[burst_time].erase(pool[burst_time].begin() + i);
							pool[burst_time - 1].push_back(s);
							
							break;
						}
					}
					if (burstID >= 0) break;
				}

				if (burstID == -1) return false;

				int bw_increase = bandwidth[burstID] - current_flow[burstID];
				G[tick].add_bw(burstID, bw_increase);
				flag = G[tick].binPacking_bestFit(tick, server_order);
			}
		}
		return true;
	}

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

	Solution ans;
    void main(int pass) {
		rng.seed(114514 + pass);
		Solver_init();
		
		int BATCHNUM = 1;
		for (int K = 0; K < BATCHNUM; K++){
			/*
			TODO: set server_order here.
			Default order is order(i) = i.
			*/
			set_server_order();

			/*
			TODO: number of high servers can be changed in each round.
			*/
			init_classification(n);

			/*
			TODO: set all servers' flow95 distribution.
			Default distribution is ratio(s) = [s \in high_server].
			*/
			set_distribution();

			int L = 0, R = 0, last = -1;
			for (int s = 0; s < n; s++)
				R = max(R, (int)((double)bandwidth[s] / server_flow95_distribution[s] + 0.5));
			while (L < R){
				cerr << "L = " << L << " R = " << R << endl;
				int mid = (L + R) / 2;
				if (check(last = mid)) R = mid;
				else L = mid + 1;

				if (!time_limit_ok(FIRST_STAGE_TIME_LIMIT)) {cerr << "TLE\n"; break;}
			}
			if (last != R) check(R);
			cerr << "pass " << pass << " round " << K << " : R = " << R << endl;
			
			Solution cur = getSolution();
			if (ans.value > cur.value) ans = std::move(cur);

			if (!time_limit_ok(FIRST_STAGE_TIME_LIMIT)) break;
		}
    }
};

// multi-thread

void set_affinity (int core) {
    cpu_set_t mask;
	CPU_ZERO(&mask);
	CPU_SET((int)core, &mask);
	assert(!pthread_setaffinity_np(pthread_self(), sizeof(mask), &mask));
}

#ifdef LOCAL
const int NThread = 1;
#else
const int NThread = 4;
#endif

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

    init();

    dispatch_threads (call_worker_main);

    global_ans.output();
}

