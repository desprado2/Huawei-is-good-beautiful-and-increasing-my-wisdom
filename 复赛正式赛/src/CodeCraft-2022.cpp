#include <bits/stdc++.h>
#include <sys/mman.h>
using namespace std;


#define cur_time() chrono::system_clock::now().time_since_epoch().count()

auto clock_start = cur_time();

#ifdef LOCAL
const int BINARY_SERACH_TIME_LIMIT = 15;
const int FIRST_STAGE_TIME_LIMIT = 45;
const int HARD_TIME_LIMIT = 40;
#else
const int BINARY_SERACH_TIME_LIMIT = 80;
const int FIRST_STAGE_TIME_LIMIT = 285;
const int HARD_TIME_LIMIT = 290;
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

const int NThread = 4;

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
const int SPECIAL = 10;

int n;      // input number of servers
int m;      // input number of customers
int t;      // input number of time slots
int q;      // qos limit
int V;      // base cost
int cnt95;  // ceil (0.95 * t)
int cnt5;  // t - cnt95

int cnt90;  // ceil (0.90 * t)
int cnt10;  // t - cnt90

auto sqr (double x) { return x * x; }
string cid[M], sid[N];
vector <string> stream_id[TS];
int demand[TS][M][S], customer_demand[TS][M], total_demand[TS];
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
void init() {
	input_data();
	cnt95 = (int)ceil(0.95 * t);
	cnt5 = t - cnt95;
	// now you tell wtf is this
	cnt10 = max(cnt5, cnt5 * 2 - 1);
	cnt90 = t - cnt10;
}

struct server_record {
	int customer, stream, flow;
	bool operator == (const server_record& a) const {
		return customer == a.customer && stream == a.stream && flow == a.flow;
	}
	bool operator < (const server_record& a) const {
		if (flow != a.flow) return flow > a.flow;
		if (customer != a.customer) return customer < a.customer;
		return stream < a.stream;
	}
};

struct customer_record {
	int server, stream, flow;
	bool operator == (const customer_record& a) const {
		return server == a.server && stream == a.stream && flow == a.flow;
	}
	bool operator < (const customer_record& a) const {
		if (flow != a.flow) return flow > a.flow;
		if (server != a.server) return server < a.server;
		return stream < a.stream;
	}
};

struct Solution {

	/* allocation[s][tick] : {<customer, stream, flow>...}*/

    vector <vector<server_record>> allocation[N];
    int flow95[N], value;
	set <int> special;
    Solution () {value = INT_MAX;}

	void init() {
		value = INT_MAX;
		for (int i = 0; i < n; i++) allocation[i].clear();
		special.clear();
	}

	void calc(){
		assert (special.size() == SPECIAL);
		double sum = 0;
		for (int s = 0; s < n; s++){
			vector <int> flows(t);
			for (int tick = 0; tick < t; tick++){
				flows[tick] = 0;
				for (auto r : allocation[s][tick]) {
					assert (qos[s][r.customer]);
					flows[tick] += r.flow;
				}
				assert (flows[tick] <= bandwidth[s]);
			}
			sort(flows.begin(), flows.end());

			flow95[s] = flows[(special.count(s) ? cnt90 : cnt95) - 1];
			if (flows[t - 1]){
				if (flow95[s] <= V) sum += V;
				else sum += flow95[s] + sqr((double)(flow95[s] - V)) / bandwidth[s];
			}
		}
		cerr << "calc result : " << sum << endl;
		value = int(sum + 0.5);
	}

    void output() {
        cerr << "value = " << value << endl;
        freopen (FILE_OUTPUT, "w", stdout);
		{
			bool first = true;
			for (auto s : special) {
				if (!first) printf(",");
				first = false;
				printf("%s", sid[s].c_str());
			}
			printf("\n");
		}
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
					for (auto stream : streamVec) {
						printf(",%s", stream_id[tick][stream].c_str());
					}
					putchar('>');
				}
				putchar('\n');
			}
        }
        fclose(stdout);
    }
};
struct FlowSeries{
	int low, high;
	multiset<int> low_flows; // size == LOW

	void init ( vector <int> used, int LOW, int HIGH) {
		low = LOW;
		high = HIGH;
		sort(used.begin(), used.end());
		for (int i = 0; i < low; i++) low_flows.insert(used[i]);
	}
	int flow95() const{return *low_flows.rbegin();}

	void modify(int oldVal, int newVal){
		assert (low_flows.size() == low);
		low_flows.erase(low_flows.find(oldVal));
		low_flows.insert(newVal);
	}
	int peek(int before, int after) {
		return max(flow95(), after);
	}
};

// Solution global_ans;
pthread_mutex_t global_mutex = PTHREAD_MUTEX_INITIALIZER;

struct Solver {
    mt19937 rng; // absolute hero
	Solution ans;

	void main (int pass) {
		ans.init();
		for (int s = 0; s < n; s++) {
			ans.allocation[s].clear();
			ans.allocation[s].resize(t);
		}
		vector <vector <vector<server_record>>> unsatisfied(t); //unsatisfied[tick][c] = {<c, stream, flow>, ...}
		vector <vector <int>> sum_demands(t);
		for (int tick = 0; tick < t; tick++) {
			unsatisfied[tick].resize(m);
			sum_demands[tick].resize(m);
			int stream_cnt = (int) stream_id[tick].size();
			for (int c = 0; c < m; c++) {
				for (int k = 0; k < stream_cnt; k++) if (demand[tick][c][k]) {
					unsatisfied[tick][c].push_back({c, k, demand[tick][c][k]});
					sum_demands[tick][c] += demand[tick][c][k];
				}
			}
		}
		vector <int> servers(n);
		vector <double> weight(n);
		vector <int> high_ticks(t);
		iota (high_ticks.begin(), high_ticks.end(), 0);
		sort(high_ticks.begin(), high_ticks.end(), [] (auto x, auto y) { return total_demand[x] > total_demand[y];});
		int high_num;
		
		switch (pass >> 1) {
			case 0: high_num = int(t * 0.3); break;	
			case 1: high_num = int(t * 0.7); break;	
			case 2: high_num = int(t * 0.5); break;	
			case 3: high_num = int(t * 0.8); break;	
			default: high_num = int(t * ((pass % 40) * 0.01 + 0.6)); break;
		}
		high_ticks.resize(high_num);

		for (auto tick : high_ticks) {
			for (int s = 0; s < n; s++) {
				int sum = 0;
				for (int c = 0; c < m; c++) if (qos[s][c]) {
					sum += sum_demands[tick][c];
				}
				weight[s] += min(bandwidth[s], sum);
			}
		}

		iota(servers.begin(), servers.end(), 0);
		sort(servers.begin(), servers.end(), [&] (auto u, auto v){
			return weight[u] > weight[v];
		});
		for (int i = 0; i < SPECIAL; i++) {
			ans.special.insert (servers[i]);
		}

		if (pass & 1) reverse(servers.begin(), servers.end());
		vector <vector <int>> server_used(n, vector <int> (t));

		auto allocate_stage1 = [&] (int tick, int s, server_record i) {
			assert ((server_used[s][tick] + i.flow) <= bandwidth[s]);
			server_used[s][tick] += i.flow;
			sum_demands[tick][i.customer] -= i.flow;
			unsatisfied[tick][i.customer].erase(find(unsatisfied[tick][i.customer].begin(), unsatisfied[tick][i.customer].end(), i));
			ans.allocation[s][tick].push_back(i);
		};
		for (int s : servers) {
			vector <int> time_slots(t), s_demand(t);
			for (int tick = 0; tick < t; tick ++) {
				for (int c = 0; c < m; c++) if (qos[s][c]) {
					for (auto [customer, stream, flow] : unsatisfied[tick][c]){
						if (flow <= bandwidth[s]) {s_demand[tick] += flow;}
					}
				}
			}
			int tot_valid_demand = 0;
			iota(time_slots.begin(), time_slots.end(), 0);
			sort(time_slots.begin(), time_slots.end(), [&] (auto u, auto v) {
				return s_demand[u] > s_demand[v];
			});
			int HIGH = ans.special.count(s) ? cnt10 : cnt5;
			time_slots.resize(HIGH);
			for (auto tick : time_slots)
				tot_valid_demand += s_demand[tick];
			for (auto tick : time_slots) {
				vector <server_record> rs;
				for (int c = 0; c < m; c++) if (qos[s][c]) {
					for (auto i : unsatisfied[tick][c]) {
						rs.push_back(i);
					}
				}
				sort(rs.begin(), rs.end());
				int flow = bandwidth[s];
				for (auto i : rs) {
					if (i.flow <= flow) {
						allocate_stage1(tick, s, i);
						flow -= i.flow;
					}
				}
			}
		}
		cerr << "burst complete" << endl; 
		vector < pair <server_record, int> > all_unsatisfied;
		int tot_flow = 0;
		for (int tick = 0; tick < t; tick ++) {
			for (int c = 0; c < m; c++) {
				for (auto i : unsatisfied[tick][c]) {
					all_unsatisfied.push_back({i, tick});
					tot_flow += i.flow;
				}
			}
		}
		cerr << "tot_flow = " << tot_flow << endl;
		cerr << "all = " << all_unsatisfied.size() << endl;
		if (!time_limit_ok()) return;
		sort(all_unsatisfied.begin(), all_unsatisfied.end());
		if (!time_limit_ok()) return;
		
		vector <FlowSeries> server_flow(n);
		for (int s = 0; s < n; s++) {
			bool is_special = ans.special.count(s);
			server_flow[s].init(server_used[s], is_special ? cnt90 : cnt95, is_special ? cnt10 : cnt5);
		}
		auto allocate_stage2 = [&] (int tick, int s, server_record i) {
			int before = server_used[s][tick];
			int after = before + i.flow;
			assert (after <= bandwidth[s]);
			server_used[s][tick] = after;
			server_flow[s].modify(before, after);
			ans.allocation[s][tick].push_back(i);
		};
		auto cost = [] (int flow, int band) {
			flow = max(flow, V);
			return 1. / band * sqr(flow - V) + flow;
		};
		auto calc_diff_cost = [&] (int tick, int s, int flow) {
			int b95 = server_flow[s].flow95();
			int a95 = server_flow[s].peek(server_used[s][tick], server_used[s][tick] + flow);
			return cost(a95, bandwidth[s]) - cost(b95, bandwidth[s]);
		};
		int cnt = 0;
		for (auto [i, tick] : all_unsatisfied) {
			if (!time_limit_ok()) return;
			if (cnt++ % 10000 == 0) cerr << "cnt = " << cnt << endl;
			vector < pair <pair <double, int>, int> > candidate;
			for (int s = 0; s < n; s++) {
				if (qos[s][i.customer] && server_used[s][tick] + i.flow <= bandwidth[s]) {
					auto inc = calc_diff_cost(tick, s, i.flow);
					candidate.push_back({{inc, server_used[s][tick] - server_flow[s].flow95()}, s});
				}
			}
			assert (candidate.size() > 0);
			allocate_stage2(tick, min_element(candidate.begin(), candidate.end())->second, i);
		}
		ans.calc();
	}
};

// multi-thread

void set_affinity (int core) {
#ifndef __APPLE__
    cpu_set_t mask;
	CPU_ZERO(&mask);
	CPU_SET((int)core, &mask);
	assert(!pthread_setaffinity_np(pthread_self(), sizeof(mask), &mask));
#endif
}

Solver worker[NThread];

void *call_worker_main(void *arg) {
    int pass = (int)(unsigned long long)arg;
	set_affinity (pass % 4);
	worker[pass % 4].main(pass);
    return NULL;
}
int history = 0;
void dispatch_threads (void* (*func)(void*)) {
    pthread_t thrd[NThread]; void *thrd_ret;
	for(int i = 0; i < NThread; i++)
		pthread_create(&thrd[i], NULL, func, (void*)(1ull*i + history));
	history += 4;
    for(int i = 0; i < NThread; i++) {
		pthread_join(thrd[i], &thrd_ret);
    }
}

int main() {
	init();
	Solution ans;
	while (time_limit_ok()) {
	    dispatch_threads (call_worker_main);
		for (int i = 0 ; i < NThread; i++) if (ans.value > worker[i].ans.value) {
			ans = worker[i].ans;
		}
	}
	ans.output();	
}


