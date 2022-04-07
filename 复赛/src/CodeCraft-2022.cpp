#include <bits/stdc++.h>
#include <sys/mman.h>
using namespace std;


#define cur_time() chrono::system_clock::now().time_since_epoch().count()
#define sqr(x) ((x) * (x))

auto clock_start = cur_time();

#ifdef LOCAL
const int BINARY_SERACH_TIME_LIMIT = 15;
const int FIRST_STAGE_TIME_LIMIT = 35;
const int HARD_TIME_LIMIT = 40;
#else
const int BINARY_SERACH_TIME_LIMIT = 80;
const int FIRST_STAGE_TIME_LIMIT = 270;
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

int flow_s, flow_c;
void init() {
	cnt95 = ceil(0.95 * t);
	cnt5 = t - cnt95;
	flow_c = 2, flow_s = m + 2;

	for (int tick = 0; tick < t; tick++)
		for (int c = 0; c < m; c++){
			int stream_num = stream_id[tick].size();
			customer_demand[tick][c] = 0;
			for (int stream = 0; stream < stream_num; stream++)
				customer_demand[tick][c] += demand[tick][c][stream];
		}
}

struct FlowGraph{
	struct edge {
        int v, nxt, f;
    };
    vector <edge> e;
    int n, S, T;
    vector <int> q, tag, head, cur;

    FlowGraph() { reset(-1, -1, -1); }
    void reset(int _n, int _S, int _T){
        n = _n;
        S = _S;
        T = _T;
        e.resize(2);
        q.resize(_n + 1);
        tag.resize(_n + 1);
        head.resize(_n + 1);
        fill(head.begin(), head.end(), 0);
        cur.resize(_n + 1);
    }
    void add(int u, int v, int f) {
        e.push_back({v, head[u], f});
        head[u] = int(e.size()) - 1;
        e.push_back({u, head[v], 0});
        head[v] = int(e.size()) - 1;
    }

    bool bfs() {
        for (int i = 0; i <= n; i++) tag[i] = 0;
        int he = 0, ta = 1;
        q[0] = S;
        tag[S] = 1;
        while (he < ta) {
            int x = q[he++];
            if (x == T) return true;
            for (int o = head[x]; o; o = e[o].nxt)
                if (e[o].f && !tag[e[o].v])
                    tag[e[o].v] = tag[x] + 1, q[ta++] = e[o].v;
        }
        return false;
    }
    int dfs(int x, int flow) {
        if (x == T) return flow;
        int used = 0;
        for (int &o = cur[x]; o; o = e[o].nxt) {
            if (e[o].f && tag[x] < tag[e[o].v]) {
                int ret = dfs(e[o].v, min(flow - used, e[o].f));
                if (ret) {
                    e[o].f -= ret; e[o ^ 1].f += ret;
                    used += ret;
                    if (used == flow) return flow;
                }
            }
        }
        return used;
    }
    int dinic() {
        int ans = 0;
        while (bfs()) {
            for (int i = 0; i <= n; i++) cur[i] = head[i];
            ans += dfs(S, INF);
        }
        return ans;
    }
	void print_graph(){
		cerr << "S : ";
		for (int o = head[S]; o; o = e[o].nxt){
			if (e[o].f) cerr << "<server " << e[o].v - flow_s << ", " << e[o].f << "> ";  
		}
		cerr << endl;
		for (int s = 0; s < n; s++){
			cerr << "server " << s << " : ";
			for (int o = head[s + flow_s]; o; o = e[o].nxt){
				int v = e[o].v;
				if (!e[o].f) continue;
				if (v == S) cerr << "<S, " << e[o].f << "> ";
				else cerr << "<customer " << v - flow_c << ", INF> "; 
			}
			cerr << endl;
		}
		for (int c = 0; c < m; c++){
			cerr << "customer " << c << " : ";
			for (int o = head[c + flow_c]; o; o = e[o].nxt){
				int v = e[o].v;
				if (!e[o].f) continue;
				if (v == T) cerr << "<T, " << e[o].f << "> ";
				else cerr << "<server " << v - flow_s << ", " << e[o].f << "> "; 
			}
			cerr << endl;
		}
		cerr << "T : ";
		for (int o = head[T]; o; o = e[o].nxt){
			if (e[o].f) cerr << "<customer " << e[o].v - flow_c << ", " << e[o].f << "> ";  
		}
		cerr << endl;
	}
};

struct FlowSolution {
    struct record {
        int c;
        int flow;
    };
    vector <vector<record>> allocation[N]; // allocation[s][tick] = {<c1, flow1>, <c2, flow2>, ...}
    int flow95[N], value;
    FlowSolution () {value = INT_MAX;}
};

struct Solution {
	struct server_record {
		int customer, stream, flow;
	};

    struct customer_record {
        int server, stream, flow;
    };
	
	/* allocation[s][tick] : {<customer, stream, flow>...}*/
    vector <vector<server_record>> allocation[N];
    int flow95[N], value;
    Solution () {value = INT_MAX;}

	void calc(){
		double sum = 0;
		for (int s = 0; s < n; s++){
			int flows[t];
			for (int tick = 0; tick < t; tick++){
				flows[tick] = 0;
				for (auto r : allocation[s][tick])
					flows[tick] += r.flow;
			}
			sort(flows, flows + t);

			flow95[s] = flows[cnt95 - 1];
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

struct FlowSeries{
	multiset<int> low_flows; // size == cnt95
	multiset<int> high_flows; // size == cnt5

	void clear(){low_flows.clear(); high_flows.clear();}
	int get_flow95()const{return *low_flows.rbegin();}
	void modify(int oldVal, int newVal){
		int flow95 = *low_flows.rbegin();
		if (oldVal <= flow95){
			low_flows.erase(low_flows.find(oldVal));
			low_flows.insert(newVal);
		}
		else{
			high_flows.erase(high_flows.find(oldVal));
			high_flows.insert(newVal);
		}
		if (!high_flows.empty() && *low_flows.rbegin() > *high_flows.begin()){
			int val1 = *low_flows.rbegin(), val2 = *high_flows.begin();
			low_flows.erase(--low_flows.end()), high_flows.erase(high_flows.begin());
			low_flows.insert(val2), high_flows.insert(val1);
		}
	}
};

struct StreamRecord{
	int tick, c, stream;
	bool operator < (const StreamRecord& other)const{
		return demand[tick][c][stream] > demand[other.tick][other.c][other.stream];
	}
};



Solution global_ans;
pthread_mutex_t global_mutex = PTHREAD_MUTEX_INITIALIZER;

struct Solver {
    mt19937 rng; // absolute hero

	FlowGraph G[TS], Gbase;

	vector< vector<int> > stream_allocated[TS]; // stream_allocated[tick][c][stream] : int
	vector<int> occupied_bandwidth[TS]; // occupied_bandwidth[tick][s] : int
	vector<int> assigned_flow[TS]; // assigned_flow[tick][c] : int
	vector<bool> server_bursted[TS]; // server_bursted[tick][s] : bool
	int tot_assigned_flow[TS];
	void Solver_init(int pass){
		rng.seed(114514 + pass);
		for (int tick = 0; tick < t; tick++){
			stream_allocated[tick].resize(m);
			for (int c = 0; c < m; c++)
				stream_allocated[tick][c].resize(stream_id[tick].size(), -1);
			
			occupied_bandwidth[tick].resize(n, 0);
			assigned_flow[tick].resize(m, 0);
			server_bursted[tick].resize(n, false);
			tot_assigned_flow[tick] = 0;
		}

		Gbase.reset(n + m + 2, 1, n + m + 2);
		for (int s = 0; s < n; s++)
			for (int c = 0; c < m; c++) {
				if (qos[s][c])
					Gbase.add(s + flow_s, c + flow_c, INF);
			}
	}

	vector<int> server_indices;
	/*
	TODO: set server order
	*/
	void set_server_order(){
		server_indices.clear();
		for (int s = 0; s < n; s++)
			server_indices.push_back(s);
		sort(server_indices.begin(), server_indices.end(), [&](auto x, auto y){
			return bandwidth[x] > bandwidth[y];
		});
	}

	struct ServerInformation {
        pair<int, int> flow[TS]; // first: flow, second: tick
        int unassigned_burst;
		int type;
    } server[N];

	void flow_in(double estimate_ratio){
        auto indices = server_indices;
        for (int tick = 0; tick < t; tick++){
            G[tick].reset(n + m + 2, 1, n + m + 2);
            for (int c = 0; c < m; c++){
                G[tick].add(c + flow_c, G[tick].T, customer_demand[tick][c]);
            }
            for (auto s : indices) {
                G[tick].add(G[tick].S, s + flow_s, int(bandwidth[s] * estimate_ratio));
                for (int c = 0; c < m; c++){
                    if (qos[s][c]) {
                        G[tick].add(s + flow_s, c + flow_c, INF);
                    }
                }
                G[tick].dinic();
            }
        }
    }

	vector<int> burst_pool[TS / 20 + 5];
	vector <int> burst_server[TS], high_server[TS];
    vector <int> sortedMoments;
    int bursted_value[TS];
	/*
	flow_in() and parse the result, set the servers' type.
	The stratagy is according to the pass.
	*/
	void init_classification(int pass){
		memset(bursted_value, 0, sizeof(bursted_value));

		if (pass == 0){
			/*
			Thread 0: set a estimate ratio for flow_in(), so that we can set how many servers are HIGH and the others are LOW
			Only HIGH servers' have possitive flow95, LOW servers' flow95 are all zero.
			Finetune the flow solution in step 4.
			*/
			flow_in(1.0);

			for (int tick = 0; tick < t; tick++){
				burst_server[tick].clear();
				high_server[tick].clear();
				for (int o = G[tick].head[G[tick].S]; o; o = G[tick].e[o].nxt){
					int s = G[tick].e[o].v - flow_s;
					int flow = G[tick].e[o ^ 1].f;
					server[s].flow[tick] = {flow, tick};
				}
			}

			for (auto s : server_indices){
				auto &si = server[s];

				sort(si.flow, si.flow + t, [&] (auto x, auto y) {
					if (x.first != y.first) return x.first < y.first;
					return total_demand[x.second] < total_demand[y.second];
				});
				si.unassigned_burst = 0;
				if (si.flow[cnt95 - 1].first > 0){
					si.type = 1; // HIGH
					for (int tick = 0; tick < cnt95; tick++){
						high_server[si.flow[tick].second].push_back(s);
					}
					for (int tick = cnt95; tick < t; tick++) {
						burst_server[si.flow[tick].second].push_back(s);
						bursted_value[si.flow[tick].second] += int(bandwidth[s] * 0.15);
					}
				}
				else {
					si.type = 2; // LOW
					si.unassigned_burst = cnt5;
				}
			}
		}
		else if (pass == 1){
			/*
			Thread 1: All the servers have possitive flow95, regardless it is HIGH or LOW.
			Finetune the flow solution in step 4.
			*/
			flow_in(1.0);

			for (int tick = 0; tick < t; tick++){
				burst_server[tick].clear();
				high_server[tick].clear();
				for (int o = G[tick].head[G[tick].S]; o; o = G[tick].e[o].nxt){
					int s = G[tick].e[o].v - m - 2;
					int flow = G[tick].e[o ^ 1].f;
					server[s].flow[tick] = {flow, tick};
				}
			}

			for (auto s : server_indices){
				auto &si = server[s];

				sort(si.flow, si.flow + t, [&] (auto x, auto y) {
					if (x.first != y.first) return x.first < y.first;
					return total_demand[x.second] < total_demand[y.second];
				});
				si.unassigned_burst = 0;
				if (si.flow[cnt95 - 1].first > 0){
					si.type = 1; // HIGH
					for (int tick = 0; tick < cnt95; tick++){
						high_server[si.flow[tick].second].push_back(s);
					}
					for (int tick = cnt95; tick < t; tick++) {
						//si.unassigned_burst++;
						//high_server[si.flow[tick].second].push_back(s);
						//totalCost[si.flow[tick].second] += si.flow[tick].first;
						burst_server[si.flow[tick].second].push_back(s);
						bursted_value[si.flow[tick].second] += int(bandwidth[s] * 0.3);
					}
				}
				else {
					si.type = 2; // LOW
					for (int tick = 0; tick < t; tick++)
						high_server[tick].push_back(s);
					si.unassigned_burst = cnt5;
				}
			}			
		}

        sortedMoments.resize(t);
        for (int tick = 0; tick < t; tick++)
            sortedMoments[tick] = tick;
        sort(sortedMoments.begin(), sortedMoments.end(), [&](int x, int y){
            return total_demand[x] - bursted_value[x] > total_demand[y] - bursted_value[y];
        });

        for (int b = 0; b <= t / 20; b++) burst_pool[b].clear();
        for (int s = 0; s < n; s++) {
            burst_pool[server[s].unassigned_burst].push_back(s);
        }
        for (int b = 1; b <= t / 20; b++) {
            shuffle(burst_pool[b].begin(), burst_pool[b].end(), rng);
        }
	}

	/*
	Reset high_server, burst_server, ServerInformation::unassigned_burst
	Recalculate the bursted_value, and sort the moments again according to the new result.
	*/
	void reclassify(int pass){
		for (int tick = 0; tick < t; tick++){
			high_server[tick].clear(), burst_server[tick].clear();
			bursted_value[tick] = 0;
		}
		// flow_in again, but consider the predetermined burst server
		/*for (int tick = 0; tick < t; tick++){
			init_graph(tick);
			for (int s : server_indices){
				G[tick].add(G[tick].S, s + flow_s, bandwidth[s] - occupied_bandwidth[tick][s]);
				G[tick].dinic();
			}
		}*/

		if (pass == 0){
			for (int s = 0; s < n; s++)
				server[s].unassigned_burst = cnt5;
			for (int tick = 0; tick < t; tick++){
				for (int s = 0; s < n; s++){
					if (server_bursted[tick][s]){
						burst_server[tick].push_back(s);
						server[s].unassigned_burst--;
						bursted_value[tick] += occupied_bandwidth[tick][s] + int((bandwidth[s] - occupied_bandwidth[tick][s]) * 0.3);
					}
					else{
						if (server[s].type == 1) // only consider HIGH server.
							high_server[tick].push_back(s);
					}
				}
			}
		}
		else if (pass == 1){
			for (int s = 0; s < n; s++)
				server[s].unassigned_burst = cnt5;
			for (int tick = 0; tick < t; tick++){
				for (int s = 0; s < n; s++){
					if (server_bursted[tick][s]){
						burst_server[tick].push_back(s);
						server[s].unassigned_burst--;
						bursted_value[tick] += occupied_bandwidth[tick][s] + int((bandwidth[s] - occupied_bandwidth[tick][s]) * 0.3);
					}
					else{
						high_server[tick].push_back(s);
					}
				}
			}
		}

		sortedMoments.resize(t);
        for (int tick = 0; tick < t; tick++)
            sortedMoments[tick] = tick;
        sort(sortedMoments.begin(), sortedMoments.end(), [&](int x, int y){
            return total_demand[x] - bursted_value[x] > total_demand[y] - bursted_value[y];
        });

        for (int b = 0; b <= t / 20; b++) burst_pool[b].clear();
        for (int s = 0; s < n; s++)
            burst_pool[server[s].unassigned_burst].push_back(s);
	}

	double server_flow95_distribution[N];
	int MID_FLOW;
	double LARGE_STREAM_RATIO;
	/*
	Set high servers' ratio
	It can be set in other ways, such as simulated annealing.
	*/
	void set_distribution(int pass){
		if (pass == 0){
			for (int s = 0; s < n; s++)
				if (server[s].type == 1) server_flow95_distribution[s] = 1.0;
				else server_flow95_distribution[s] = 0;
		}
		else if (pass == 1){
			for (int s = 0; s < n; s++)
				server_flow95_distribution[s] = 1.0;
		}
	}
	/*
	Add the normal edges
	*/
	void init_graph(int tick){
        G[tick] = Gbase;
        for (int c = 0; c < m; c++)
            G[tick].add(c + flow_c, G[tick].T, customer_demand[tick][c] - assigned_flow[tick][c]);
    }
	vector<bool> burst_plan[TS]; // burst_plan[tick][s] : bool
	/*
	Check whether the midFlow is large enough for customer demand.
	*/
	bool check(int midFlow){
		//cerr << "check " << midFlow << " ticks = " << sortedMoments.size() <<endl;
        bool ok = true;
        vector<int> pool[t / 20 + 1];
        for (int b = 1; b <= t / 20; b++) {
            pool[b] = burst_pool[b];
        }
        int remain_burst[n];
        for (int s = 0; s < n; s++) {
            remain_burst[s] = server[s].unassigned_burst;
        }

        for (auto tick : sortedMoments) {
            init_graph(tick);

            burst_plan[tick].resize(n);
            int current_flow[n];
            memset (current_flow, -1, sizeof(current_flow));
			fill (burst_plan[tick].begin(), burst_plan[tick].end(), false);

            for (int s : burst_server[tick]){
                current_flow[s] = bandwidth[s] - occupied_bandwidth[tick][s];
                burst_plan[tick][s] = true;
                G[tick].add(G[tick].S, s + flow_s, current_flow[s]);
            }

            int ans = G[tick].dinic();

            for (int s : high_server[tick]){
                current_flow[s] = min(int(midFlow * server_flow95_distribution[s]), bandwidth[s] - occupied_bandwidth[tick][s]);
                G[tick].add(G[tick].S, s + flow_s, current_flow[s]);
            }
            ans += G[tick].dinic();

			for (int s = 0; s < n; s++) if (current_flow[s] == -1){
				current_flow[s] = 0;
				G[tick].add(G[tick].S, s + flow_s, current_flow[s]);
			}

			vector<int> refund_burst; //RNM, 退钱！！！
            while (ans < total_demand[tick] - tot_assigned_flow[tick]){
                int burstID = -1;
                for (int b = t / 20; b; b--){
                    for (size_t i = 0; i < pool[b].size(); i++){
                        int s = pool[b][i];
                        if (!burst_plan[tick][s]) {
                            burstID = s;
                            burst_plan[tick][s] = true;
                            remain_burst[s]--;

                            pool[b].erase(pool[b].begin() + i);
                            pool[b - 1].push_back(s);
                            break;
                        }
                    }
                    if (burstID >= 0) break;
                }
                if (burstID < 0) {
                    ok = false;
                    break;
                }

                int flow_increase = bandwidth[burstID] - occupied_bandwidth[tick][burstID] - current_flow[burstID];
                current_flow[burstID] += flow_increase;
                for (int o = G[tick].head[G[tick].S]; o; o = G[tick].e[o].nxt) if (G[tick].e[o].v == burstID + flow_s) {
                    G[tick].e[o].f += flow_increase;
                    break;
                }

                int flow_utility = G[tick].dinic();
                ans += flow_utility;
                if (flow_utility == 0) {
                    pool[remain_burst[burstID]++].pop_back();
                    pool[remain_burst[burstID]].insert (pool[remain_burst[burstID]].begin(), burstID);
					refund_burst.push_back(burstID);
                    for (int o = G[tick].head[G[tick].S]; o; o = G[tick].e[o].nxt) if (G[tick].e[o].v == burstID + flow_s) {
                        G[tick].e[o].f -= flow_increase;
                        break;
                    }
                }
				//else cerr << "flow_utility of server " << burstID << " at tick " << tick << " is " << (double)flow_utility / (bandwidth[burstID] - occupied_bandwidth[tick][burstID]) << endl, this_thread::sleep_for(chrono::milliseconds(500));
            }
			for (int s : refund_burst)
				burst_plan[tick][s] = false;
			//cerr << "tick = " << tick << "ok = " << ok << endl;

            if (!ok) {
                return false;
            }
        }
        return ok;
	}

	FlowSolution get_flow_solution () {
        FlowSolution sol;
        {   // calculate value of solution, and validate bandwidth limit
            priority_queue<int, vector <int>, greater <int> > top_traffic[n];
            for (int s = 0; s < n; s++) top_traffic[s].push(0);
            size_t time5 = t - cnt95 + 1;
            for (int tick = 0; tick < t; tick++){
                for (int o = G[tick].head[G[tick].S]; o; o = G[tick].e[o].nxt){
                    int s = G[tick].e[o].v - flow_s;
                    int flow = G[tick].e[o ^ 1].f + occupied_bandwidth[tick][s];
                    assert (flow <= bandwidth[s]); // validatation
                    top_traffic[s].push(flow);
                    if (top_traffic[s].size() > time5) top_traffic[s].pop();
                }
            }
            for (int s = 0; s < n; s++){
                sol.flow95[s] = top_traffic[s].top();
            }
            sol.value = 0;
            for (int s = 0; s < n; s++) {
                sol.value += sol.flow95[s];
            }
        }
        {   // collect allocations, and check feasibility
            for (int tick = 0; tick < t; tick++){
                int provided[m] = {0};
                for (int s = 0; s < n; s++){
                    sol.allocation[s].push_back({});
                    for (int o = G[tick].head[s + flow_s]; o; o = G[tick].e[o].nxt){
                        int c = G[tick].e[o].v - 2;
                        int flow = G[tick].e[o ^ 1].f;
                        if (c >= 0 && c < m && flow) {
                            provided[c] += flow;
                            sol.allocation[s].back().push_back({c, flow});
                        }
                    }
                }
                for (int c = 0; c < m; c++) {
                    assert (provided[c] == customer_demand[tick][c] - assigned_flow[tick][c]); // validation
                }
            }
        }
        return sol;
    }

	static const int DEFAULT_LOW_REFUND = 20;
    static const int DEFAULT_RANDOM_RANGE = 100;
    FlowSolution finetune_with_error (
        const FlowSolution& ans,
        int time_limit = FIRST_STAGE_TIME_LIMIT,
        int random_range = DEFAULT_RANDOM_RANGE,
        int low_refund = DEFAULT_LOW_REFUND)
    {
        auto st = cur_time();
        FlowSolution cur = ans;
        vector <int> order = server_indices;
        reverse(order.begin(), order.end());
        int error[n], cnt_burst[n], extra_burst[n];
        auto count_bursts = [&] () {
            memset(cnt_burst, 0, sizeof(int) * n);
            for (int tick = 0; tick < t; tick ++) {
                for (int s = 0; s < n; s++) {
                    int flow = 0;
                    for (auto r : cur.allocation[s][tick]) {
                        flow += r.flow;
                    }
					flow += occupied_bandwidth[tick][s];
                    if (flow > cur.flow95[s]) {
                        cnt_burst[s] ++;
                    }
                }
            }
        };
        count_bursts();
        for (auto rs : order) if (cur.flow95[rs] > low_refund) {
			// Try to reduce the flow95 of server rs while allow other servers' flow95 increase by a random number in [0, random_range)
            memset(extra_burst, 0, sizeof(extra_burst));
            for (int s = 0; s < n; s++) {
                error[s] = min(bandwidth[s], cur.flow95[s] + int(rng() % random_range)) - cur.flow95[s];
            }
            int min_refund = INT_MAX;
            for (int tick = 0; tick < t && min_refund > low_refund; tick++) {
                int allow[n];
                vector <int> burst, normal;
                for (int i = 0; i < n; i++) {
                    int flow = 0;
                    for (auto r : cur.allocation[i][tick]) {
                        flow += r.flow;
                    }
					flow += occupied_bandwidth[tick][i];
                    if (flow > cur.flow95[i]) {
                        burst.push_back(i);
                        allow[i] = bandwidth[i] - occupied_bandwidth[tick][i];
                    } else {
                        normal.push_back(i);
                        allow[i] = cur.flow95[i] + error[i] - occupied_bandwidth[tick][i];
                    }
                }
                init_graph(tick);

                for (auto s : normal)
                    if (s == rs || allow[s])
                        G[tick].add(G[tick].S, s + flow_s, s == rs ? 0 : allow[s]);
                for (auto s : burst)
                    if (s == rs || allow[s])
                        G[tick].add(G[tick].S, s + flow_s, s == rs ? 0 : allow[s]);
                int tot = G[tick].dinic();
                if (allow[rs] != bandwidth[rs]) {
                    int cur_refund = allow[rs] - error[rs] - (total_demand[tick] - tot_assigned_flow[tick] - tot);
                    if (cur_refund <= low_refund) {
                        // if rs can burst, then take it as a burst
                        if (cnt_burst[rs] + extra_burst[rs] < cnt5) {
                            extra_burst[rs] ++;
                            allow[rs] = bandwidth[rs] - occupied_bandwidth[tick][rs];
                            cur_refund = cur.flow95[rs];
                        }
                        else { // otherwise we burst other nodes
                            for (auto s : normal) if (cnt_burst[s] + extra_burst[s] < cnt5) {
                                extra_burst[s] ++;
                                if (allow[s]) {
                                    int diff = (bandwidth[s] - occupied_bandwidth[tick][s]) - allow[s];
                                    allow[s] += diff;
                                    for (int o = G[tick].head[G[tick].S]; o; o = G[tick].e[o].nxt)
                                        if (G[tick].e[o].v == s + flow_s) {
                                            G[tick].e[o].f += diff;
                                            break;
                                        }
                                } else {
                                    G[tick].add(G[tick].S, s + flow_s, allow[s] = bandwidth[s] - occupied_bandwidth[tick][s]);
                                }
                                cur_refund += G[tick].dinic();
                                if (cur_refund > low_refund) break;
                            }
                        }
                    }
                    min_refund = min(min_refund, cur_refund);
                }
                // unlock blocked flows
                for (int o = G[tick].head[G[tick].S]; o; o = G[tick].e[o].nxt)
                    if (G[tick].e[o].v == rs + flow_s) {
                        G[tick].e[o].f = allow[rs];
                        break;
                    }
            }
            if (min_refund > low_refund) {
                // flow all blocked edges
                for (int tick = 0; tick < t; tick ++) {
                    G[tick].dinic();
                }
                auto nxt = get_flow_solution();
                if (nxt.value < cur.value) {
                    cur = std::move(nxt);
                    count_bursts();
                }
            }
            if (!time_limit_ok(time_limit)) break;
        }
        cerr << "[finetune] " << (cur_time() - st) / 1e9 << "s" << endl;
        return cur;
    }

    inline FlowSolution finetune (const FlowSolution& ans, int time_limit = HARD_TIME_LIMIT) {
        return finetune_with_error (ans, time_limit, 1, 0);
    }

    FlowSolution disperse (const FlowSolution& ans, int time_limit = HARD_TIME_LIMIT) {
        auto st = cur_time();
        FlowSolution cur = ans;
        vector <int> order = server_indices;
        sort(order.begin(), order.end(), [&] (auto u, auto v) {
            return ans.flow95[u] < ans.flow95[v];
        });
        int error[n];
        memset(error, 0, sizeof(error));
        for (int iu = 0; iu < n; iu++) {
            int u = order[iu];
            for (int iw = n - 1; iw > iu && cur.flow95[u] < bandwidth[u]; iw--) if (rng() & 1) { // randomized to utilize cpu
                int w = order[iw];
                if (cur.flow95[w] > 0) {
                    if (!time_limit_ok(time_limit)) break;
                    int vv = min(bandwidth[u] - cur.flow95[u], cur.flow95[w]);
                    error[u] = +vv;
                    error[w] = -vv;
                    bool ok = true;
                    int tick = -1;
                    vector <int> last_to_check;
                    for (tick = 0; tick < t; tick++) {
                        int allow[n];
                        bool no_check = false;
                        int total_allow = 0;
                        for (int i = 0; i < n; i++) {
                            int flow = 0;
                            for (auto r : cur.allocation[i][tick]) {
                                flow += r.flow;
                            }
							flow += occupied_bandwidth[tick][i];
                            if (flow > cur.flow95[i]) {
                                if (i == w) {
                                    last_to_check.push_back(tick);
                                    no_check = true;
                                }
                                allow[i] = bandwidth[i] - occupied_bandwidth[tick][i];
                            } else {
                                allow[i] = cur.flow95[i] + error[i] - occupied_bandwidth[tick][i];
                            }
                            total_allow += allow[i];
                        }
                        if (total_allow < total_demand[tick] - tot_assigned_flow[tick]) {
                            ok = false;
                            break;
                        }
                        init_graph(tick);
                        for (int s = 0; s < n; s++) if (allow[s])
                            G[tick].add(G[tick].S, s + flow_s, allow[s]);
                        if (no_check) continue;
                        int tot = G[tick].dinic();
                        if (tot != total_demand[tick] - tot_assigned_flow[tick]) {
                            ok = false;
                            break;
                        }
                    }
                    if (ok) {
                        for (auto tick : last_to_check) G[tick].dinic();
                        auto nxt = get_flow_solution();
                        assert (nxt.value <= cur.value);
                        cur = std::move(nxt);
                    }
                    error[u] = 0;
                    error[w] = 0;
                }
            }
        }
        cerr << "[disperse] " << fixed << setprecision(2) << (cur_time() - st) / 1e9 << "s" << endl;
        return cur;
    }

	void finetune_flow_solution(FlowSolution& sol){
		do{
        	shuffle(server_indices.begin(), server_indices.end(), rng);
			cerr << "value = " << sol.value << ", finetune.\n";
         	sol = finetune(finetune_with_error(sol, FIRST_STAGE_TIME_LIMIT), FIRST_STAGE_TIME_LIMIT);
        	sol = disperse(sol, FIRST_STAGE_TIME_LIMIT);
			cerr << "after finetune, value = " << sol.value << endl;
		} while (time_limit_ok(FIRST_STAGE_TIME_LIMIT));
	}

	Solution get_solution(const FlowSolution& flow_sol){
		// first round, best fit
		const double CAPACITY_RATIO = 1.0;
		for (int tick = 0; tick < t; tick++){
			// parse the flow_sol
			vector<pair<int, int>> bins[m]; // first: capacity, second: server index
			for (int s = 0; s < n; s++){
				for (auto r : flow_sol.allocation[s][tick])
					bins[r.c].push_back({int(r.flow * CAPACITY_RATIO), s});
			}

			// bin packing.
			int stream_num = stream_id[tick].size();
			for (int c = 0; c < m; c++){
				vector<StreamRecord> streams;
				for (int stream = 0; stream < stream_num; stream++)
					if (demand[tick][c][stream] && stream_allocated[tick][c][stream] == -1)
						streams.push_back({tick, c, stream});
				sort(streams.begin(), streams.end());

				set<pair<int, int>> remain_capacity; // first: remain capacity, second: server index
				for (auto [capa, s] : bins[c])
					remain_capacity.insert({min(capa, bandwidth[s] - occupied_bandwidth[tick][s]), s});

				for (auto [tick, c, stream] : streams){
					int F = demand[tick][c][stream];
					auto it = remain_capacity.lower_bound({F, -1});
					if (it == remain_capacity.end()) continue;

					int capa = it->first, s = it->second;
					occupied_bandwidth[tick][s] += F;
					stream_allocated[tick][c][stream] = s;

					remain_capacity.erase(it);
					remain_capacity.insert({capa - F, s});
				}
			}
		}
		cerr << "first round finished.\n";

		// before second round, record current server_flow and fucked streams.
		FlowSeries flow_series[n];
		for (int s = 0; s < n; s++){
			vector<int> flows(t);
			for (int tick = 0; tick < t; tick++)
				flows[tick] = occupied_bandwidth[tick][s];
			sort(flows.begin(), flows.end());
			
			flow_series[s].clear();
			for (int tick = 0; tick < cnt95; tick++)
				flow_series[s].low_flows.insert(flows[tick]);
			for (int tick = cnt95; tick < t; tick++)
				flow_series[s].high_flows.insert(flows[tick]);
		}

		vector<StreamRecord> remain_streams;
		for (int tick = 0; tick < t; tick++)
			for (int c = 0; c < m; c++)
				for (int stream = 0; stream < stream_id[tick].size(); stream++)
					if (demand[tick][c][stream] && stream_allocated[tick][c][stream] == -1)
						remain_streams.push_back({tick, c, stream});
		sort(remain_streams.begin(), remain_streams.end(), [&](auto x, auto y){
			return demand[x.tick][x.c][x.stream] < demand[y.tick][y.c][y.stream];
		});

		/* second round. if a stream is fucked because it's too large
		find a server whose cost will increase the least when allocate the stream to it. */
		for (auto [tick, c, stream] : remain_streams){
			int F = demand[tick][c][stream];
			int best_s = -1;
			double min_cost_increase = 1e10;
			for (int s = 0; s < n; s++) if (qos[s][c]){
				if (bandwidth[s] - occupied_bandwidth[tick][s] < F) continue;
				
				int old_flow95 = flow_series[s].get_flow95();
				flow_series[s].modify(occupied_bandwidth[tick][s], occupied_bandwidth[tick][s] + F);
				int new_flow95 = flow_series[s].get_flow95();
				double old_cost = sqr((double)(old_flow95 - V)) / bandwidth[s];
				double new_cost = sqr((double)(new_flow95 - V)) / bandwidth[s];
				if (new_cost - old_cost < min_cost_increase){
					min_cost_increase = new_cost - old_cost;
					best_s = s;
				}

				flow_series[s].modify(occupied_bandwidth[tick][s] + F, occupied_bandwidth[tick][s]);
			}
			
			if (best_s == -1) while(1); // TLE

			stream_allocated[tick][c][stream] = best_s;
			flow_series[best_s].modify(occupied_bandwidth[tick][best_s], occupied_bandwidth[tick][best_s] + F);
			occupied_bandwidth[tick][best_s] += F;
		}

		Solution ret;
		for (int s = 0; s < n; s++)
			ret.allocation[s].resize(t);
		for (int tick = 0; tick < t; tick++)
			for (int c = 0; c < m; c++)
				for (int stream = 0; stream < stream_id[tick].size(); stream++)
					if (demand[tick][c][stream]){
						assert(stream_allocated[tick][c][stream] != -1);
						int s = stream_allocated[tick][c][stream];
						ret.allocation[s][tick].push_back({c, stream, demand[tick][c][stream]});
					}

		ret.calc();
		return ret;
	}

	Solution ans;
    void main(int pass) {
		Solver_init(pass);
		
		/*
		Step 1: flow_in, set server type (HIGH/LOW/ZERO)
		*/ 
		set_server_order();
		init_classification(pass);

		/*
		Step 2: binary search. Get midFlow.
		*/
		set_distribution(pass);
		int L = 0, R = 0;
		for (int s = 0; s < n; s++)
			R = max(R, bandwidth[s]);
		while (L < R){
			int mid = (L + R) / 2;
			if (check(mid)) R = mid;
			else L = mid + 1;
		}
		MID_FLOW = R;
		LARGE_STREAM_RATIO = 1.0;
		cerr <<"MID_FLOW = " << MID_FLOW << endl;

		/*
		Step 3: determine where the "large streams" should be allocated.
		*/
		// get all the streams
		vector<StreamRecord> stream_rec;
		for (int tick = 0; tick < t; tick++){
			for (int c = 0; c < m; c++)
				for (int stream = 0; stream < stream_id[tick].size(); stream++)
					if (demand[tick][c][stream])
						stream_rec.push_back({tick, c, stream});
		}
		sort(stream_rec.begin(), stream_rec.end());
		// reset unassigned_burst
		for (int s = 0; s < n; s++)
			server[s].unassigned_burst = cnt5;
		// determine the "large streams"
		for (auto [tick, c, stream] : stream_rec){
			int F = demand[tick][c][stream];
			if (F <= LARGE_STREAM_RATIO * MID_FLOW) continue; // this stream isn't "large"
			
			int best_s = -1;
			// if the stream can be allocated to a server that has already bursted, allocate it.
			for (int s = 0; s < n; s++)
				if (qos[s][c] && server_bursted[tick][s] && bandwidth[s] - occupied_bandwidth[tick][s] >= F){
					if (best_s == -1 || bandwidth[best_s] - occupied_bandwidth[tick][best_s] > bandwidth[s] - occupied_bandwidth[tick][s])
						best_s = s;
				}
			if (best_s >= 0){
				occupied_bandwidth[tick][best_s] += F;
				assigned_flow[tick][c] += F;
				tot_assigned_flow[tick] += F;
				stream_allocated[tick][c][stream] = best_s;
				continue;
			}
			// else, if the stream can be allocated to a server that has bursted in former check(), let the server burst and allocate it.
			for (int s = 0; s < n; s++)
				if (qos[s][c] && !server_bursted[tick][s] && burst_plan[tick][s] && bandwidth[s] >= F){
					if (best_s == -1 || bandwidth[best_s] > bandwidth[s])
						best_s = s;
				}
			if (best_s >= 0){
				server[best_s].unassigned_burst--;	
				server_bursted[tick][best_s] = true;

				occupied_bandwidth[tick][best_s] += F;
				assigned_flow[tick][c] += F;
				tot_assigned_flow[tick] += F;
				stream_allocated[tick][c][stream] = best_s;
				continue;
			}
			// else, if the stream can be allocated to some unburst server, let the server burst and allocate it.
			for (int s = 0; s < n; s++)
				if (qos[s][c] && !server_bursted[tick][s] && server[s].unassigned_burst > 0 && bandwidth[s] >= F){
					if (best_s == -1 ||
					    make_pair(server[best_s].unassigned_burst, -bandwidth[best_s]) < make_pair(server[s].unassigned_burst, -bandwidth[s]))
						best_s = s;
				}
			if (best_s >= 0){
				server[best_s].unassigned_burst--;	
				server_bursted[tick][best_s] = true;

				occupied_bandwidth[tick][best_s] += F;
				assigned_flow[tick][c] += F;
				tot_assigned_flow[tick] += F;
				stream_allocated[tick][c][stream] = best_s;
				continue;
			}
			// fail to allocate this stream
			cerr << "fail to allocate large stream: tick = " << tick << " c = " << c << " stream = " << stream << " flow =  " << demand[tick][c][stream] << endl;
		}

		/*
		Step 4: binary seach again, calculate a flow_solution for the "small streams"
		*/
		FlowSolution flow_ans;
		reclassify(pass);
		//for (int b = t / 20; b; b--)
		//	cerr << "burst_time = " << b << " server number = " << burst_pool[b].size() << endl;
		const int BATCHNUM = 1000;
		for (int K = 0; K < BATCHNUM; K++){
			for (int b = 1; b <= t / 20; b++)
				shuffle(burst_pool[b].begin(), burst_pool[b].end(), rng);
			
			set_distribution(pass);

			int L = 0, R = 0, last = -1;
			for (int s = 0; s < n; s++)
				if (server_flow95_distribution[s] > 0)
					R = max(R, (int)((double)bandwidth[s] / server_flow95_distribution[s] + 0.5));
			while (L < R){
				//cerr << "L = " << L << " R = " << R << endl;
				int mid = (L + R) / 2;
				if (check(last = mid)) R = mid;
				else L = mid + 1;

				if (!time_limit_ok(FIRST_STAGE_TIME_LIMIT)) {cerr << "TLE\n"; break;}
			}
			if (last != R) check(R);
			cerr << "pass " << pass << " round " << K << " : R = " << R << endl;
			
			FlowSolution flow_sol = get_flow_solution();
			
			if (pass == 0 || pass == 1){
				flow_sol = finetune(flow_sol);
				if (flow_sol.value < flow_ans.value)
					flow_ans = std::move(flow_sol);
				cerr << "this flow_sol's value = " << flow_sol.value << endl;
				if (!time_limit_ok(BINARY_SERACH_TIME_LIMIT)) break;
			}
			else{

			}
		}
		if (pass == 0 || pass == 1){
			finetune_flow_solution(flow_ans);
			cerr << "flow_ans value = " << flow_ans.value << endl;

			ans = get_solution(flow_ans);
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
const int NThread = 2;
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

