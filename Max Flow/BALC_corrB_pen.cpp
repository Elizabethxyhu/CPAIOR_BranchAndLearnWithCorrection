#include <iostream>
#include <vector>
#include <climits>
#include <algorithm>
#include <vector>
#include <math.h>
#include <stack>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <string>
#include <sstream>
#include <time.h>

#include "PiecewiseLinearFunction.h"

using namespace std;

const int maxnodes = 100;
const double penaltyCons = 0;

int nodes = maxnodes, src, dest;
int dist[maxnodes], q[maxnodes], work[maxnodes];
int indexTemp;  // store the minimum edge index

struct Edge {
	int from;
	int to;
	double cap;
	PiecewiseLinearFunction flow;
	PiecewiseLinearFunction remain;                 // cap - flow
	int rev;
	double preCap;
	bool use = false;
	bool preUse = false;
};

vector<Edge> g[maxnodes];       // store all edges, including the reverse edges (2D matrix)
vector<Edge> path;              // store each path from s to t
vector<int> indexPath;              // help to store path
vector<PiecewiseLinearFunction> remainPath;    // store the remain cap on the path
vector<double> capPath;
vector<double> preCapPath;
int totalPath = 0;
int infeaPath = 0;
int feaPath = 0;
int violate = 0;
int totalP = 0;
int totalFea = 0;
int totalVio = 0;
int PR = 0;
int TR = 0;
int PW = 0;

// Adds bidirectional edge
void addEdge(int s, int t, double c, double trueC, double start, double end, double a, double b) {
	Edge e1, e2;

	e1.from = s;
	e1.to = t;
	e1.preCap = c;
	e1.cap = trueC;
	e1.flow.assign(start, end, 0, 0, 0);
	e1.remain.assign(start, end, a, b, c);  // id stores the real capacity
	e1.rev = g[t].size();

	e2.from = t;
	e2.to = s;
	e2.preCap = 0;
	e2.cap = 0;
	e2.flow.assign(start, end, 0, 0, 0);
	e2.remain.assign(start, end, 0, 0, 0);
	e2.rev = g[s].size();

	g[s].push_back(e1);
	g[t].push_back(e2);
}

void Pre_addEdge(int s, int t, double realC, double preC) {
	Edge e1, e2;

	e1.from = s;
	e1.to = t;
	e1.cap = realC;
	e1.preCap = preC;
	e1.rev = g[t].size();

	e2.from = t;
	e2.to = s;
	e2.cap = 0;
	e2.preCap = 0;
	e2.rev = g[s].size();

	g[s].push_back(e1);
	g[t].push_back(e2);
}

bool dinic_bfs() {
	fill(dist, dist + nodes, -1);
	dist[src] = 0;
	int qt = 0;
	q[qt++] = src;
	for (int qh = 0; qh < qt; qh++) {
		int u = q[qh];
		for (int j = 0; j < (int)g[u].size(); j++) {
			Edge& e = g[u][j];
			int v = e.to;
			if (dist[v] < 0 && (e.remain.a[0] != 0 || e.remain.b[0] != 0)) {
				dist[v] = dist[u] + 1;
				q[qt++] = v;
			}
		}
	}
	return dist[dest] >= 0;
}

bool dinic_dfs(int u) {
	if (u == dest)
		return 1;
	for (int& i = work[u]; i < (int)g[u].size(); i++) {
		Edge& e = g[u][i];
		if (e.remain.a[0] == 0 && e.remain.b[0] == 0)
			continue;
		int v = e.to;
		if (dist[v] == dist[u] + 1) {
			bool d = dinic_dfs(v);
			if (d > 0) {
				path.push_back(e);      // record the edge
				indexPath.push_back(i);
				return d;
			}
		}
	}
	return 0;
}

PiecewiseLinearFunction maxFlow(int _src, int _dest) {
	src = _src;
	dest = _dest;
	PiecewiseLinearFunction result;

	while (dinic_bfs()) {
		fill(work, work + nodes, 0);
		while (dinic_dfs(src)) {
			for (int j = 0; j < path.size(); j++) {
				remainPath.push_back(path[j].remain);
			}
			PiecewiseLinearFunction df = MinCompare(remainPath);
			//df.merge();
			int range = df.a.size();

			for (int m = 0; m < range; m++) {

				vector<Edge> gtemp[maxnodes];
				for (int k = 0; k < nodes; k++) {
					gtemp[k] = g[k];
				}

				int i = 0;
				while (i < path.size()) {
					g[path[i].from][indexPath[i]].remain = MinusSpecial(g[path[i].from][indexPath[i]].remain, df, m);
					g[path[i].to][path[i].rev].remain = PlusSpecial(g[path[i].to][path[i].rev].remain, df, m);
					i++;
				}

				for (int k = 0; k < nodes; k++) {
					for (int l = 0; l < g[k].size(); l++) {
						g[k][l].remain.start_x[0] = df.start_x[m];
						g[k][l].remain.end_x[0] = df.end_x[m];
					}
				}
				result.assign(df.start_x[m], df.end_x[m], df.a[m], df.b[m], df.id[m]);
				vector<Edge> pathT = path;
				path.clear();
				vector<int> indexT = indexPath;
				indexPath.clear();
				remainPath.clear();
				//result.output();
				//std::cout << "----------------------------------------------" << endl;

				PiecewiseLinearFunction df2 = maxFlow(src, dest);
				if (df2.a.size() != 0) {
					result = Plus(result, df2);
				}

				if (df2.a.size() == 0 && result.a.size() != 0) {
					for (int k = 0; k < nodes; k++) {
						for (int l = 0; l < g[k].size(); l++) {
							for (int m = 0; m < g[k][l].remain.id.size(); m++) {
								/*if (g[k][l].remain.id[m] < 0) {
									// std::cout << "penalty" << endl;
									int tail = result.id.size();
									result.id[tail - 1] = result.id[tail - 1] + g[k][l].remain.id[m] * penalty;
								}*/
							}
						}
					}
				}
				//result.output();
				//std::cout << "----------------------------------------------" << endl;

				if (m < (range - 1)) {
					for (int k = 0; k < nodes; k++) {
						for (int l = 0; l < g[k].size(); l++) {
							g[k][l].remain = gtemp[k][l].remain;
						}
					}
				}

				path = pathT;
				indexPath = indexT;
			}
			path.clear();
			indexPath.clear();
		}
	}
	//result.merge();
	return result;
}

bool De_bfs() {
	fill(dist, dist + nodes, -1);
	dist[src] = 0;
	int qt = 0;
	q[qt++] = src;
	for (int qh = 0; qh < qt; qh++) {
		int u = q[qh];
		for (int j = 0; j < (int)g[u].size(); j++) {
			Edge& e = g[u][j];
			int v = e.to;
			if (dist[v] < 0 && e.cap > 0) {
				dist[v] = dist[u] + 1;
				q[qt++] = v;
			}
		}
	}
	return dist[dest] >= 0;
}

bool De_dfs(int u) {
	if (u == dest)
		return 1;
	for (int& i = work[u]; i < (int)g[u].size(); i++) {
		Edge& e = g[u][i];
		if (e.cap <= 0)
			continue;
		int v = e.to;
		if (dist[v] == dist[u] + 1) {
			bool d = De_dfs(v);
			if (d > 0) {
				path.push_back(e);      // record the edge
				indexPath.push_back(i);
				return d;
			}
		}
	}
	return 0;
}

double dinic_min(vector<double> p) {
	double mini = p[0];
	indexTemp = 0;
	for (int i = 1; i < p.size(); i++) {
		if (mini > p[i]) {
			mini = p[i];
			indexTemp = i;
		}
	}
	return mini;
}

double De_maxFlow(int _src, int _dest) {
	src = _src;
	dest = _dest;
	double result = 0;
	while (De_bfs()) {
		fill(work, work + nodes, 0);
		while (De_dfs(src)) {
			for (int j = 0; j < path.size(); j++) {
				capPath.push_back(path[j].cap);
			}
			double df = dinic_min(capPath);
			for (int i = 0; i < path.size(); i++) {
				g[path[i].from][indexPath[i]].cap -= df;
				g[path[i].from][indexPath[i]].use = true;
				g[path[i].to][path[i].rev].cap += df;
			}
			result += df;
			path.clear();
			indexPath.clear();
			capPath.clear();
		}
	}
	return result;
}

bool Pre_bfs() {
	fill(dist, dist + nodes, -1);
	dist[src] = 0;
	int qt = 0;
	q[qt++] = src;
	for (int qh = 0; qh < qt; qh++) {
		int u = q[qh];
		for (int j = 0; j < (int)g[u].size(); j++) {
			Edge& e = g[u][j];
			int v = e.to;
			if (dist[v] < 0 && e.preCap > 0) {
				dist[v] = dist[u] + 1;
				q[qt++] = v;
			}
		}
	}
	return dist[dest] >= 0;
}

bool Pre_dfs(int u) {
	if (u == dest)
		return 1;
	for (int& i = work[u]; i < (int)g[u].size(); i++) {
		Edge& e = g[u][i];
		if (e.preCap <= 0)
			continue;
		int v = e.to;
		if (dist[v] == dist[u] + 1) {
			bool d = Pre_dfs(v);
			if (d > 0) {
				path.push_back(e);      // record the edge
				indexPath.push_back(i);
				return d;
			}
		}
	}
	return 0;
}

double Pre_error_maxFlow(int _src, int _dest) {
	src = _src;
	dest = _dest;
	double result = 0;
	while (Pre_bfs()) {
		fill(work, work + nodes, 0);
		while (Pre_dfs(src)) {
			totalPath++;
			for (int j = 0; j < path.size(); j++) {
				capPath.push_back(path[j].cap);
				preCapPath.push_back(path[j].preCap);
			}
			double preDf = dinic_min(preCapPath);
			double df = capPath.at(indexTemp);
			for (int i = 0; i < path.size(); i++) {
				g[path[i].from][indexPath[i]].cap -= df;
				g[path[i].to][path[i].rev].cap += df;
				g[path[i].from][indexPath[i]].preUse = true;
				g[path[i].from][indexPath[i]].preCap -= preDf;
				g[path[i].to][path[i].rev].preCap += preDf;
			}
			vector<double>::iterator minimum = min_element(begin(capPath), end(capPath));
			if (df > * minimum) {
				infeaPath++;
			}
			/*else if (df <= *minimum) {     // error or feasible error
				result += df;
			}*/
			result += df;
			path.clear();
			indexPath.clear();
			capPath.clear();
			preCapPath.clear();
		}
	}
	return result;
}

double Pre_feaE_maxFlow(int _src, int _dest) {
	src = _src;
	dest = _dest;
	double result = 0;
	while (Pre_bfs()) {
		fill(work, work + nodes, 0);
		while (Pre_dfs(src)) {
			totalPath++;
			for (int j = 0; j < path.size(); j++) {
				capPath.push_back(path[j].cap);
				preCapPath.push_back(path[j].preCap);
			}
			/*for (int k = path.size() - 1; k >= 0; k--) {
				std::cout << path[k].from + 1 << " -> ";
				if (k == 0) {
					std::cout << path[k].from + 1 << " -> " << path[k].to + 1 << endl;
				}
			}*/
			double preDf = dinic_min(preCapPath);
			double df = capPath.at(indexTemp);
			for (int i = 0; i < path.size(); i++) {
				g[path[i].from][indexPath[i]].cap -= df;
				g[path[i].to][path[i].rev].cap += df;
				g[path[i].from][indexPath[i]].preUse = true;
				g[path[i].from][indexPath[i]].preCap -= preDf;
				g[path[i].to][path[i].rev].preCap += preDf;
			}
			vector<double>::iterator minimum = min_element(begin(capPath), end(capPath));
			if (df > * minimum) {
				infeaPath++;
			}
			else if (df <= *minimum) {     // error or feasible error
				result += df;
			}
			//result += df;
			path.clear();
			indexPath.clear();
			capPath.clear();
			preCapPath.clear();
		}
	}
	return result;
}

double Pre_2S_maxFlow(int _src, int _dest) {
	src = _src;
	dest = _dest;
	double result = 0;
	while (Pre_bfs()) {
		fill(work, work + nodes, 0);
		while (Pre_dfs(src)) {
			totalPath++;
			for (int j = 0; j < path.size(); j++) {
				capPath.push_back(path[j].cap);
				preCapPath.push_back(path[j].preCap);
			}
			/*for (int k = path.size() - 1; k >= 0; k--) {
				std::cout << path[k].from + 1 << " -> ";
				if (k == 0) {
					std::cout << path[k].from + 1 << " -> " << path[k].to + 1 << endl;
				}
			}*/
			double preDf = dinic_min(preCapPath);
			vector<double>::iterator minimum = min_element(begin(capPath), end(capPath));
			double df = *minimum;
			//double df = capPath.at(indexTemp);
			for (int i = 0; i < path.size(); i++) {
				g[path[i].from][indexPath[i]].cap -= df;
				g[path[i].to][path[i].rev].cap += df;
				g[path[i].from][indexPath[i]].preUse = true;
				g[path[i].from][indexPath[i]].preCap -= preDf;
				g[path[i].to][path[i].rev].preCap += preDf;
			}
      if (df <= 0){
        result = result - penaltyCons;
      }
			/*if (df > * minimum) {
				infeaPath++;
			}
			else if (df <= *minimum) {     // error or feasible error
				result += df;
			}*/
			result += df;
			path.clear();
			indexPath.clear();
			capPath.clear();
			preCapPath.clear();
		}
	}
	return result;
}

double Preprocessing(vector<PiecewiseLinearFunction> v) {
	double start_x_new = -INF;
	if (v[0].a[0] != 0) {
		vector<double> vTemp;
		int size = v.size();
		vTemp.resize(size);
		for (int i = 0; i < size; i++) {
			double start = -INF;
			if (v[i].a[0] != 0) {
				start = -(v[i].b[0] / v[i].a[0]);
			}
			vTemp[i] = start;
			//            std::cout << vTemp[i] << " ";
		}
		vector<double>::iterator newStart = min_element(begin(vTemp), end(vTemp));  // at least one edge is positive
		if (*newStart > start_x_new) {
			start_x_new = *newStart;
		}
	}
	return start_x_new;
}


int main(int argc, char* argv[]) {

	if (argc != 9)
		printf("usage: ./BALC_corrB_pen [graph file] [train data] [test data] [prediction file] [runtime file] [iteration] [result file] [Ridge para]");

	// from to realCap fea1 fea2
	// exp 1
	int n = 12;
	nodes = n;
	int nodeNum = n;
	int linkNum = 18;
	int featureNum = 8;
	int benchmarkNum = 610;
	vector<double> alpha(featureNum, 1);                   // coefficients
	double regret = INF;
	vector<PiecewiseLinearFunction> ans;
	PiecewiseLinearFunction finalRes;

	vector<vector<double>> data(linkNum, vector<double>(featureNum, 1));    // [linkNum, featureNum]
	vector<double> realCap(linkNum);
	double benchmarkId;

	int testNum = 179;
	double testmarkId;
	vector<vector<double>> testData(linkNum, vector<double>(featureNum, 0));
	vector<double> predictCap(linkNum);

    ifstream RidgeParaFile;
    RidgeParaFile.open(argv[8]);
    for (int i = 0; i < featureNum; i++) {
        RidgeParaFile >> alpha[i];
    }
    RidgeParaFile.close();

	int source, sink;
	vector<vector<int>> graph(linkNum, vector<int>(2, 0));      // graph information (s,t)
	/* read graph */
	ostringstream oss1;
	oss1 << argv[1];
	string graphfileName = oss1.str();
	ifstream graphfile;
	graphfile.open(graphfileName);
	graphfile >> source;
	graphfile >> sink;
	source = 0;
	sink = nodeNum - 1;
  //sink = 31;
	for (int i = 0; i < linkNum; i++) {
		for (int j = 0; j < 2; j++) {
			graphfile >> graph[i][j];
		}
	}
	graphfile.close();


	ostringstream oss2;
	oss2 << argv[2];
	string trainfileName = oss2.str();
	ostringstream oss3;
	oss3 << argv[3];
	string testfileName = oss3.str();
	ostringstream oss4;
	oss4 << argv[4];
	string prefileName = oss4.str();
	ostringstream oss5;
	oss5 << argv[5];
	string runtimefileName = oss5.str();

	time_t startTime, endTime;
	startTime = time(NULL);

	vector<double> regretRec(featureNum, 0.1);
	int iteration = atoi(argv[6]);
	std::cout << "==================================================== training ====================================================" << endl;
	for (int T = 0; T < iteration; T++) {
		for (int k = 0; k < featureNum; k++) {
			ifstream trainfile;
			trainfile.open(trainfileName);
			for (int benN = 0; benN < benchmarkNum; benN++) {

				/* read data */
				for (int i = 0; i < linkNum; i++) {
					trainfile >> benchmarkId;
					for (int j = 0; j < featureNum; j++) {
						trainfile >> data[i][j];
					}
					trainfile >> realCap[i];
				}

				vector<double> A(linkNum);
				vector<double> B(linkNum);
				for (int i = 0; i < linkNum; i++) {
					A[i] = data[i][k];
					B[i] = 0;
					for (int d = 0; d < featureNum; d++) {
						B[i] = B[i] + data[i][d] * alpha[d];
					}
					B[i] = B[i] - data[i][k] * alpha[k];
				}

				vector<PiecewiseLinearFunction> temp;
				temp.resize(linkNum);
				for (int i = 0; i < linkNum; i++) {
					temp[i].assign(-INF, INF, A[i], B[i], 0);
				}
				double newSt = Preprocessing(temp);
				temp.clear();

				for (int i = 0; i < linkNum; i++) {
					addEdge(graph[i][0], graph[i][1], realCap[i], realCap[i], newSt, INF, A[i], B[i]);
				}

				double realValue = De_maxFlow(source, sink);
//				std::cout << "realValue: " << realValue << endl;

				PiecewiseLinearFunction predictValue = maxFlow(source, sink);
				if (newSt != -INF) {
					predictValue.start_x.insert(predictValue.start_x.begin(), -INF);
					predictValue.end_x.insert(predictValue.end_x.begin(), newSt);
					predictValue.a.insert(predictValue.a.begin(), 0);
					predictValue.b.insert(predictValue.b.begin(), 0);
					predictValue.id.insert(predictValue.id.begin(), 0);
				}
//				std::cout << "origin: " << endl;
				predictValue.merge();
//				predictValue.output();
				for (int predictValueCnt = 0; predictValueCnt < predictValue.a.size(); predictValueCnt++) {
					vector<double> alphaTemp(featureNum, 1);
					for (int h = 0; h < nodes; h++) {
						g[h].clear();
					}
					for (int h = 0; h < featureNum; h++) {
						alphaTemp[h] = alpha[h];
					}
					if (predictValueCnt == 0) {
						alphaTemp[k] = predictValue.end_x[predictValueCnt] - 1;
					}
					else if (predictValueCnt == 1 && predictValue.a.size() > 2) {
						alphaTemp[k] = predictValue.end_x[predictValueCnt];
					}
					else if (predictValueCnt == predictValue.start_x.size() - 1) {
						alphaTemp[k] = predictValue.start_x[predictValueCnt] + 1;
					}
					else {
						alphaTemp[k] = (predictValue.start_x[predictValueCnt] + predictValue.end_x[predictValueCnt]) / 2;
					}
					//                    std::cout << alphaTemp[k] << endl;
					for (int linkCnt = 0; linkCnt < linkNum; linkCnt++) {
						double predValTemp = 0;
						for (int j = 0; j < featureNum; j++) {
							predValTemp = predValTemp + alphaTemp[j] * data[linkCnt][j];
						}
						//                        std::cout << predValTemp << endl;
						addEdge(graph[linkCnt][0], graph[linkCnt][1], predValTemp, realCap[linkCnt], predictValue.start_x[predictValueCnt], predictValue.end_x[predictValueCnt], A[linkCnt], B[linkCnt]);
						//                        std::cout << graph[linkCnt][0]  << " " << graph[linkCnt][1] << " " << realCap[linkCnt] << endl;
					}
					double preObjTemp = Pre_2S_maxFlow(source, sink);
					//                    for (int h = 0; h < nodes; h++){
					//                        for (int hI = 0; hI < g[h].size(); hI++){
					//                            if (g[h][hI].from < g[h][hI].to){
					//                                std::cout << g[h][hI].from << " -> " << g[h][hI].to << " trueCap: " << g[h][hI].trueCap << " ";
					//                                g[h][hI].flow.output();
					//                            }
					//                        }
					//                    }
					//                    std::cout << preObjTemp << endl;
					predictValue.id[predictValueCnt] = preObjTemp;
				}
				//std::cout << "after: " << endl;
				//predictValue.output();
				predictValue.computeRegret(realValue);
				predictValue.merge();
				ans.push_back(predictValue);


				for (int h = 0; h < nodes; h++) {
					g[h].clear();
				}
			}

			PiecewiseLinearFunction res;
			res.assign(-INF, INF, 0, 0, 0);
			int g = 0;
			while (g < ans.size()) {
				res = Plus(res, ans[g]);
				g++;
			}
			std::cout << "=============================== the " << k << "th alpha done================================" << endl;
			//res.output();
			res.merge();

			vector<double>::iterator smallest = min_element(begin(res.id), end(res.id));
			int index = distance(begin(res.id), smallest);
			double Rtemp = *smallest;
			std::cout << "minimum regret = " << Rtemp << " at position " << index << endl;

			if (Rtemp < regret) {
				regret = Rtemp;
				if (index == 0) {
					alpha[k] = res.end_x[index] - 1;
				}
				else if (index == res.start_x.size() - 1) {
					alpha[k] = res.start_x[index] + 1;
				}
				else {
					alpha[k] = (res.start_x[index] + res.end_x[index]) / 2;
				}
				std::cout << "alpha[" << k << "]: " << alpha[k] << endl;
			}
			else {
				std::cout << "alpha[" << k << "]: " << alpha[k] << endl;
			}
			std::cout << "minimum regret = " << regret << endl;
			regretRec[k] = regret;

			if (k == (featureNum - 1)) {
				finalRes = res;
			}


			ans.clear();
			trainfile.close();
		}
		if (regretRec[0] - regretRec[featureNum - 1] < 1)
			break;
	}

	//std::cout << "================================================= training done ==================================================" << endl;
	vector<double>::iterator finalSmallest = min_element(begin(finalRes.id), end(finalRes.id));
	double totalRegret = *finalSmallest / benchmarkNum;
	//std::cout << "trainRegret = " << totalRegret << endl;
	endTime = time(NULL);
	double diff = difftime(endTime, startTime);
	cout << "Time = " << diff << endl;
	ofstream runtimeFile;
	runtimeFile.open(runtimefileName);
	runtimeFile << "Time = " << diff << endl;
	runtimeFile.close();

	/* testing */
	//std::cout << endl;
	//std::cout << "==================================================== testing ====================================================" << endl;
	ifstream testFile;
	//    testFile.open("test_GEANT.txt");
	testFile.open(testfileName);
	ofstream preFile;
	//    prefile.open(argv[4]);
	preFile.open(prefileName);
	for (int testT = 0; testT < testNum; testT++) {
		/* read data */
		for (int i = 0; i < linkNum; i++) {
			testFile >> testmarkId;
			for (int j = 0; j < featureNum; j++) {
				testFile >> testData[i][j];
			}
			testFile >> realCap[i];
		}

		for (int i = 0; i < linkNum; i++) {
			double temp2 = 0;
			for (int j = 0; j < featureNum; j++) {
				temp2 = temp2 + alpha[j] * testData[i][j];
			}
			predictCap[i] = temp2;
			preFile << testT << " " << realCap[i] << " " << predictCap[i] << endl;
		}
	}
	testFile.close();
	preFile.close();

	//std::cout << "====================================================  2S ====================================================" << endl;

	testFile.open(testfileName);
	ofstream resFile;
	ostringstream oss7;
	oss7 << argv[7];
	string resfileName = oss7.str();
	resFile.open(resfileName);

	double testRegret = 0;
	double totalFea = 0;
	double totalP = 0;
	double TOV = 0;
	for (int benN = 0; benN < testNum; benN++) {

		totalPath = 0;
		infeaPath = 0;
		violate = 0;
		TR = 0;
		PR = 0;
		PW = 0;
		/* read data */
		for (int i = 0; i < linkNum; i++) {
			testFile >> benchmarkId;
			for (int j = 0; j < featureNum; j++) {
				testFile >> testData[i][j];
			}
			testFile >> realCap[i];
		}

		for (int i = 0; i < linkNum; i++) {
			double temp2 = 0;
			for (int j = 0; j < featureNum; j++) {
				temp2 = temp2 + alpha[j] * testData[i][j];
			}
			predictCap[i] = temp2;
		}


		for (int i = 0; i < linkNum; i++) {
			Pre_addEdge(graph[i][0], graph[i][1], realCap[i], predictCap[i]);
		}
		vector<Edge> gtemp[maxnodes];
		for (int k1 = 0; k1 < nodes; k1++) {
			gtemp[k1] = g[k1];
		}
		double realT = De_maxFlow(source, sink);
		for (int k2 = 0; k2 < nodes; k2++) {
			for (int l = 0; l < g[k2].size(); l++) {
				g[k2][l].cap = gtemp[k2][l].cap;
			}
		}
		double preT = Pre_2S_maxFlow(source, sink);
		for (int h = 0; h < nodes; h++) {
			g[h].clear();
		}

		testRegret = testRegret + abs(realT - preT);
    
		TOV = TOV + realT;
		feaPath = totalPath - infeaPath;
		totalP = totalP + totalPath;
		totalFea = totalFea + feaPath;
	}
	//std::std::cout << "================================================= 2S done ==================================================" << endl;
    std::cout << testRegret / testNum << " " << TOV / testNum << endl;
	resFile << testRegret / testNum << " " << TOV / testNum << endl;
	//acceptRatio = double(totalFea / double(totalP));
	//std::cout << "totalP: " << totalP << " totalFea: " << totalFea << endl;
	//std::cout << "acceptRatio: " << acceptRatio << endl;
	testFile.close();
	//}

}


