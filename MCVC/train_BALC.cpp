// MinWeightVertexCover.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>  
#include <time.h>
#include <set>
#include "PiecewiseLinearFunction.h"

using namespace std;

class MWVC {
public:
	int edgeNum, nodeNum;
	vector<double> realCost;
	vector<PiecewiseLinearFunction> cost;
	vector<vector<bool>> edge;

	MWVC(int eN, int nN, vector<double> rC, vector<PiecewiseLinearFunction> c, vector<vector<bool>> e) {
		edgeNum = eN;
		nodeNum = nN;
		realCost = rC;
		cost = c;
		edge = e;
	}

};

class Result {
public:
	double totalWeight;
	vector<int> finalSelected;
};

Result De_MWVC(MWVC G, int cnt, vector<int> selected, int haveToCover) {

	Result totalResult;
	cnt = cnt + 1;
	if (cnt == G.nodeNum || haveToCover == 0) {
		double totalWeight = 0;

		/*for (int i = 0; i < G.nodeNum; i++) {
			for (int j = 0; j < G.nodeNum; j++) {
				if (G.edge[i][j] == 1) {
					totalWeight = INFINITY;
				}
			}
		}*/
		if (haveToCover > 0) {
			totalWeight = INFINITY;
		}

		if (totalWeight == 0) {
			for (int i = 0; i < selected.size(); i++) {
				totalWeight = totalWeight + G.realCost[selected[i]];
			}
		}

		totalResult.totalWeight = totalWeight;
		totalResult.finalSelected = selected;

		/*for (int i = 0; i < selected.size(); i++) {
			cout << selected[i] << " ";
		}
		cout << ": " << totalWeight << endl;*/

		return totalResult;
	}

	else {
		vector<int> selected1 = selected;
		selected1.push_back(cnt);
		MWVC G1 = G;
		int haveToCover1 = haveToCover;
		for (int i = 0; i < G1.nodeNum; i++) {
			if (G1.edge[i][cnt] == 1) {
				G1.edge[i][cnt] = 0;
				G1.edge[cnt][i] = 0;
				haveToCover1 = haveToCover1 - 1;
			}
		}

		vector<int> selected2 = selected;
		MWVC G2 = G;
		int haveToCover2 = haveToCover;

		Result result1 = De_MWVC(G1, cnt, selected1, haveToCover1);
		Result result2 = De_MWVC(G2, cnt, selected2, haveToCover2);
		Result result;
		//double weight = min(result1.totalWeight, result2.totalWeight);
		//double weight = min(De_MWVC(G1, cnt, selected1), De_MWVC(G2, cnt, selected2));
		if (result1.totalWeight > result2.totalWeight) {
			result = result2;
		}
		else {
			result = result1;
		}

		return result;
	}
}

PiecewiseLinearFunction CorrectionalCostTrain_MWVC(MWVC G, int cnt, vector<int> selected, int src, int des, bool useCorrectional, int haveToCover) {

	PiecewiseLinearFunction totalResult;
	cnt = cnt + 1;
	if (cnt == G.nodeNum || haveToCover == 0) {
		PiecewiseLinearFunction totalWeight;
		totalWeight.assign(-INF, INF, 0, 0, 0);

		/*for (int i = 0; i < G.nodeNum; i++) {
			for (int j = 0; j < G.nodeNum; j++) {
				if (G.edge[i][j] == 1) {
					totalWeight.a[0] = INF;
					totalWeight.b[0] = INF;
					totalWeight.id[0] = INF;
				}
			}
		}*/
		if (haveToCover > 0) {
			totalWeight.a[0] = INF;
			totalWeight.b[0] = INF;
			totalWeight.id[0] = INF;
		}

		if (totalWeight.a[0] == 0) {
			for (int i = 0; i < selected.size(); i++) {
				totalWeight = Plus(totalWeight, G.cost[selected[i]]);
			}

			if (useCorrectional) {
				bool checkSrc = count(selected.begin(), selected.end(), src);
				bool checkDes = count(selected.begin(), selected.end(), des);
				if (checkSrc == 0 && checkDes == 0) {
					//totalWeight = Plus(totalWeight, G.cost[src]);
					//totalWeight = Plus(totalWeight, G.cost[des]);
					totalWeight.id[0] = totalWeight.id[0] + G.cost[src].id[0] + G.cost[des].id[0];
				}
			}
		}

		totalResult = totalWeight;

		/*if (totalWeight.id[0] != INF) {
			for (int i = 0; i < selected.size(); i++) {
				cout << selected[i] << " ";
			}
			cout << ": " << totalWeight.id[0] << endl;
			//totalWeight.output();
		}*/


		return totalResult;
	}

	else {
		vector<int> selected1 = selected;
		selected1.push_back(cnt);
		MWVC G1 = G;
		int haveToCover1 = haveToCover;
		for (int i = 0; i < G1.nodeNum; i++) {
			if (G1.edge[i][cnt] == 1) {
				G1.edge[i][cnt] = 0;
				G1.edge[cnt][i] = 0;
				haveToCover1 = haveToCover1 - 1;
			}
		}

		vector<int> selected2 = selected;
		MWVC G2 = G;
		int haveToCover2 = haveToCover;

		PiecewiseLinearFunction result1 = CorrectionalCostTrain_MWVC(G1, cnt, selected1, src, des, useCorrectional, haveToCover1);
		PiecewiseLinearFunction result2 = CorrectionalCostTrain_MWVC(G2, cnt, selected2, src, des, useCorrectional, haveToCover2);
		PiecewiseLinearFunction result;
		//double weight = min(result1.totalWeight, result2.totalWeight);
		//double weight = min(CostTrain_MWVC(G1, cnt, selected1), CostTrain_MWVC(G2, cnt, selected2));
		if (result1.id[0] != INF && result2.id[0] != INF) {
			result = InfCombination(result1, result2);
			result.merge();
			/*cout << "result" << endl;
			result.output();
			cout << endl;*/
		}
		else if (result1.id[0] == INF) {
			result = result2;
		}
		else {
			result = result1;
		}

		return result;
	}
}

PiecewiseLinearFunction CostTrain_MWVC(MWVC G, int cnt, vector<int> selected) {

	PiecewiseLinearFunction totalResult;
	cnt = cnt + 1;
	if (cnt == G.nodeNum) {
		PiecewiseLinearFunction totalWeight;
		totalWeight.assign(-INF, INF, 0, 0, 0);

		for (int i = 0; i < G.nodeNum; i++) {
			for (int j = 0; j < G.nodeNum; j++) {
				if (G.edge[i][j] == 1) {
					totalWeight.a[0] = INF;
					totalWeight.b[0] = INF;
					totalWeight.id[0] = INF;
				}
			}
		}

		if (totalWeight.a[0] == 0) {
			for (int i = 0; i < selected.size(); i++) {
				totalWeight = Plus(totalWeight, G.cost[selected[i]]);
			}
		}

		totalResult = totalWeight;

		/*for (int i = 0; i < selected.size(); i++) {
			cout << selected[i] << " ";
		}
		cout << ": " << totalWeight.id[0] << endl;*/
		//totalWeight.output();

		return totalResult;
	}

	else {
		vector<int> selected1 = selected;
		selected1.push_back(cnt);
		MWVC G1 = G;
		for (int i = 0; i < G1.nodeNum; i++) {
			G1.edge[i][cnt] = 0;
			G1.edge[cnt][i] = 0;
		}

		vector<int> selected2 = selected;
		MWVC G2 = G;

		PiecewiseLinearFunction result1 = CostTrain_MWVC(G1, cnt, selected1);
		PiecewiseLinearFunction result2 = CostTrain_MWVC(G2, cnt, selected2);
		PiecewiseLinearFunction result;
		//double weight = min(result1.totalWeight, result2.totalWeight);
		//double weight = min(CostTrain_MWVC(G1, cnt, selected1), CostTrain_MWVC(G2, cnt, selected2));
		if (result1.id[0] != INF && result2.id[0] != INF) {
			result = InfCombination(result1, result2);
			result.merge();
			/*cout << "result" << endl;
			result.output();
			cout << endl;*/
		}
		else if (result1.id[0] == INF) {
			result = result2;
		}
		else {
			result = result1;
		}

		return result;
	}
}

double CostTest_MWVC(vector<int> selected, vector<double> realCost) {
	double totalWeight = 0;
	for (int i = 0; i < selected.size(); i++) {
		totalWeight = totalWeight + realCost[selected[i]];
	}
	return totalWeight;
}

double Preprocessing(vector<PiecewiseLinearFunction> v) {
	double start_x_new = -INF;
	//if (v[0].a[0] != 0 || v[1].a[0] != 0) {
	vector<double> vTemp;
	int size = v.size();
	vTemp.resize(size);
	for (int i = 0; i < size; i++) {
		double start = -INF;
		if (v[i].a[0] != 0) {
			start = -(v[i].b[0] / v[i].a[0]);
		}
		vTemp[i] = start;
	}
	vector<double>::iterator newStart = max_element(begin(vTemp), end(vTemp));
	start_x_new = *newStart;
	//}
	return start_x_new;
}

int main(int argc, char* argv[])
{
	if (argc != 10) {
		std::cerr << "usage: ./train_BALC [graph] [train edge data] [test edge data] [train cost data] [test cost data] [pre edge file] [pre cost file] [iteration] [runtime file]" << std::endl;
		exit(0);
	}

	int nodeNum = 11, edgeNum = 34;
	int featureNum = 8;
	vector<double> alpha(featureNum, 1);
	vector<double> beta(featureNum, 1);
	vector<double> alphaRegretRec(featureNum, 0);
	vector<double> betaRegretRec(featureNum, 0);
	double alphaRegret = INF;
	double betaRegret = INF;

	int trainNum = 70;
	int testNum = 30;

	ifstream infile;
	ifstream infile2;
	double benchmarkId;

    string graphType = argv[1];
    if(graphType=="PDH"){
        infile.open("edge_PDH.txt");
        nodeNum = 11, edgeNum = 34;
    }
    else if(graphType=="ABILENE"){
        infile.open("edge_ABILENE.txt");
        nodeNum = 12, edgeNum = 15;
    }
	vector<vector<int>> edgeTemp(edgeNum, vector<int>(2, 0));
	for (int i = 0; i < edgeNum; i++) {
		infile >> edgeTemp[i][0];
		infile >> edgeTemp[i][1];
	}
	infile.close();

	std::cout << "==================================================== training ====================================================" << endl;
	time_t Tstart, Tend;
	Tstart = time(NULL);

	int iteration = atoi(argv[8]);
	int convergenIter = 0;
	for (int T = 0; T < iteration; T++) {
		convergenIter = convergenIter + 1;
		vector<PiecewiseLinearFunction> ans;
		PiecewiseLinearFunction finalRes;
		int cntAB = 0;

		for (; cntAB < featureNum; cntAB++) {
			infile.open(argv[2]);
			infile2.open(argv[4]);
			for (int benN = 0; benN < trainNum; benN++) {

				//cout << benN << " ";
				/* predict cost */
				/* compute the cost value */
				vector<vector<double>> trainCostFeature(nodeNum, vector<double>(featureNum, 0));
				vector<double> realCost(nodeNum, 0);
				vector<double> preCost(nodeNum, 0);
				for (int j = 0; j < nodeNum; j++) {
					double benchmarkId;
					infile2 >> benchmarkId;
					double predictCost = 0.0;
					for (int k = 0; k < featureNum; k++) {
						infile2 >> trainCostFeature[j][k];
						predictCost += trainCostFeature[j][k] * beta[k];
					}
					preCost[j] = predictCost;
					// std::cout << std::endl;
					infile2 >> realCost[j];
					//cout << j << " " << realCost[j] << " " << preCost[j] << std::endl;
				}

				vector<PiecewiseLinearFunction> costTemp;
				costTemp.resize(nodeNum);
				for (int i = 0; i < nodeNum; i++) {
					costTemp[i].assign(-INF, INF, 0, 0, 0);
				}

				/* train edge */
				vector<PiecewiseLinearFunction> flow;
				flow.resize(edgeNum);
				vector<vector<double>> trainFlowFeature(edgeNum, vector<double>(featureNum, 0));
				vector<double> trainFlow(edgeNum, 0);

				for (int i = 0; i < edgeNum; i++) {
					infile >> benchmarkId;
					for (int j = 0; j < featureNum; j++) {
						infile >> trainFlowFeature[i][j];
					}
					infile >> trainFlow[i];
				}

				vector<double>::iterator smallest = min_element(begin(trainFlow), end(trainFlow));
				int realIndex = distance(begin(trainFlow), smallest);
				//cout << "real index: " << realIndex << endl;

				vector<vector<bool>> realEdge(nodeNum, vector<bool>(nodeNum, 0));
				for (int i = 0; i < edgeNum; i++) {
					if (i != realIndex) {
						int src = edgeTemp[i][0];
						int des = edgeTemp[i][1];
						realEdge[src][des] = 1;
						realEdge[des][src] = 1;
					}
				}

				MWVC GTrue(edgeNum, nodeNum, realCost, costTemp, realEdge);
				int HTCTrue = edgeNum - 1;
				vector<int> STrue;
				Result realValue = De_MWVC(GTrue, -1, STrue, HTCTrue);
				STrue.clear();

				/* compute the piecewise function */
				PiecewiseLinearFunction predictValue;

				vector<double> A(edgeNum, 0);
				vector<double> B(edgeNum, 0);
				for (int i = 0; i < edgeNum; i++) {
					A[i] = trainFlowFeature[i][cntAB];
					B[i] = 0;
					for (int d = 0; d < featureNum; d++) {
						B[i] = B[i] + trainFlowFeature[i][d] * alpha[d];
					}
					B[i] = B[i] - trainFlowFeature[i][cntAB] * alpha[cntAB];
				}

				for (int i = 0; i < edgeNum; i++) {
					flow[i].assign(-INF, INF, A[i], B[i], i);
					//flow[i].output();
				}

				PiecewiseLinearFunction infOfFlow = MinCompare(flow);
				//infOfFlow.output();

				double newSt = Preprocessing(flow);
				//cout << "newSt: " << newSt << endl;
				set<double> transPointTemp;
				transPointTemp.insert(newSt);
				for (int i = 0; i < infOfFlow.start_x.size(); i++) {
					transPointTemp.insert(infOfFlow.start_x[i]);
				}

				vector<double> transPoint;
				for (set<double>::iterator it = transPointTemp.begin(); it != transPointTemp.end(); ++it) {
					if (*it >= newSt) {
						//cout << *it << " ";
						double pt = *it;
						transPoint.push_back(pt);
					}
				}
				//cout << endl;
				if (transPoint[0] > -INF) {
					transPoint.insert(transPoint.begin(), -INF);
				}
				else if (transPoint[0] < -INF) {
					transPoint.erase(transPoint.begin());
				}

				for (int i = 0; i < transPoint.size(); i++) {
					//cout << transPoint[i] << " ";
					double start_x, end_x;
					double currentPt;
					if (transPoint.size() == 1) {
						currentPt = 1;
						start_x = -INF;
						end_x = INF;
					}
					else {
						if (i == 0) {
							currentPt = transPoint[i + 1] - 1;
							start_x = -INF;
							end_x = transPoint[i + 1];
							predictValue.assign(start_x, end_x, 0, 0, INF);
							continue;
						}
						else if (i == transPoint.size() - 1) {
							currentPt = transPoint[i] + 1;
							start_x = transPoint[i];
							end_x = INF;
						}
						else {
							currentPt = (transPoint[i + 1] + transPoint[i]) / 2;
							start_x = transPoint[i];
							end_x = transPoint[i + 1];
						}
					}
					//cout << "currentPt: " << currentPt << endl;

					vector<double> flowTemp(edgeNum, 0);
					for (int i = 0; i < edgeNum; i++) {
						flowTemp[i] = A[i] * currentPt + B[i];
						//cout << flowTemp[i] << endl;
					}
					//cout << endl;

					vector<double>::iterator smallest = min_element(begin(flowTemp), end(flowTemp));
					int preIndex = distance(begin(flowTemp), smallest);
					//cout << "pre index: " << preIndex << endl;

					vector<vector<bool>> edge(nodeNum, vector<bool>(nodeNum, 0));
					for (int i = 0; i < edgeNum; i++) {
						if (i != preIndex) {
							int src = edgeTemp[i][0];
							int des = edgeTemp[i][1];
							edge[src][des] = 1;
							edge[des][src] = 1;
						}
					}

					MWVC G(edgeNum, nodeNum, preCost, costTemp, edge);
					int HTC = edgeNum - 1;
					vector<int> S;
					Result weight = De_MWVC(G, -1, S, HTC);
					double preValue = CostTest_MWVC(weight.finalSelected, realCost);

					/*for (int i = 0; i < weight.finalSelected.size(); i++) {
						cout << weight.finalSelected[i] << " ";
					}
					cout << endl;*/

					if (realIndex != preIndex) {
						bool checkSrc = count(weight.finalSelected.begin(), weight.finalSelected.end(), edgeTemp[preIndex][0]);
						bool checkDes = count(weight.finalSelected.begin(), weight.finalSelected.end(), edgeTemp[preIndex][1]);
						if (checkSrc == 0 && checkDes == 0) {
							preValue = preValue + realCost[edgeTemp[preIndex][0]] + realCost[edgeTemp[preIndex][1]];
						}
					}

					predictValue.assign(start_x, end_x, 0, 0, preValue);
				}

				predictValue.merge();
				//predictValue.output();
				predictValue.computeRegret(realValue.totalWeight);
//				std::cout << realValue.totalWeight << " ";
				//predictValue.output();
				//cout << endl;
				ans.push_back(predictValue);
			}

			PiecewiseLinearFunction res;
			res.assign(-INF, INF, 0, 0, 0);
			int g = 0;
			while (g < ans.size()) {
				res = Plus(res, ans[g]);
				g++;
			}
			std::cout << endl;
			std::cout << "=============================== the " << cntAB << "th alpha done================================" << endl;
			res.merge();
			//res.output();

			vector<double>::iterator smallest = min_element(begin(res.id), end(res.id));
			int index = distance(begin(res.id), smallest);
			double Rtemp = *smallest;
			std::cout << "minimum regret = " << Rtemp << " at position " << index << endl;

			if (Rtemp < alphaRegret) {
				alphaRegret = Rtemp;
				if (index == 0) {
					alpha[cntAB] = res.end_x[index] - 1;
				}
				else if (index == res.start_x.size() - 1) {
					alpha[cntAB] = res.start_x[index] + 1;
				}
				else {
					alpha[cntAB] = (res.start_x[index] + res.end_x[index]) / 2;
				}
				std::cout << "alpha[" << cntAB << "]: " << alpha[cntAB] << endl;
			}
			else {
				std::cout << "alpha[" << cntAB << "]: " << alpha[cntAB] << endl;
			}
			std::cout << "minimum regret = " << alphaRegret << endl;
			alphaRegretRec[cntAB] = alphaRegret;

			if (cntAB == (featureNum - 1)) {
				finalRes = res;
			}

			ans.clear();
			infile.close();
			infile2.close();

		}
		double Tmid = time(NULL);
		double diffMid = difftime(Tmid, Tstart);
		std::cout << "Time = " << diffMid << endl;


		ans.clear();
		finalRes.start_x.clear();
		finalRes.end_x.clear();
		finalRes.a.clear();
		finalRes.b.clear();
		finalRes.id.clear();

		for (; cntAB < 2 * featureNum; cntAB++) {

			infile.open(argv[4]);
			infile2.open(argv[2]);
			for (int benN = 0; benN < trainNum; benN++) {

				/* predict edge */
				vector<vector<double>> trainFlowFeature(edgeNum, vector<double>(featureNum, 0));
				vector<double> realFlow(edgeNum, 0);
				vector<double> preFlow(edgeNum, 0);
				for (int j = 0; j < edgeNum; j++) {
					double benchmarkId;
					infile2 >> benchmarkId;
					double predictflow = 0.0;
					for (int k = 0; k < featureNum; k++) {
						infile2 >> trainFlowFeature[j][k];
						predictflow += trainFlowFeature[j][k] * alpha[k];
					}
					preFlow[j] = predictflow;
					// std::cout << std::endl;
					infile2 >> realFlow[j];
					//cout << j << " " << realFlow[j] << " " << preFlow[j] << std::endl;
				}

				vector<double>::iterator smallest = min_element(begin(realFlow), end(realFlow));
				int realIndex = distance(begin(realFlow), smallest);
				//cout << "real index: " << realIndex << endl;

				vector<vector<bool>> realEdge(nodeNum, vector<bool>(nodeNum, 0));
				for (int i = 0; i < edgeNum; i++) {
					if (i != realIndex) {
						int src = edgeTemp[i][0];
						int des = edgeTemp[i][1];
						realEdge[src][des] = 1;
						realEdge[des][src] = 1;
					}
				}

				vector<double>::iterator preSmallest = min_element(begin(preFlow), end(preFlow));
				int preIndex = distance(begin(preFlow), preSmallest);
				//cout << "pre index: " << preIndex << endl;

				int delSrc = edgeTemp[preIndex][0];
				int delDes = edgeTemp[preIndex][1];
				bool useCorrectional = false;
				if (realIndex != preIndex) {
					useCorrectional = true;
				}

				vector<vector<bool>> preEdge(nodeNum, vector<bool>(nodeNum, 0));
				for (int i = 0; i < edgeNum; i++) {
					if (i != preIndex) {
						int src = edgeTemp[i][0];
						int des = edgeTemp[i][1];
						preEdge[src][des] = 1;
						preEdge[des][src] = 1;
					}
				}

				/* train cost */
				vector<PiecewiseLinearFunction> cost;
				cost.resize(nodeNum);
				vector<vector<double>> trainFeature(nodeNum, vector<double>(featureNum, 0));
				vector<double> trainCost(nodeNum, 0);

				for (int i = 0; i < nodeNum; i++) {
					infile >> benchmarkId;
					for (int j = 0; j < featureNum; j++) {
						infile >> trainFeature[i][j];
					}
					infile >> trainCost[i];
				}

				vector<double> A(nodeNum, 0);
				vector<double> B(nodeNum, 0);
				for (int i = 0; i < nodeNum; i++) {
					A[i] = trainFeature[i][cntAB - featureNum];
					B[i] = 0;
					for (int d = 0; d < featureNum; d++) {
						B[i] = B[i] + trainFeature[i][d] * beta[d];
					}
					B[i] = B[i] - trainFeature[i][cntAB - featureNum] * beta[cntAB - featureNum];
				}

				for (int i = 0; i < nodeNum; i++) {
					cost[i].assign(-INF, INF, A[i], B[i], trainCost[i]);
					//cost[i].output();
				}

				MWVC GTrue(edgeNum, nodeNum, trainCost, cost, realEdge);
				int HTCtrue = edgeNum - 1;
				vector<int> S;
				Result realValue = De_MWVC(GTrue, -1, S, HTCtrue);

				S.clear();
				MWVC Gpre(edgeNum, nodeNum, trainCost, cost, preEdge);
				int HTCpre = edgeNum - 1;
				PiecewiseLinearFunction predictValue = CorrectionalCostTrain_MWVC(Gpre, -1, S, delSrc, delDes, useCorrectional, HTCpre);

				//cout << "realWeight: " << realValue.totalWeight << endl;
//				std::cout << realValue.totalWeight << " ";
				/*cout << "finalSelected: ";
				for (int i = 0; i < realWeight.finalSelected.size(); i++) {
					cout << realWeight.finalSelected[i] << " ";
				}
				cout << endl;
				preWeight.output();
				cout << endl;*/
				//predictValue.output();
				//cout << endl;
				predictValue.computeRegret(realValue.totalWeight);
				//predictValue.output();
				predictValue.merge();
				ans.push_back(predictValue);
			}

			PiecewiseLinearFunction res;
			res.assign(-INF, INF, 0, 0, 0);
			int g = 0;
			while (g < ans.size()) {
				res = Plus(res, ans[g]);
				g++;
			}
			std::cout << endl;
			std::cout << "=============================== the " << cntAB - featureNum << "th beta done================================" << endl;
			//res.output();
			res.merge();

			vector<double>::iterator smallest = min_element(begin(res.id), end(res.id));
			int index = distance(begin(res.id), smallest);
			double Rtemp = *smallest;
			std::cout << "minimum regret = " << Rtemp << " at position " << index << endl;

			if (Rtemp < betaRegret) {
				betaRegret = Rtemp;
				if (index == 0) {
					beta[cntAB - featureNum] = res.end_x[index] - 1;
				}
				else if (index == res.start_x.size() - 1) {
					beta[cntAB - featureNum] = res.start_x[index] + 1;
				}
				else {
					beta[cntAB - featureNum] = (res.start_x[index] + res.end_x[index]) / 2;
				}
				std::cout << "beta[" << cntAB - featureNum << "]: " << beta[cntAB - featureNum] << endl;
			}
			else {
				std::cout << "beta[" << cntAB - featureNum << "]: " << beta[cntAB - featureNum] << endl;
			}
			std::cout << "minimum regret = " << betaRegret << endl;
			betaRegretRec[cntAB - featureNum] = betaRegret;

			if (cntAB == (2 * featureNum - 1)) {
				finalRes = res;
			}

			ans.clear();
			infile.close();
			infile2.close();
		}


		if (alphaRegretRec[featureNum - 1] == alphaRegretRec[0] && betaRegretRec[featureNum - 1] == betaRegretRec[0])
			break;

	}

	std::cout << "================================================= training done ==================================================" << endl;
	//vector<double>::iterator finalSmallest = min_element(begin(finalRes.id), end(finalRes.id));
	//double totalRegret = *finalSmallest / trainNum;
	//cout << "trainRegret = " << totalRegret << endl;

	Tend = time(NULL);
	double diff = difftime(Tend, Tstart);
	std::cout << "Time = " << diff << endl;
	ofstream runfile;
	runfile.open(argv[9]);
	runfile << convergenIter << " " << diff << "s" << endl;


	// read test data 
	// int testbenchmarkNum, testfeatureNum; 
	infile.open(argv[3]);
	infile2.open(argv[5]);
	double testFlowFeatures[testNum][edgeNum][featureNum];
	double testRealFlow[testNum][edgeNum];
	double testCostFeatures[testNum][nodeNum][featureNum];
	double testRealCost[testNum][nodeNum];
	ofstream ofile;
	ofstream ofile2;
	ofile.open(argv[6]);
	ofile2.open(argv[7]);
	for (int i = 0; i < testNum; i++) {
		for (int j = 0; j < edgeNum; j++) {
			double benchmarkId;
			infile >> benchmarkId;
			double predictedFlow = 0.0;
			for (int k = 0; k < featureNum; k++) {
				infile >> testFlowFeatures[i][j][k];
				predictedFlow += testFlowFeatures[i][j][k] * alpha[k];
			}
			infile >> testRealFlow[i][j];
			ofile << i << " " << testRealFlow[i][j] << " " << predictedFlow << std::endl;
		}

		for (int j = 0; j < nodeNum; j++) {
			double benchmarkId;
			infile2 >> benchmarkId;
			double predictedCost = 0.0;
			for (int k = 0; k < featureNum; k++) {
				infile2 >> testCostFeatures[i][j][k];
				predictedCost += testCostFeatures[i][j][k] * beta[k];
			}
			// std::cout << std::endl;
			infile2 >> testRealCost[i][j];
			ofile2 << i << " " << testRealCost[i][j] << " " << predictedCost << std::endl;
		}
	}
	//ofile << "Time = " << diff << endl;
	ofile.close();
	ofile2.close();
	infile.close();
	infile2.close();
}

