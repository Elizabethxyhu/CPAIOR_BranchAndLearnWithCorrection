#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>  
#include <time.h>
#include <set>
#include "PiecewiseLinearFunction.h"

using namespace std;

class MCVC {
public:
	int edgeNum, nodeNum;
	vector<double> realCost;
	vector<PiecewiseLinearFunction> cost;
	vector<vector<bool>> edge;

	MCVC(int eN, int nN, vector<double> rC, vector<PiecewiseLinearFunction> c, vector<vector<bool>> e) {
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

Result De_MCVC(MCVC G, int cnt, vector<int> selected) {

	Result totalResult;
	cnt = cnt + 1;
	if (cnt == G.nodeNum) {
		double totalWeight = 0;

		for (int i = 0; i < G.nodeNum; i++) {
			for (int j = 0; j < G.nodeNum; j++) {
				if (G.edge[i][j] == 1) {
					totalWeight = INFINITY;
				}
			}
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
		MCVC G1 = G;
		for (int i = 0; i < G1.nodeNum; i++) {
			G1.edge[i][cnt] = 0;
			G1.edge[cnt][i] = 0;
		}

		vector<int> selected2 = selected;
		MCVC G2 = G;

		Result result1 = De_MCVC(G1, cnt, selected1);
		Result result2 = De_MCVC(G2, cnt, selected2);
		Result result;
		//double weight = min(result1.totalWeight, result2.totalWeight);
		//double weight = min(De_MCVC(G1, cnt, selected1), De_MCVC(G2, cnt, selected2));
		if (result1.totalWeight > result2.totalWeight) {
			result = result2;
		}
		else {
			result = result1;
		}

		return result;
	}
}

PiecewiseLinearFunction CorrectionalCostTrain_MCVC(MCVC G, int cnt, vector<int> selected, int src, int des, bool useCorrectional) {

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
			
			if (useCorrectional) {
				bool checkSrc = count(selected.begin(), selected.end(), src);
				bool checkDes = count(selected.begin(), selected.end(), des);
				if (checkSrc == 0 && checkDes == 0) {
					totalWeight = Plus(totalWeight, G.cost[src]);
					totalWeight = Plus(totalWeight, G.cost[des]);
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
		MCVC G1 = G;
		for (int i = 0; i < G1.nodeNum; i++) {
			G1.edge[i][cnt] = 0;
			G1.edge[cnt][i] = 0;
		}

		vector<int> selected2 = selected;
		MCVC G2 = G;

		PiecewiseLinearFunction result1 = CorrectionalCostTrain_MCVC(G1, cnt, selected1, src, des, useCorrectional);
		PiecewiseLinearFunction result2 = CorrectionalCostTrain_MCVC(G2, cnt, selected2, src, des, useCorrectional);
		PiecewiseLinearFunction result;
		//double weight = min(result1.totalWeight, result2.totalWeight);
		//double weight = min(CostTrain_MCVC(G1, cnt, selected1), CostTrain_MCVC(G2, cnt, selected2));
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

PiecewiseLinearFunction CostTrain_MCVC(MCVC G, int cnt, vector<int> selected) {

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
		MCVC G1 = G;
		for (int i = 0; i < G1.nodeNum; i++) {
			G1.edge[i][cnt] = 0;
			G1.edge[cnt][i] = 0;
		}

		vector<int> selected2 = selected;
		MCVC G2 = G;

		PiecewiseLinearFunction result1 = CostTrain_MCVC(G1, cnt, selected1);
		PiecewiseLinearFunction result2 = CostTrain_MCVC(G2, cnt, selected2);
		PiecewiseLinearFunction result;
		//double weight = min(result1.totalWeight, result2.totalWeight);
		//double weight = min(CostTrain_MCVC(G1, cnt, selected1), CostTrain_MCVC(G2, cnt, selected2));
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

double CostTest_MCVC(vector<int> selected, vector<double> realCost) {
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

int main(int argc, char *argv[])
{
	if (argc != 5) {
		std::cerr << "usage: ./test [graph] [pre edge data] [pre cost data] [result file]" << std::endl;
		exit(0);
	}

	int nodeNum = 11, edgeNum = 34;
//    int nodeNum = 12, edgeNum = 15;
	
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

	// read test data 
	infile.open(argv[2]);
	infile2.open(argv[3]);
	double totalRegret = 0;
	double TOV = 0;
	for (int benN = 0; benN < testNum; benN++) {

		vector<double> preCost(nodeNum);
		vector<double> realCost(nodeNum);
		vector<double> preFlowTemp(edgeNum);
		vector<double> realFlowTemp(edgeNum);

		for (int j = 0; j < edgeNum; j++) {
			infile >> benchmarkId;
			infile >> realFlowTemp[j];
			infile >> preFlowTemp[j];
		}

		for (int j = 0; j < nodeNum; j++) {
			infile2 >> benchmarkId;
			infile2 >> realCost[j];
			infile2 >> preCost[j];
		}

		vector<PiecewiseLinearFunction> costTemp;
		costTemp.resize(nodeNum);
		for (int i = 0; i < nodeNum; i++) {
			costTemp[i].assign(-INF, INF, 0, 0, 0);
		}

		vector<double>::iterator realSmallest = min_element(begin(realFlowTemp), end(realFlowTemp));
		int realIndex = distance(begin(realFlowTemp), realSmallest);

		vector<vector<bool>> realEdge(nodeNum, vector<bool>(nodeNum, 0));
		for (int i = 0; i < edgeNum; i++) {
			if (i != realIndex) {
				int src = edgeTemp[i][0];
				int des = edgeTemp[i][1];
				realEdge[src][des] = 1;
				realEdge[des][src] = 1;
			}
		}

		vector<double>::iterator preSmallest = min_element(begin(preFlowTemp), end(preFlowTemp));
		int preIndex = distance(begin(preFlowTemp), preSmallest);

		vector<vector<bool>> preEdge(nodeNum, vector<bool>(nodeNum, 0));
		for (int i = 0; i < edgeNum; i++) {
			if (i != preIndex) {
				int src = edgeTemp[i][0];
				int des = edgeTemp[i][1];
				preEdge[src][des] = 1;
				preEdge[des][src] = 1;
			}
		}

		MCVC GTrue(edgeNum, nodeNum, realCost, costTemp, realEdge);
		vector<int> S;
		Result realValue = De_MCVC(GTrue, -1, S);
//		cout << realValue.totalWeight << " ";
		TOV = TOV + realValue.totalWeight;

		MCVC GPre(edgeNum, nodeNum, preCost, costTemp, preEdge);
		S.clear();
		Result preSelected = De_MCVC(GPre, -1, S);
		double preValue = CostTest_MCVC(preSelected.finalSelected, realCost);

		if (realIndex != preIndex) {
			bool checkSrc = count(preSelected.finalSelected.begin(), preSelected.finalSelected.end(), edgeTemp[preIndex][0]);
			bool checkDes = count(preSelected.finalSelected.begin(), preSelected.finalSelected.end(), edgeTemp[preIndex][1]);
			if (checkSrc == 0 && checkDes == 0) {
				preValue = preValue + realCost[edgeTemp[preIndex][0]] + realCost[edgeTemp[preIndex][1]];
			}
		}
//		cout << preValue << endl;

		double regret = abs(preValue - realValue.totalWeight);
		totalRegret = totalRegret + regret;
	}
	//cout << endl;
	totalRegret = totalRegret / testNum;
	cout << "avgTotalRegret: " << totalRegret << endl;
	ofstream ofile;
	ofile.open(argv[4]);
	double avgTOV = TOV / testNum;
	//ofile << avgTOV << endl;
	ofile << totalRegret << endl;
	ofile.close();


}

