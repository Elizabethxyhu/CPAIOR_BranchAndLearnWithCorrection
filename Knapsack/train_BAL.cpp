// KnapsackWithWeight.cpp : This file contains the 'main' function. Program executiocntbegins and ends there.
//
#include <iostream>
#include <algorithm>
#include <numeric>
#include <vector>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include "PiecewiseLinearFunction.h"

using namespace std;
double remainCap = 0;

class KS {
public:
	int itemNum;
	double remainCap;
	PiecewiseLinearFunction preRemainCap;
	vector<double> realWeight;
	vector<double> value;
	vector<PiecewiseLinearFunction> weight;

	KS(int itemN, double Cap, vector<double> rW, vector<PiecewiseLinearFunction> w, vector<double> v) {
		itemNum = itemN;
		remainCap = Cap;
		realWeight = rW;
		value = v;
		weight = w;
		preRemainCap.assign(-INF, INF, 0, remainCap, remainCap);
	}

};

class Result {
public:
	double totalVal;
	vector<int> finalSelected;
};

Result De_KS(KS P, int cnt, vector<int> selected) {

	Result totalResult;
	if (cnt == -1 || P.remainCap <= 0) {
		double totalVal = 0;
		for (int i = 0; i < selected.size(); i++) {
			totalVal = totalVal + P.value[selected[i]];
		}
		totalResult.totalVal = totalVal;
		totalResult.finalSelected = selected;

		/*for (int i = 0; i < selected.size(); i++) {
			std::cout << selected[i] << " ";
		}
		std::cout << ": " << totalVal << endl;*/

		return totalResult;
	}

	if (P.realWeight[cnt] > P.remainCap) {
		return De_KS(P, cnt - 1, selected);
	}

	else {
		vector<int> selected1 = selected;
		selected1.push_back(cnt);
		KS P1 = P;
		P1.remainCap = P1.remainCap - P1.realWeight[cnt];

		vector<int> selected2 = selected;
		KS P2 = P;

		Result result1 = De_KS(P1, cnt - 1, selected1);
		Result result2 = De_KS(P2, cnt - 1, selected2);
		Result result;

		if (result1.totalVal < result2.totalVal) {
			result = result2;
		}
		else {
			result = result1;
		}

		return result;
	}
}

PiecewiseLinearFunction WeightTrain_KS(KS P, int cnt, vector<int> selected) {

	PiecewiseLinearFunction totalResult;
	//std::cout << "cnt: " << cnt << endl;

	if (cnt == -1) {
		totalResult.assign(P.preRemainCap.start_x[0], P.preRemainCap.end_x[0], 0, 0, 0);

		double totalVal = 0;
		for (int i = 0; i < selected.size(); i++) {
			totalVal = totalVal + P.value[selected[i]];
		}
		totalResult.id[0] = totalVal;
		totalResult.selected.push_back(selected);

		//totalResult.output();
		//return totalResult;
	}

	else {
		PiecewiseLinearFunction canSelectOrNot = InfCombination(P.weight[cnt], P.preRemainCap);
		/*std::cout << "canSelectOrNot: " << endl;
		canSelectOrNot.output();
		std::cout << endl;*/

		int range = canSelectOrNot.a.size();
		for (int cntRange = 0; cntRange < range; cntRange++) {
			//std::cout << "cntRange: " << cntRange << endl;
			int cntTemp = cnt;

			if (canSelectOrNot.a[cntRange] == P.weight[cnt].a[0] && canSelectOrNot.b[cntRange] == P.weight[cnt].b[0]) {
				vector<int> selected1 = selected;
				selected1.push_back(cnt);
				KS P1 = P;
				/*update coefficient range*/
				for (int i = 0; i < P1.itemNum; i++) {
					P1.weight[i].start_x[0] = canSelectOrNot.start_x[cntRange];
					P1.weight[i].end_x[0] = canSelectOrNot.end_x[cntRange];
				}
				P1.preRemainCap = MinusSpecial(P1.preRemainCap, P1.weight[cnt], 0);
				//P1.remainCap = P1.remainCap - P1.realWeight[cnt];

				vector<int> selected2 = selected;
				KS P2 = P;
				/*update coefficient range*/
				for (int i = 0; i < P2.itemNum; i++) {
					P2.weight[i].start_x[0] = canSelectOrNot.start_x[cntRange];
					P2.weight[i].end_x[0] = canSelectOrNot.end_x[cntRange];
				}
				P2.preRemainCap.start_x[0] = canSelectOrNot.start_x[cntRange];
				P2.preRemainCap.end_x[0] = canSelectOrNot.end_x[cntRange];

				cnt = cnt - 1;
				PiecewiseLinearFunction result1 = WeightTrain_KS(P1, cnt, selected1);
				PiecewiseLinearFunction result2 = WeightTrain_KS(P2, cnt, selected2);
				result1.mergeWithSelect();
				result2.mergeWithSelect();

				/*std::cout << "result 1" << endl;
				result1.output();
				std::cout << "result 2" << endl;
				result2.output();
				std::cout << endl;*/

				PiecewiseLinearFunction result;
				result = SupValCombination(result1, result2);
				result.mergeWithSelect();

				/*std::cout << "result" << endl;
				result.output();
				std::cout << endl;*/

				for (int i = 0; i < result.start_x.size(); i++) {
					totalResult.assign(result.start_x[i], result.end_x[i], 0, 0, result.id[i]);
					totalResult.selected.push_back(result.selected[i]);
				}

				if (cntRange != range - 1) {
					cnt = cntTemp;
				}
				//return result;
			}

			else {
				vector<int> selected2 = selected;
				KS P2 = P;
				/*update coefficient range*/
				for (int i = 0; i < P2.itemNum; i++) {
					P2.weight[i].start_x[0] = canSelectOrNot.start_x[cntRange];
					P2.weight[i].end_x[0] = canSelectOrNot.end_x[cntRange];
				}
				P2.preRemainCap.start_x[0] = canSelectOrNot.start_x[cntRange];
				P2.preRemainCap.end_x[0] = canSelectOrNot.end_x[cntRange];

				cnt = cnt - 1;
				PiecewiseLinearFunction result2 = WeightTrain_KS(P2, cnt, selected2);
				result2.mergeWithSelect();

				/*std::cout << "result 2:" << endl;
				result2.output();
				std::cout << endl;*/

				for (int i = 0; i < result2.start_x.size(); i++) {
					totalResult.assign(result2.start_x[i], result2.end_x[i], 0, 0, result2.id[i]);
					totalResult.selected.push_back(result2.selected[i]);
				}

				if (cntRange != range - 1) {
					cnt = cntTemp;
				}

				//return result2;
			}
		}
	}

	/*std::cout << "totalResult:" << endl;
	totalResult.output();
	std::cout << endl;*/
	return totalResult;
}

double WeightTest_KS(KS P, vector<int> selected) {
	double totalVal = 0;
	for (int i = 0; i < selected.size(); i++) {
		totalVal = totalVal + P.value[selected[i]];
	}
	return totalVal;
}


int main(int argc, char* argv[])
{
	if (argc != 9) {
		std::cerr << "usage: ./train_BAL [RidgePara file] [iteration] [train weight] [train value] [runtime file] [test weight] [predict weight] [cap]" << std::endl;
		exit(0);
	}

	int itemNum = 10;
	int featureNum = 8;
	vector<double> alpha(featureNum, 1);

	vector<double> regretRec(featureNum, 0);
	double regret = INF;
	vector<PiecewiseLinearFunction> ans;
	PiecewiseLinearFunction finalRes;
	int trainNum = 210;
	int testNum = 90;
    remainCap = atoi(argv[8]);

	ifstream trainWeightFile;
	ifstream trainValFile;
	ifstream RidgeParaFile;
	//ifstream flowfile;
	double benchmarkId;

	RidgeParaFile.open(argv[1]);
	for (int i = 0; i < featureNum; i++) {
		RidgeParaFile >> alpha[i];
	}

	std::cout << "==================================================== training ====================================================" << endl;
	time_t Tstart, Tend;
	Tstart = time(NULL);

	int iteration = atoi(argv[2]);
	for (int T = 0; T < iteration; T++) {
		for (int k = 0; k < featureNum; k++) {
			trainWeightFile.open(argv[3]);
			trainValFile.open(argv[4]);

			for (int benN = 0; benN < trainNum; benN++) {
				vector<PiecewiseLinearFunction> weight;
				weight.resize(itemNum);
				vector<double> trainValue;
				trainValue.resize(itemNum);
				vector<vector<double>> trainFeature(itemNum, vector<double>(featureNum, 0));
				vector<double> trainWeight(itemNum, 0);

				for (int i = 0; i < itemNum; i++) {
					trainWeightFile >> benchmarkId;
					for (int j = 0; j < featureNum; j++) {
						trainWeightFile >> trainFeature[i][j];
					}
					trainWeightFile >> trainWeight[i];
					trainValFile >> trainValue[i];
				}

				vector<double> A(itemNum, 0);
				vector<double> B(itemNum, 0);
				for (int i = 0; i < itemNum; i++) {
					A[i] = trainFeature[i][k];
					B[i] = 0;
					for (int d = 0; d < featureNum; d++) {
						B[i] = B[i] + trainFeature[i][d] * alpha[d];
					}
					B[i] = B[i] - trainFeature[i][k] * alpha[k];
				}

				for (int i = 0; i < itemNum; i++) {
					weight[i].assign(-INF, INF, A[i], B[i], trainWeight[i]);
					//cost[i].output();
				}

				KS P(itemNum, remainCap, trainWeight, weight, trainValue);

				vector<int> S;
				Result realValue = De_KS(P, itemNum - 1, S);

				/*std::cout << "totalVal: " << realValue.totalVal << endl;
				std::cout << "finalSelected: ";
				for (int i = 0; i < realValue.finalSelected.size(); i++) {
					std::cout << realValue.finalSelected[i] << " ";
				}
				std::cout << endl;*/

				S.clear();
				PiecewiseLinearFunction predictValue = WeightTrain_KS(P, itemNum - 1, S);
				//predictValue.outputSelect();

				//cout << "=============== after correction ===============" << endl;
				//PiecewiseLinearFunction corrValue = Correction1ForTrain(predictValue, P, removalFactor, PenaltyVORW);
        //PiecewiseLinearFunction corrValue = Correction4ForTrain(predictValue, P, removalFactor);
				//corrValue.outputSelect();

				predictValue.computeRegret(realValue.totalVal);
				//predictValue.merge();
				ans.push_back(predictValue);
			}

			PiecewiseLinearFunction res;
			res.assign(-INF, INF, 0, 0, 0);
			int g = 0;
			while (g < ans.size()) {
				res = Plus(res, ans[g]);
				g++;
			}
			cout << endl;
			cout << "=============================== the " << k << "th alpha done================================" << endl;
			//res.output();
			res.merge();

			vector<double>::iterator smallest = min_element(begin(res.id), end(res.id));
			int index = distance(begin(res.id), smallest);
			double Rtemp = *smallest;
			cout << "minimum regret = " << Rtemp / trainNum << " at position " << index << endl;

			if (Rtemp / trainNum < regret) {
				regret = Rtemp / trainNum;
				if (index == 0) {
					alpha[k] = res.end_x[index] - 1;
				}
				else if (index == res.start_x.size() - 1) {
					alpha[k] = res.start_x[index] + 1;
				}
				else {
					alpha[k] = (res.start_x[index] + res.end_x[index]) / 2;
				}
				cout << "alpha[" << k << "]: " << alpha[k] << endl;
			}
			else {
				cout << "alpha[" << k << "]: " << alpha[k] << endl;
			}
			cout << "minimum regret = " << regret << endl;
			regretRec[k] = regret;

			if (k == (featureNum - 1)) {
				finalRes = res;
			}

			ans.clear();
			trainWeightFile.close();
			trainValFile.close();
		}
		if (regretRec[featureNum - 1] == regretRec[0])
			break;
	}
	cout << "================================================= training done ==================================================" << endl;
	vector<double>::iterator finalSmallest = min_element(begin(finalRes.id), end(finalRes.id));
	double totalRegret = *finalSmallest / trainNum;
	//cout << "trainRegret = " << totalRegret << endl;
	Tend = time(NULL);
	double diff = difftime(Tend, Tstart);
	cout << "Time = " << diff << endl;
	ofstream runtimefile;
	runtimefile.open(argv[5]);
	runtimefile << diff << endl;
	runtimefile.close();


	// read test data 
	// int testbenchmarkNum, testfeatureNum; 
	ifstream testWeightFile;
	testWeightFile.open(argv[6]);
	ofstream preWeightFile;
	preWeightFile.open(argv[7]);
	vector<vector<vector<double>>> testFeatures(testNum, vector<vector<double>>(itemNum, vector<double>(featureNum)));
	vector<vector<double>> testRealWeight(testNum, vector<double>(itemNum));

	for (int i = 0; i < testNum; i++) {
		for (int j = 0; j < itemNum; j++) {
			double benchmarkId;
			testWeightFile >> benchmarkId;
			double predictedWeight = 0.0;
			for (int k = 0; k < featureNum; k++) {
				testWeightFile >> testFeatures[i][j][k];
				predictedWeight += testFeatures[i][j][k] * alpha[k];
			}
			// std::cout << std::endl;
			testWeightFile >> testRealWeight[i][j];
			preWeightFile << i << " " << testRealWeight[i][j] << " " << predictedWeight << std::endl;
		}
	}
	//ofile << "Time = " << diff << endl;
	testWeightFile.close();
	preWeightFile.close();
}
