# CPAIOR_BranchAndLearnWithCorrection

This repository is the official implementation of the paper: Branch & Learn with Post-hoc Correction for Predict+Optimize with Unknown Parameters in Constraints.

## “Max Flow” Folder
This is the package for maximum flow with unknown edge capacities. For Correction Function A, the training process and the test process are divided in two files. To compare the solution quality with the state-of-the-art approximation method (IntOpt-C), we use the same python test file as provided. For Correction Function B and Penalty Function I, the training process and the test process are in one file. Please use the following commands:

Conduct training with B&L: 
1.	Compile the program: g++ -std=c++11 BAL.cpp -o BAL
2.	Conduct training: ./BAL [graph file] [train data] [test data] [prediction file] [runtime file] [iteration] [result file] [RidgePara file]  \
E.g.
```
./BAL data/graph_POLSKA.txt data/train_POLSKA/train_POLSKA\(0\).txt data/test_POLSKA/test_POLSKA\(0\).txt data/BAL_POLSKA/BAL_POLSKA\(0\).txt runtime.txt 5 result.txt data/Ridge_para_POLSKA/Ridge_para_POLSKA\(0\).txt
```

Conduct training with B&L-C:
1.	Choose the correction function and penalty function, and compile the corresponding program. 
- 	Train with Correction Function A and no penalty function: g++ -std=c++11 BALC_corrA.cpp -o BALC_corrA
-   Train with Correction Function B and penalty function: g++ -std=c++11 BALC_corrB_pen.cpp -o BALC_corrB_pen
2.	Conduct training:
-	Train with Correction Function A and no penalty function: ./BALC_corrA [graph file] [train data] [test data] [prediction file] [runtime file] [iteration] [RidgePara file]  \
E.g.
```
./BALC_corrA data/graph_POLSKA.txt data/train_POLSKA/train_POLSKA\(0\).txt data/test_POLSKA/test_POLSKA\(0\).txt data/BALC_POLSKA_corrA/BALC_POLSKA_corrA\(0\).txt runtime.txt 5 data/Ridge_para_POLSKA/Ridge_para_POLSKA\(0\).txt
```
-	Train with Correction Function B and penalty function: ./BALC_corrB_pen [graph file] [train data] [test data] [prediction file] [runtime file] [iteration] [result file] [RidgePara file]  \
E.g.
```
./BALC_corrB_pen data/graph_POLSKA.txt data/train_POLSKA/train_POLSKA\(0\).txt data/test_POLSKA/test_POLSKA\(0\).txt data/BALC_POLSKA_corrB/BALC_POLSKA_corrB\(0\).txt runtime.txt 5 result.txt data/Ridge_para_POLSKA/Ridge_para_POLSKA\(0\).txt
```

Conduct testing:
1.	Test with Correction Function A and no penalty function: python3 test_corrA.py [start_simulation_id] [end_simulation_id]  \
E.g.
```
python3 test_corrA.py 0 1
```

## “Knapsack” Folder
This is the package for 0-1 knapsack with unknown weights. The training process and the test process are divided in two files. Please use the following commands:

Conduct training with B&L: 
1.	Compile the program: g++ -std=c++11 train_BAL.cpp -o train_BAL
2.	Conduct training: ./train_BAL [RidgePara file] [iteration] [train weight] [train value] [runtime file] [test weight] [predict weight] [cap]
E.g. ./train_BAL data/ridge_para/ridge_para\(0\).txt 3 data/train_weight/train_weight\(0\).txt data/train_weakly_correlated_value/train_weakly_correlated_value\(0\).txt data/BAL_runtime_100_weakly_correlated/runtime\(0\).txt data/test_weight/test_weight\(0\).txt data/BAL_100_weakly_correlated/BAL_100_weakly_correlated\(0\).txt 100

Conduct training with B&L-C:
1.	Compile the program: g++ -std=c++11 train_BALC.cpp -o train_BALC 
2.	Choose the correction function and penalty function, conduct training: ./train [RidgePara file] [iteration] [train weight] [train value] [runtime file] [test weight] [predict weight] [cap] [correction func] [penalty func]
where 
[correction func] = 1,2,3 corresponding to training with Correction Function A, B, C
[penalty func] = 0, 1 corresponding to training with Penalty Function I, II
E.g. ./train_BALC data/ridge_para/ridge_para\(0\).txt 3 data/train_weight/train_weight\(0\).txt data/train_weakly_correlated_value/train_weakly_correlated_value\(0\).txt data/BALC_runtime_100_weakly_correlated/runtime\(0\).txt data/test_weight/test_weight\(0\).txt data/BALC_100_weakly_correlated/BALC_100_weakly_correlated\(0\).txt 100 1 0

Conduct testing:
1.	Compile the program: g++ -std=c++11 test.cpp -o test
2.	Choose the correction function and penalty function, conduct testing: ./test [test value] [predict weight] [result file] [cap] [correction func] [penalty func]
where 
[correction func] = 1,2,3 corresponding to training with Correction Function A, B, C
[penalty func] = 0, 1 corresponding to training with Penalty Function I, II
E.g. ./test data/test_weakly_correlated_value/test_weakly_correlated_value\(0\).txt data/BALC_100_weakly_correlated/BALC_100_weakly_correlated\(0\).txt data/result/BALCResult\(0\).txt 100 1 0

 
## “MCVC” Folder
This is the package for minimum cost vertex cover with unknown costs and edge values. The training process and the test process are divided in two files. Please use the following commands:

Conduct training with B&L: 
1.	Compile the program: g++ -std=c++11 train_BAL.cpp -o train_BAL
2.	Conduct training: ./train_BAL [graph] [train edge data] [test edge data] [train cost data] [test cost data] [pre edge file] [pre cost file] [iteration] [runtime file]
E.g. ./train_BAL ABILENE data/ABILENE/100/edge/train_ABILENE_100\(0\).txt data/ABILENE/100/edge/test_ABILENE_100\(0\).txt data/ABILENE/100/cost/train_ABILENE_100\(0\).txt data/ABILENE/100/cost/test_ABILENE_100\(0\).txt data/ABILENE/100/BAL_edge/BAL_edge\(0\).txt data/ABILENE/100/BAL_cost/BAL_cost\(0\).txt 1 runtime.txt

Conduct training with B&L-C: 
1.	Compile the program: g++ -std=c++11 train_BALC.cpp -o train_BALC
2.	Conduct training: ./train_BAL [graph] [train edge data] [test edge data] [train cost data] [test cost data] [pre edge file] [pre cost file] [iteration] [runtime file]
E.g. ./train_BALC ABILENE data/ABILENE/100/edge/train_ABILENE_100\(0\).txt data/ABILENE/100/edge/test_ABILENE_100\(0\).txt data/ABILENE/100/cost/train_ABILENE_100\(0\).txt data/ABILENE/100/cost/test_ABILENE_100\(0\).txt data/ABILENE/100/BALC_edge/BALC_edge\(0\).txt data/ABILENE/100/BALC_cost/BALC_cost\(0\).txt 1 runtime.txt
Conduct testing:
1.	Compile the program: g++ -std=c++11 test.cpp -o test
2.	Conduct testing: ./test [graph] [pre edge data] [pre cost data] [result file]
E.g. ./test ABILENE data/ABILENE/100/BAL_edge/BAL_edge\(0\).txt data/ABILENE/100/BAL_cost/BAL_cost\(0\).txt result.txt
