# Python program for implementation
# of Ford Fulkerson algorithm
import sys
from collections import defaultdict
import numpy as np
import gurobipy as gp
from gurobipy import GRB
from sklearn.metrics import mean_squared_error
from numpy import inf

# This class represents a directed graph
# using adjacency matrix representation
class Graph:

    def __init__(self, arc, realCap, preCap):
        self.arc = arc
        self.realCap = realCap
        self.preCap = preCap
        self.edgeID = np.zeros(edgeNum)
        
        for i in range(edgeNum):
            self.edgeID[i] = 10 * arc[i][0] + arc[i][1]

        self.realGraph = np.zeros((nodeNum, nodeNum))
        for i in range(edgeNum):
            self.realGraph[arc[i][0], arc[i][1]] = realCap[i]

        self.preGraph = np.zeros((nodeNum, nodeNum))
        for i in range(edgeNum):
            self.preGraph[arc[i][0], arc[i][1]] = preCap[i]

        self.ROW = len(self.realGraph)

    # self.COL = len(gr[0])

    '''Returns true if there is a path from source 's' to sink 't' in
    residual graph. Also fills parent[] to store the path '''

    def BFS(self, s, t, parent, graph):

        # Mark all the vertices as not visited
        visited = [False] * (self.ROW)
        # Create a queue for BFS
        queue = []
        # Mark the source node as visited and enqueue it
        queue.append(s)
        visited[s] = True
        # Standard BFS Loop
        while queue:
            # Dequeue a vertex from queue and print it
            u = queue.pop(0)
            # Get all adjacent vertices of the dequeued vertex u
            # If a adjacent has not been visited, then mark it
            # visited and enqueue it
            for ind, val in enumerate(graph[u]):
                if visited[ind] == False and val > 0:
                    # If we find a connection to the sink node,
                    # then there is no point in BFS anymore
                    # We just have to set its parent and can return true
                    queue.append(ind)
                    visited[ind] = True
                    parent[ind] = u
                    if ind == t:
                        return True
        # We didn't reach sink in BFS starting
        # from source, so return false
        return False

    # Returns tne maximum flow from s to t in the given graph
    def maxFlow(self, source, sink):
        x_sol = np.zeros(allPathNum)
        # graphTemp = self.realGraph.copy()
        # print(self.realGraph)
        # This array is filled by BFS and to store path
        parent = [-1] * (self.ROW)
        max_flow = 0  # There is no flow initially
        # Augment the flow while there is path from source to sink
        while self.BFS(source, sink, parent, self.realGraph):
            # Find minimum residual capacity of the edges along the
            # path filled by BFS. Or we can say find the maximum flow
            # through the path found.
            path_flow = float("Inf")
            s = sink
            while (s != source):
                path_flow = min(path_flow, self.realGraph[parent[s]][s])
                s = parent[s]
            # Add path flow to overall flow
            max_flow += path_flow
            # update residual capacities of the edges and reverse edges
            # along the path
            v = sink
            currPath = []
#            print(v, end=" ")
            while (v != source):
                u = parent[v]
                currEdgeID = 10*u+v
                currPath.append(currEdgeID)
#                print(u, end=" ")
                self.realGraph[u][v] -= path_flow
                self.realGraph[v][u] += path_flow
                v = parent[v]
#            print("")
            currPath.reverse()
            currPath_str = [str(e) for e in currPath]
            currPathVal = " ".join(currPath_str)
#            print(currPathVal)
            currPathID = pathDic_new[currPathVal]
            x_sol[currPathID] = path_flow
#            print(currPathID, end=" ")
#            print(currPath, path_flow)

        # reset realGraph
        # for i in range(nodeNum):
        #     for j in range(nodeNum):
        #         self.realGraph[i][j] = graphTemp[i][j]
        # print(self.realGraph)
#        max_flow = sum(x_sol)
#        print(x_sol)
        return max_flow
        
    def corr_maxFlow(self, source, sink):
        graphTemp = self.preGraph.copy()
        realGraphTemp = self.realGraph.copy()
        preFlow = np.zeros([nodeNum, nodeNum])
        # print(self.realGraph)
        # This array is filled by BFS and to store path
        parent = [-1] * (self.ROW)
        max_flow = 0  # There is no flow initially
        # Augment the flow while there is path from source to sink
        while self.BFS(source, sink, parent, self.preGraph):
            # Find minimum residual capacity of the edges along the
            # path filled by BFS. Or we can say find the maximum flow
            # through the path found.
            path_flow = float("Inf")
            path_flow_index = sink
            s = sink
            while (s != source):
                # path_flow = min(path_flow, self.preGraph[parent[s]][s])
                # print(parent[s], "->", end="")
                if path_flow > self.preGraph[parent[s]][s]:
                    path_flow = self.preGraph[parent[s]][s]
                    path_flow_index = s
                s = parent[s]

            preDf = path_flow # path_flow = self.preGraph[parent[path_flow_index]][path_flow_index]
            df = self.realGraph[parent[path_flow_index]][path_flow_index]
            # Add path flow to overall flow
            max_flow += df
            # print(df)
            # print(preDf, df)
            # update residual capacities of the edges and reverse edges
            # along the path
            v = sink
#            print(v, end=" ")
            while (v != source):
                u = parent[v]
#                print(u, end=" ")
                # self.preGraph[u][v] -= path_flow
                # self.preGraph[v][u] += path_flow
                self.preGraph[u][v] -= preDf
                self.realGraph[u][v] -= df
                preFlow[u][v] = preFlow[u][v] + df
                self.preGraph[v][u] += preDf
                self.realGraph[v][u] += df
                v = parent[v]
#            print("\n")

        # print(realGraphTemp)
        # np.savetxt("realGraphTemp.txt", realGraphTemp, fmt="%.2f")
        # print(self.realGraph)
        # np.savetxt("realGraph.txt", self.realGraph, fmt="%.2f")
        # print(preFlow)
        # np.savetxt("preFlow.txt", preFlow, fmt="%.2f")
        # check whether violate constraint and find tau
        tauTemp = []
        for i in range(nodeNum):
            for j in range(nodeNum):
                if preFlow[i][j] > realGraphTemp[i][j]:
                    newT = realGraphTemp[i][j]/preFlow[i][j]
#                    print(i, j, preFlow[i][j], realGraphTemp[i][j], newT)
                    tauTemp.append(newT)
                # if self.realGraph[i][j] < 0 and realGraphTemp[i][j] > 0:
                #     newT = realGraphTemp[i][j]/(-self.realGraph[i][j]+realGraphTemp[i][j])
                #     tauTemp.append(newT)
        # print(tauTemp)
        if len(tauTemp) > 0:
            tau = min(tauTemp)
            max_flow = tau * max_flow

        # reset realGraph
        for i in range(nodeNum):
            for j in range(nodeNum):
                self.preGraph[i][j] = graphTemp[i][j]
        # print(self.realGraph)
        return max_flow
        
    
    def sol_for_forprop(self, source, sink):
        x_sol = np.zeros(allPathNum)
        parent = [-1] * (self.ROW)
        max_flow = 0  # There is no flow initially
        # Augment the flow while there is path from source to sink
        while self.BFS(source, sink, parent, self.preGraph):
            path_flow = float("Inf")
            path_flow_index = sink
            s = sink
            while (s != source):
                # path_flow = min(path_flow, self.preGraph[parent[s]][s])
                if path_flow > self.preGraph[parent[s]][s]:
                    path_flow = self.preGraph[parent[s]][s]
                    path_flow_index = s
                s = parent[s]

            preDf = path_flow # path_flow = self.preGraph[parent[path_flow_index]][path_flow_index]
            df = self.realGraph[parent[path_flow_index]][path_flow_index]
            # Add path flow to overall flow
            max_flow += df
            # print(preDf, df)
            # update residual capacities of the edges and reverse edges
            # along the path
            v = sink
            currPath = []
            while (v != source):
                u = parent[v]
                currEdgeID = 10*u+v
                currPath.append(currEdgeID)
                # self.preGraph[u][v] -= path_flow
                # self.preGraph[v][u] += path_flow
                self.preGraph[u][v] -= preDf
                self.realGraph[u][v] -= df
                self.preGraph[v][u] += preDf
                self.realGraph[v][u] += df
                v = parent[v]
                
            currPath.reverse()
            currPath_str = [str(e) for e in currPath]
            currPathVal = " ".join(currPath_str)
#            print(currPathVal)
            currPathID = pathDic_new[currPathVal]
            x_sol[currPathID] = df
        # reset realGraph
        # for i in range(nodeNum):
        #     for j in range(nodeNum):
        #         self.preGraph[i][j] = graphTemp[i][j]
        # print(self.realGraph)
        return x_sol
        

    def maxFlow_model(self, source, sink):
        # graphTemp = self.preGraph.copy()
        # print(self.realGraph)
        # This array is filled by BFS and to store path
        parent = [-1] * (self.ROW)
        max_flow = 0  # There is no flow initially
        # Augment the flow while there is path from source to sink
        while self.BFS(source, sink, parent, self.preGraph):
            # Find minimum residual capacity of the edges along the
            # path filled by BFS. Or we can say find the maximum flow
            # through the path found.
            path_flow = float("Inf")
            path_flow_index = sink
            s = sink
            while (s != source):
                # path_flow = min(path_flow, self.preGraph[parent[s]][s])
                if path_flow > self.preGraph[parent[s]][s]:
                    path_flow = self.preGraph[parent[s]][s]
                    path_flow_index = s
                s = parent[s]

            preDf = path_flow # path_flow = self.preGraph[parent[path_flow_index]][path_flow_index]
            df = self.realGraph[parent[path_flow_index]][path_flow_index]
            # Add path flow to overall flow
            max_flow += df
            # print(preDf, df)
            # update residual capacities of the edges and reverse edges
            # along the path
            v = sink
            while (v != source):
                u = parent[v]
                # self.preGraph[u][v] -= path_flow
                # self.preGraph[v][u] += path_flow
                self.preGraph[u][v] -= preDf
                self.realGraph[u][v] -= df
                self.preGraph[v][u] += preDf
                self.realGraph[v][u] += df
                v = parent[v]
        # reset realGraph
        # for i in range(nodeNum):
        #     for j in range(nodeNum):
        #         self.preGraph[i][j] = graphTemp[i][j]
        # print(self.realGraph)
        return max_flow


#def actual_obj(source, sink, arc, capT, n_instance):
#    obj_list = []
#    for num in range(n_instance):
#        cap = np.zeros(edgeNum)
#        cnt = num * edgeNum
#        for i in range(edgeNum):
#            cap[i] = capT[cnt]
#            cnt = cnt + 1
#        g = Graph(arc, cap, cap)
#        objective = g.maxFlow(source, sink)
#        obj_list.append(objective)
#    return np.array(obj_list)

def actual_obj(c, G, hTemp, n_instance):
    obj_list = []
    for num in range(n_instance):
        h = np.zeros(edgeNum)
        cnt = num * edgeNum
        for i in range(edgeNum):
            h[i] = hTemp[cnt]
            cnt = cnt + 1

        m = gp.Model()
        m.setParam('OutputFlag', 0)
        x = m.addVars(allPathNum, vtype=GRB.CONTINUOUS, name='x')
        m.setObjective(x.prod(c), GRB.MAXIMIZE)
#        m.addConstr((x.prod(G)) <= h)
        for i in range(edgeNum):
            m.addConstr((x.prod(G[i])) <= h[i])

        m.optimize()
        sol = []
        for i in range(allPathNum):
            sol.append(x[i].x)
#        for i in range(allPathNum):
#            if sol[i] != 0:
#                print(i, end=" ")
        objective = m.objVal
        obj_list.append(objective)
#        print("")
    return np.array(obj_list)

def correction_obj(c, G, realhTemp, prehTemp, n_instance):
    obj_list = []
    for num in range(n_instance):
        realh = np.zeros(edgeNum)
        preh = np.zeros(edgeNum)
        cnt = num * edgeNum
        for i in range(edgeNum):
            realh[i] = realhTemp[cnt]
            preh[i] = prehTemp[cnt]
            cnt = cnt + 1

        m = gp.Model()
        m.setParam('OutputFlag', 0)
        x = m.addVars(allPathNum, vtype=GRB.CONTINUOUS, name='x')
        m.setObjective(x.prod(c), GRB.MAXIMIZE)
#        m.addConstr((x.prod(G)) <= h)
        for i in range(edgeNum):
            m.addConstr((x.prod(G[i])) <= preh[i])

        m.optimize()
        sol = []
        for i in range(allPathNum):
            sol.append(x[i].x)
#        for i in range(allPathNum):
#            if sol[i] != 0:
#                print(i, end=" ")
        objective = m.objVal
        
        #correction
        Gx = np.dot(G, sol)
#        print(Gx, objective)
        tauTemp = []
        for i in range(edgeNum):
            if Gx[i] > realh[i]:
                newT = realh[i]/Gx[i]
                tauTemp.append(newT)
        if len(tauTemp) > 0:
            tau = min(tauTemp)
            objective = tau * objective
            sol = np.multiply(sol, tau)

        obj_list.append(objective)
#        print(np.dot(G, sol), objective)
#        print("")
    return np.array(obj_list)
    
def correction_single_obj(c, G, realhTemp, prehTemp):
    realh = np.zeros(edgeNum)
    preh = np.zeros(edgeNum)
    for i in range(edgeNum):
        realh[i] = realhTemp[i]
        preh[i] = prehTemp[i]
    if min(preh) >= 0:
        m = gp.Model()
        m.setParam('OutputFlag', 0)
        x = m.addVars(allPathNum, vtype=GRB.CONTINUOUS, name='x')
        m.setObjective(x.prod(c), GRB.MAXIMIZE)
    #        m.addConstr((x.prod(G)) <= h)
        for i in range(edgeNum):
            m.addConstr((x.prod(G[i])) <= preh[i])

        m.optimize()
        sol = []
    #    print(x)
        for i in range(allPathNum):
            sol.append(x[i].x)
    #        for i in range(allPathNum):
    #            if sol[i] != 0:
    #                print(i, end=" ")
        objective = m.objVal
        
        #correction
        Gx = np.dot(G, sol)
    #        print(Gx, objective)
        tauTemp = []
        for i in range(edgeNum):
            if Gx[i] > realh[i]:
                newT = realh[i]/Gx[i]
                tauTemp.append(newT)
        if len(tauTemp) > 0:
            tau = min(tauTemp)
            objective = tau * objective
            sol = np.multiply(sol, tau)
            
    else:
        objective = 0

#        print(np.dot(G, sol), objective)
#        print("")
    return objective
    

def correction2_single_obj(c, G, realhTemp, prehTemp):
    realh = np.zeros(edgeNum)
    preh = np.zeros(edgeNum)
    for i in range(edgeNum):
        realh[i] = realhTemp[i]
        preh[i] = prehTemp[i]
    if min(preh) >= 0:
        m = gp.Model()
        m.setParam('OutputFlag', 0)
        x = m.addVars(allPathNum, vtype=GRB.CONTINUOUS, name='x')
        m.setObjective(x.prod(c), GRB.MAXIMIZE)
    #        m.addConstr((x.prod(G)) <= h)
        for i in range(edgeNum):
            m.addConstr((x.prod(G[i])) <= preh[i])

        m.optimize()
        sol = []
    #    print(x)
        for i in range(allPathNum):
            sol.append(x[i].x)
    #        for i in range(allPathNum):
    #            if sol[i] != 0:
    #                print(i, end=" ")
        objective = m.objVal
        
        #correction
        Gx = np.dot(G, sol)
    #        print(Gx, objective)
        tauTemp = []
        for i in range(edgeNum):
            if realh[i] != 0:
                newT = Gx[i]/realh[i]
            else:
                newT = 0
            tauTemp.append(newT)
        if len(tauTemp) > 0:
            tau = max(tauTemp)
            if tau == 0:
                objective = 0
            else:
                objective = objective / tau
            sol = np.multiply(sol, tau)
            
    else:
        objective = 0

#        print(np.dot(G, sol), objective)
#        print("")
    return objective
    
def correction3_single_obj(c, G, realhTemp, prehTemp):
    realh = np.zeros(edgeNum)
    preh = np.zeros(edgeNum)
    for i in range(edgeNum):
        realh[i] = realhTemp[i]
        preh[i] = prehTemp[i]
    if min(preh) >= 0:
        m = gp.Model()
        m.setParam('OutputFlag', 0)
        x = m.addVars(allPathNum, vtype=GRB.CONTINUOUS, name='x')
        m.setObjective(x.prod(c), GRB.MAXIMIZE)
    #        m.addConstr((x.prod(G)) <= h)
        for i in range(edgeNum):
            m.addConstr((x.prod(G[i])) <= preh[i])

        m.optimize()
        sol = []
    #    print(x)
        for i in range(allPathNum):
            sol.append(x[i].x)
    #        for i in range(allPathNum):
    #            if sol[i] != 0:
    #                print(i, end=" ")
        objective = m.objVal
        
        #correction
        Gx = np.dot(G, sol)
    #        print(Gx, objective)
        tauTemp = []
        tauTemp.append(1)
        for i in range(edgeNum):
            if realh[i] != 0:
                newT = Gx[i]/realh[i]
            else:
                newT = float('inf')
            tauTemp.append(newT)
        if len(tauTemp) > 0:
            tau = max(tauTemp)
            if tau is inf:
                objective = 0
            else:
                objective = objective / tau
#            sol = np.multiply(sol, tau)
            
    else:
        objective = 0

#        print(np.dot(G, sol), objective)
#        print("")
    return objective

    
def pred_single_obj(c, G, prehTemp):
    preh = np.zeros(edgeNum)
    for i in range(edgeNum):
        preh[i] = prehTemp[i]

    m = gp.Model()
    m.setParam('OutputFlag', 0)
    x = m.addVars(allPathNum, vtype=GRB.CONTINUOUS, name='x')
    m.setObjective(x.prod(c), GRB.MAXIMIZE)
#        m.addConstr((x.prod(G)) <= h)
    for i in range(edgeNum):
        m.addConstr((x.prod(G[i])) <= preh[i])

    m.optimize()
    sol = []
    for i in range(allPathNum):
        sol.append(x[i].x)
    objective = m.objVal
    
    return objective



# Create a graph given in the above diagram

nodeNum = 12
edgeNum = 18
testmarkNum = 179
allPathNum = 13

arc_temp = np.loadtxt('./data/graph_POLSKA.txt', dtype=int)
arc = arc_temp[1:,:]
src = 0
dist = nodeNum - 1
c_data = np.loadtxt('./data/POLSKA_C011.txt')
G_data = np.loadtxt('./data/POLSKA_G011.txt')
c_data = c_data.tolist()
G_data = G_data.tolist()

startmark = int(sys.argv[1])
endmark = int(sys.argv[2])

for testmark in range(startmark, endmark):
    capT = np.loadtxt('./data/BALC_POLSKA_corrA/BALC_POLSKA_corrA(' + str(testmark) + ').txt')
#        capT = np.loadtxt('BALC_POLSKA.txt')
    realCapT = capT[:, 1]
    preCapT = capT[:, 2]


    real_obj = actual_obj(c_data, G_data, realCapT, n_instance=testmarkNum)


    corr_obj_list = []
    pre_obj_list = []
    relativeErr_list = []
#        print(real_obj)
    for testNum in range(testmarkNum):
#        print(testNum)
        realCap = {}
        predCap = {}
        for i in range(edgeNum):
            realCap[i] = realCapT[i+testNum*edgeNum]
            predCap[i] = preCapT[i+testNum*edgeNum]
        g = Graph(arc, realCap, predCap)
        corrrlst = correction_single_obj(c_data, G_data, realCap, predCap)

        predrlst = g.maxFlow_model(src, dist)
        corr_obj_list.append(corrrlst)
        pre_obj_list.append(predrlst)

        if real_obj[testNum] != 0:
            relativeErr = abs(real_obj[testNum] - corrrlst)/real_obj[testNum]
            relativeErr_list.append(relativeErr)

    print("corrReg: ", sum(abs(real_obj - np.array(corr_obj_list)))/testmarkNum, "TOV: ", sum(real_obj)/testmarkNum)
print("\n")
#    print("TOV: ", sum(real_obj))
        
