import sys
import time
import numpy as np
from itertools import combinations, permutations
import gurobipy as gp
from gurobipy import GRB, tupledict
import igraph as ig
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

version = "-TS"
outputFlag =1
solutionCheckingFlag = 0
f = 2.8
instanceList = []
numFiles = 10

# for name in ['Champaign']:

for name in ['Cook', 'Winnebago', 'Champaign', 'LaSalle', 'Adams', 'Fulton', 'Jefferson', 'Johnson', 'Cumberland']:

    for fileID in range(0, numFiles):
        instanceList.append((name, fileID + 1))

# for B in [3]:
for B in [1, 2, 3, 4, 5]:
    for instance in instanceList:

        start_time = time.time()
        name = instance[0]
        fileID = instance[1]
        print("TS-RCT2-CAVAD"+ "_B_" + str(B)+"_County_" + str(name)+"_ID_" + str(fileID) )

        filename = "./Test_Instances/Urban_Rural_Instances/" + str(name) + "_" + str(fileID) + ".csv"
        filename_Driving = "./Test_Instances/Urban_Rural_Instances/" + str(name) + "_" + str(fileID) + "_Driving.csv"
        filename_Walking = "./Test_Instances/Urban_Rural_Instances/" + str(name) + "_" + str(fileID) + "_Walking.csv"
        loc_pd = pd.read_csv(filename)
        n = loc_pd.shape[0]
        cn = n - 1
        D_pd = pd.read_csv(filename_Driving)
        W_pd = pd.read_csv(filename_Walking)
        Driving_list = D_pd[['key', 'Time']].to_numpy()
        drivingTime = {}
        for d in Driving_list:
            drivingTime[tuple(map(int, d[0][1:-1].split(', ')))] = d[1]


        Walking_list = W_pd[['key', 'Time']].to_numpy()
        walkingTime = {}
        for w in Walking_list:
            walkingTime[tuple(map(int, w[0][1:-1].split(', ')))] = w[1]

        drivingAlwaysFasterFlag = 1

        "Build resitricted walking graph (rwG)"
        rwG = ig.Graph(directed=True)
        rwG.add_vertices(n)
        rwG.vs["name_id"]=list(range(n))
        eWweights = []
        eDweights = []
        for i, j in walkingTime.keys():
            if (i != 0) and (j != 0) and (i != j) and ((f + drivingTime[i, j]) > walkingTime[i, j]):
                rwG.add_edge(i, j)
                eWweights.append(walkingTime[i, j])
                eDweights.append(drivingTime[i, j])
            if drivingAlwaysFasterFlag==1 and (drivingTime[i, j] > walkingTime[i, j]):
                drivingAlwaysFasterFlag = 0
                print("Driving is not always faster")
        rwG.es["wTime"] = eWweights
        rwG.es["dTime"] = eDweights

        # if drivingAlwaysFasterFlag == 0:
        #     continue

        "TS Formulation"
        TS = gp.Model("CAVADP_TS")

        x = gp.tupledict()
        for i in range(n):
            for j in range(n):
                x[i, j] = TS.addVar(obj=drivingTime[i, j], vtype=GRB.BINARY, name='x[%d,%d]' % (i, j))
        xInd = dict((i, x[i].index) for i in x.keys())

        v = gp.tupledict()
        for i in range(n):
            for j in range(1,n):
                if i != j:
                    v[i, j] = TS.addVar(obj=0, vtype=GRB.INTEGER, name='v[%d,%d]' % (i, j))
        vInd = dict((i, v[i].index) for i in v.keys())

        "Generating the set of y variables"
        y = gp.tupledict()  # initialize a set of y variables
        P = dict()  # paths
        Pdetail = dict()
        reverseP = dict()  # contains the set of loops that contains node i in Q
        for i in range(n):
            reverseP[i] = set()
        Psize = np.zeros((n, n), dtype=int)
        QpSize = dict()
        numQ =0
        if drivingAlwaysFasterFlag == 1:
            "claim 2"
            for i in range(1, n):
                "Lemma 2"
                Q = [i]
                P[i, i, Psize[i, i]] = Q
                Pdetail[i, i, Psize[i, i]] = []
                cpValue = f
                y[i, i, Psize[i, i]] = TS.addVar(obj=cpValue, vtype=GRB.BINARY,
                                                     name='y[%d,%d,%d]' % (i, i, Psize[i, i]))
                QpSize[i, i, Psize[i, i]] = len(Q)
                for j in Q:
                    reverseP[j].add((i, i, Psize[i, i]))
                Psize[i, i] += 1

            for i in range(1, n):
                for j in range(1, n):
                    simple_path = rwG.get_all_simple_paths(i, to=j, cutoff=int(B - 1), mode='out')
                    for path in simple_path:
                        Q = path
                        w_hatTime = 0
                        for ii in range(len(path) - 1):
                            w_hatTime += walkingTime[path[ii], path[ii + 1]]
                        P[i, j, Psize[i, j]] = Q
                        Pdetail[i, j, Psize[i, j]] = [(path[node_id], path[node_id + 1]) for node_id in
                                                      range(len(path) - 1)]
                        cpValue = f + w_hatTime
                        y[i, j, Psize[i, j]] = TS.addVar(obj=cpValue, vtype=GRB.BINARY,
                                                             name='y[%d,%d,%d]' % (i, j, Psize[i, j]))
                        QpSize[i, j, Psize[i, j]] = len(Q)
                        for k in Q:
                            reverseP[k].add((i, j, Psize[i, j]))
                        Psize[i, j] += 1
        else:
            for b in range(1, B + 1):
                for Q in combinations(range(1, n), b):
                    c1c2Flag = 1
                    "c1c2 assumption checking"
                    for i in Q:
                        for k in range(1, n):
                            if drivingTime[i, k] <= walkingTime[i, k] and drivingTime[k, i] <= walkingTime[k, i]:
                                c1c2Flag = 1
                            else:
                                c1c2Flag = 0
                                break
                        if c1c2Flag == 0:
                            break
                    if c1c2Flag == 1:
                        "c1c2 assumption are satisfied"
                        for path in permutations(Q):
                            pathFlag = 1
                            w_hatTime = 0
                            for ii in range(len(path) - 1):
                                if (f + drivingTime[path[ii], path[ii + 1]]) > walkingTime[path[ii], path[ii + 1]]:
                                    w_hatTime += walkingTime[path[ii], path[ii + 1]]
                                else:
                                    pathFlag = 0
                                    break
                            if pathFlag == 1:
                                i = path[0]
                                j = path[-1]
                                P[i, j, Psize[i, j]] = path
                                cpValue = f + max(w_hatTime, drivingTime[i, j])
                                y[i, j, Psize[i, j]] = TS.addVar(obj=cpValue, vtype=GRB.BINARY,
                                                                     name='y[%d,%d,%d]' % (i, j, Psize[i, j]))
                                Pdetail[i, j, Psize[i, j]] = [(path[node_id], path[node_id + 1]) for node_id in
                                                              range(len(path) - 1)]
                                QpSize[i, j, Psize[i, j]] = len(Q)
                                for k in Q:
                                    reverseP[k].add((i, j, Psize[i, j]))
                                Psize[i, j] += 1
                                if (f + drivingTime[j, i]) > walkingTime[j, i]:
                                    "path that starts and ends at the same node i"
                                    cpValue = f + w_hatTime + walkingTime[j, i]
                                    P[i, i, Psize[i, i]] = path
                                    y[i, i, Psize[i, i]] = TS.addVar(obj=cpValue, vtype=GRB.BINARY,
                                                                         name='y[%d,%d,%d]' % (i, i, Psize[i, i]))
                                    Pdetail[i, i, Psize[i, i]] = [(path[node_id], path[node_id + 1]) for node_id in
                                                                  range(len(path) - 1)]
                                    Pdetail[i, i, Psize[i, i]].append((j, i))
                                    QpSize[i, i, Psize[i, i]] = len(Q)
                                    for k in Q:
                                        reverseP[k].add((i, i, Psize[i, i]))
                                    Psize[i, i] += 1

                    else:
                        "c1 and c2 are not satisfied"

                        "both i and j are in Q. Thus, C3 and C4 cannot be satisfied"
                        for path in permutations(Q):
                            i = path[0]
                            j = path[-1]
                            w_hatTime = 0
                            for ii in range(len(path) - 1):
                                w_hatTime += walkingTime[path[ii], path[ii + 1]]
                            P[i, j, Psize[i, j]] = path
                            cpValue = f + max(w_hatTime, drivingTime[i, j])
                            y[i, j, Psize[i, j]] = TS.addVar(obj=cpValue, vtype=GRB.BINARY,
                                                                 name='y[%d,%d,%d]' % (i, j, Psize[i, j]))
                            Pdetail[i, j, Psize[i, j]] = [(path[node_id], path[node_id + 1]) for node_id in
                                                          range(len(path) - 1)]
                            QpSize[i, j, Psize[i, j]] = len(Q)
                            for k in Q:
                                reverseP[k].add((i, j, Psize[i, j]))
                            Psize[i, j] += 1

                            "path that starts and ends at the same node i"
                            cpValue = f + w_hatTime + walkingTime[j, i]
                            P[i, i, Psize[i, i]] = path
                            y[i, i, Psize[i, i]] = TS.addVar(obj=cpValue, vtype=GRB.BINARY,
                                                                 name='y[%d,%d,%d]' % (i, i, Psize[i, i]))
                            Pdetail[i, i, Psize[i, i]] = [(path[node_id], path[node_id + 1]) for node_id in
                                                          range(len(path) - 1)]
                            Pdetail[i, i, Psize[i, i]].append((j, i))
                            QpSize[i, i, Psize[i, i]] = len(Q)
                            for k in Q:
                                reverseP[k].add((i, i, Psize[i, i]))
                            Psize[i, i] += 1

                        "At least one of i and j is not in Q. Thus, C3 and C4 might be satisfied"
                        for i in range(1, n):
                            for j in range(1, n):
                                if i in Q and j in Q:
                                    continue
                                else:
                                    Qlist = list(Q)
                                    c3Flag, c4Flag, QFlag, iInQFlag, jInQFlag = 0, 0, 1, 0, 0
                                    if i not in Qlist:
                                        "need to check claim 3"
                                        c3Flag = 1
                                    else:
                                        iInQFlag = 1

                                    if j not in Qlist:
                                        "need to check claim 4"
                                        c4Flag = 1
                                    else:
                                        jInQFlag = 1

                                    if c3Flag == 1 and QFlag == 1:
                                        min_what = np.inf
                                        min_path = []
                                        if iInQFlag == 0 and jInQFlag == 0:
                                            for path in permutations(Qlist):
                                                w_hatTime = 0
                                                for ii in range(len(path) - 1):
                                                    w_hatTime += walkingTime[path[ii], path[ii + 1]]
                                                if w_hatTime < min_what:
                                                    min_what = w_hatTime
                                                    min_path = list(path)
                                        else:
                                            if jInQFlag == 1:
                                                Qlist.remove(j)
                                            if len(Qlist) >= 1:
                                                for path in permutations(Qlist):
                                                    w_hatTime = walkingTime[path[-1], j]
                                                    for ii in range(len(path) - 1):
                                                        w_hatTime += walkingTime[path[ii], path[ii + 1]]
                                                    if w_hatTime < min_what:
                                                        min_what = w_hatTime
                                                        min_path = list(path)
                                            else:
                                                min_what = 0
                                            min_path.append(j)

                                        if walkingTime[i, min_path[0]] - drivingTime[i, min_path[0]] >= max(0,
                                                                                                            drivingTime[
                                                                                                                min_path[
                                                                                                                    0], j] - (
                                                                                                                    min_what +
                                                                                                                    walkingTime[
                                                                                                                        min_path[
                                                                                                                            -1], j])):
                                            QFlag = 0

                                    if c4Flag == 1 and QFlag == 1:
                                        min_what = np.inf
                                        min_path = []
                                        if iInQFlag == 0 and jInQFlag == 0:
                                            for path in permutations(Qlist):
                                                w_hatTime = 0
                                                for ii in range(len(path) - 1):
                                                    w_hatTime += walkingTime[path[ii], path[ii + 1]]
                                                if w_hatTime < min_what:
                                                    min_what = w_hatTime
                                                    min_path = list(path)
                                        else:
                                            if iInQFlag == 1:
                                                Qlist.remove(i)
                                            if (len(Qlist)) >= 1:
                                                for path in permutations(Qlist):
                                                    w_hatTime = walkingTime[i, path[0]]
                                                    for ii in range(len(path) - 1):
                                                        w_hatTime += walkingTime[path[ii], path[ii + 1]]
                                                    if w_hatTime < min_what:
                                                        min_what = w_hatTime
                                                        min_path = list(path)
                                            else:
                                                min_what = 0
                                            min_path.insert(0, i)
                                        if walkingTime[min_path[-1], j] - drivingTime[min_path[-1], j] >= max(0,
                                                                                                              drivingTime[
                                                                                                                  i,
                                                                                                                  min_path[
                                                                                                                      -1]] - (
                                                                                                                      min_what +
                                                                                                                      walkingTime[
                                                                                                                          i,
                                                                                                                          min_path[
                                                                                                                              0]])):
                                            QFlag = 0

                                    if QFlag == 1:
                                        if (len(Qlist)) >= 1:
                                            for path in permutations(Qlist):
                                                w_hatTime = walkingTime[i, path[0]] + walkingTime[path[-1], j]
                                                for ii in range(len(path) - 1):
                                                    w_hatTime += walkingTime[path[ii], path[ii + 1]]

                                                P[i, j, Psize[i, j]] = Q
                                                cpValue = f + max(w_hatTime, drivingTime[i, j])
                                                y[i, j, Psize[i, j]] = TS.addVar(obj=cpValue, vtype=GRB.BINARY,
                                                                                     name='y[%d,%d,%d]' % (
                                                                                         i, j, Psize[i, j]))

                                                Pdetail[i, j, Psize[i, j]] = [(path[node_id], path[node_id + 1]) for
                                                                              node_id in
                                                                              range(len(path) - 1)]
                                                if j != path[-1]:
                                                    Pdetail[i, j, Psize[i, j]].append((path[-1], j))
                                                if i != path[0]:
                                                    Pdetail[i, j, Psize[i, j]].insert(0,(i, path[0]))
                                                QpSize[i, j, Psize[i, j]] = len(Q)
                                                for k in Q:
                                                    reverseP[k].add((i, j, Psize[i, j]))
                                                Psize[i, j] += 1
                                        else:
                                            w_hatTime = walkingTime[i, j]
                                            P[i, j, Psize[i, j]] = Q
                                            cpValue = f + max(w_hatTime, drivingTime[i, j])
                                            y[i, j, Psize[i, j]] = TS.addVar(obj=cpValue, vtype=GRB.BINARY,
                                                                                 name='y[%d,%d,%d]' % (
                                                                                     i, j, Psize[i, j]))
                                            Pdetail[i, j, Psize[i, j]] = [(path[node_id], path[node_id + 1]) for
                                                                          node_id in
                                                                          range(len(path) - 1)]
                                            QpSize[i, j, Psize[i, j]] = len(Q)
                                            for k in Q:
                                                reverseP[k].add((i, j, Psize[i, j]))
                                            Psize[i, j] += 1


        yInd = dict((i, y[i].index) for i in y.keys())
        modelingYVariablesTime = time.time()
        print("The number of y variables: "+str(len(yInd)))
        print("Model Building Time: " + str(modelingYVariablesTime - start_time))

        TS.update()
        # adding constraints
        "Con1"
        conSet1 = TS.addConstr(gp.quicksum(x[i, 0] for i in range(1,n)) == 1)
        "Con2"
        conSet2 = TS.addConstr(gp.quicksum(x[0, i] for i in range(1,n))== 1)
        "Con3"
        conSet3 = TS.addConstrs(gp.quicksum(y[p] for p in reverseP[i]) == 1 for i in range(1, n))
        "Con4"
        conSet4 = TS.addConstrs(y.sum(i, '*', '*')== x.sum('*', i) for i in range(1, n))
        "Con5"
        conSet5 = TS.addConstrs(y.sum('*', i, '*') == x.sum(i,'*') for i in range(1, n))
        "Con6"
        conSet6 = TS.addConstr(gp.quicksum(v[0,i] for i in range(1,n)) == cn)
        "Con7"
        conSet7 = TS.addConstrs(v[0,i] <= cn*x[0,i] for i in range(1, n))
        "Con8"
        conSet8 = TS.addConstrs(v[i,k] <= cn * (x[i,k]+y.sum(i,k,'*')) for i in range(1, n) for k in range(1, n) if i!=k)
        "Con9"
        conSet9 = TS.addConstrs( v.sum('*',i) - v.sum(i,'*') == x.sum(i,'*') + gp.quicksum(gp.quicksum((QpSize[i, k, pid]-1)*y[i, k, pid] for pid in range(Psize[i,k])) for k in range(1,n)) for i in range(1,n))
        if drivingAlwaysFasterFlag ==1:
            "add claim 5 cuts"
            conSetC5 = TS.addConstrs(gp.quicksum(y[p] for p in reverseP[i] if p in reverseP[k])+ x[i,k] <= 1 for i in range(1, n) for k in range(1, n) if i!=k)
            "add claim 6 cuts"
            conSetC6 = TS.addConstrs(y.sum(k, i, '*') + x[k,i] +y.sum(i,k, '*') + x[i,k]<=1 for i in range(1, n) for k in range(i+1,n))
            "variable fixing"
            for i in range(1, n):
                customerTooFarFlag = 1
                for j in range(1, n):
                    if i != j:
                        if drivingTime[i, j] + f > walkingTime[i, j] or drivingTime[j, i] + f > walkingTime[j, i]:
                            customerTooFarFlag = 0
                            break
                if customerTooFarFlag == 1:
                    TS.addConstr(x.sum(i, '*') == 1)
                    TS.addConstr(x.sum('*', i) == 1)
                    print("fix customer " + str(i))


        TS.update()
        modelingFinishedTime = time.time()
        TS.optimize()
        FinishingTime = time.time()
        print("Model Building Time: " + str(modelingFinishedTime - start_time))
        print("Model Solving Time: " + str(FinishingTime - modelingFinishedTime))
        print("Total Running Time: " + str(FinishingTime - start_time))
        if outputFlag ==1:
            ObjectiveValue = TS.objVal
            file = open("TS-V"+str(version)+"-RCT2-CAVAD_County_" + str(name) + "_ModelBuildingTime.log", 'a')
            file.write("\n" + str(modelingFinishedTime - start_time))
            file.close()
            file = open("TS-V"+str(version)+"-RCT2-CAVAD_County_" + str(name) + "_ModelSolvingTime.log", 'a')
            file.write("\n" + str(FinishingTime - modelingFinishedTime))
            file.close()
            file = open("TS-V"+str(version)+"-RCT2-CAVAD_County_" + str(name)  + "_TotalRunningTime.log", 'a')
            file.write("\n" + str(FinishingTime - start_time))
            file.close()
            file = open("TS-V"+str(version)+"-RCT2-CAVAD_County_" + str(name)  + "_ObjectiveValue.log", 'a')
            file.write("\n" + str(ObjectiveValue))
            file.close()
            file = open("TS-V"+str(version)+"-RCT2-CAVAD_County_" + str(name) + "_NumberOfyVariables.log", 'a')
            file.write("\n" + str(len(yInd)))
            file.close()



