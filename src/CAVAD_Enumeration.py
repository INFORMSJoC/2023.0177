import sys
import time
import math
import random
import numpy as np
from itertools import combinations, permutations
import gurobipy as gp
from gurobipy import GRB, tupledict
import networkx as nx
import igraph as ig
import graph_tool.all as gt
from copy import deepcopy, copy
import pandas as pd
import matplotlib.pyplot as plt
import multiprocessing as mp
from functools import partial


"this is create for IJOC"
"The codes are updated for Gurobi 11: Model.getAttr and Model.setAttr now raise an exception if the modeling objects passed to them do not belong to the model. Previously, the attribute values would be returned or set for a variable on the model with the same index as the modeling object which is typically not the desired behavior."
"full formulation with enumerating y variable"


def subtourElim_Vehicle(model, where):
    if where == GRB.Callback.MIPNODE and userCutsFlag == 1:
        node_count = model.cbGet(GRB.Callback.MIPNODE_NODCNT)
        if node_count == 0 and model.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.OPTIMAL:
            foundNewRow = 0
            x_val = model.cbGetNodeRel(model._x)
            z_val = model.cbGetNodeRel(model._z)
            t_val = model.cbGetNodeRel(model._t)
            Gs = gt.Graph()
            Gs.add_vertex(n + n)
            arc_weight = Gs.new_edge_property("double")

            "depot node 0 only has no arcs"
            for i in range(1, n):
                for j in range(1, n):
                    if i != j and (x_val[i, j] + z_val[i, j]) > 0.001:
                        arc = Gs.add_edge(i, j)
                        arc_weight[arc] = max(1 - x_val[i, j] - z_val[i, j], 0)
            "node n to 2n-1 only has incoming arcs and depot node 0 cannot go to these nodes directly"
            for i in range(1, n):
                for j in range(n + 1, 2 * n):
                    if i != j - n and (x_val[i, j - n] + z_val[i, j - n]) > 0.001:
                        arc = Gs.add_edge(i, j)
                        arc_weight[arc] = max(1 - x_val[i, j - n] - z_val[i, j - n], 0)
                        # arc_weight[arc] = 1 - x_val[i, j - n] - (ySum[i] + ySum[j - n]) / 2 +1
            # nodeVisitedSet = set()
            for i in range(n):
                # if i not in nodeVisitedSet:
                dist, path_map = gt.shortest_distance(Gs, source=Gs.vertex(i), target=Gs.vertex(i + n),
                                                      weights=arc_weight, pred_map=True, negative_weights=True)
                path = gt.shortest_path(Gs, source=Gs.vertex(i), target=Gs.vertex(i + n), pred_map=path_map)
                cycle = list()

                if dist < 1 - epsilon:
                    for node in range(len(path[0]) - 1):
                        cycle.append(Gs.vertex_index[path[0][node]])
                    # nodeVisitedSet.union(set(cycle))
                    foundNewRow = 1
                    for k in cycle:
                        # numNewRow += 1
                        Sk = cycle.copy()
                        Sk.remove(k)
                        model.cbCut(gp.quicksum(
                            model._x[j, l] + model._z[j, l] for j in cycle for l in cycle if
                            j != l) <= gp.quicksum(model._t[i] for i in Sk))

            "this is the exact separation procedure for fractional solutions"
            if foundNewRow == 0:
                sep.setObjective(
                    gp.quicksum((x_val[e] + z_val[e]) * u_sep[e] for e in u_sep.keys()) - gp.quicksum(
                        (t_val[i]) * v_sep[i] for i in v_sep.keys()) + gp.quicksum(
                        (t_val[i]) * r_sep[i] for i in r_sep.keys()), GRB.MAXIMIZE)
                sep.update()
                sep.optimize()

                if sep.objVal > epsilon:
                    v_sep_val = sep.getAttr('x', v_sep)
                    r_sep_val = sep.getAttr('x', r_sep)
                    Sk = list()
                    cycle = list()
                    for i in range(1, n):
                        if v_sep_val[i] == 1:
                            cycle.append(i)
                            Sk.append(i)
                            if r_sep_val[i] == 1:
                                Sk.remove(i)
                    model.cbCut(gp.quicksum(
                        model._x[j, l] + model._z[j, l] for j in cycle for l in cycle if
                        j != l) <= gp.quicksum(model._t[i] for i in Sk))
    elif where == GRB.Callback.MIPSOL and lazyConsFlag == 1:
        # make a list of edges selected in the solution
        x_val = model.cbGetSolution(model._x)
        z_val = model.cbGetSolution(model._z)
        t_val = model.cbGetSolution(model._t)


        selectedX = list((i, j) for i, j in model._x.keys() if x_val[i, j] > 0.5)
        selectedZ = list((i, j) for i, j in model._z.keys() if z_val[i, j] > 0.5)

        Gs = gt.Graph()
        Gs.add_edge_list(selectedX)
        Gs.add_edge_list(selectedZ)
        findRowFlag = 0
        for cycle in gt.all_circuits(Gs):
            cycle=cycle.tolist()
            # if len(cycle) < n - Qsum:
            if 0 not in cycle:
                findRowFlag = 1
                "CAVAD SSE"
                for k in cycle:
                    Sk = cycle.copy()
                    Sk.remove(k)
                    model.cbLazy(gp.quicksum((model._x[j, i]+model._z[j, i]) for j in cycle for i in cycle if j != i) <= gp.quicksum(model._t[i] for i in Sk))

        if findRowFlag == 0:
            print("Lazy Constraint found a feasible solution")

if __name__ == '__main__':
    version = "-Enumeration"
    outputFlag = 1
    solutionCheckingFlag = 1
    disaggregateCon8Flag = 0
    userCutsFlag = 1
    lazyConsFlag = 0
    terminateGap = 0.05
    epsilon = 0.01
    LOG = 1
    printData = 0
    printRouteFlag = 0
    B = 1
    f = 2.8
    instanceList = []
    numFiles = 10
    # for name in ['Cook']:

    for name in ['Cook','Winnebago','Champaign','LaSalle','Adams','Fulton','Jefferson','Johnson','Cumberland']:
        for fileID in range(0, numFiles):
            instanceList.append((name, fileID+1))

    # for B in [2]:
    for B in [1, 2, 3, 4, 5]:

        for instance in instanceList:
            sepTime, pricingTime, enumeratingTime, TimeBE, TimeAE, BE_UB, Global_y_numberBE = 0, 0, 0, 0, 0, 0, 0
            name = instance[0]
            fileID = instance[1]
            print("BCP-CAVAD" + "_B_" + str(B) + "_County_" + str(name) + "_ID_" + str(fileID))

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
            "Checking drivingAlwaysFasterFlag"
            for i, j in walkingTime.keys():
                if drivingTime[i, j] > walkingTime[i, j]:
                    drivingAlwaysFasterFlag = 0
                    print("driving is NOT always faster!!")
                    break

            "Building restricted walking graph"
            rwG_ind = np.zeros((n,n))
            rwG = ig.Graph(directed=True)
            rwG.add_vertices(n)
            for i, j in walkingTime.keys():
                if drivingAlwaysFasterFlag == 1:
                    if i != 0 and j != 0 and i != j and f + drivingTime[i, j] > walkingTime[i, j]:
                        rwG.add_edge(i, j)
                        rwG_ind[i,j] = 1
                else:
                    if i != 0 and j != 0 and i != j:
                        rwG.add_edge(i, j)
                        rwG_ind[i, j] = 1

            ijPath = dict()
            for i in range(n):
                ijPath[i]=set()


            for i in range(1,n):
                ijPath[i].add(i)
                for j in range(1,n):
                    if i!=j:
                        results = rwG.get_shortest_paths(i, to=j, output="vpath")
                        if drivingAlwaysFasterFlag == 0:
                            if results[0] != []:
                                ijPath[i].add(j)
                        else:
                            if len(results[0]) <= B:
                                ijPath[i].add(j)

            start_time = time.time()
            "this part is for the exact separation problem"
            sep = gp.Model("sepEx")
            u_sep = gp.tupledict()
            for j in range(1, n):
                for k in range(1, n):
                    if j != k:
                        u_sep[j, k] = sep.addVar(vtype=GRB.BINARY, name='u[%d,%d]' % (j, k))

            v_sep = gp.tupledict()
            for k in range(1, n):
                v_sep[k] = sep.addVar(vtype=GRB.BINARY, name='v[%d]' % k)

            r_sep = gp.tupledict()
            for k in range(1, n):
                r_sep[k] = sep.addVar(vtype=GRB.BINARY, name='gamma[%d]' % k)

            sep.addConstrs(u_sep[e] <= v_sep[e[0]] for e in u_sep.keys())
            sep.addConstrs(u_sep[e] <= v_sep[e[1]] for e in u_sep.keys())
            sep.addConstrs(r_sep[i] <= v_sep[i] for i in v_sep.keys())
            sep.addConstr(gp.quicksum(v_sep[i] for i in v_sep.keys()) >= 2)
            sep.addConstr(gp.quicksum(r_sep[i] for i in r_sep.keys()) == 1)
            sep.Params.OutputFlag = 0


            Global_x_val = dict()
            Global_y_val = dict()
            Global_z_val = dict()
            Global_t_val = dict()
            Global_y_costCoeff = dict()
            Global_P = dict()
            Global_Pdetail = dict()
            Global_LB = -np.inf
            Global_UB = np.inf
            master = gp.Model("CAVADP_master")

            t = master.addVars(range(1, n), obj=0, vtype=GRB.BINARY, name='t')
            tInd = dict((i, t[i].index) for i in t.keys())

            x = gp.tupledict()
            for i in range(n):
                for j in range(n):
                    if i!=j:
                        x[i, j] = master.addVar(obj=drivingTime[i, j], vtype=GRB.BINARY, name='x[%d,%d]' % (i, j))
            xInd = dict((i, x[i].index) for i in x.keys())

            z = gp.tupledict()
            for i in range(1,n):
                for j in range(1,n):
                    if i != j:
                        z[i, j] = master.addVar(obj= 0, vtype=GRB.BINARY, name='z[%d,%d]' % (i, j))
            zInd = dict((i, z[i].index) for i in z.keys())

            v = gp.tupledict()
            for i in range(n):
                for j in range(1, n):
                    if i != j:
                        v[i, j] = master.addVar(obj=0, vtype=GRB.INTEGER, name='v[%d,%d]' % (i, j))
            vInd = dict((i, v[i].index) for i in v.keys())

            "Initialized a set of y variables"
            y = gp.tupledict()  # initialize a set of y variables
            P = dict()  # paths
            Pdetail = dict()
            reverseP = dict()  # contains the set of loops that contains node i in Q
            for i in range(n):
                reverseP[i] = set()
            Psize = np.zeros((n, n), dtype=int)
            QpSize = dict()
            numQ = 0
            if drivingAlwaysFasterFlag == 1:
                "claim 2"
                for i in range(1, n):
                    "Lemma 2"
                    Q = [i]
                    P[i, i, Psize[i, i]] = Q
                    Pdetail[i, i, Psize[i, i]] = []
                    cpValue = f
                    y[i, i, Psize[i, i]] = master.addVar(obj=cpValue, vtype=GRB.BINARY,
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
                            Pdetail[i, j, Psize[i, j]] = [(path[node_id], path[node_id+1]) for node_id in range(len(path)-1) ]
                            cpValue = f + w_hatTime
                            y[i, j, Psize[i, j]] = master.addVar(obj=cpValue, vtype=GRB.BINARY,
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
                                    y[i, j, Psize[i, j]] = master.addVar(obj=cpValue, vtype=GRB.BINARY,
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
                                        y[i, i, Psize[i, i]] = master.addVar(obj=cpValue, vtype=GRB.BINARY,
                                                                         name='y[%d,%d,%d]' % (i, i, Psize[i, i]))
                                        Pdetail[i, i, Psize[i, i]] = [(path[node_id], path[node_id + 1]) for node_id in
                                                                      range(len(path) - 1)]
                                        Pdetail[i, i, Psize[i, i]].append((j,i))
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
                                y[i, j, Psize[i, j]] = master.addVar(obj=cpValue, vtype=GRB.BINARY,
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
                                y[i, i, Psize[i, i]] = master.addVar(obj=cpValue, vtype=GRB.BINARY,
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
                                                    y[i, j, Psize[i, j]] = master.addVar(obj=cpValue, vtype=GRB.BINARY,
                                                                                     name='y[%d,%d,%d]' % (
                                                                                     i, j, Psize[i, j]))
                                                    Pdetail[i, j, Psize[i, j]] = [(path[node_id], path[node_id + 1]) for
                                                                                  node_id in
                                                                                  range(len(path) - 1)]
                                                    if j != path[-1]:
                                                        Pdetail[i, j, Psize[i, j]].append((path[-1], j))
                                                    if i != path[0]:
                                                        Pdetail[i, j, Psize[i, j]].insert(0, (i, path[0]))
                                                    QpSize[i, j, Psize[i, j]] = len(Q)
                                                    for k in Q:
                                                        reverseP[k].add((i, j, Psize[i, j]))
                                                    Psize[i, j] += 1
                                            else:
                                                w_hatTime = walkingTime[i, j]
                                                P[i, j, Psize[i, j]] = Q
                                                cpValue = f + max(w_hatTime, drivingTime[i, j])
                                                y[i, j, Psize[i, j]] = master.addVar(obj=cpValue, vtype=GRB.BINARY,
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

            master.update()
            # adding constraints
            "Con1"
            conSet1 = master.addConstrs(x.sum('*', i)== 1 for i in [0])
            conSet1Ind = dict((i, conSet1[i].index) for i in conSet1.keys())
            "Con2"
            conSet2 = master.addConstrs(x.sum(i, '*') == 1 for i in [0])
            conSet2Ind = dict((i, conSet2[i].index) for i in conSet2.keys())

            "Con5"
            conSet5 = master.addConstrs(x.sum('*', i) + z.sum('*', i) == x.sum(i, '*') + z.sum(i, '*') for i in range(1, n))
            conSet5Ind = dict((i, conSet5[i].index) for i in conSet5.keys())
            "Con6"
            conSet6 = master.addConstrs(gp.quicksum(y[p] for p in reverseP[i]) == 1 for i in range(1, n))
            conSet6Ind = dict((i, conSet6[i].index) for i in conSet6.keys())

            "Con9"
            conSet9 = master.addConstrs(
                gp.quicksum(y[i, j, p] for p in range(Psize[i, j])) == z[i, j] for i, j in z.keys())
            conSet9Ind = dict((i, conSet9[i].index) for i in conSet9.keys())

            "SCF"
            master.addConstrs(v[0, i] <= cn * x[0, i] for i in range(1, n))
            master.addConstrs(
                v[i, k] <= cn * (x[i, k] + z[i, k]) for i in range(1, n) for k in range(1, n) if i != k)
            conSet12 = master.addConstrs(v.sum('*', i) - v.sum(i, '*') == gp.quicksum(
                y[i, j, p] * QpSize[i, j, p] for j in range(1, n) for p in range(Psize[i, j])) for i in range(1, n))
            conSet12Ind = dict()

            if drivingAlwaysFasterFlag == 0:
                "Con8"
                if disaggregateCon8Flag == 1:
                    conSet8 = master.addConstrs(y[i, i, p] <=(x.sum('*', i) + z.sum('*', i)) for i in range(1, n) for p in range(Psize[i, i]) )
                else:
                    conSet8 = master.addConstrs(gp.quicksum(y[i, i, p] for p in range(Psize[i, i]))<=cn*(x.sum('*', i) + z.sum('*', i)) for i in range(1, n)  )
                conSet8Ind = dict((i, conSet8[i].index) for i in conSet8.keys())
                conSet11 = dict()
                conSet11Ind = dict()
            else:
                "Con8"
                conSet8 = master.addConstrs(gp.quicksum(y[i, i, p] for p in range(Psize[i, i]))<=x.sum('*', i) + z.sum('*', i) for i in range(1, n))
                conSet8Ind = dict((i, conSet8[i].index) for i in conSet8.keys())
                "Con3"
                conSet3 = master.addConstrs(x[i, j] + z[i, j] + x[j, i] + z[j, i] - t[k] <= 0 for i in range(1, n) for j in range(i + 1, n) for k in [i, j])
                conSet3Ind = dict((i, conSet3[i].index) for i in conSet3.keys())
                conSet10 = master.addConstrs(t[i] == x.sum('*', i) + z.sum('*', i) for i in range(1, n))
                conSet10Ind = dict((i, conSet10[i].index) for i in conSet10.keys())
                conSet11 = master.addConstrs(x.sum('*', i) <= 1 for i in range(1, n))
                conSet11Ind = dict((i, conSet11[i].index) for i in conSet11.keys())

                "variable fixing"
                for i in range(1,n):
                    customerTooFarFlag =1
                    for j in range(1,n):
                        if i!=j:
                            if drivingTime[i,j] +f > walkingTime[i,j] or drivingTime[j,i] +f > walkingTime[j,i]:
                                customerTooFarFlag =0
                                break
                    if customerTooFarFlag ==1:
                        master.addConstr(x.sum(i, '*') == 1)
                        master.addConstr(x.sum( '*' ,i ) == 1)
                        print("fix customer "+str(i))


            master.update()
            master._x = x
            master._y = y
            master._z = z
            master._t = t
            master.Params.OutputFlag = 0  # silent mode
            if LOG:
                master.Params.OutputFlag = 1  # verbose mode
            relax = master.relax()
            relax.Params.OutputFlag = 0

            if drivingAlwaysFasterFlag == 0:
                master.optimize()
            else:
                master.Params.PreCrush = 1
                master.optimize(subtourElim_Vehicle)



            if master.objVal < Global_UB:

                Global_x_val = master.getAttr('x', x)
                Global_y_val = master.getAttr('x', master._y)
                Global_z_val = master.getAttr('x', z)
                Global_y_costCoeff = master.getAttr('obj', master._y)
                Global_P = P
                Global_Pdetail = Pdetail
                Global_UB = master.objVal
                Global_Gap = (Global_UB - Global_LB) / Global_UB
                print('The number of y variables: ' + str(len(yInd)))

            x_sol = gp.tupledict()
            y_sol = gp.tupledict()
            z_sol = gp.tupledict()


            print('After enumerating, CAVAD Optimal Cost of the Master: %g' % Global_UB)
            print("Global UB: " + str(Global_UB))
            print("--- Running Time %s seconds ---" % (time.time() - start_time))


            if outputFlag ==1:
                file = open("CAVAD-Full-V" + str(version) + "_CAVAD_County_" + str(name) + "_N" + str(n)  + "_TGap" + str(terminateGap) + "_UB.log", 'a')
                file.write("\n" + str(Global_UB))
                file.close()
                file = open("CAVAD-Full-V" + str(version) + "_CAVAD_County_" + str(name) + "_N" + str(n)  + "_TGap" + str(terminateGap) + "_LB.log", 'a')
                file.write("\n" + str(Global_LB))
                file.close()
                file = open("CAVAD-Full-V" + str(version) + "_CAVAD_County_" + str(name) + "_N" + str(n)  + "_TGap" + str(terminateGap) + "_GlobalGap.log", 'a')
                file.write("\n" + str(Global_Gap))
                file.close()
                file = open("CAVAD-Full-V" + str(version) + "_CAVAD_County_" + str(name) + "_N" + str(n) + "_TGap" + str(
                    terminateGap) + "_TotalRunningTime.log", 'a')
                file.write("\n" + str(time.time() - start_time))
                file.close()
                file = open(
                    "CAVAD-Full-V" + str(version) + "_CAVAD_County_" + str(name) + "_N" + str(n) + "_TGap" + str(terminateGap) + "_NumberOfyVariables.log", 'a')
                file.write("\n" + str(len(Global_y_val)))
                file.close()
                file = open(
                    "CAVAD-Full-V" + str(version) + "_CAVAD_County_" + str(name) + "_N" + str(n) + "_TGap" + str(terminateGap) + "_SeparationTime.log", 'a')
                file.write("\n" + str(sepTime))
                file.close()
                file.close()
