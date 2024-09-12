import sys
import time
import math
import random
import numpy as np
from itertools import combinations
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


"this is create for IJOC R1"
"The codes are updated for Gurobi 11: Model.getAttr and Model.setAttr now raise an exception if the modeling objects passed to them do not belong to the model. Previously, the attribute values would be returned or set for a variable on the model with the same index as the modeling object which is typically not the desired behavior."

class Label:
    def __init__(self, load=0, visited_flag=None, reduced_cost=0.0, semi_reduced_cost=0.0, path=None, serviced = None, what =0, dual_cost_hat = 0):
        # self.time = time
        self.load = load
        self.visited_flag = visited_flag if visited_flag is not None else []
        self.reduced_cost = reduced_cost
        self.semi_reduced_cost = semi_reduced_cost
        self.path = path if path is not None else []
        self.serviced = serviced if serviced is not None else set()
        self.what = what
        self.dual_cost_hat = dual_cost_hat

def initialize_label(origin, n_nodes):
    visited_flag = [0]*(n_nodes)
    visited_flag[origin] =1
    return Label(0, visited_flag, np.inf, np.inf, [origin], set(), 0 ,0)

def should_rollback(dual_cost, p, n, cn, B, walkingTime, drivingTime, f):
    if len(p.path) < 3:
        return False
    v_i, v_k, v_j = p.path[-3], p.path[-2], p.path[-1]
    reduced_cost_p = walkingTime[v_i, v_k] + walkingTime[v_k, v_j] - dual_cost[v_k] - dual_cost[v_j]
    reduced_cost_pp = walkingTime[v_i, v_j]- dual_cost[v_j]
    return reduced_cost_pp <= reduced_cost_p

def isbounded(dual_val, p, primal_bounds, lower_bounds, dc_lower_bounds, root, destination, n, cn, B, conSet8Ind, conSet9Ind, conSet12Ind, walkingTime, drivingTime, f, threshold, collect_flag, dual_cost, orgin_served):
    v_i = p.path[-1]

    bound = np.inf
    bound = min(bound, primal_bounds[root].reduced_cost)


    k = B-p.load
    if lower_bounds[v_i, k].semi_reduced_cost < np.inf:
        if v_i == destination-cn:
            offset = (dual_val[conSet8Ind[v_i]]+dual_val[conSet12Ind[v_i]])
        else:
            offset = (dual_val[conSet9Ind[v_i,destination-cn]]+dual_val[conSet12Ind[v_i]])

        if orgin_served == True:
            offset += dual_cost[v_i]

        if p.what + lower_bounds[v_i, k].what >= drivingTime[root, destination-cn]:
            return f+ p.what + lower_bounds[v_i, k].what -(p.dual_cost_hat + lower_bounds[v_i, k].dual_cost_hat-offset) > bound
        else:
            return f + drivingTime[root, destination - cn] - (
                        p.dual_cost_hat + dc_lower_bounds[v_i, k].dual_cost_hat - offset) > bound
    else:
        return False

def pulse_procedure(dual_val, rwG_E, dual_cost, origin, destination, p, primal_bounds, lower_bounds, dc_lower_bounds, below_threshold_routes, init_capacity, root, destination_served, threshold, collect_flag, n, cn, B, conSet8Ind, conSet9Ind, conSet12Ind, walkingTime, drivingTime, f, orgin_served):
    "the last node of the current path"
    v_i = p.path[-1]
    if v_i == destination:
        p.reduced_cost = f+max(p.what, drivingTime[origin, destination - cn]) - p.dual_cost_hat
        p.semi_reduced_cost = p.what - p.dual_cost_hat
        p.path[-1] = destination -cn
        if p.reduced_cost < primal_bounds[root].reduced_cost:
            primal_bounds[root] = p

        if p.semi_reduced_cost < lower_bounds[root, init_capacity].semi_reduced_cost :
            lower_bounds[root, init_capacity] = p

        if p.dual_cost_hat > dc_lower_bounds[root, init_capacity].dual_cost_hat :
            dc_lower_bounds[root, init_capacity] = p

        if collect_flag == True:
            if p.reduced_cost<= threshold:
                below_threshold_routes.append(p)
            "The place changes for collecting different size of enumerating paths."

            # if 0.9*threshold <=p.reduced_cost<= threshold:
            #     below_threshold_routes.append(p)

        return

    if isbounded(dual_val, p, primal_bounds, lower_bounds, dc_lower_bounds, root, destination, n, cn, B, conSet8Ind, conSet9Ind, conSet12Ind, walkingTime, drivingTime, f, threshold, collect_flag, dual_cost, orgin_served):
        return


    if should_rollback(dual_cost, p, n, cn, B, walkingTime, drivingTime, f):
        return
    if p.load < B:
        for node in rwG_E[v_i]:
            v_j = node[0]
            if v_j >= n and v_j !=destination:
                continue
            if p.visited_flag[v_j] == 0:
                new_load = p.load +1
                # if new_load <= B:
                if v_j >=n:
                    new_what = p.what + walkingTime[v_i, v_j-cn]

                    if origin != destination -cn:
                        new_dual_cost_hat = p.dual_cost_hat + dual_cost[v_j-cn]
                    else:
                        new_dual_cost_hat = p.dual_cost_hat
                else:
                    new_what = p.what + walkingTime[v_i, v_j]
                    new_dual_cost_hat = p.dual_cost_hat + dual_cost[v_j]
                new_reduced_cost = np.inf
                # f+max(new_what, drivingTime[origin, destination-cn]) - new_dual_cost_hat
                pp = Label(new_load, copy(p.visited_flag), new_reduced_cost, np.inf, copy(p.path), copy(p.serviced), new_what, new_dual_cost_hat)
                pp.path.append(v_j)
                if new_dual_cost_hat != p.dual_cost_hat:
                    if v_j == destination:
                        pp.serviced.add(v_j - cn)
                    else:
                        pp.serviced.add(v_j)
                pp.visited_flag[v_j] = 1
                pulse_procedure(dual_val, rwG_E, dual_cost, origin, destination, pp, primal_bounds, lower_bounds, dc_lower_bounds, below_threshold_routes, init_capacity, root, destination_served, threshold, collect_flag, n, cn, B, conSet8Ind, conSet9Ind, conSet12Ind, walkingTime, drivingTime, f, orgin_served)

    if p.load == B and origin == destination-cn:
        v_j = destination
        new_load = p.load
        new_what = p.what + walkingTime[v_i, v_j - cn]
        new_dual_cost_hat = p.dual_cost_hat
        new_reduced_cost = np.inf
        #
        # new_reduced_cost = f+max(new_what, drivingTime[origin, destination - cn]) - new_dual_cost_hat
        pp = Label(new_load, copy(p.visited_flag), new_reduced_cost, np.inf, copy(p.path), copy(p.serviced), new_what, new_dual_cost_hat)
        pp.path.append(v_j)
        pp.visited_flag[v_j] = 1
        pulse_procedure(dual_val, rwG_E, dual_cost, origin, destination, pp, primal_bounds, lower_bounds, dc_lower_bounds, below_threshold_routes,init_capacity, root, destination_served, threshold, collect_flag, n, cn, B, conSet8Ind, conSet9Ind, conSet12Ind, walkingTime, drivingTime, f, orgin_served)

    if destination_served == False:
        if p.load == B and origin != destination-cn and destination in rwG_E[v_i]:
            v_j = destination
            new_load = p.load
            new_what = p.what + walkingTime[v_i, v_j - cn]
            new_dual_cost_hat = p.dual_cost_hat
            new_reduced_cost = np.inf

            # new_reduced_cost = f+max(new_what, drivingTime[origin, destination - cn]) - new_dual_cost_hat
            pp = Label(new_load, deepcopy(p.visited_flag), new_reduced_cost, np.inf, deepcopy(p.path), deepcopy(p.serviced), new_what, new_dual_cost_hat)
            pp.path.append(v_j)
            pp.visited_flag[v_j] = 1
            pulse_procedure(dual_val, rwG_E, dual_cost, origin, destination, pp, primal_bounds, lower_bounds, dc_lower_bounds, below_threshold_routes, init_capacity, root, destination_served, threshold, collect_flag, n, cn, B, conSet8Ind, conSet9Ind, conSet12Ind, walkingTime, drivingTime, f, orgin_served)

def pulse_with_fixed_destination(destination, dual_cost, dual_val, rwG_E, threshold, collect_flag, orgin_served, destination_served, n, cn, B, conSet8Ind, conSet9Ind, conSet12Ind, walkingTime, drivingTime, f):
    "store the best labels from i to destination"
    below_threshold_routes = []
    most_neg_cost_routes = []
    primal_bounds = {}
    for i in range(1, n):
        primal_bounds[i] = initialize_label(i, n + cn)

    "store the best labels from i to destination when remaining capacity is k"
    lower_bounds = {}
    for i in range(1, n):
        for k in range(B + 1):
            lower_bounds[i, k] = initialize_label(i, n + cn)

    dc_lower_bounds = {}
    for i in range(1, n):
        for k in range(B + 1):
            dc_lower_bounds[i, k] = initialize_label(i, n + cn)

    "bounding procedure for remaing capacity form 0 to B"
    "bounding procedure serving node v_i: 0, ...., B-1"
    "bounding procedure not serving node v_i: 0, ...., B"
    if orgin_served == True:
        remainingCapacity = B
    else:
        remainingCapacity = B +1
    for k in range(remainingCapacity):
        for v_i in range(1, n):
            p = initialize_label(v_i, n + cn)
            p.load = B - k
            if orgin_served == True:
                p.dual_cost_hat = dual_cost[v_i]
                p.serviced.add(v_i)

            if v_i == destination - cn:
                p.dual_cost_hat += (dual_val[conSet8Ind[destination - cn]]+dual_val[conSet12Ind[v_i]])
            else:
                p.dual_cost_hat += (dual_val[conSet9Ind[v_i, destination - cn]]+dual_val[conSet12Ind[v_i]])


            pulse_procedure(dual_val, rwG_E, dual_cost, v_i, destination, p, primal_bounds, lower_bounds, dc_lower_bounds,
                            below_threshold_routes, k, v_i, destination_served, threshold, collect_flag, n, cn, B, conSet8Ind, conSet9Ind, conSet12Ind, walkingTime, drivingTime, f, orgin_served)


    if collect_flag == True:
        return below_threshold_routes
    else:
        for i in range(1, n):
            if primal_bounds[i].reduced_cost < 0:
                most_neg_cost_routes.append(primal_bounds[i])
        return most_neg_cost_routes

def CAVAD_Pulse_Pricing(master, relax, yInd, P, Pdetail, QpSize, reverseP, Psize, conSet8Ind, conSet9Ind, conSet12Ind):
    dual_val = relax.getAttr("Pi")
    # y = master._y
    foundNewColumn = 0
    numNewColumn = 0
    "dual cost in the objective value except the beta and gamma term "
    dual_cost = {}
    for l in range(1,n):
        if drivingAlwaysFasterFlag == 1:
            dual_cost[l] = (dual_val[conSet6Ind[l]] + dual_val[conSet11Ind[l]])
        else:
            dual_cost[l] = (dual_val[conSet6Ind[l]] )

    rwG_E = dict()
    for i in range(1,n):
        rwG_E[i]=list()
        if drivingAlwaysFasterFlag == 1:
            for j in rwG.successors(i):
                rwG_E[i].append((j,  walkingTime[i,j]+dual_cost[j]))
                rwG_E[i].append((j+cn, walkingTime[i,j]+dual_cost[j]))
        else:
            for j in range(1,n):
                if i!=j:
                    rwG_E[i].append((j,  walkingTime[i,j]+dual_cost[j]))
                    rwG_E[i].append((j+cn, walkingTime[i,j]+dual_cost[j]))
        rwG_E[i].sort(key=lambda x: x[1])

    threshold = 0
    collect_flag = False

    destination_arg = range(1+cn, n+cn)
    pool = mp.Pool(mp.cpu_count())
    # pool = mp.Pool(2)
    if drivingAlwaysFasterFlag == 1:
        orgin_served = True
        destination_served = True
        "both start and end nodes are served"
        results = pool.map(partial(pulse_with_fixed_destination, dual_cost=dual_cost, dual_val=dual_val, rwG_E=rwG_E, threshold=threshold,
                                         collect_flag=collect_flag, orgin_served=orgin_served, destination_served=destination_served, n=n, cn=cn, B=B, conSet8Ind=conSet8Ind, conSet9Ind=conSet9Ind, conSet12Ind=conSet12Ind, walkingTime= walkingTime, drivingTime = drivingTime, f=f), destination_arg)

    else:
        destination_served = False # if it is False, end node is not served

        "start node is served and end node might not be served"
        orgin_served = True
        results1 = pool.map(partial(pulse_with_fixed_destination, dual_cost=dual_cost, dual_val=dual_val, rwG_E=rwG_E,
                         threshold=threshold,
                         collect_flag=collect_flag,
                         orgin_served=orgin_served, destination_served=destination_served, n=n, cn=cn, B=B, conSet8Ind=conSet8Ind, conSet9Ind=conSet9Ind, conSet12Ind=conSet12Ind, walkingTime= walkingTime, drivingTime = drivingTime, f=f), destination_arg)

        "start node is not served and end node might not be served"
        orgin_served = False
        results2 = pool.map(partial(pulse_with_fixed_destination, dual_cost=dual_cost, dual_val=dual_val, rwG_E=rwG_E,
                         threshold=threshold,
                         collect_flag=collect_flag,
                         orgin_served=orgin_served, destination_served=destination_served, n=n, cn=cn, B=B, conSet8Ind=conSet8Ind, conSet9Ind=conSet9Ind, conSet12Ind=conSet12Ind, walkingTime= walkingTime, drivingTime = drivingTime, f=f), destination_arg)

        results = results1 + results2
    for most_neg_cost_routes in results:
        for p in most_neg_cost_routes:
            i = p.path[0]
            j = p.path[-1]

            Q = p.serviced
            old_path_flag = 0
            for pid in range(Psize[i, j]):
                if P[i, j, pid] == Q:
                    old_path_flag = 1
                    break
            if old_path_flag == 1:
                continue
            cp_val = f + max(p.what, drivingTime[i, j])
            col = gp.Column()
            for k in Q:
                col.addTerms(1, master.getConstrs()[conSet6Ind[k]])
                if drivingAlwaysFasterFlag == 1:
                    if k != i and k !=j:
                        col.addTerms(1, master.getConstrs()[conSet11Ind[k]])

            if i != j:
                col.addTerms(1, master.getConstrs()[conSet9Ind[i, j]])
            else:
                col.addTerms(1, master.getConstrs()[conSet8Ind[i]])

            col.addTerms(-len(Q), master.getConstrs()[conSet12Ind[i]])
            numNewColumn += 1
            foundNewColumn = 1
            master._y[i, j, Psize[i,j]] = master.addVar(obj=cp_val, vtype=GRB.BINARY, name='y[%d,%d,%d]' % (i, j, Psize[i,j]), column=col)
            yInd[i, j, Psize[i,j]] = master._y[i, j, Psize[i, j]].index

            master.update()
            P[i, j, Psize[i,j]] = Q
            Qdetail = [(p.path[node_id], p.path[node_id+1]) for node_id in range(len(p.path)-1) ]
            Pdetail[i, j, Psize[i, j]] = Qdetail
            QpSize[i, j, Psize[i,j]] = len(Q)
            for k in Q:
                reverseP[k].add((i, j, Psize[i,j]))
            Psize[i,j] += 1


    print("number of new columns: "+str(numNewColumn))
    return foundNewColumn, master, relax, yInd, P, Pdetail, QpSize, reverseP, Psize, numNewColumn

def CAVAD_Pulse_Enumerating(master, relax, yInd, P, Pdetail, QpSize, reverseP, Psize, conSet8Ind, conSet9Ind, conSet12Ind, Threshold):
    dual_val = relax.getAttr("Pi")
    foundNewColumn = 0
    numNewColumn = 0
    "dual cost in the objective value except the beta and gamma term "
    dual_cost = {}
    for l in range(1,n):
        if drivingAlwaysFasterFlag == 1:
            dual_cost[l] = (dual_val[conSet6Ind[l]] + dual_val[conSet11Ind[l]])
        else:
            dual_cost[l] = (dual_val[conSet6Ind[l]] )

    rwG_E = dict()
    for i in range(1,n):
        rwG_E[i]=list()
        if drivingAlwaysFasterFlag == 1:
            for j in rwG.successors(i):
                rwG_E[i].append((j, walkingTime[i, j] + dual_cost[j]))
                rwG_E[i].append((j + cn, walkingTime[i, j] + dual_cost[j]))
        else:
            for j in range(1, n):
                if i != j:
                    rwG_E[i].append((j, walkingTime[i, j] + dual_cost[j]))
                    rwG_E[i].append((j + cn, walkingTime[i, j] + dual_cost[j]))
        rwG_E[i].sort(key=lambda x: x[1])

    threshold = Threshold
    collect_flag = True

    destination_arg = range(1+cn, n+cn)
    pool = mp.Pool(mp.cpu_count())
    # pool = mp.Pool(2)
    if drivingAlwaysFasterFlag == 1:
        orgin_served = True
        destination_served = True
        "both start and end nodes are served"
        results = pool.map(partial(pulse_with_fixed_destination, dual_cost=dual_cost, dual_val=dual_val, rwG_E=rwG_E, threshold=threshold,
                                         collect_flag=collect_flag, orgin_served=orgin_served, destination_served=destination_served, n=n, cn=cn, B=B, conSet8Ind=conSet8Ind, conSet9Ind=conSet9Ind, conSet12Ind=conSet12Ind, walkingTime= walkingTime, drivingTime = drivingTime, f=f), destination_arg)


    else:
        destination_served = False # if it is False, end node is not served

        "start node is served and end node might not be served"
        orgin_served = True
        results1 = pool.map(partial(pulse_with_fixed_destination, dual_cost=dual_cost, dual_val=dual_val, rwG_E=rwG_E,
                         threshold=threshold,
                         collect_flag=collect_flag,
                         orgin_served=orgin_served, destination_served=destination_served, n=n, cn=cn, B=B, conSet8Ind=conSet8Ind, conSet9Ind=conSet9Ind, conSet12Ind=conSet12Ind, walkingTime= walkingTime, drivingTime = drivingTime, f=f), destination_arg)

        "start node is not served and end node might not be served"
        orgin_served = False
        results2 = pool.map(partial(pulse_with_fixed_destination, dual_cost=dual_cost, dual_val=dual_val, rwG_E=rwG_E,
                         threshold=threshold,
                         collect_flag=collect_flag,
                         orgin_served=orgin_served, destination_served=destination_served, n=n, cn=cn, B=B, conSet8Ind=conSet8Ind, conSet9Ind=conSet9Ind, conSet12Ind=conSet12Ind, walkingTime= walkingTime, drivingTime = drivingTime, f=f), destination_arg)

        results = results1 + results2
    for most_neg_cost_routes in results:
        for p in most_neg_cost_routes:
            i = p.path[0]
            j = p.path[-1]

            Q = p.serviced
            old_path_flag = 0
            for pid in range(Psize[i, j]):
                if P[i, j, pid] == Q:
                    old_path_flag = 1
                    break
            if old_path_flag == 1:
                continue
            cp_val = f + max(p.what, drivingTime[i, j])
            col = gp.Column()
            for k in Q:
                col.addTerms(1, master.getConstrs()[conSet6Ind[k]])
                if drivingAlwaysFasterFlag == 1:
                    if k != i and k !=j:
                        col.addTerms(1, master.getConstrs()[conSet11Ind[k]])


            if i != j:
                col.addTerms(1, master.getConstrs()[conSet9Ind[i, j]])
            else:
                col.addTerms(1, master.getConstrs()[conSet8Ind[i]])
            col.addTerms(-len(Q), master.getConstrs()[conSet12Ind[i]])
            numNewColumn += 1
            foundNewColumn = 1
            master._y[i, j, Psize[i, j]] = master.addVar(obj=cp_val, vtype=GRB.BINARY,
                                                         name='y[%d,%d,%d]' % (i, j, Psize[i, j]), column=col)
            yInd[i, j, Psize[i, j]] = master._y[i, j, Psize[i, j]].index
            master.update()
            P[i, j, Psize[i,j]] = Q
            Qdetail = [(p.path[node_id], p.path[node_id+1]) for node_id in range(len(p.path)-1) ]
            Pdetail[i, j, Psize[i, j]] = Qdetail
            QpSize[i, j, Psize[i,j]] = len(Q)
            for k in Q:
                reverseP[k].add((i, j, Psize[i,j]))
            Psize[i,j] += 1


    print("number of new columns: "+str(numNewColumn))
    return foundNewColumn, master, relax, yInd, P, Pdetail, QpSize, reverseP, Psize, numNewColumn

def checkSolutionStaus(model):
    if model.status == GRB.OPTIMAL:
        print('Optimal objective: %g' % model.objVal)
    elif model.status == GRB.INF_OR_UNBD:
        print('Model is infeasible or unbounded')
        sys.exit(0)
    elif model.status == GRB.INFEASIBLE:
        print('Model is infeasible')
        sys.exit(0)
    elif model.status == GRB.UNBOUNDED:
        print('Model is unbounded')
        sys.exit(0)
    else:
        print('Optimization ended with status %d' % model.status)
        sys.exit(0)

# Calculate the total travel time
def totalTourTravelTimeCalculator(tour, travelTime):
    if tour[len(tour)-1] > tour[0]:
        totalTravelTime = travelTime[tour[len(tour)-1], tour[0]]
    else:
        totalTravelTime = travelTime[tour[0], tour[len(tour)-1]]
    for i in range(len(tour)-1):
        if tour[i] > tour[i + 1]:
            totalTravelTime += travelTime[tour[i], tour[i + 1]]
        else:
            totalTravelTime += travelTime[tour[i + 1], tour[i]]
    return totalTravelTime

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
            for i in range(n):
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
            if 0 not in cycle:
                findRowFlag = 1
                "CAVAD SSE"
                for k in cycle:
                    Sk = cycle.copy()
                    Sk.remove(k)
                    model.cbLazy(gp.quicksum((model._x[j, i]+model._z[j, i]) for j in cycle for i in cycle if j != i) <= gp.quicksum(model._t[i] for i in Sk))

        if findRowFlag == 0:
            print("Lazy Constraint found a feasible solution")

class Node:
    def __init__(self, master, relax, conSet1Ind,conSet2Ind, conSet3Ind, conSet4Ind,conSet5Ind, conSet6Ind, conSet7Ind, conSet8Ind,conSet9Ind,conSet10Ind, conSet11Ind, conSet12Ind, xInd, yInd, zInd, tInd, P, Pdetail, reverseP, Psize, QpSize, parentNode):
        self.master = master
        self.relax = relax
        self.conSet1Ind = conSet1Ind
        self.conSet2Ind = conSet2Ind
        self.conSet3Ind = conSet3Ind
        self.conSet4Ind = conSet4Ind
        self.conSet5Ind = conSet5Ind
        self.conSet6Ind = conSet6Ind
        self.conSet7Ind = conSet7Ind
        self.conSet8Ind = conSet8Ind
        self.conSet9Ind = conSet9Ind
        self.conSet10Ind = conSet10Ind
        self.conSet11Ind = conSet11Ind
        self.conSet12Ind = conSet12Ind
        # self.hInd = hInd
        self.xInd = xInd
        self.yInd = yInd
        self.zInd = zInd
        self.tInd = tInd
        self.P = P
        self.Pdetail = Pdetail
        self.reverseP = reverseP
        self.Psize = Psize
        self.QpSize = QpSize
        self.parentNode = parentNode
        self.LB = -np.infty

    def branching(self):
        t_val = self.relax.getAttr('x', t)
        "constraint branching t_i  = 1 or t_i  = 0"
        conBranchingID = 0
        fracAmonMost = 1
        for i in range(1,n):
            fracAmon = abs(t_val[i]  - 0.5)
            if fracAmon < fracAmonMost:
                fracAmonMost = fracAmon
                conBranchingID = i
        return  conBranchingID

def CAVAD_Row_Generation(master, relax, sep):
    sol_val = relax.getAttr('x')

    numNewRow = 0
    foundNewRow = 0
    integralityFlag = 1
    for p in yInd.keys():
        if sol_val[yInd[p]] > 0.001:
            # Qsum += QpSize[p] * y_val[p]
            if integralityFlag == 1 and abs(sol_val[yInd[p]] - 1) > 1e-5:
                integralityFlag = 0

    for i in range(n):
        for j in range(n):
            if i != j and sol_val[xInd[i, j]] > 0.001:
                if integralityFlag == 1 and abs(sol_val[xInd[i, j]] - 1) > 1e-5:
                    integralityFlag = 0

    for i in range(1,n):
        for j in range(1,n):
            if i != j and sol_val[zInd[i, j]] > 0.001:
                if integralityFlag == 1 and abs(sol_val[zInd[i, j]] - 1) > 1e-5:
                    integralityFlag = 0
    # FeasibleSolutionFlag = 0
    if integralityFlag == 1:
        selectedX = list((i, j) for i, j in xInd.keys() if sol_val[xInd[i, j]] > 0.5)
        selectedZ = list((i, j) for i, j in zInd.keys() if sol_val[zInd[i, j]] > 0.5)
        Gs = gt.Graph()
        Gs.add_edge_list(selectedX)
        Gs.add_edge_list(selectedZ)
        for cycle in gt.all_circuits(Gs):
            cycle=cycle.tolist()
            # if len(cycle) < n - Qsum:
            if 0 not in cycle:
                foundNewRow = 1

                for k in cycle:
                    numNewRow += 1
                    Sk = cycle.copy()
                    Sk.remove(k)
                    # print("Sk is "+str(Sk))
                    master.addConstr(gp.quicksum((master.getVars()[xInd[j,l]]+master.getVars()[zInd[j,l]]) for j in cycle for l in cycle if j != l) <= gp.quicksum(master.getVars()[tInd[i]] for i in Sk ))
            else:
                # FeasibleSolutionFlag = 1
                print("found a cycle with depot")
                # break
    else:
        Gs = gt.Graph()
        Gs.add_vertex(n + n)
        arc_weight = Gs.new_edge_property("double")

        "depot node 0 only has no arcs"
        for i in range(1,n):
            for j in range(1,n):
                if i != j and (sol_val[xInd[i, j]]+ sol_val[zInd[i, j]]) > 0.001:
                    arc = Gs.add_edge(i, j)
                    arc_weight[arc] = max(1 - sol_val[xInd[i, j]]-sol_val[zInd[i, j]],0)
        "node n to 2n-1 only has incoming arcs and depot node 0 cannot go to these nodes directly"
        for i in range(1,n):
            for j in range(n+1, 2 * n):
                if i != j - n and (sol_val[xInd[i, j - n]]+sol_val[zInd[i, j - n]]) > 0.001:
                    arc = Gs.add_edge(i, j)
                    arc_weight[arc] = max(1 - sol_val[xInd[i, j - n]]-sol_val[zInd[i, j - n]],0)
                    # arc_weight[arc] = 1 - x_val[i, j - n] - (ySum[i] + ySum[j - n]) / 2 +1
        # nodeVisitedSet = set()
        for i in range(n):
            # if i not in nodeVisitedSet:
            dist, path_map = gt.shortest_distance(Gs, source=Gs.vertex(i), target=Gs.vertex(i + n), weights=arc_weight, pred_map=True, negative_weights=True)
            path = gt.shortest_path(Gs, source=Gs.vertex(i), target=Gs.vertex(i + n), pred_map=path_map)
            cycle = list()

            if dist < 1-epsilon:
                for node in range(len(path[0]) - 1):
                    cycle.append(Gs.vertex_index[path[0][node]])
                # nodeVisitedSet.union(set(cycle))
                foundNewRow = 1
                for k in cycle:
                    numNewRow += 1
                    Sk = cycle.copy()
                    Sk.remove(k)
                    master.addConstr(gp.quicksum(master.getVars()[xInd[j,l]]+master.getVars()[zInd[j,l]] for j in cycle for l in cycle if j != l) <= gp.quicksum(master.getVars()[tInd[i]] for i in Sk ))

        "this is the exact separation procedure for fractional solutions"
        if foundNewRow ==0:
            sep.setObjective( gp.quicksum((sol_val[xInd[e]]+sol_val[zInd[e]]) * u_sep[e] for e in u_sep.keys())-gp.quicksum((sol_val[tInd[i]])*v_sep[i] for i in v_sep.keys())+gp.quicksum((sol_val[tInd[i]])*r_sep[i] for i in r_sep.keys()), GRB.MAXIMIZE)
            sep.update()
            sep.optimize()

            if sep.objVal > epsilon:
                print("found exact row separation")
                foundNewRow = 1
                numNewRow += 1
                v_sep_val = sep.getAttr('x', v_sep)
                r_sep_val = sep.getAttr('x', r_sep)
                Sk = list()
                cycle = list()
                for i in range(1,n):
                    if v_sep_val[i] ==1:
                        cycle.append(i)
                        Sk.append(i)
                        if r_sep_val[i] ==1:
                            Sk.remove(i)
                master.addConstr(gp.quicksum(master.getVars()[xInd[j,l]]+master.getVars()[zInd[j,l]] for j in cycle for l in cycle if j != l) <= gp.quicksum(master.getVars()[tInd[i]] for i in Sk))

    print("number of new rows: " + str(numNewRow))
    return foundNewRow, master

def CP_RoutineNode(node):
    global Global_LB, Global_UB, x_sol, y_sol, Global_Gap, Global_x_val, Global_y_val, Global_z_val, Global_r_val, Global_t_val, Global_y_costCoeff, Global_P, Global_Pdetail, sepTime, pricingTime, enumeratingTime, TimeBE, TimeAE, start_time, BE_UB, Global_y_numberBE, Global_y_numberAE
    master= node.master
    relax = node.relax
    conSet1Ind = node.conSet1Ind
    conSet2Ind = node.conSet2Ind
    conSet3Ind = node.conSet3Ind
    conSet4Ind = node.conSet4Ind
    conSet5Ind = node.conSet5Ind
    conSet6Ind = node.conSet6Ind
    conSet7Ind = node.conSet7Ind
    conSet8Ind = node.conSet8Ind
    conSet9Ind = node.conSet9Ind
    conSet10Ind = node.conSet10Ind
    conSet11Ind = node.conSet11Ind
    conSet12Ind = node.conSet12Ind
    # hInd = node.hInd
    xInd = node.xInd
    yInd = node.yInd
    zInd = node.zInd
    tInd = node.tInd
    P = deepcopy(node.P)
    Pdetail = deepcopy(node.Pdetail)
    reverseP = deepcopy(node.reverseP)
    Psize = deepcopy(node.Psize)
    QpSize = deepcopy(node.QpSize)
    parentNode = node.parentNode
    LB = -np.infty

    if drivingAlwaysFasterFlag == 1 and BPCFlag == 1:
        repeatFlag = 1
        repeatIter = 0
        numNewColumn = n
        rootNodeSP_Flag = 1
        rootNodeSP_Done = 0

        while repeatFlag:
            repeatIter += 1
            print("Repeat Full Iter: " + str(repeatIter))
            print("Global UB: " + str(Global_UB))
            print("Global LB: " + str(Global_LB))
            repeatFlag = 0
            # "flag for controlling row generation"
            "Column generation procedure"
            relax = master.relax()
            relax.optimize()
            iterCol = 0
            foundNewColumn = 1
            time0 = time.time()
            while foundNewColumn:
                # spExecuteFlag = 0
                foundNewColumn = 0
                iterCol += 1
                print("Repeat Full Iter: " + str(repeatIter) + " Column Generation Iter: " + str(iterCol))
                # master.write("MP_F"+str(repeatIter)+"_CG" + str(iter) + ".lp")
                "pricing procedure. Column Generation."
                if node.parentNode == -1:
                    foundNewColumn, master, relax, yInd, P, Pdetail, QpSize, reverseP, Psize, numNewColumn = CAVAD_Pulse_Pricing(
                        master, relax, yInd, P, Pdetail, QpSize, reverseP, Psize, conSet8Ind, conSet9Ind, conSet12Ind)

                else:
                    foundNewColumn, master, relax, yInd, P, Pdetail, QpSize, reverseP, Psize, numNewColumn = CAVAD_Pulse_Pricing(
                        master, relax, yInd, P, Pdetail, QpSize, reverseP, Psize, conSet8Ind, conSet9Ind, conSet12Ind)

                if foundNewColumn == 1:
                    relax = master.relax()
                    relax.optimize()
                    checkSolutionStaus(relax)
                    repeatFlag = 1
                else:
                    print("===No new columns!")
                    LB = relax.objVal
            pricingTime += (time.time() - time0)

            "Row generation procedure"
            foundNewRow = 1
            iterRow = 0
            time0 = time.time()
            while foundNewRow:
                iterRow += 1
                print("Repeat Full Iter: " + str(repeatIter) + " Row Generation Iter: " + str(iterRow))
                master.update()
                # master.write("MP_F"+str(repeatIter)+"_CG" + str(iter) + ".lp")
                relax = master.relax()
                relax.optimize()
                checkSolutionStaus(relax)
                # print('Relaxation Obj: %g' % relax.objVal)
                # foundNewRow = 0 #need to remove after debugging
                "row procedure. Row Generation."
                foundNewRow, master = CAVAD_Row_Generation(master, relax, sep)
                if foundNewRow == 1:
                    repeatFlag = 1
            sepTime += (time.time() - time0)
    else:
        print("Global UB: " + str(Global_UB))
        print("Global LB: " + str(Global_LB))
        repeatFlag = 0
        "Column generation procedure"
        relax = master.relax()
        relax.optimize()
        iterCol = 0
        foundNewColumn = 1
        time0 = time.time()
        while foundNewColumn:
            # spExecuteFlag = 0
            foundNewColumn = 0
            iterCol += 1
            print(" Column Generation Iter: " + str(iterCol))
            # master.write("MP_F"+str(repeatIter)+"_CG" + str(iter) + ".lp")
            "pricing procedure. Column Generation."
            if node.parentNode == -1:
                foundNewColumn, master, relax, yInd, P, Pdetail, QpSize, reverseP, Psize, numNewColumn = CAVAD_Pulse_Pricing(
                    master, relax, yInd, P, Pdetail, QpSize, reverseP, Psize, conSet8Ind, conSet9Ind, conSet12Ind)

            else:
                foundNewColumn, master, relax, yInd, P, Pdetail, QpSize, reverseP, Psize, numNewColumn = CAVAD_Pulse_Pricing(
                    master, relax, yInd, P, Pdetail, QpSize, reverseP, Psize, conSet8Ind, conSet9Ind, conSet12Ind)

            if foundNewColumn == 1:
                relax = master.relax()
                relax.optimize()
                checkSolutionStaus(relax)
                repeatFlag = 1
            else:
                print("===No new columns!")
                LB = relax.objVal
        pricingTime += (time.time() - time0)

        Global_LB = relax.objVal

    node.master = master
    node.relax = relax
    node.conSet1Ind = conSet1Ind
    node.conSet2Ind = conSet2Ind
    node.conSet3Ind = conSet3Ind
    node.conSet4Ind = conSet4Ind
    node.conSet5Ind = conSet5Ind
    node.conSet6Ind = conSet6Ind
    node.conSet7Ind = conSet7Ind
    node.conSet8Ind = conSet8Ind
    node.conSet9Ind = conSet9Ind
    node.conSet10Ind = conSet10Ind
    node.conSet11Ind = conSet11Ind
    node.conSet12Ind = conSet12Ind
    # node.hInd = hInd
    node.xInd = xInd
    node.yInd = yInd
    node.zInd = zInd
    node.tInd = tInd
    node.P = P
    node.Pdetail = Pdetail
    node.reverseP = reverseP
    node.Psize = Psize
    node.QpSize = QpSize
    node.parentNode = parentNode
    node.LB = LB

    if parentNode <= nodeUBFlag:
        if parentNode == -1:
            print("Start Solving Root Node Master")
            Global_LB = relax.objVal
        else:
            print("Start Solving for an UB")

        master.Params.OutputFlag = 0
        master.Params.TimeLimit = 3600
        if LOG:
            master.Params.OutputFlag = 1  # verbose mode

        if drivingAlwaysFasterFlag == 0:
            master.optimize()
        else:
            master.Params.PreCrush = 1
            master.optimize(subtourElim_Vehicle)

        TimeBE = time.time() - start_time
        BE_UB = master.objVal
        print('CAVAD BCP Node Optimal Cost of the Master: %g' % master.objVal)
        if master.objVal < Global_UB:

            Global_x_val = master.getAttr('x', x)
            Global_y_val = master.getAttr('x', master._y)
            Global_z_val = master.getAttr('x', z)
            Global_t_val = master.getAttr('x', t)
            Global_y_costCoeff = master.getAttr('obj', master._y)
            Global_P = P
            Global_Pdetail = Pdetail
            Global_UB = master.objVal
            Global_Gap = (Global_UB - Global_LB) / Global_UB
            Global_y_numberBE = len(yInd)
            print("Global LB: " + str(Global_LB))
            print('CAVAD Global Gap: ' + str(Global_Gap))
            print('The number of y variables: ' + str(len(yInd)))

    Threshold = Global_UB - Global_LB
    time0 = time.time()
    foundNewColumn, master, relax, yInd, P, Pdetail, QpSize, reverseP, Psize, numNewColumn = CAVAD_Pulse_Enumerating(master, relax, yInd, P, Pdetail, QpSize, reverseP, Psize, conSet8Ind, conSet9Ind, conSet12Ind, Threshold)
    enumeratingTime = (time.time() - time0)

    if drivingAlwaysFasterFlag == 0:
        master.optimize()
    else:
        master.Params.PreCrush = 1
        master.optimize(subtourElim_Vehicle)
    TimeAE = time.time() - time0
    print('After enumerating, CAVAD BCP Node Optimal Cost of the Master: %g' % master.objVal)
    Global_y_numberAE = len(yInd)
    print('After enumerating, The number of y variables: ' + str(Global_y_numberAE))

    if master.objVal < Global_UB:

        Global_x_val = master.getAttr('x', x)
        Global_y_val = master.getAttr('x', master._y)
        Global_z_val = master.getAttr('x', z)
        Global_t_val = master.getAttr('x', t)
        Global_y_costCoeff = master.getAttr('obj', master._y)
        Global_P = P
        Global_Pdetail = Pdetail
        Global_UB = master.objVal
        Global_Gap = (Global_UB - Global_LB) / Global_UB
        print("Global LB: " + str(Global_LB))
        print('CAVAD Global Gap: ' + str(Global_Gap))



    return  node

if __name__ == '__main__':
    version = "-BP"
    outputFlag = 1
    BPCFlag = 0
    userCutsFlag = 1
    lazyConsFlag = 0
    solutionCheckingFlag = 0
    disaggregateCon8Flag = 0
    terminateGap = 0.05
    epsilon = 0.01
    LOG = 0
    printData = 0
    # maximumNode = 640000000
    maximumNode = -1
    nodeUBFlag = -1 # -1 only root node. maximumNode all nodes
    # nodeUBFlag = maximumNode
    printRouteFlag = 0
    f = 2.8
    instanceList = []
    numFiles = 10
    # for name in ['Cook']:
    for name in ['Cook','Winnebago','Champaign','LaSalle','Adams','Fulton','Jefferson','Johnson','Cumberland']:

        for fileID in range(0, numFiles):
            instanceList.append((name, fileID+1))
    # for B in [3]:
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
            Global_r_val = dict()
            Global_t_val = dict()
            Global_y_costCoeff = dict()
            Global_P = dict()
            Global_Pdetail = dict()
            Global_LB = -np.inf
            Global_UB = np.inf
            master = gp.Model("CAVADP_master")  # master LP problem

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
            y = gp.tupledict() #initialize a set of y variables
            P = dict() #paths
            Pdetail = dict()
            reverseP = dict() #contains the set of loops that contains node i in Q
            for i in range(n):
                reverseP[i]=set()
            Psize = np.zeros((n,n), dtype=int)
            QpSize = dict()
            for i in range(1, n):  # (1,2,1), (2,3,2),......(n-2,n-1, n-2), the last one (n-1, 1, n-1)
                Q = set() #the set of nodes serviced in this path
                if i == n - 1:
                    j = 1
                else:
                    j = i + 1
                Q.add(j)

                P[i,i, Psize[i,i]] = Q
                Qdetail = [(i, j), (j, i)]
                Pdetail[i,i,Psize[i,i]]=Qdetail
                cpValue = f  + walkingTime[i,j]+walkingTime[j,i]
                y[i, i, Psize[i,i]] = master.addVar(obj=cpValue, vtype=GRB.BINARY, name='y[%d,%d,%d]' % (i, i, Psize[i, i]))
                QpSize[i, i, Psize[i, i]] = len(Q)
                for j in Q:
                    reverseP[j].add((i,i, Psize[i, i]))
                Psize[i,i] += 1
            for i in range(1,n): #(1), (2), (3), ......(n-2) , (n-1)
                Q = set()
                Q.add(i)
                P[i, i, Psize[i, i]] = Q
                Qdetail = []
                Pdetail[i, i, Psize[i, i]] = Qdetail
                cpValue = f
                y[i, i, Psize[i, i]] = master.addVar(obj=cpValue, vtype=GRB.BINARY, name='y[%d,%d,%d]' % (i, i, Psize[i, i]))
                QpSize[i, i, Psize[i, i]] = len(Q)
                for j in Q:
                    reverseP[j].add((i, i, Psize[i, i]))
                Psize[i, i] += 1
                if B >= 2: # (1,2), (2,3),......(n-2,n-1), the last one (n-1, 1)
                    for j in rwG.successors(i):
                        Q = set()
                        Q.add(i)
                        Q.add(j)
                        P[i, j, Psize[i, j]] = Q
                        Qdetail = [(i, j)]
                        Pdetail[i, j, Psize[i, j]] = Qdetail
                        cpValue = f + max(walkingTime[i, j], drivingTime[i,j])
                        y[i, j, Psize[i, j]] = master.addVar(obj=cpValue, vtype=GRB.BINARY, name='y[%d,%d,%d]' % (i, j, Psize[i, j]))
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
            conSet5 = master.addConstrs(
                x.sum('*', i) + z.sum('*', i) == x.sum(i, '*') + z.sum(i, '*') for i in range(1, n))
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
            conSet12 = master.addConstrs(v.sum('*', i) - v.sum(i, '*') - gp.quicksum(
                y[i, j, p] * QpSize[i, j, p] for j in range(1, n) for p in range(Psize[i, j])) == 0 for i in
                                         range(1, n))
            conSet12Ind = dict((i, conSet12[i].index) for i in conSet12.keys())

            if drivingAlwaysFasterFlag == 0:
                "Con8"
                if disaggregateCon8Flag == 1:
                    conSet8 = master.addConstrs(y[i, i, p] <=(x.sum('*', i) + z.sum('*', i)) for i in range(1, n) for p in range(Psize[i, i]) )
                else:
                    conSet8 = master.addConstrs(gp.quicksum(y[i, i, p] for p in range(Psize[i, i]))<=cn*(x.sum('*', i) + z.sum('*', i)) for i in range(1, n)  )
                conSet8Ind = dict((i, conSet8[i].index) for i in conSet8.keys())
                conSet3 = dict()
                conSet3Ind = dict()
                conSet4 = dict()
                conSet4Ind = dict()
                conSet7 = dict()
                conSet7Ind = dict()
                conSet10 = dict()
                conSet10Ind = dict()
                conSet11 = dict()
                conSet11Ind = dict()
            else:
                conSet4 = dict()
                conSet4Ind = dict()
                conSet7 = dict()
                conSet7Ind = dict()
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
            # if lazyConsFlag ==1:
            #     master.Params.lazyConstraints = 1
            # if userCutsFlag == 1:
            #     master.Params.PreCrush = 1
            relax = master.relax()
            relax.Params.OutputFlag = 0
            x_sol = gp.tupledict()
            y_sol = gp.tupledict()
            z_sol = gp.tupledict()

            "branch-price-cut full routine"
            tree = dict()
            nodeCount = 0
            tree[nodeCount] = Node(master, relax, conSet1Ind, conSet2Ind, conSet3Ind, conSet4Ind, conSet5Ind, conSet6Ind, conSet7Ind, conSet8Ind, conSet9Ind, conSet10Ind, conSet11Ind, conSet12Ind, xInd, yInd,
                                   zInd, tInd, P, Pdetail, reverseP, Psize, QpSize, -1)
            nextLevel = [nodeCount]
            levelNum = 0
            while nextLevel != []:
                print("This is level: " + str(levelNum))
                print("Global UB: " + str(Global_UB))
                print("Global LB: " + str(Global_LB))
                Global_Gap = (Global_UB - Global_LB) / Global_UB
                print('Global Gap: ' + str(Global_Gap))

                currentLevel = nextLevel.copy()
                nextLevel = []
                for nodeID in currentLevel:
                    print("This is node: " + str(nodeID))
                    CP_RoutineNode(tree[nodeID])
                "update global LB"
                currentLevelLB = np.infty
                for nodeID in currentLevel:
                    if tree[nodeID].LB <= currentLevelLB:
                        currentLevelLB = tree[nodeID].LB
                Global_LB = currentLevelLB
                "branching"
                if nodeCount < maximumNode:
                    for nodeID in currentLevel:
                        if (Global_UB - tree[nodeID].LB) / Global_UB > terminateGap:
                            conBranchingID = tree[nodeID].branching()
                            rhs = 1
                            if conBranchingID != 0:
                                for i in range(2):
                                    nextMaster = tree[nodeID].master.copy()
                                    hvar = nextMaster.getVarByName("t[" + str(conBranchingID) + "]")
                                    nextMaster.addConstr(hvar  == rhs)
                                    nextMaster.update()
                                    nextRelax = nextMaster.relax()
                                    nodeCount += 1
                                    tree[nodeCount] = Node(nextMaster, nextRelax, tree[nodeID].conSet1Ind, tree[nodeID].conSet2Ind,
                                                           tree[nodeID].conSet3Ind, tree[nodeID].conSet4Ind, tree[nodeID].conSet5Ind,
                                                           tree[nodeID].conSet6Ind, tree[nodeID].conSet7Ind,
                                                           tree[nodeID].conSet8Ind, tree[nodeID].conSet9Ind, tree[nodeID].conSet10Ind,
                                                           tree[nodeID].conSet11Ind, tree[nodeID].conSet12Ind, tree[nodeID].xInd, tree[nodeID].yInd, tree[nodeID].zInd,
                                                           tree[nodeID].tInd,
                                                           tree[nodeID].P, tree[nodeID].Pdetail, deepcopy(tree[nodeID].reverseP), tree[nodeID].Psize,
                                                           tree[nodeID].QpSize,  nodeID)
                                    nextLevel.append(nodeCount)
                                    rhs = 0
                    levelNum += 1
                else:
                    nextLevel = []

            print("=========== Finishing BCP =================")
            print('CAVAD BCP Node Optimal Cost of the Master Before Enumerating: %g' % BE_UB)
            print('After enumerating, CAVAD BCP Node Optimal Cost of the Master: %g' % Global_UB)
            print("Global UB: " + str(Global_UB))
            print("Global LB: " + str(Global_LB))
            Global_Gap = (Global_UB - Global_LB) / Global_UB
            print('CAVAD BCP Global Gap: ' + str(Global_Gap))
            print("---BCP Running Time Before Enumerating %s seconds ---" % (TimeBE))
            print("---BCP Running Time After Enumerating %s seconds ---" % (TimeAE))
            print(f"---BCP Pricing Time {pricingTime} seconds. It is {pricingTime/(time.time() - start_time)*100}%. ---")
            print(f"---BCP Enumerating Time {enumeratingTime} seconds. It is {enumeratingTime/(time.time() - start_time)*100}%. ---")

            if outputFlag ==1:
                file = open("Plus_H_V"+str(version)+"_CAVAD_County_" + str(name) + "_N" + str(n) + "_TGap" + str(terminateGap) + "_UB_BeforeEnumerating.log", 'a')
                file.write("\n" + str(BE_UB))
                file.close()
                file = open("Plus_H_V"+str(version)+"_CAVAD_County_" + str(name) + "_N" + str(n)  + "_TGap" + str(terminateGap) + "_UB.log", 'a')
                file.write("\n" + str(Global_UB))
                file.close()
                file = open("Plus_H_V"+str(version)+"_CAVAD_County_" + str(name) + "_N" + str(n)  + "_TGap" + str(terminateGap) + "_LB.log", 'a')
                file.write("\n" + str(Global_LB))
                file.close()
                file = open("Plus_H_V"+str(version)+"_CAVAD_County_" + str(name) + "_N" + str(n)  + "_TGap" + str(terminateGap) + "_GlobalGap.log", 'a')
                file.write("\n" + str(Global_Gap))
                file.close()
                file = open("Plus_H_V"+str(version)+"_CAVAD_County_" + str(name) + "_N" + str(n)  + "_TGap" + str(terminateGap) + "_RunningTimeBeforeEnumerating.log", 'a')
                file.write("\n" + str(TimeBE))
                file.close()
                file = open("Plus_H_V"+str(version)+"_CAVAD_County_" + str(name) + "_N" + str(n) + "_TGap" + str(
                    terminateGap) + "_RunningTimeAfterEnumerating.log", 'a')
                file.write("\n" + str(TimeAE))
                file.close()
                file = open("Plus_H_V"+str(version)+"_CAVAD_County_" + str(name) + "_N" + str(n) + "_TGap" + str(
                    terminateGap) + "_TotalRunningTime.log", 'a')
                file.write("\n" + str(TimeAE+TimeBE))
                file.close()
                file = open("Plus_H_V"+str(version)+"_CAVAD_County_" + str(name) + "_N" + str(n)  + "_TGap" + str(terminateGap) + "_NodesCount.log", 'a')
                file.write("\n" + str(nodeCount))
                file.close()
                file = open(
                    "Plus_H_V"+str(version)+"_CAVAD_County_" + str(name) + "_N" + str(n) + "_TGap" + str(
                        terminateGap) + "_NumberOfyVariablesBeforeEnumerating.log", 'a')
                file.write("\n" + str(Global_y_numberBE))
                file.close()
                file = open(
                    "Plus_H_V"+str(version)+"_CAVAD_County_" + str(name) + "_N" + str(n) + "_TGap" + str(terminateGap) + "_NumberOfyVariables.log", 'a')
                file.write("\n" + str(Global_y_numberAE))
                file.close()
                # file = open(
                #     "Plus_H_V"+str(version)+"_CAVAD_County_" + str(name) + "_N" + str(n) + "_TGap" + str(terminateGap) + "_SeparationTime.log", 'a')
                # file.write("\n" + str(sepTime))
                # file.close()
                file.close()
                file = open("Plus_H_V"+str(version)+"_CAVAD_County_" + str(name) + "_N" + str(n) + "_TGap" + str(terminateGap) + "_PricingTime.log",
                            'a')
                file.write("\n" + str(pricingTime))
                file.close()

                file = open(
                    "Plus_H_V"+str(version)+"_CAVAD_County_" + str(name) + "_N" + str(n) + "_TGap" + str(terminateGap) + "_EnumeratingTime.log",
                    'a')
                file.write("\n" + str(enumeratingTime))
                file.close()


