

"""
A modular of solvers and algorithms
"""

import numpy as np
import networkx as nx
import os
import json
import gurobipy as gp
from gurobipy import GRB
import utils
from copy import deepcopy
from ortools.constraint_solver import pywrapcp
from ortools.constraint_solver import routing_enums_pb2
import sys
from termcolor import colored




"""
Grid-like graph generator
"""
def grid_like_graph_generator(x_range, y_range, n, I):
    
    grid = utils.av_grid_pts_generator(x_range, y_range, n)
    
    zone_dict = utils.node_zoning(I.coordinates, grid)
    
    G = utils.add_customer_n_depot_to_grid(grid, I.coordinates, zone_dict)
    
    
    #TODO: The logic is not correct 
    # for i, coor_i in enumerate(I.coordinates):
    
    #     for j, coor_j in enumerate(I.coordinates):

    #         node_i = tuple(coor_i)
    #         node_j = tuple(coor_j)
            
    #         no_shortest_path = (np.linalg.norm(np.array(coor_i) - np.array(coor_j), ord=1) < nx.shortest_path_length(G, node_i, node_j, weight = 'length'))
    #         if j > i and no_shortest_path:
    #             new_node = (node_i[0], node_j[1], 'a')
    #             G.add_node(new_node)
    #             # G.add_edge(node_i, new_node, length = abs(node_i[1] - new_node[1]), is_av = 0)
    #             # G.add_edge(new_node, node_i, length = abs(node_i[1] - new_node[1]), is_av = 0)
    #             # G.add_edge(node_j, new_node, length = abs(node_j[0] - new_node[0]), is_av = 0)
    #             # G.add_edge(new_node, node_j, length = abs(node_j[0] - new_node[0]), is_av = 0)
    
    
    G = utils.mark_customer_n_depot(G, I.coordinates)
    
    
    
    return G

"""
Grid-like graph generator (dense)
"""
def grid_like_graph_generator_dense(x_range, y_range, n, I):
    
    
    grid = utils.av_grid_pts_generator(x_range, y_range, n)
    grid_pts = [[node[0], node[1]] for node in grid.nodes]
    
    
    G = utils.generate_grid_by_points(I.coordinates + grid_pts)
    
    
    # Mark all edges aligned with grid as av-enabled roads
    x_step = grid.graph['x_step']
    y_step = grid.graph['y_step']
    
    vertical = np.arange(x_range[0], x_range[0] + (n + 1) * x_step, x_step)
    horizontal = np.arange(y_range[0], y_range[0] + (n + 1) * y_step, y_step)
    
    for e in G.edges:
        if e[0][0] == e[1][0] and e[0][0] in vertical:
            G.edges[e]['is_av'] = 1
        if e[0][1] == e[1][1] and e[0][1] in horizontal:
            G.edges[e]['is_av'] = 1
    
    
    
    # Mark all customers and depot
    G = utils.mark_customer_n_depot(G, I.coordinates)
    
    return G






"""
An efficient generic VRP solver (by Gurobi)
"""

def generic_vrp_efficient_solver(G, customer, demand, depot, capacity, vehicle_cost, dir):
    # Compute all shortest paths among clients
    K_D = depot + customer
    cost_non_av = np.zeros((len(customer) + 1, len(customer) + 1))
    cost_av = np.zeros((len(customer) + 1, len(customer) + 1))
    path_non_av = {}
    path_av = {}

    for i in range(len(K_D)):
        for j in range(len(K_D)):
            if i!=j:
                cost_non_av[i][j] = nx.shortest_path_length(G, K_D[i], K_D[j], weight='non_av_cost')
                cost_av[i][j] = nx.shortest_path_length(G, K_D[i], K_D[j], weight='av_cost')

                path_non_av[str(K_D[i])+'-'+str(K_D[j])] = nx.shortest_path(G, K_D[i], K_D[j], weight='non_av_cost')
                path_av[str(K_D[i])+'-'+str(K_D[j])] = nx.shortest_path(G, K_D[i], K_D[j], weight='av_cost')
         
    # Create a new model
    model = gp.Model('generic_vrp')

    # Set the file path for logging
    path = dir + '/generic_vrp_log.txt'
    if os.path.exists(path):
        os.remove(path)
    model.Params.LogFile = path

    # Parameters tuning for the MIP solver
    model.Params.MIPFocus = 2 # Focusing on the phase 3 (proving optimality)
    # model.Params.MIPGap = 0 # Specify the terminating duality gap
    
    # Create variables
    x_a = model.addMVar((len(K_D), len(K_D)), vtype=GRB.BINARY)
    x_h = model.addMVar((len(K_D), len(K_D)), vtype=GRB.BINARY)
    h = model.addMVar(len(customer), vtype=GRB.INTEGER, ub=np.repeat(capacity, len(customer)), lb=demand)
    
    
    # Constraints
    model.addConstrs(sum(x_h[:, i]) == sum(x_h[i, :]) for i in range(len(K_D)))
    model.addConstrs(sum(x_a[:, i]) == sum(x_a[i, :]) for i in range(len(K_D)))
    model.addConstrs(sum(x_h[i, :]) + sum(x_a[i, :]) == 1 for i in range(1, len(K_D)))    
    
    for i in range(len(K_D)):
        model.addConstr(x_h[i, i] == 0)
        model.addConstr(x_a[i, i] == 0)
            
    for i in range(len(customer)):
        for j in range(len(customer)):
            if i != j:
                model.addConstr(h[j] >= h[i] + capacity * (x_h[i+1, j+1] + x_a[i+1, j+1] - 1) + demand[j])


    # Set the bbjective
    model.setObjective(sum(x_a[i, j] * cost_av[i][j] for i in range(len(K_D)) for j in range(len(K_D))) +
                        sum(x_h[i, j] * cost_non_av[i][j] for i in range(len(K_D)) for j in range(len(K_D))) +
                        x_a[0, :] @ np.repeat(vehicle_cost, len(K_D)) + 
                        x_h[0, :] @ np.repeat(vehicle_cost, len(K_D)),
                        GRB.MINIMIZE)
    
    # Solve the VRP
    model.optimize()
    
    
    # Write the gurobi solution to a file
    # model.write(dir+'/solutiom.sol')
    
   
    np.savetxt(dir + '/x_h.txt', x_h.X)
    np.savetxt(dir + '/x_a.txt', x_a.X)
    
    
    # Extract the routes
    # routes = []
    # for m in fleet:
    #     if y[0, m].X == 0:
    #         routes.append([])
    #         continue
    #     route_m = []
    #     start = 0
    #     next = np.where(x[0, :, m].X > .5)[0][0]
    #     while next != 0:
    #         if m < num_av:
    #             route_m = route_m + path_av[str(K_D[start])+'-'+str(K_D[next])][1:]
    #         else:
    #             route_m = route_m + path_non_av[str(K_D[start])+'-'+str(K_D[next])][1:]
    #         start = next
    #         next = np.where(x[next, :, m].X > .5)[0][0]
    #     if m < num_av:
    #         route_m = route_m + path_av[str(K_D[start])+'-'+str(K_D[0])][1:]
    #     else:
    #         route_m = route_m + path_non_av[str(K_D[start])+'-'+str(K_D[0])][1:]
    #     route_m.insert(0, K_D[0])
    #     routes.append(route_m)
    
       



"""
A generic VRP solver (by Gurobi)
"""


def generic_vrp_solver(G, customer, demand, depot, fleet, num_av, capacity, vehicle_cost, dir, t_max, gran, vehicle_range=None):
    
    
    # Compute all shortest paths among clients
    K_D = depot + customer
    cost_non_av = np.zeros((len(customer) + 1, len(customer) + 1))
    cost_av = np.zeros((len(customer) + 1, len(customer) + 1))
    path_non_av = {}
    path_av = {}

    for i in range(len(K_D)):
        for j in range(len(K_D)):
            if i!=j:
                cost_non_av[i][j] = nx.shortest_path_length(G, K_D[i], K_D[j], weight='non_av_cost')
                cost_av[i][j] = nx.shortest_path_length(G, K_D[i], K_D[j], weight='av_cost')

                path_non_av[str(K_D[i])+'-'+str(K_D[j])] = nx.shortest_path(G, K_D[i], K_D[j], weight='non_av_cost')
                path_av[str(K_D[i])+'-'+str(K_D[j])] = nx.shortest_path(G, K_D[i], K_D[j], weight='av_cost')
            


    # Create a new model
    model = gp.Model('generic_vrp')

    # Set the file path for logging
    path = dir + '/generic_vrp_log.txt'
    if os.path.exists(path):
        os.remove(path)
    model.Params.LogFile = path

    # Parameters tuning for the MIP solver
    model.Params.MIPFocus = 3 # Focusing on the phase 2 (moving bounds)
    # model.Params.MIPGap = 0 # Specify the terminating duality gap

    # Create variables
    x = model.addMVar((len(K_D), len(K_D), len(fleet)), vtype=GRB.BINARY)
    y = model.addMVar((len(K_D), len(fleet)), vtype=GRB.BINARY)
    h = model.addMVar((len(customer), len(fleet)), vtype=GRB.CONTINUOUS, ub=sum(demand))

    # Specify an initial feasible solution
    # x.setAttr("Start", x_init)
    # y.setAttr("Start", y_init)
    # h.setAttr("Start", h_init)

    # Constraints
    for i in range(len(K_D)):
        model.addConstrs(sum(x[:, i, m]) == sum(x[i, :, m]) for m in fleet)
        model.addConstrs(x[i, i, m] == 0 for m in fleet)
        model.addConstrs(sum(x[:, i, m]) == y[i, m] for m in fleet)

        #! We should not add the same constraint twice
        model.addConstrs(sum(y[i, :]) == 1 for i in range(1, len(K_D)))

        model.addConstr(sum(y[0, :]) <= len(fleet))

        if vehicle_range != None:
            model.addConstrs(sum(x[i, j, m] * cost_non_av[i][j] for i in range(len(K_D)) for j in range(len(K_D))) <= range[m] 
                        for m in fleet)

        model.addConstrs(sum(demand[i] * y[i+1, m] for i in range(len(customer))) <= capacity
                    for m in fleet)

    for i in range(len(customer)):
        for j in range(len(customer)):
            model.addConstrs(h[j, m] >= h[i, m] + capacity * (x[i+1, j+1, m] - 1) + demand[i] for m in fleet)
    model.addConstrs(h[i, m] >= demand[i] for m in fleet for i in range(len(customer)))

    # Set the bbjective
    model.setObjective(sum(x[i, j, m] * cost_av[i][j] for i in range(len(K_D)) for j in range(len(K_D)) for m in fleet[: num_av]) +
                        sum(x[i, j, m] * cost_non_av[i][j] for i in range(len(K_D)) for j in range(len(K_D)) for m in fleet[num_av :]) +
                        y[0, :] @ np.repeat(vehicle_cost, len(fleet)), GRB.MINIMIZE)

    # Solve the VRP
    model.optimize()
    
    # Extract the routes
    routes = []
    for m in fleet:
        if y[0, m].X == 0:
            routes.append([])
            continue
        route_m = []
        start = 0
        next = np.where(x[0, :, m].X > .5)[0][0]
        while next != 0:
            if m < num_av:
                route_m = route_m + path_av[str(K_D[start])+'-'+str(K_D[next])][1:]
            else:
                route_m = route_m + path_non_av[str(K_D[start])+'-'+str(K_D[next])][1:]
            start = next
            next = np.where(x[next, :, m].X > .5)[0][0]
        if m < num_av:
            route_m = route_m + path_av[str(K_D[start])+'-'+str(K_D[0])][1:]
        else:
            route_m = route_m + path_non_av[str(K_D[start])+'-'+str(K_D[0])][1:]
        route_m.insert(0, K_D[0])
        routes.append(route_m)
        

    # Write the routes to a json file
    dict_routes = dict(zip(range(len(routes)), routes))
    with open(dir + '/routes.json', 'w') as outfile:
        json.dump(dict_routes, outfile)
    
    
    # Write the time stamps to a json file
    time_stamps = utils.time_tracker_all_routes(G, routes)
    dict_time_stamps = dict(zip(range(len(time_stamps)), time_stamps))
    with open(dir + '/all_time_stamps.json', 'w') as outfile:
        json.dump(dict_time_stamps, outfile)
        
        
    # Write the control vector to a json file
    control = utils.control_tracker(G, routes, fleet, num_av, t_max, gran)
    num_control = {'control_vector': list(np.sum(control, axis=0))}
    with open(dir + '/control.json', 'w') as outfile:
        json.dump(num_control, outfile)
      
      
    # Output which vehicle serve which customers
    clusters = []
    for m in range(y.shape[1]):
        if sum(y[:, m].X) == 0:
            continue
        # customers = np.where(y[:, m].X == 1)[0].tolist()
        
        
        K_D_subset = []
        for ind in np.where(y[:, m].X == 1)[0].tolist():
            K_D_subset.append(K_D[ind])

        clusters.append(K_D_subset)

    dict_clusters = dict(zip(range(len(clusters)), clusters))
    with open(dir + '/clusters.json', 'w') as outfile:
        json.dump(dict_clusters, outfile)
            
    
    # Visualize routes  
    # utils.visualize_separate_routes(G, routes, num_av, dir)
    
    


"""
Re-scheduling FRP
"""

#TODO: The MILP is different from the paper: travel time is not flexible in ordinary roads. May need to change it later. 

def frp_reschedule(G, routes, all_time_stamps, t_max, budget, dir, gamma, obj="makespan"):
    
    num_job = utils.convert_route_to_job(G, routes)
    all_job_start, all_job_end = utils.extract_start_n_end_time(G, all_time_stamps)
    
    # Write donw the old schedule
    dict_all_job_start = dict(zip(range(len(all_job_start)), all_job_start))
    dict_all_job_end = dict(zip(range(len(all_job_end)), all_job_end))
    with open(dir + '/old_schedule_start.json', 'w') as outfile:
        json.dump(dict_all_job_start, outfile)
    with open(dir + '/old_schedule_end.json', 'w') as outfile:
        json.dump(dict_all_job_end, outfile)
    
    
    
    # max_job_end = [np.max(sublist) for sublist in all_job_end]
    all_job_start = np.array([job for sublist in all_job_start for job in sublist])
    all_job_end = np.array([job for sublist in all_job_end for job in sublist])
    duration = all_job_end - all_job_start
    
    
    # A helper function determineing whether job i and job j are in the same route
    def in_the_same_route(i, j, num_job):
    
        b0 = False
        
        cum_num = np.cumsum(num_job)
        cum_num = np.insert(cum_num, 0, 0)
        for ind in range(1, len(cum_num)):
            b1 = cum_num[ind-1] <= i and i < cum_num[ind]
            b2 = cum_num[ind-1] <= j and j < cum_num[ind]
            if b1 and b2:
                b0 = True
                
        return b0
    
    # Create a new model
    model = gp.Model('rescheduling')

    # Set the file path for logging and disable the console output
    path = dir + '/re_scheduling_frp_log.log'
    if os.path.exists(path):
        os.remove(path)
    model.Params.LogFile = path
    model.Params.LogToConsole = 0
    
    # strategies
    # model.Params.MIPFocus = 1 # Focusing on phase 1 (getting the first feasible solution)
    model.Params.MIPFocus = 3 # Focusing on phase 2 (moving the best upper bound)
    model.Params.ImproveStartTime = 1 # Give up proving optimality after 2 min
    
    # Termination    
    model.Params.TimeLimit = 240 # Terminate the MILP solver after 5 min
    
    

    # Define variables
    z = model.addMVar((sum(num_job), budget), vtype=GRB.BINARY)
    c = model.addMVar(sum(num_job), vtype=GRB.CONTINUOUS, lb = all_job_start + duration, ub = np.repeat(t_max, sum(num_job)))
    f = model.addMVar((sum(num_job), sum(num_job)), vtype=GRB.BINARY)
    
    if obj == "makespan":
        c_max = model.addVar(lb = max(all_job_end), ub = t_max, vtype=GRB.CONTINUOUS)


    # Constraints

    # Unique assignment of job
    model.addConstrs(sum(z[i, :]) == 1 for i in range(sum(num_job)))

    # Maximum delay for each pair of job
    model.addConstrs(c[j] - duration[j] - c[i] <= gamma[1] * (all_job_start[j] - all_job_start[i] - duration[i]) 
                    for i in range(sum(num_job)) for j in range(sum(num_job)) if i == j-1 and j not in np.cumsum(num_job))

    # The precedence constraint
    model.addConstrs(c[j] - duration[j] - c[i] >= (all_job_start[j] - all_job_start[i] - duration[i])
                    for i in range(sum(num_job)) for j in range(sum(num_job)) if i == j-1 and j not in np.cumsum(num_job))


    # Non-overlapping jobs in a single machine
    model.addConstrs(c[i] <= t_max * (3 - z[i, b] - z[j, b] - f[i, j]) + c[j] - duration[j]
                    for i in range(sum(num_job)) for j in range(sum(num_job)) for b in range(budget)
                    if not in_the_same_route(i, j, num_job))

    model.addConstrs(c[j] <= t_max * (2 - z[i, b] - z[j, b] + f[i, j]) + c[i] - duration[i]
                    for i in range(sum(num_job)) for j in range(sum(num_job)) for b in range(budget)
                    if not in_the_same_route(i, j, num_job))

    # f_ij = 1 iff the start time of job i is ahead of the start time of job j
    model.addConstrs(c[i] <= c[j] + t_max * (1 - f[i, j])
                    for i in range(sum(num_job)) for j in range(sum(num_job)) if not in_the_same_route(i, j, num_job))
    model.addConstrs(c[j] <= c[i] + t_max * f[i, j]
                    for i in range(sum(num_job)) for j in range(sum(num_job)) if not in_the_same_route(i, j, num_job))

    # c <= c_max for all c
    if obj == "makespan":
        model.addConstrs(c[i] <= c_max for i in range(sum(num_job)))


    # print(c[np.cumsum(num_job)-1])
    
    # Set the objective
    if obj == "makespan":
        model.setObjective(c_max, GRB.MINIMIZE)
    if obj == "average_completion":
        weight = np.repeat(1/len(num_job), len(num_job))
        model.setObjective(c[np.cumsum(num_job)-1] @ weight, GRB.MINIMIZE)

    # Solve the scheduling feasibility recovering problem
    model.optimize()    
    
    # Return 0 if the re-scheduling FRP is not feasible
    
    
    if model.Status == 3:
        print(colored("The re-scheduling MILP is infeasible.", 'yellow'))
        return 0
    
    if model.Status == 9:
        print(colored("No feasible solution found within the time limit.", 'yellow'))
        return -1
    
    
    print(colored("Optimal solution found for the re-scheduling MILP!", 'green'))
    
    
    #* [Comment out this block if output is not needed]
    # Write down the new schedule
    new_start = c.X - duration
    
    cum_sum_job = np.cumsum(num_job)    
    new_schedule_end = [list(c[: cum_sum_job[0]].X)] + [list(c[cum_sum_job[i]: cum_sum_job[i+1]].X) for i in range(len(cum_sum_job)-1)]
    new_schedule_start = [list(new_start[: cum_sum_job[0]])] + [list(new_start[cum_sum_job[i]: cum_sum_job[i+1]]) for i in range(len(cum_sum_job)-1)]
    
    dict_new_schedule_end = dict(zip(range(len(new_schedule_end)), new_schedule_end))
    dict_new_schedule_start = dict(zip(range(len(new_schedule_start)), new_schedule_start))
    
    with open(dir + '/new_schedule_end.json', 'w') as outfile:
        json.dump(dict_new_schedule_end, outfile)
    with open(dir + '/new_schedule_start.json', 'w') as outfile:
        json.dump(dict_new_schedule_start, outfile)
        
    return 1
        
        


"""
Graph pre-pruning algorithm
"""

def graph_pre_pruning(G, depot, customer):

    E = list(G.edges)
    V = list(G.nodes)
    
    for e in E:
        if G.edges[e]['is_av'] == 1:
            G.edges[e]['control_time'] = 0
        else:
            G.edges[e]['control_time'] = G.edges[e]['length']

    remaining_edge = []
    remaining_node = []
    for i in depot + customer:
        for j in depot + customer:
            if i == j:
                continue
            p1 = nx.shortest_path(G, i, j, weight = 'length')
            p2 = nx.shortest_path(G, i, j, weight = 'control_time')
            remaining_node += p1
            remaining_node += p2
            for k in range(len(p1)-1):
                remaining_edge.append((p1[k], p1[k+1]))
            for k in range(len(p2)-1):
                remaining_edge.append((p2[k], p2[k+1]))
    
    remaining_edge = set(remaining_edge)
    remaining_node = set(remaining_node)
    edge_to_delete = set(E) - remaining_edge

    # Remove a set of nodes
    #* Comment this line to visualize how many nodes are deleted 
    G = G.subgraph(remaining_node)

    # Make a copy to unfreeze the graph
    G_s = G.copy()
    # Remove a set of edges 
    G_s.remove_edges_from(edge_to_delete)
    # G_s = utils.generate_av_roads_primary(G_s)

    E = list(G_s.edges)
    V = list(G_s.nodes)
        
    
    for v in list(set(V) - set(customer) - set(depot)):

        succ = list(G_s.successors(v))
        pred = list(G_s.predecessors(v))

        # Merge all one way roads
        if len(succ) == 1 and len(pred) == 1:
            # print(v)
            # print(G_s[pred[0]][v])
            # print("============")

            # Road types
            rt1 = G_s[pred[0]][v]['is_av']
            rt2 = G_s[v][succ[0]]['is_av']

            # Road length
            rl1 = G_s[pred[0]][v]['length']
            rl2 = G_s[v][succ[0]]['length']

            if rt1 == 1 and rt2 == 1:
                G_s.add_edge(pred[0], succ[0], is_av = rt1, is_artificial = 0, length = rl1 + rl2)
                G_s.remove_edge(pred[0], v)
                G_s.remove_edge(v, succ[0])

                #* Comment this line to visualize how many nodes are deleted 
                G_s.remove_node(v)
            
            if rt1 == 0 and rt2 == 0:
                G_s.add_edge(pred[0], succ[0], is_av = rt1, is_artificial = 0, length = rl1 + rl2)
                G_s.remove_edge(pred[0], v)
                G_s.remove_edge(v, succ[0])
                
                #* Comment this line to visualize how many nodes are deleted 
                G_s.remove_node(v)


        elif len(pred) == 2 and set(pred) == set(succ):
            
            # Road types
            rt1 = G_s[pred[0]][v]['is_av']
            rt2 = G_s[v][pred[1]]['is_av']

            # Road length
            rl1 = G_s[pred[0]][v]['length']
            rl2 = G_s[v][pred[1]]['length']

            if rt1 == 1 and rt2 == 1:
                G_s.add_edge(pred[0], pred[1], is_av = rt1, is_artificial = 0, length = rl1 + rl2)
                G_s.add_edge(pred[1], pred[0], is_av = rt1, is_artificial = 0, length = rl1 + rl2)
                G_s.remove_edge(pred[0], v)
                G_s.remove_edge(v, pred[0])
                G_s.remove_edge(pred[1], v)
                G_s.remove_edge(v, pred[1])
                G_s.remove_node(v)

            if rt1 == 0 and rt2 == 0:
                G_s.add_edge(pred[0], pred[1], is_av = rt1, is_artificial = 0, length = rl1 + rl2)
                G_s.add_edge(pred[1], pred[0], is_av = rt1, is_artificial = 0, length = rl1 + rl2)
                G_s.remove_edge(pred[0], v)
                G_s.remove_edge(v, pred[0])
                G_s.remove_edge(pred[1], v)
                G_s.remove_edge(v, pred[1])
                G_s.remove_node(v)

        else:
            pass
    
    return G_s



"""
Re-routing FRP
"""

class ReRoutingFRPSolver():
    
    
    def __init__(self, G, I, routes, all_time_stamps, clusters, t_max, budget, gamma, num_layer, dir) -> None:
        
        # Members from input
        self.G = G
        for e in self.G.edges:
            self.G[e[0]][e[1]]['is_artificial'] = 0
        self.I = I
        
        self.routes = routes
        self.all_time_stamps = all_time_stamps
        self.clusters = clusters
        self.demand = I.demands
        self.capacity = I.capacity
        self.budget = budget
        self.t_max = t_max
        self.gamma = gamma
        self.num_layer = num_layer
        self.dir = dir
        
        
        # Members for output (dynamically updated)
        self.new_routes = {}
        self.new_all_time_stamps = {}
        self.G_adjusted = G # Need to be output as well

        

        
        
    def budget_violation_check(self):
        
        all_job_start, all_job_end = utils.extract_start_n_end_time(self.G_adjusted, self.new_all_time_stamps)
        new_all_time_stamps_list = list(set([element for innerList in all_job_start + all_job_end for element in innerList]))
        new_all_time_stamps_list.sort()
        new_all_time_stamps_list.append(self.t_max)
        
        Q = np.arange(len(new_all_time_stamps_list)-1)
        

        for i in range(len(Q)):
            ss = new_all_time_stamps_list[i]
            tt = new_all_time_stamps_list[i+1]
            
            control_q = 0
            for v in range(len(self.new_all_time_stamps)):
                bb = False
                sss = all_job_start[v]
                ttt = all_job_end[v]
                for j in range(len(sss)):
                    if ss >= sss[j] and tt <= ttt[j]:
                        bb = True
                if bb == True:
                    control_q += 1
            
            if control_q > self.budget:
                return True

        return False
                
    
    
    def constrained_tsp_solver(self, v_2_re_route, budget_indicator, new_all_time_stamps_list):
        
        
        Q = np.arange(len(budget_indicator))
        K = self.clusters[v_2_re_route][:1]
        D = self.clusters[v_2_re_route][1:]
        
        demand_subset = []
        for cus in D:
            cus_index = int(cus[-1][1:])
            demand_subset.append(self.demand[cus_index])
          
        G_s = graph_pre_pruning(self.G, K, D)
        
        # Add merged edges in the pruned graph to the original graph for time stamps extraction in next iteration
        self.G_adjusted = utils.graph_edge_adjustment(G_s, self.G_adjusted)
    
        V = list(G_s.nodes)
        _, new_to_old = utils.re_index_nodes(V)
        
        G_e, departure, arrival, D_e, D_e_2d, demand_dict = utils.contruct_time_expanded_graph(G_s, self.num_layer, K, D, demand_subset)
    
            
        #TODO: The discount_factor and inflated_factor for AV_enabled roads are needed here to conduct sensitivity analysis
        utils.assign_cost_n_travel_time(G_e)
        
        E_e = list(G_e.edges)    
        V_e = list(G_e.nodes)
        edge_to_no, _ = utils.re_index_edges(E_e)
        

        # Create a new model
        model = gp.Model('re-routing-vehicle-'+str(v_2_re_route))

        # Model params
        
        # logging
        path = self.dir + '/re_routing_log_vehicle_'+str(v_2_re_route)+'.log'
        if os.path.exists(path):
            os.remove(path)
        model.Params.LogFile = path
        model.Params.LogToConsole = 0

        # strategies
        # model.Params.MIPFocus = 1 # Focusing on phase 1 (getting the first feasible solution)
        model.Params.MIPFocus = 3 # Focusing on phase 2 (moving the best upper bound)
        model.Params.ImproveStartTime = 120 # Give up proving optimality after 2 min
        
        # Termination
        model.Params.BestObjStop = 1.05 * self.original_route_cost(v_2_re_route)  # The original route cost is a natural lower bound. 
                                                                             # A re-route is good if the cost is < 1.05 original route cost
        
        model.Params.TimeLimit = 240 # Terminate the MILP solver after 5 min
        
        #TODO: Replacing AV cost with HDV cost naively is incorrect. We need to re-solve the TSP
        model.Params.Cutoff = self.HDV_route_cost(v_2_re_route) # cutoff the solution with even worse cost than dispatching an HDV
        

        # Create variables
        #* The lower bounds for all variables are zero by default
        x = model.addMVar(len(E_e), vtype=GRB.BINARY)
        y = model.addMVar(len(V_e), vtype=GRB.BINARY)
        t = model.addMVar(len(V_e), ub = np.repeat(self.t_max, len(V_e)), vtype=GRB.CONTINUOUS)

        alpha = model.addMVar((len(E_e), len(Q)), vtype=GRB.BINARY)
        beta = model.addMVar((len(E_e), len(Q)), vtype=GRB.BINARY)

        #// makespan = model.addVar(ub = t_max, vtype=GRB.CONTINUOUS)

        # Add constraints

        ## Flow conservation at all nodes (except the depot)
        for i in list(set(V_e) - set([departure, arrival])):
            out_no = [edge_to_no[e] for e in list(G_e.out_edges(i))]
            in_no = [edge_to_no[e] for e in list(G_e.in_edges(i))]
            model.addConstr(sum(x[out_no]) == sum(x[in_no]))
            model.addConstr(sum(x[out_no]) <= 1)
        
        ## Flow conservation at the depot
        departure_out = [edge_to_no[e] for e in list(G_e.out_edges(departure))]
        arrival_in = [edge_to_no[e] for e in list(G_e.in_edges(arrival))]

        model.addConstr(sum(x[departure_out]) == y[departure])
        model.addConstr(sum(x[arrival_in]) == y[arrival])
        model.addConstr(y[departure] == y[arrival])

        ## Prohibit vehicles from passing through the depot
        departure_in = [edge_to_no[e] for e in list(G_e.in_edges(departure))]
        model.addConstr(x[departure_in] == 0)
        
        ## Allow a vehicle to pass thru a customer without serving it 
        for i in D_e:
            i_in = [edge_to_no[e] for e in list(G_e.in_edges(i))]
            model.addConstr(sum(x[i_in]) >= y[i])

        ## Each customer is only served by one vehicle
        model.addConstrs(sum(y[i] for i in D_e_2d[:, h]) == 1 for h in range(len(D)))

        ## Transition
        for i in D_e:
            model.addConstr(sum(x[edge_to_no[(i, j)]] for j in D_e if G_e.has_edge(i, j) and G_e[i][j]['is_artificial'] == 1) <= y[i])

        ## Capacity limit
        model.addConstr(sum(demand_dict[i] * y[i] for i in D_e) <= self.capacity)


        ## Subtour elimination + time tracking
        for i in V_e:
            for j in V_e:
                if G_e.has_edge(i, j):
                    model.addConstr(t[j] >= t[i] + self.gamma[0] * G_e[i][j]['travel_time'] + self.t_max * (x[edge_to_no[(i, j)]] - 1))
                    model.addConstr(t[j] <= t[i] + self.gamma[1] * G_e[i][j]['travel_time'] + self.t_max * (1 - x[edge_to_no[(i, j)]]))
                    
        # The following constraints can be captured by the previous constraint
        #// for e in E_e:
        #//     i, j = e[:2]
        #//     # We need one more set of important constraints, which equals the time stamp of the start node to that of the end node in any aritificial edges.
        #//     if G_e[i][j]['is_artificial'] == 1:
        #//         model.addConstr(t[j] <= t[i] + t_max * (1 - x[edge_to_no[(i, j)]]))
            
        #//     if G_e[i][j]['is_av'] == 0:
        #//         model.addConstr(t[j] <= t[i] + G_e[i][j]['travel_time'] + t_max * (1 - x[edge_to_no[(i, j)]]))
        
        
        ## Do not violate the budget
        #* This is simplied constraint compared to the original MILP
        for e in E_e:
            i, j = e[:2]
            if G_e[i][j]['is_av'] == 0:
                model.addConstrs(x[edge_to_no[(i, j)]] == alpha[edge_to_no[(i, j)], q] + beta[edge_to_no[(i, j)], q]
                                for q in Q if budget_indicator[q] == 0)
                model.addConstrs(t[j] <= new_all_time_stamps_list[q] + self.t_max * (1 - alpha[edge_to_no[(i, j)], q])
                                for q in Q if budget_indicator[q] == 0)
                model.addConstrs(t[i] >= new_all_time_stamps_list[q+1] + self.t_max * (beta[edge_to_no[(i, j)], q] - 1)
                            for q in Q if budget_indicator[q] == 0)

        # Set the makespan
        #// model.addConstrs(t[i] <= makespan for i in range(len(V_e)))

        
        # Add the objective
        model.setObjective(sum(x[edge_to_no[(i, j)]] * G_e[i][j]['av_cost'] for i in V_e for j in V_e if G_e.has_edge(i, j)),
                        GRB.MINIMIZE)

        # Solve the TSP
        model.optimize()
        
        
        # Infeasible (status code = 3), or objective is worse than the cutoff value (status code = 6), or no feasible solution found within the time limit (status code = 9)
        if model.Status == 3:
            print(colored("Re-routing is infeasible. Unserved customers will be served by HDVs.", 'yellow'))
            return False, 0, _, _
        
        if model.Status == 6:
            print(colored("Re-routing is more expensive than disptaching HDVs. Unserved customers will be served by HDVs.", 'yellow'))
            return False, 0, _, _
        
        if model.Status == 9 and model.SolCount == 0:
            print(colored("Re-routing fails within the computational time limit. Unserved customers will be served by HDVs.", 'yellow'))
            return False, 0, _, _
        
        
        
        

        # Extract the routes and time stamps
        new_route, new_time_stamps = utils.extract_single_route_n_times(G_e, len(V), x, t, new_to_old, edge_to_no, departure, arrival)
        
    
        #* Returned values are: re_routing success, new routing cost, new route, new timestamps 
        return True, model.ObjVal, new_route, new_time_stamps
    
    
    def original_route_cost(self, vehicle):
        cost = 0
        route = self.routes[vehicle]
        for i, j in zip(route[:-1], route[1:]):
                cost += self.G[i][j]['av_cost']
        
        return cost
    
    
    #TODO: This is not the correct cutoff value. We should re-solve the TSP given the customers originally served by the vehicle
    def HDV_route_cost(self, vehicle):
        cost = 0
        route = self.routes[vehicle]
        for i, j in zip(route[:-1], route[1:]):
                cost += self.G[i][j]['non_av_cost']
        
        return cost

            
    def sort_routes(self, priority = 'min_av_time'):
        
        if priority == 'min_av_time':
            discounted_factor = 0.5
            inflated_factor = 1.2
            
            all_job_start, all_job_end = utils.extract_start_n_end_time(self.G, self.all_time_stamps)
            
            # The amount of time spent on ordinary roads
            control_time_list = []
            for i, j in zip(all_job_start, all_job_end):
                control_time = sum(np.array(j) - np.array(i))
                control_time_list.append(control_time)
                
            # The amount of increased cost if replaced by an HDV
            increased_cost_list = []
            t_max_original_list = [sequence[-1][-1] for sequence in self.all_time_stamps.values()]
            for control_time, total_time in zip(control_time_list, t_max_original_list):
                AV_time = total_time - control_time
                increased_cost_list.append( (1 - discounted_factor) * AV_time/ discounted_factor + (1 - inflated_factor) * control_time / inflated_factor)
            
            order = np.argsort(increased_cost_list)
        
            
        elif priority == 'random':
            order = np.arange(len(self.all_time_stamps))
            np.random.shuffle(order)
            
        else:
            pass
        
        return order
    
    
            
    def compute_budget_indicator(self):
        
        all_job_start, all_job_end = utils.extract_start_n_end_time(self.G_adjusted, self.new_all_time_stamps)
        new_all_time_stamps_list = list(set([element for innerList in all_job_start + all_job_end for element in innerList]))
        new_all_time_stamps_list.sort()
        
        Q = np.arange(len(new_all_time_stamps_list))
        new_all_time_stamps_list.append(self.t_max)

        budget_indicator = []
        
        for i in range(len(Q)):
            ss = new_all_time_stamps_list[i]
            tt = new_all_time_stamps_list[i+1]
            
            control_q = 0
            for v in range(len(self.new_all_time_stamps)):
                bb = False
                sss = all_job_start[v]
                ttt = all_job_end[v]
                for j in range(len(sss)):
                    if ss >= sss[j] and tt <= ttt[j]:
                        bb = True
                if bb == True:
                    control_q += 1 
            
            if control_q >= self.budget:
                budget_indicator.append(0)
            else:
                budget_indicator.append(1)
            
        return budget_indicator, new_all_time_stamps_list
            
            
            
    def solve(self):
        
        cost = 0
        D_bar = []
        
        # Determine the priority of re-routing
        order = self.sort_routes(priority = 'random')
        print("The re-routing order is: " + str(order))
        
        print(self.routes.keys())
        
        for v_2_re_route in order:
            
            
            self.new_routes[v_2_re_route] = self.routes[v_2_re_route]
            self.new_all_time_stamps[v_2_re_route] = self.all_time_stamps[v_2_re_route]
            
            violated = self.budget_violation_check()
            
            if not violated: 
                cost +=  self.original_route_cost(v_2_re_route)
                print("Route-"+str(v_2_re_route)+" does not violate the budget. Continue.")
                continue
            
            # Temperarily delete it since it violates the budget
            del self.new_routes[v_2_re_route]
            del self.new_all_time_stamps[v_2_re_route]
            print("Route-"+str(v_2_re_route)+" violates the budget. Try to re-route it.")
            
            
            # Compute the the constraints for the tsp below
            budget_indicator, new_all_time_stamps_list = self.compute_budget_indicator()
            
            # Call the constrained TSP solver
            success, route_cost, new_route, new_time_stamps = self.constrained_tsp_solver(v_2_re_route, budget_indicator, new_all_time_stamps_list)
            
            if success:
                cost += route_cost
                self.new_routes[v_2_re_route] = new_route
                self.new_all_time_stamps[v_2_re_route] = new_time_stamps
                print(colored("Re-routing succeeds!", 'green'))
            
            else:
                D_bar += [int(cus[-1][1:]) for cus in self.clusters[v_2_re_route][1:]]
        
        
        return cost, self.new_routes, self.new_all_time_stamps, D_bar, self.G_adjusted
    



# class ReRoutingFRPSolver_RandomOrder():
    
#     def __init__(self, G, I, routes, all_time_stamps, clusters, t_max, budget, gamma, num_layer, dir) -> None:
        
#         # Members from input
#         self.G = G
#         self.routes = routes
#         self.all_time_stamps = all_time_stamps
#         self.clusters = clusters
#         self.num_routes = len(clusters.keys())
#         self.demand = I.demands
#         self.capacity = I.capacity
#         self.budget = budget
#         self.t_max = t_max
#         self.gamma = gamma
#         self.num_layer = num_layer
#         self.dir = dir
        
#         # private members for inner use (dynamically updated)
#         self.D_bar = []
#         self.discarded_routes = [] # Need to be output as well
#         self.G_adjusted = G # Need to be output as well
#         self.budget_indicator = None 
#         self.other_time_stamps = None

        
#     def budget_violation_check(self, v_2_re_route):
        
#         all_job_start, all_job_end = utils.extract_start_n_end_time(self.G_adjusted, self.all_time_stamps)
        
#         # if v_2_re_route == 0:
#         #     print(np.array(all_job_start))
#         # all_job_start_dict = dict(zip(range(len(clusters)), all_job_start))
#         # all_job_end_dict = dict(zip(range(len(clusters)), all_job_end)) 
#         dispatched_vehicles = [key for key in self.clusters if self.clusters[key] != [] and key not in self.discarded_routes]
                
#         other_job_start = [r for i, r in enumerate(all_job_start) if i != v_2_re_route and i not in self.discarded_routes]
#         other_job_end = [r for i, r in enumerate(all_job_end) if i != v_2_re_route and i not in self.discarded_routes]
        
#         self.other_time_stamps = list(set([element for innerList in other_job_start + other_job_end for element in innerList]))
#         self.other_time_stamps.sort()
#         Q = np.arange(len(self.other_time_stamps))
#         self.other_time_stamps.append(self.t_max)
        
#         remaining_v = deepcopy(dispatched_vehicles)
#         remaining_v.remove(v_2_re_route)

#         self.budget_indicator = []
        
#         for i in range(len(Q)):
#             ss = self.other_time_stamps[i]
#             tt = self.other_time_stamps[i+1]
            
#             control_q = 0
#             for v in remaining_v:
#                 bb = False
#                 sss = all_job_start[v]
#                 ttt = all_job_end[v]
#                 for j in range(len(sss)):
#                     if ss >= sss[j] and tt <= ttt[j]:
#                         bb = True
#                 if bb == True:
#                     control_q += 1
            
#             if control_q >= self.budget:
#                 self.budget_indicator.append(0)
#             else:
#                 self.budget_indicator.append(1)
                
#         # Check whether the current route violates the budget
#         violated = False
#         for i in range(len(Q)):
            
#             if self.budget_indicator[i] == 1:
#                 continue
            
#             ss = self.other_time_stamps[i]
#             tt = self.other_time_stamps[i+1]
#             sss = all_job_start[v_2_re_route]
#             ttt = all_job_end[v_2_re_route]
#             for j in range(len(sss)):
#                 if ss >= sss[j] and tt <= ttt[j]:
#                     violated = True
                    
#         return violated
    
    
#     def re_route_one_vehicle(self, v_2_re_route):
        
#         # if v_2_re_route != 0:
#         #     return 0, 0,0 ,0
        
#         # print("===========")
#         # print(len(other_time_stamps))
#         # print(len(budget_indicator))
        
#         Q = np.arange(len(self.other_time_stamps)-1)
        
        
#         K = self.clusters[v_2_re_route][:1]
#         D = self.clusters[v_2_re_route][1:]
        
#         demand_subset = []
#         for cus in D:
#             cus_index = int(cus[-1][1:])
#             demand_subset.append(self.demand[cus_index])


#         for e in self.G.edges:
#             self.G[e[0]][e[1]]['is_artificial'] = 0
  
    
        
#         G_s = graph_pre_pruning(self.G, K, D)
        
    
#         # Add merged edges in the pruned graph to the original graph for time stamps extraction in next iteration
#         self.G_adjusted = utils.graph_edge_adjustment(G_s, self.G_adjusted)
    
#         V = list(G_s.nodes)
#         _, new_to_old = utils.re_index_nodes(V)
        
#         G_e, departure, arrival, D_e, D_e_2d, demand_dict = utils.contruct_time_expanded_graph(G_s, self.num_layer, K, D, demand_subset)
    
        
#         # if v_2_re_route == 0:
#         #     print("G_e information ================ ")
#         #     print(len(G_e.edges))
#         #     for e in G_e.edges:
#         #         print(G_e[e[0]][e[1]])
            
#         #TODO: The discount_factor and inflated_factor for AV_enabled roads are needed here to conduct sensitivity analysis
#         utils.assign_cost_n_travel_time(G_e)
        
#         E_e = list(G_e.edges)    
#         V_e = list(G_e.nodes)
#         edge_to_no, _ = utils.re_index_edges(E_e)
        

#         # Create a new model
#         model = gp.Model('re-routing-vehicle-'+str(v_2_re_route))

#         # Model params
        
#         # logging
#         path = self.dir + '/re_routing_log_vehicle_'+str(v_2_re_route)+'.txt'
#         if os.path.exists(path):
#             os.remove(path)
#         model.Params.LogFile = path
#         model.Params.LogToConsole = 0

#         # strategies
#         # model.Params.MIPFocus = 1 # Focusing on phase 1 (getting the first feasible solution)
#         model.Params.MIPFocus = 3 # Focusing on phase 2 (moving the best upper bound)
#         model.Params.ImproveStartTime = 120 # Give up proving optimality after 2 min
        
#         # Termination
#         model.Params.BestObjStop = 1.05 * self.original_route_cost(v_2_re_route)  # The original route cost is a natural lower bound. 
#                                                                              # A re-route is good if the cost is < 1.05 original route cost
        
#         model.Params.TimeLimit = 240 # Terminate the MILP solver after 5 min
#         model.Params.Cutoff = self.HDV_route_cost(v_2_re_route) # cutoff the solution with even worse cost than dispatching an HDV
        

#         # Create variables
#         #* The lower bounds for all variables are zero by default
#         x = model.addMVar(len(E_e), vtype=GRB.BINARY)
#         y = model.addMVar(len(V_e), vtype=GRB.BINARY)
#         t = model.addMVar(len(V_e), ub = np.repeat(self.t_max, len(V_e)), vtype=GRB.CONTINUOUS)

#         alpha = model.addMVar((len(E_e), len(Q)), vtype=GRB.BINARY)
#         beta = model.addMVar((len(E_e), len(Q)), vtype=GRB.BINARY)

#         #// makespan = model.addVar(ub = t_max, vtype=GRB.CONTINUOUS)

#         # Add constraints

#         ## Flow conservation at all nodes (except the depot)
#         for i in list(set(V_e) - set([departure, arrival])):
#             out_no = [edge_to_no[e] for e in list(G_e.out_edges(i))]
#             in_no = [edge_to_no[e] for e in list(G_e.in_edges(i))]
#             model.addConstr(sum(x[out_no]) == sum(x[in_no]))
#             model.addConstr(sum(x[out_no]) <= 1)
        
#         ## Flow conservation at the depot
#         departure_out = [edge_to_no[e] for e in list(G_e.out_edges(departure))]
#         arrival_in = [edge_to_no[e] for e in list(G_e.in_edges(arrival))]

#         model.addConstr(sum(x[departure_out]) == y[departure])
#         model.addConstr(sum(x[arrival_in]) == y[arrival])
#         model.addConstr(y[departure] == y[arrival])

#         ## Prohibit vehicles from passing through the depot
#         departure_in = [edge_to_no[e] for e in list(G_e.in_edges(departure))]
#         model.addConstr(x[departure_in] == 0)
        
#         ## Allow a vehicle to pass thru a customer without serving it 
#         for i in D_e:
#             i_in = [edge_to_no[e] for e in list(G_e.in_edges(i))]
#             model.addConstr(sum(x[i_in]) >= y[i])

#         ## Each customer is only served by one vehicle
#         model.addConstrs(sum(y[i] for i in D_e_2d[:, h]) == 1 for h in range(len(D)))

#         ## Transition
#         for i in D_e:
#             model.addConstr(sum(x[edge_to_no[(i, j)]] for j in D_e if G_e.has_edge(i, j) and G_e[i][j]['is_artificial'] == 1) <= y[i])

#         ## Capacity limit
#         model.addConstr(sum(demand_dict[i] * y[i] for i in D_e) <= self.capacity)


#         ## Subtour elimination + time tracking
#         for i in V_e:
#             for j in V_e:
#                 if G_e.has_edge(i, j):
#                     model.addConstr(t[j] >= t[i] + self.gamma[0] * G_e[i][j]['travel_time'] + self.t_max * (x[edge_to_no[(i, j)]] - 1))
#                     model.addConstr(t[j] <= t[i] + self.gamma[1] * G_e[i][j]['travel_time'] + self.t_max * (1 - x[edge_to_no[(i, j)]]))
                    
#         # The following constraints can be captured by the previous constraint
#         #// for e in E_e:
#         #//     i, j = e[:2]
#         #//     # We need one more set of important constraints, which equals the time stamp of the start node to that of the end node in any aritificial edges.
#         #//     if G_e[i][j]['is_artificial'] == 1:
#         #//         model.addConstr(t[j] <= t[i] + t_max * (1 - x[edge_to_no[(i, j)]]))
            
#         #//     if G_e[i][j]['is_av'] == 0:
#         #//         model.addConstr(t[j] <= t[i] + G_e[i][j]['travel_time'] + t_max * (1 - x[edge_to_no[(i, j)]]))
        
        
#         ## Do not violate the budget
#         #* This is simplied constraint compared to the original MILP
#         for e in E_e:
#             i, j = e[:2]
#             if G_e[i][j]['is_av'] == 0:
#                 model.addConstrs(x[edge_to_no[(i, j)]] == alpha[edge_to_no[(i, j)], q] + beta[edge_to_no[(i, j)], q]
#                                 for q in Q if self.budget_indicator[q] == 0)
#                 model.addConstrs(t[j] <= self.other_time_stamps[q] + self.t_max * (1 - alpha[edge_to_no[(i, j)], q])
#                                 for q in Q if self.budget_indicator[q] == 0)
#                 model.addConstrs(t[i] >= self.other_time_stamps[q+1] + self.t_max * (beta[edge_to_no[(i, j)], q] - 1)
#                             for q in Q if self.budget_indicator[q] == 0)

#         # Set the makespan
#         #// model.addConstrs(t[i] <= makespan for i in range(len(V_e)))

        
#         # Add the objective
#         model.setObjective(sum(x[edge_to_no[(i, j)]] * G_e[i][j]['av_cost'] for i in V_e for j in V_e if G_e.has_edge(i, j)),
#                         GRB.MINIMIZE)

#         # Solve the TSP
#         model.optimize()
        
        
#         # Infeasible (status code = 3), or objective is worse than the cutoff value (status code = 6), or no feasible solution found within the time limit (status code = 9)
#         if model.Status == 3 or model.Status == 6 or (model.Status == 9 and model.SolCount == 0):
#             print(model.Status)
#             return False, D, 0
        

#         # Extract the routes and time stamps
#         new_route, new_time_stamps = utils.extract_single_route_n_times(G_e, len(V), x, t, new_to_old, edge_to_no, departure, arrival)
        
#         self.routes[v_2_re_route] = new_route
#         self.all_time_stamps[v_2_re_route] = new_time_stamps 
    
#         #* Returned values are: re_routing success, unserved customers, obj value, and G_adjusted
#         return True, None, model.ObjVal
    
    
#     def original_route_cost(self, vehicle):
#         cost = 0
#         route = self.routes[vehicle]
#         for i, j in zip(route[:-1], route[1:]):
#                 cost += self.G[i][j]['av_cost']
        
#         return cost
    
    
#     #TODO: This is not the correct cutoff value. We should re-solve the TSP given the customers originally served by the vehicle
#     def HDV_route_cost(self, vehicle):
#         cost = 0
#         route = self.routes[vehicle]
#         for i, j in zip(route[:-1], route[1:]):
#                 cost += self.G[i][j]['non_av_cost']
        
#         return cost

            
#     def solve(self):
        
#         cost = 0
        
#         for v_2_re_route in range(self.num_routes):
            
#             violated = self.budget_violation_check(v_2_re_route)
            
#             if not violated: 
#                 cost +=  self.original_route_cost(v_2_re_route)
#                 print("Route-"+str(v_2_re_route)+" does not violate the budget. Continue.")
#                 continue
            
#             # Call the constrained TSP solver
#             print("Route-"+str(v_2_re_route)+" violates the budget. Try to re-route it.")
#             success, unserved, route_cost = self.re_route_one_vehicle(v_2_re_route)
            
#             cost += route_cost
            
#             if success:
#                 print(colored("Re-routing succeeds!", 'green'))
            
#             else:
#                 print(colored("Re-routing fails. Unserved customers will be served by HDVs.", 'yellow'))
#                 self.D_bar += [int(cus[-1][1:]) for cus in unserved]
#                 self.discarded_routes.append(v_2_re_route)
#                 continue
    
#         #// updated_routes = [routes[i] for i in range(num_routes) if i not in discarded_routes]
#         #// updated_all_time_stamps = [all_time_stamps[i] for i in range(num_routes) if i not in discarded_routes]
        
#         # Delete the routes that serve D_bar
#         for i in range(self.num_routes):
#             if i in self.discarded_routes:
#                 del self.routes[i]
#                 del self.all_time_stamps[i]
        
        
#         return cost, self.routes, self.all_time_stamps, self.D_bar, self.discarded_routes, self.G_adjusted
    



# def frp_rerouting(G, routes, all_time_stamps, clusters, t_max, demand, capacity, budget, gamma, num_layer, dir):
    
    

#     def budget_violation_check(v_2_re_route, discarded_routes, G_adjusted):
        
#         # Compute the budget indicator for the subsequent re-routing
        
#         all_job_start, all_job_end = utils.extract_start_n_end_time(G_adjusted, all_time_stamps)
#         # all_job_start_dict = dict(zip(range(len(clusters)), all_job_start))
#         # all_job_end_dict = dict(zip(range(len(clusters)), all_job_end)) 
#         dispatched_vehicles = [key for key in clusters if clusters[key] != [] and key not in discarded_routes]
                
#         other_job_start = [r for i, r in enumerate(all_job_start) if i != v_2_re_route and i not in discarded_routes]
#         other_job_end = [r for i, r in enumerate(all_job_end) if i != v_2_re_route and i not in discarded_routes]
        
#         other_time_stamps = list(set([element for innerList in other_job_start + other_job_end for element in innerList]))
#         other_time_stamps.sort()
#         Q = np.arange(len(other_time_stamps))
#         other_time_stamps.append(t_max)
        
#         remaining_v = deepcopy(dispatched_vehicles)
#         remaining_v.remove(v_2_re_route)

#         budget_indicator = []
        
#         for i in range(len(Q)):
#             ss = other_time_stamps[i]
#             tt = other_time_stamps[i+1]
            
#             control_q = 0
#             for v in remaining_v:
#                 bb = False
#                 sss = all_job_start[v]
#                 ttt = all_job_end[v]
#                 for j in range(len(sss)):
#                     if ss >= sss[j] and tt <= ttt[j]:
#                         bb = True
#                 if bb == True:
#                     control_q += 1
            
#             if control_q >= budget:
#                 budget_indicator.append(0)
#             else:
#                 budget_indicator.append(1)
                
#         # Check whether the current route violates the budget
#         violated = False
#         for i in range(len(Q)):
            
#             if budget_indicator[i] == 1:
#                 continue
            
#             ss = other_time_stamps[i]
#             tt = other_time_stamps[i+1]
#             sss = all_job_start[v_2_re_route]
#             ttt = all_job_end[v_2_re_route]
#             for j in range(len(sss)):
#                 if ss >= sss[j] and tt <= ttt[j]:
#                     violated = True
                    
#         return violated, budget_indicator, other_time_stamps
            
    
    
#     def re_route_one_vehicle(v_2_re_route, budget_indicator, other_time_stamps, G_adjusted):
        
#         # if v_2_re_route != 0:
#         #     return 0, 0,0 ,0
        
#         # print("===========")
#         # print(len(other_time_stamps))
#         # print(len(budget_indicator))
        
#         Q = np.arange(len(other_time_stamps)-1)
        
        
#         K = clusters[v_2_re_route][:1]
#         D = clusters[v_2_re_route][1:]
        
#         demand_subset = []
#         for cus in D:
#             cus_index = int(cus[-1][1:])
#             demand_subset.append(demand[cus_index])


#         for e in G.edges:
#             G[e[0]][e[1]]['is_artificial'] = 0
  
    
        
#         G_s = graph_pre_pruning(G, K, D)
        
    
#         # Add merged edges in the pruned graph to the original graph for time stamps extraction in next iteration
#         G_adjusted = utils.graph_edge_adjustment(G_s, G_adjusted)
    
#         V = list(G_s.nodes)
#         _, new_to_old = utils.re_index_nodes(V)
        
#         G_e, departure, arrival, D_e, D_e_2d, demand_dict = utils.contruct_time_expanded_graph(G_s, num_layer, K, D, demand_subset)
    
        
#         # if v_2_re_route == 0:
#         #     print("G_e information ================ ")
#         #     print(len(G_e.edges))
#         #     for e in G_e.edges:
#         #         print(G_e[e[0]][e[1]])
            
#         #TODO: The discount_factor and inflated_factor for AV_enabled roads are needed here to conduct sensitivity analysis
#         utils.assign_cost_n_travel_time(G_e)
        
#         E_e = list(G_e.edges)    
#         V_e = list(G_e.nodes)
#         edge_to_no, _ = utils.re_index_edges(E_e)
        

#         # Create a new model
#         model = gp.Model('re-routing-vehicle-'+str(v_2_re_route))

#         # Model params
        
#         # logging
#         path = dir + '/re_routing_log_vehicle_'+str(v_2_re_route)+'.txt'
#         if os.path.exists(path):
#             os.remove(path)
#         model.Params.LogFile = path
#         model.Params.LogToConsole = 0

#         # strategies
#         # model.Params.MIPFocus = 1 # Focusing on phase 1 (getting the first feasible solution)
#         model.Params.MIPFocus = 3 # Focusing on phase 2 (moving the best upper bound)
#         model.Params.ImproveStartTime = 120 # Give up proving optimality after 2 min
        
#         # Termination
#         model.Params.BestObjStop = 1.05 * original_route_cost(v_2_re_route)  # The original route cost is a natural lower bound. 
#                                                                              # A re-route is good if the cost is < 1.05 original route cost
        
#         model.Params.TimeLimit = 5 # Terminate the MILP solver after 5 min
#         model.Params.Cutoff = HDV_route_cost(v_2_re_route) # cutoff the solution with even worse cost than dispatching an HDV
        

#         # Create variables
#         #* The lower bounds for all variables are zero by default
#         x = model.addMVar(len(E_e), vtype=GRB.BINARY)
#         y = model.addMVar(len(V_e), vtype=GRB.BINARY)
#         t = model.addMVar(len(V_e), ub = np.repeat(t_max, len(V_e)), vtype=GRB.CONTINUOUS)

#         alpha = model.addMVar((len(E_e), len(Q)), vtype=GRB.BINARY)
#         beta = model.addMVar((len(E_e), len(Q)), vtype=GRB.BINARY)

#         #// makespan = model.addVar(ub = t_max, vtype=GRB.CONTINUOUS)

#         # Add constraints

#         ## Flow conservation at all nodes (except the depot)
#         for i in list(set(V_e) - set([departure, arrival])):
#             out_no = [edge_to_no[e] for e in list(G_e.out_edges(i))]
#             in_no = [edge_to_no[e] for e in list(G_e.in_edges(i))]
#             model.addConstr(sum(x[out_no]) == sum(x[in_no]))
#             model.addConstr(sum(x[out_no]) <= 1)
        
#         ## Flow conservation at the depot
#         departure_out = [edge_to_no[e] for e in list(G_e.out_edges(departure))]
#         arrival_in = [edge_to_no[e] for e in list(G_e.in_edges(arrival))]

#         model.addConstr(sum(x[departure_out]) == y[departure])
#         model.addConstr(sum(x[arrival_in]) == y[arrival])
#         model.addConstr(y[departure] == y[arrival])

#         ## Prohibit vehicles from passing through the depot
#         departure_in = [edge_to_no[e] for e in list(G_e.in_edges(departure))]
#         model.addConstr(x[departure_in] == 0)
        
#         ## Allow a vehicle to pass thru a customer without serving it 
#         for i in D_e:
#             i_in = [edge_to_no[e] for e in list(G_e.in_edges(i))]
#             model.addConstr(sum(x[i_in]) >= y[i])

#         ## Each customer is only served by one vehicle
#         model.addConstrs(sum(y[i] for i in D_e_2d[:, h]) == 1 for h in range(len(D)))

#         ## Transition
#         for i in D_e:
#             model.addConstr(sum(x[edge_to_no[(i, j)]] for j in D_e if G_e.has_edge(i, j) and G_e[i][j]['is_artificial'] == 1) <= y[i])

#         ## Capacity limit
#         model.addConstr(sum(demand_dict[i] * y[i] for i in D_e) <= capacity)


#         ## Subtour elimination + time tracking
#         for i in V_e:
#             for j in V_e:
#                 if G_e.has_edge(i, j):
#                     model.addConstr(t[j] >= t[i] + gamma[0] * G_e[i][j]['travel_time'] + t_max * (x[edge_to_no[(i, j)]] - 1))
#                     model.addConstr(t[j] <= t[i] + gamma[1] * G_e[i][j]['travel_time'] + t_max * (1 - x[edge_to_no[(i, j)]]))
                    
#         # The following constraints can be captured by the previous constraint
#         #// for e in E_e:
#         #//     i, j = e[:2]
#         #//     # We need one more set of important constraints, which equals the time stamp of the start node to that of the end node in any aritificial edges.
#         #//     if G_e[i][j]['is_artificial'] == 1:
#         #//         model.addConstr(t[j] <= t[i] + t_max * (1 - x[edge_to_no[(i, j)]]))
            
#         #//     if G_e[i][j]['is_av'] == 0:
#         #//         model.addConstr(t[j] <= t[i] + G_e[i][j]['travel_time'] + t_max * (1 - x[edge_to_no[(i, j)]]))
        
        
#         ## Do not violate the budget
#         #* This is simplied constraint compared to the original MILP
#         for e in E_e:
#             i, j = e[:2]
#             if G_e[i][j]['is_av'] == 0:
#                 model.addConstrs(x[edge_to_no[(i, j)]] == alpha[edge_to_no[(i, j)], q] + beta[edge_to_no[(i, j)], q]
#                                 for q in Q if budget_indicator[q] == 0)
#                 model.addConstrs(t[j] <= other_time_stamps[q] + t_max * (1 - alpha[edge_to_no[(i, j)], q])
#                                 for q in Q if budget_indicator[q] == 0)
#                 model.addConstrs(t[i] >= other_time_stamps[q+1] + t_max * (beta[edge_to_no[(i, j)], q] - 1)
#                             for q in Q if budget_indicator[q] == 0)

#         # Set the makespan
#         #// model.addConstrs(t[i] <= makespan for i in range(len(V_e)))

        
#         # Add the objective
#         model.setObjective(sum(x[edge_to_no[(i, j)]] * G_e[i][j]['av_cost'] for i in V_e for j in V_e if G_e.has_edge(i, j)),
#                         GRB.MINIMIZE)

#         # Solve the TSP
#         model.optimize()
        
        
#         # Infeasible (status code = 3), or objective is worse than the cutoff value (status code = 6), or no feasible solution found within the time limit (status code = 9)
#         if model.Status == 3 or model.Status == 6 or (model.Status == 9 and model.SolCount == 0):
#             return False, D, 0, G_adjusted
        

#         # Extract the routes and time stamps
#         new_route, new_time_stamps = utils.extract_single_route_n_times(G_e, len(V), x, t, new_to_old, edge_to_no, departure, arrival)
        
#         routes[v_2_re_route] = new_route
#         all_time_stamps[v_2_re_route] = new_time_stamps 
    
#         #* Returned values are: re_routing success, unserved customers, obj value, and G_adjusted
#         return True, None, model.ObjVal, G_adjusted
        
    
#     def original_route_cost(vehicle):
#         cost = 0
#         route = routes[vehicle]
#         for i, j in zip(route[:-1], route[1:]):
#                 cost += G[i][j]['av_cost']
        
#         return cost
    
#     def HDV_route_cost(vehicle):
#         cost = 0
#         route = routes[vehicle]
#         for i, j in zip(route[:-1], route[1:]):
#                 cost += G[i][j]['non_av_cost']
        
#         return cost
        
#     def main():
        
#         D_bar = []
#         cost = 0
#         G_adjusted = G
#         num_routes = len(clusters.keys())
#         discarded_routes = []
        
#         for v_2_re_route in range(num_routes):
            
#             # print(v_2_re_route)
#             violated, budget_indicator, other_time_stamps = budget_violation_check(v_2_re_route, discarded_routes, G_adjusted)
            
#             # print(budget_indicator)
#             if not violated: 
#                 cost +=  original_route_cost(v_2_re_route)
#                 print("Route-"+str(v_2_re_route)+" does not violate the budget. Continue.")
#                 continue
            
#             # Call the constrained TSP solver
#             print("Route-"+str(v_2_re_route)+" violates the budget. Try to re-route it.")
#             success, unserved, route_cost, G_adjusted = re_route_one_vehicle(v_2_re_route, budget_indicator, other_time_stamps, G_adjusted)
            
#             cost += route_cost
            
#             if success:
#                 print(colored("Re-routing succeeds!", 'green'))
            
#             else:
#                 print(colored("Re-routing fails. Unserved customers will be served by HDVs.", 'yellow'))
#                 D_bar += [int(cus[-1][1:]) for cus in unserved]
#                 discarded_routes.append(v_2_re_route)
#                 continue
    
#         #// updated_routes = [routes[i] for i in range(num_routes) if i not in discarded_routes]
#         #// updated_all_time_stamps = [all_time_stamps[i] for i in range(num_routes) if i not in discarded_routes]
        
#         # Delete the routes that serve D_bar
#         for i in range(num_routes):
#             if i in discarded_routes:
#                 del routes[i]
#                 del all_time_stamps[i]
        
        
#         return cost, routes, all_time_stamps, D_bar, discarded_routes, G_adjusted
    
#     #* The main entry of the frp-rerouting
#     return main()


class UHGS_CVRPSolver():
    
    def __init__(self, instance_path, dir, G, I, t = 3):
    
        self.tmp_path = dir + '/tmp_file_4_hgs.vrp'
        
        with open(instance_path, "r") as source:
            contents = source.read()
            
        with open(self.tmp_path, "w") as destination:
            destination.write(contents)
    
        self.dir = dir
        self.G = G
        self.I = I
        self.t = t
        
        self.result = None
    
        
        
    def add_explicit_dist_mtx(self, dist_matrix):
        
        with open(self.tmp_path, "a") as file:
            
            file.write("DIST_SECTION\n")
            for i in range(len(self.I.coordinates)):
                file.write(str(i+1))
                for j in range(len(self.I.coordinates)):
                    file.write(" " + str(dist_matrix[i][j]))
                file.write("\n")
                
                
    def solve(self, output_file):
        
        #* Call the C++ UHGS-CVRP solver to solve the instance
        arg_1 = 'HGS-CVRP/build/hgs'
        arg_2 = self.tmp_path
        arg_3 = self.dir + '/' + output_file
        self.result = arg_3
        
        arg_4 = '-t ' + str(self.t)
        arg_5 = '-round ' + str(0)
        arg_6 = '-seed ' + str(0)
        
        args_list = [arg_1, arg_2, arg_3, arg_4, arg_5, arg_6]
        
        os.system(' '.join(args_list))
        os.remove(self.tmp_path)
        
        
    def extract_results(self, path_matrix):
        
        routes_dict = {}
        routes_list = []
        clusters_dict = {}
        
        with open(self.result, 'r') as file:
            for line in file:
                if line.startswith('Route'):
                    
                    s1, s2 = line.split(': ')
                    route_num = int(s1.split('#')[-1]) - 1
                    
                    
                    s2 = [0] + [int(node_index) for node_index in s2.split(' ')] + [0]
                    route = []
                    cluster = []
                    
                    for i in range(len(s2)-1):
                        
                        x, y = self.I.coordinates[s2[i]]
                        x_next, y_next = self.I.coordinates[s2[i+1]]
                        node = (x, y, 'c'+str(s2[i]))
                        node_next = (x_next, y_next, 'c'+str(s2[i+1]))
                        path = path_matrix[str(node)+'-'+str(node_next)]
                        if i < len(s2) - 2:
                            route += path[:-1]
                        else:
                            route += path
                        cluster.append(node)
                        
                    
                    routes_dict[route_num] = route
                    routes_list.append(route)      
                    clusters_dict[route_num] = cluster
                    
                elif line.startswith('Cost'):
                    cost = float(line.split(' ')[-1])
        
        time_stamps = utils.time_tracker_all_routes(self.G, routes_list)
        time_stamps_dict = dict(zip(range(len(time_stamps)), time_stamps))
        
        return cost, routes_dict, time_stamps_dict, clusters_dict
        
        
        


        
    

#! There may be bugs in the solver for two reasons: (1) The result is not consistent with the HGS-CVRP solver. 
        #!(2) Some feasible solutions are sovled to be infeasible.
class GoogleORCVRPSolver():
    
    def __init__(self, G, I, path_matrix, dist_matrix, customers_subset = None, num_vehicles = None):
        
        self.G = G
        self.demands = I.demands
        self.coordinates = I.coordinates
        self.num_vehicles = int(I.name.split('-')[-1][1:])
        self.capacity = np.repeat(I.capacity, self.num_vehicles)
        self.depot = 0
        self.path_matrix = path_matrix
        self.dist_matrix = dist_matrix
        self.customers_subset = customers_subset
        
        
        if customers_subset != None:

            #* Only use a subset of customers to facilitate the last step of the algorithm
            
            self.num_vehicles= num_vehicles
            self.capacity = np.repeat(I.capacity, self.num_vehicles)
            
            
            
            self.customers_subset.append(0)
            self.customers_subset.sort()
            
            self.dist_matrix = self.dist_matrix[self.customers_subset][:, self.customers_subset]
        
            
            self.demands = [I.demands[cus] for cus in range(len(I.demands)) if cus in self.customers_subset]
            self.coordinates = [I.coordinates[cus] for cus in range(len(I.coordinates)) if cus in self.customers_subset]

            
            #// print(self.demands, self.capacity, self.dist_matrix, self.coordinates)
          
            
    def print_solution(self, manager, routing, solution):
        
        """Prints solution on console."""
        print(f'Objective: {solution.ObjectiveValue()}')
        
        total_distance = 0
        total_load = 0
        for vehicle_id in range(self.num_vehicles):
            index = routing.Start(vehicle_id)
            plan_output = 'Route for vehicle {}:\n'.format(vehicle_id)
            route_distance = 0
            route_load = 0
            while not routing.IsEnd(index):
                node_index = manager.IndexToNode(index)
                route_load += self.demands[node_index]
                plan_output += ' {0} Load({1}) -> '.format(node_index, route_load)
                previous_index = index
                index = solution.Value(routing.NextVar(index))
                route_distance += routing.GetArcCostForVehicle(
                    previous_index, index, vehicle_id)
            plan_output += ' {0} Load({1})\n'.format(manager.IndexToNode(index),
                                                    route_load)
            plan_output += 'Distance of the route: {}m\n'.format(route_distance)
            plan_output += 'Load of the route: {}\n'.format(route_load)
            print(plan_output)
            total_distance += route_distance
            total_load += route_load
        print('Total distance of all routes: {}m'.format(total_distance))
        print('Total load of all routes: {}'.format(total_load))
        
    
    def output_results(self, manager, routing, solution):
        
        routes_dict = {}
        routes_list = []
        clusters_dict = {}
                
        print("==========")
        
        
        for vehicle_id in range(self.num_vehicles):
            
            cluster = []
            route = []
            index = routing.Start(vehicle_id)
            
            while not routing.IsEnd(index):
                
                # Extract next node
                next_index = solution.Value(routing.NextVar(index))
                
                # This node
                node_index = manager.IndexToNode(index)
                x_coor, y_coor = self.coordinates[node_index]
                
                if self.customers_subset != None:
                    node_index = self.customers_subset[node_index]
                    
                node = (x_coor, y_coor, 'c'+str(node_index))
                cluster.append(node)
                
                # Next node
                next_node_index = manager.IndexToNode(next_index)
                next_x_coor, next_y_coor = self.coordinates[next_node_index]
                
                if self.customers_subset != None:
                    next_node_index = self.customers_subset[next_node_index]
                next_node = (next_x_coor, next_y_coor, 'c'+str(next_node_index))
                

                
                # Add a sub-route to the route
                path = self.path_matrix[str(node) + '-' + str(next_node)]
                if routing.IsEnd(next_index):
                    route += path
                else:
                    route += path[:-1]
                
                # Move to next node
                index = next_index
            
            routes_dict[vehicle_id] = route
            clusters_dict[vehicle_id] = cluster
            routes_list.append(route)
        
        time_stamps = utils.time_tracker_all_routes(self.G, routes_list)
        time_stamps_dict = dict(zip(range(len(time_stamps)), time_stamps))
        
        return routes_dict, time_stamps_dict, clusters_dict  
    
    def solve(self):
        
        """Solve the CVRP problem."""

        # Create the routing index manager.
        manager = pywrapcp.RoutingIndexManager(len(self.demands), self.num_vehicles, self.depot)

        # Create Routing Model.
        routing = pywrapcp.RoutingModel(manager)


        # Create and register a transit callback.
        def distance_callback(from_index, to_index):
            """Returns the distance between the two nodes."""
            # Convert from routing variable Index to distance matrix NodeIndex.
            from_node = manager.IndexToNode(from_index)
            to_node = manager.IndexToNode(to_index)
            return self.dist_matrix[from_node][to_node]

        transit_callback_index = routing.RegisterTransitCallback(distance_callback)

        # Define cost of each arc.
        routing.SetArcCostEvaluatorOfAllVehicles(transit_callback_index)


        # Add Capacity constraint.
        def demand_callback(from_index):
            
            """Returns the demand of the node."""
            # Convert from routing variable Index to demands NodeIndex.
            from_node = manager.IndexToNode(from_index)
            return self.demands[from_node]

        demand_callback_index = routing.RegisterUnaryTransitCallback(
            demand_callback)
        routing.AddDimensionWithVehicleCapacity(
            demand_callback_index,
            0,  # null capacity slack
            self.capacity,  # vehicle maximum capacities
            True,  # start cumul to zero
            'Capacity')

        # Setting first solution heuristic.
        search_parameters = pywrapcp.DefaultRoutingSearchParameters()
        search_parameters.first_solution_strategy = (
            routing_enums_pb2.FirstSolutionStrategy.AUTOMATIC)
        search_parameters.local_search_metaheuristic = (
            routing_enums_pb2.LocalSearchMetaheuristic.GUIDED_LOCAL_SEARCH)
        search_parameters.time_limit.seconds = 3

        # Solve the problem.
        solution = routing.SolveWithParameters(search_parameters)
                
        if solution:
            
            # Print solution on console.
            if self.customers_subset == None:
                self.print_solution(manager, routing, solution)
                
            routes, all_time_stamps, clusters = self.output_results(manager, routing, solution)
            
            return solution.ObjectiveValue(), routes, all_time_stamps, clusters
        else:
            print("The CVRP is not feasible")
            sys.exit()
        
   
#TODO: The generic TSP solver
class GoogleORTSPSolver():
    ...










  

"""
Deprocated
#! Do not use any functions below. They are depricated.
"""


def frp_re_routing_depricate(G, control, all_time_stamps, t_max, budget, vehicle_cost, clusters, demand, capacity, dir, num_layer=2):
#! The strategy that simultaneously remove more than one routes does NOT work
    
    if max(control['control_vector']) - budget <= 0:
        return (None, None)

    
    all_job_start, all_job_end = utils.extract_start_n_end_time(G, all_time_stamps)
    
    dispatched_vehicles = [int(key) for key in clusters if clusters[key] != []]
    vehicles_to_re_route = np.random.choice(dispatched_vehicles, int(max(control['control_vector'])
 - budget))
    remaining_dispacthed_vehicles = list(set(dispatched_vehicles) - set(vehicles_to_re_route))
    
    K = None
    D = []
    for key in vehicles_to_re_route:
        customer_subset = clusters[str(key)]
        if K == None:
            K = customer_subset[:1]
        D += customer_subset[1:]


    old_routes_start = [all_job_start[v] for v in remaining_dispacthed_vehicles]
    old_routes_end = [all_job_end[v] for v in remaining_dispacthed_vehicles]
    
    
    tmp = list(set([item for sublist in old_routes_start for item in sublist] + [item for sublist in old_routes_end for item in sublist]))
    tmp.sort()

    new_control_list = []
    budget_list = []
    for i in range(len(tmp)-1):
        s = tmp[i]
        t = tmp[i+1]
        
        new_control = 0
        for j in range(len(old_routes_start)):
            bb = False
            ss = old_routes_start[j]
            tt = old_routes_end[j]
            for k in range(len(ss)):
                if s >= ss[k] and t <= tt[k]:
                    bb = True
            if bb == True:
                new_control += 1
                
        budget_list.append(max(budget - new_control, 0))
        new_control_list.append(new_control)
        
    num_interval = len(tmp)
    Q = np.arange(num_interval)
    a = tmp
    b = tmp[1:] + [t_max * 3600]
    new_control_list.append(0)
    budget_list.append(budget)
    fleet = np.arange(len(vehicles_to_re_route))
    num_av = len(vehicles_to_re_route)
    

    # Pre-prune the graph before expanding it

    G_s = graph_pre_pruning(G, K, D)
    V = list(G_s.nodes)
    _, new_to_old = utils.re_index_nodes(V)
    utils.visualize_graph(G_s, filepath = dir + "/graph_pruned.png", customer = D, depot = K)
    
    G_e, departure, arrival, D_e, D_e_2d, demand_dict = utils.contruct_time_expanded_graph(G_s, num_layer, K, D, demand)
    utils.assign_cost_n_travel_time(G_e, 0.5, 1)
    
    E_e = list(G_e.edges())    
    V_e = list(G_e.nodes)
    edge_to_no, _ = utils.re_index_edges(E_e)





    # Create a new model
    model = gp.Model('re-routing')

    # Set the file path for logging
    path = dir + '/re_routing_log.txt'
    if os.path.exists(path):
        os.remove(path)
    model.Params.LogFile = path
    

    # Parameters tuning for the MIP solver
    # model.Params.MIPFocus = 1 # Focusing on phase 1 (getting the first feasible solution)
    model.Params.MIPFocus = 2 # Focusing on phase 3 (proving optimality) 
    # model.Params.Symmetry = 0
    # model.Params.MIPGap = 0.05 # Specify the terminating duality gap


    # Create variables
    #* The lower bounds for all variables are zero by default
    x = model.addMVar((len(E_e), len(fleet)), vtype=GRB.BINARY)
    y = model.addMVar((len(V_e), len(fleet)), vtype=GRB.BINARY)

    #* The following way to create variables save some space (helpful but not that much)
    # y = model.addMVar(num_layer * (len(D) + 1), num_truck, vtype=GRB.BINARY)


    t = model.addMVar((len(V_e), len(fleet)), vtype=GRB.CONTINUOUS)
    u = model.addMVar((num_interval, num_av), vtype=GRB.BINARY)
    alpha = model.addMVar((num_interval, len(E_e), num_av), vtype=GRB.BINARY)
    beta = model.addMVar((num_interval, len(E_e), num_av), vtype=GRB.BINARY)

    #? Do we really need the makespan variable? Because the routing cost has positive correlation with the makespan
    #! We do not need the makespan variable.
    makespan = model.addVar(ub = t_max * 3600, vtype=GRB.CONTINUOUS)


    # Add constraints

    ## Flow conservation at all nodes (except the depot)
    for i in list(set(V_e) - set([departure, arrival])):
        out_no = [edge_to_no[e] for e in list(G_e.out_edges(i))]
        in_no = [edge_to_no[e] for e in list(G_e.in_edges(i))]
        model.addConstrs(sum(x[out_no, m]) == sum(x[in_no, m]) for m in fleet)
        model.addConstrs(sum(x[out_no, m]) <= 1 for m in fleet)

    ## Flow conservation at the depot
    departure_out = [edge_to_no[e] for e in list(G_e.out_edges(departure))]
    arrival_in = [edge_to_no[e] for e in list(G_e.in_edges(arrival))]

    model.addConstrs(sum(x[departure_out, m]) == y[departure, m]
                    for m in fleet)
    model.addConstrs(sum(x[arrival_in, m]) == y[arrival, m]
                    for m in fleet)
    model.addConstrs(y[departure, m] == y[arrival, m]
                    for m in fleet)

    ## The # of truck used
    model.addConstr(sum(y[departure, :]) <= len(fleet))
    #? Why can't we set it to <= ? Anything wrong in the model? 
    #! The bug comes from the fact that we did not set flow conservation at the depot 
        #! while allowing vehicles to pass throught it without coming back
    model.addConstr(sum(y[departure, :]) <= len(fleet))
    #* We can prohibit vehicles from passing through the depot
    departure_in = [edge_to_no[e] for e in list(G_e.in_edges(departure))]
    model.addConstrs(x[departure_in, m] == 0 for m in fleet)

    ## Allow a vehicle to pass thru a customer without serving it 
    for i in D_e:
        i_in = [edge_to_no[e] for e in list(G_e.in_edges(i))]
        model.addConstrs(sum(x[i_in, m]) >= y[i, m] for m in fleet)

    ## Each customer is only served by one vehicle
    model.addConstrs(sum(y[i, m] for i in D_e_2d[:, h] for m in fleet) == 1 for h in range(len(D)))


    ## Transition
    for i in D_e:
        model.addConstrs(sum(x[edge_to_no[(i, j)], m] for j in D_e if G_e.has_edge(i, j) and G_e[i][j][0]['highway'] == 'artificial') <= y[i, m] 
                        for m in fleet)

    ## Range limit
    # model.addConstrs(sum(G_e[i][j][0]['length'] * x[edge_to_no[(i, j)], m] for i in V_e for j in V_e if G_e.has_edge(i, j)) <= vehicle_range[m]
    #                 for m in fleet)

    ## Capacity limit
    model.addConstrs(sum(demand_dict[i] * y[i, m] for i in D_e) <= capacity
                    for m in fleet)
                
    ## Subtour elimination + time tracking
    for i in V_e:
        for j in V_e:
            if G_e.has_edge(i, j):
                model.addConstrs(t[j, m] >= t[i, m] + G_e[i][j][0]['travel_time'] + t_max * 3600 * (x[edge_to_no[(i, j)], m] - 1) for m in fleet)
    for e in E_e:
        i, j = e[:2]
        
        #* We need one more set of important constraints, which equals the time stamp of the start node to that of the end node in any aritificial edges.
        if G_e[i][j][0]['highway'] == 'artificial':
            model.addConstrs(t[j, m] <= t[i, m] + t_max * 3600 * (1 - x[edge_to_no[(i, j)], m]) for m in fleet)
        
        if G_e[i][j][0]['is_av'] == 0:
            model.addConstrs(t[j, m] <= t[i, m] + G_e[i][j][0]['travel_time'] + t_max * 3600 * (1 - x[edge_to_no[(i, j)], m]) for m in fleet)
        
            
            
    #? To make the LP relaxation tighter, we can apply two more constraints: 
    #? (1) remove the flexibility of slowing down AVs in ordinary roads
    #? (2) Force all fake time stamps to be zero
    # for i in V_e:
    #     # if i == departure:
    #     #     continue
    #     i_out = [edge_to_no[e] for e in list(G_e.out_edges(i))]
    #     model.addConstrs(t_max * 3600 * sum(x[i_out, m]) >= t[i, m] for m in fleet)



    ### Budget of remote control
    model.addConstrs(sum(u[q, :]) <= budget_list[q] for q in Q)

    ### Time consistency of a truck being remotely controlled
    for i in V_e:
        for j in V_e:
            if G_e.has_edge(i, j) and G_e[i][j][0]['is_av'] == 0:
                model.addConstrs(3 * u[q, m] >= alpha[q, edge_to_no[(i, j)], m] + beta[q, edge_to_no[(i, j)], m] + x[edge_to_no[(i, j)], m] - 2
                                for q in Q
                                for m in fleet[: num_av])
                model.addConstrs(t[j, m]  <= a[q] + t_max * 3600 * (alpha[q, edge_to_no[(i, j)], m])
                                for q in Q
                                for m in fleet[: num_av])
                model.addConstrs(t[j, m] >= a[q] + t_max * 3600 * (alpha[q, edge_to_no[(i, j)], m] - 1)
                                for q in Q
                                for m in fleet[: num_av])
                model.addConstrs(t[i, m] >= b[q] - t_max * 3600 * (beta[q, edge_to_no[(i, j)], m])
                                for q in Q
                                for m in fleet[: num_av])
                model.addConstrs(t[i, m] <= b[q] - t_max * 3600 * (beta[q, edge_to_no[(i, j)], m] - 1)
                                for q in Q
                                for m in fleet[: num_av])
    ## Set the makespan 
    model.addConstrs(t[i, m] <= makespan for i in range(len(V_e)) for m in fleet)

    # Add the objective
    model.setObjective(sum(x[edge_to_no[(i, j)], m] * G_e[i][j][0]['av_cost'] for i in V_e for j in V_e for m in fleet[: num_av] if G_e.has_edge(i, j)) + 
                    sum(x[edge_to_no[(i, j)], m] * G_e[i][j][0]['non_av_cost'] for i in V_e for j in V_e for m in fleet[num_av :] if G_e.has_edge(i, j)) +
                    y[departure, :] @ np.repeat(vehicle_cost, len(fleet)) + makespan,
    #                 #    sum(u[q, :] @ control_cost for q in Q),
    #                    #* The last term by itself accounts for the extra cost of remote control
    #                    #* The last term can also avoid unneccesary positive u_{qm}
    #                    #* Adding this term affect the computational efficiency
                    GRB.MINIMIZE)


    # Solve the VRP
    model.optimize()
    
    # Extract the routes and time stamps

    new_routes, new_time_stamps = utils.extract_routes_n_times(G_e, len(V), x, t, new_to_old, edge_to_no, fleet, departure, arrival)
    new_routes = dict(zip(range(len(new_routes)), new_routes))
    new_time_stamps = dict(zip(range(len(new_time_stamps)), new_time_stamps))
    
    with open(dir + '/new_routes.json', 'w') as outfile:
        json.dump(new_routes, outfile)
    with open(dir + '/new_time_stamps.json', 'w') as outfile:
        json.dump(new_time_stamps, outfile)
    
    return (vehicles_to_re_route, G_s)
    


def frp_rerouting_depricated_2(G, routes, all_time_stamps, clusters, t_max, demand, capacity, budget, num_layer, dir):        
    
    v_to_re_route = 0
    G_adjusted = G
    while v_to_re_route != -1:
        v_to_re_route, G_adjusted, routes, all_time_stamps, clusters = re_route_one_vehicle_depricated(G, G_adjusted, routes, all_time_stamps, clusters, t_max, demand, capacity, budget, num_layer, dir)
        
    return routes, all_time_stamps


def re_route_one_vehicle_depricated(G, G_adjusted, routes, all_time_stamps, clusters, t_max, demand, capacity, budget, num_layer, dir):
    
    # f1 = open(dir + '/routes.json')
    # f2 = open(dir + '/all_time_stamps.json')
    # f3 = open(dir + '/clusters.json')

    # routes = json.load(f1)
    # all_time_stamps = json.load(f2)
    # clusters = json.load(f3)
  
    all_job_start, all_job_end = utils.extract_start_n_end_time(G_adjusted, all_time_stamps)
    
    print(all_job_start)

    # if vvv == 2:
    #     print(all_job_start)
    #     print(all_job_end)
    #     print("============")

    dispatched_vehicles = [int(key) for key in clusters if clusters[key] != []]

    all_time_stamps_list = list(set([element for innerList in all_job_start + all_job_end for element in innerList]))
    all_time_stamps_list.sort()

    Q = np.arange(len(all_time_stamps_list))
    all_time_stamps_list.append(t_max)
    # # Rank the vehicles by the occupancy of control
    # occupancy_time = []
    # for v in dispatched_vehicles:
    #     occ = sum(np.array(all_job_end[v]) - np.array(all_job_start[v]))
    #     occupancy_time.append(occ)

    
    # if vvv == 2:
    #     print(all_time_stamps_list)
    #     print("============")

    # extract violated v
    violated_v = []

    print(dispatched_vehicles)
    for i in range(len(Q)):
        ss = all_time_stamps_list[i]
        tt = all_time_stamps_list[i+1]
        
        control_q = 0
        tmp_v = []
        for v in dispatched_vehicles:
            bb= False
            sss = all_job_start[v]
            ttt = all_job_end[v]
            for j in range(len(sss)):
                if ss >= sss[j] and tt <= ttt[j]:
                    bb = True
            if bb == True:
                control_q += 1
                tmp_v.append(v)
                
            # if vvv == 2:
            #     print(bb, v, control_q)
            #     print(sss, ttt)
            #     print(ss, tt)
            #     print("================")
        
        if control_q > budget:
            violated_v += tmp_v

    violated_v = list(set(violated_v))
    
    # if vvv == 2:
    #     print(violated_v)
    #     print("================")
    
    #* Set v_to_re_route = -1 to indicate the completion of the rerouting FRP
    if violated_v == []:
        return (-1, G_adjusted, routes, all_time_stamps, clusters)
    
    # pick a violated vehicle
    # TODO: Debugging is needed. Why the algo keeps selecting vehicle 0 and vehicle 1?
    v_to_re_route = violated_v[0]

    # Depot, subset of customers
    K = clusters[v_to_re_route][:1]
    D = clusters[v_to_re_route][1:]
    
    # Demand of the subset of customers
    demand_subset = []
    for cus in D:
        demand_subset.append(demand[int(cus[-1][1:])])
        
        
    remaining_v = deepcopy(dispatched_vehicles)
    remaining_v.remove(v_to_re_route)


    budget_indicator = []
    for i in range(len(Q)):
        ss = all_time_stamps_list[i]
        tt = all_time_stamps_list[i+1]
        
        control_q = 0
        for v in remaining_v:
            bb = False
            sss = all_job_start[v]
            ttt = all_job_end[v]
            for j in range(len(sss)):
                if ss >= sss[j] and tt <= ttt[j]:
                    bb = True
            if bb == True:
                control_q += 1
        
        if control_q >= budget:
            budget_indicator.append(0)
        else:
            budget_indicator.append(1)
                
    
    for e in G.edges:
        G[e[0]][e[1]]['is_artificial'] = 0
    
    G_s = graph_pre_pruning(G, K, D)
    
    # Add merged edges in the pruned graph to the original graph for time stamps extraction in next iteration
    G_adjusted = utils.graph_edge_adjustment(G_s, G_adjusted)
    
    V = list(G_s.nodes)
    _, new_to_old = utils.re_index_nodes(V)
        
    G_e, departure, arrival, D_e, D_e_2d, demand_dict = utils.contruct_time_expanded_graph(G_s, num_layer, K, D, demand_subset)
    
    #TODO: The discount_factor and inflated_factor for AV_enabled roads are needed here to conduct sensitivity analysis
    utils.assign_cost_n_travel_time(G_e)
        
    E_e = list(G_e.edges)    
    V_e = list(G_e.nodes)
    edge_to_no, _ = utils.re_index_edges(E_e)

    # Create a new model
    model = gp.Model('re-routing-vehicle-'+str(v_to_re_route))

    # Set the file path for logging
    path = dir + '/re_routing_log_vehicle_'+str(v_to_re_route)+'.txt'
    if os.path.exists(path):
        os.remove(path)
    model.Params.LogFile = path
    model.Params.LogToConsole = 0

    model.Params.MIPFocus = 1 # Focusing on phase 1 (getting the first feasible solution)

    # Create variables
    #* The lower bounds for all variables are zero by default
    x = model.addMVar(len(E_e), vtype=GRB.BINARY)
    y = model.addMVar(len(V_e), vtype=GRB.BINARY)
    t = model.addMVar(len(V_e), ub = np.repeat(t_max, len(V_e)), vtype=GRB.CONTINUOUS)

    alpha = model.addMVar((len(E_e), len(Q)), vtype=GRB.BINARY)
    beta = model.addMVar((len(E_e), len(Q)), vtype=GRB.BINARY)

    makespan = model.addVar(ub = t_max, vtype=GRB.CONTINUOUS)

    # Add constraints

    ## Flow conservation at all nodes (except the depot)
    for i in list(set(V_e) - set([departure, arrival])):
        out_no = [edge_to_no[e] for e in list(G_e.out_edges(i))]
        in_no = [edge_to_no[e] for e in list(G_e.in_edges(i))]
        model.addConstr(sum(x[out_no]) == sum(x[in_no]))
        model.addConstr(sum(x[out_no]) <= 1)
        
    ## Flow conservation at the depot
    departure_out = [edge_to_no[e] for e in list(G_e.out_edges(departure))]
    arrival_in = [edge_to_no[e] for e in list(G_e.in_edges(arrival))]

    model.addConstr(sum(x[departure_out]) == y[departure])
    model.addConstr(sum(x[arrival_in]) == y[arrival])
    model.addConstr(y[departure] == y[arrival])

    ## Prohibit vehicles from passing through the depot
    departure_in = [edge_to_no[e] for e in list(G_e.in_edges(departure))]
    model.addConstr(x[departure_in] == 0)
        
    ## Allow a vehicle to pass thru a customer without serving it 
    for i in D_e:
        i_in = [edge_to_no[e] for e in list(G_e.in_edges(i))]
        model.addConstr(sum(x[i_in]) >= y[i])

    ## Each customer is only served by one vehicle
    model.addConstrs(sum(y[i] for i in D_e_2d[:, h]) == 1 for h in range(len(D)))


    ## Transition
    for i in D_e:
        model.addConstr(sum(x[edge_to_no[(i, j)]] for j in D_e if G_e.has_edge(i, j) and G_e[i][j]['is_artificial'] == 1) <= y[i])


    ## Capacity limit
    model.addConstr(sum(demand_dict[i] * y[i] for i in D_e) <= capacity)


    ## Subtour elimination + time tracking
    for i in V_e:
        for j in V_e:
            if G_e.has_edge(i, j):
                model.addConstr(t[j] >= t[i] + G_e[i][j]['travel_time'] + t_max * (x[edge_to_no[(i, j)]] - 1))
    for e in E_e:
        i, j = e[:2]
        
        #* We need one more set of important constraints, which equals the time stamp of the start node to that of the end node in any aritificial edges.
        if G_e[i][j]['is_artificial'] == 1:
            model.addConstr(t[j] <= t[i] + t_max * (1 - x[edge_to_no[(i, j)]]))
        
        if G_e[i][j]['is_av'] == 0:
            model.addConstr(t[j] <= t[i] + G_e[i][j]['travel_time'] + t_max * (1 - x[edge_to_no[(i, j)]]))
        
        
    ## Do not violate the budget
    #* This is simplied constraint compared to the original MILP
    for e in E_e:
        i, j = e[:2]
        if G_e[i][j]['is_av'] == 0:
            model.addConstrs(x[edge_to_no[(i, j)]] == alpha[edge_to_no[(i, j)], q] + beta[edge_to_no[(i, j)], q]
                            for q in Q if budget_indicator[q] == 0)
            model.addConstrs(t[j] <= all_time_stamps_list[q] + t_max * (1 - alpha[edge_to_no[(i, j)], q])
                            for q in Q if budget_indicator[q] == 0)
            model.addConstrs(t[i] >= all_time_stamps_list[q+1] + t_max * (beta[edge_to_no[(i, j)], q] - 1)
                        for q in Q if budget_indicator[q] == 0)

    ## Set the makespan
    model.addConstrs(t[i] <= makespan for i in range(len(V_e)))

        
    # Add the objective
    model.setObjective(sum(x[edge_to_no[(i, j)]] * G_e[i][j]['av_cost'] for i in V_e for j in V_e if G_e.has_edge(i, j)),
                    GRB.MINIMIZE)

    # Solve the VRP
    model.optimize()
    
    # Set v_to_re_route = -2 to indicate infeasibility
    if model.Status == 3:
        
        del routes[v_to_re_route]
        del all_time_stamps[v_to_re_route]
        del clusters[v_to_re_route]
        
        return (-2, G_adjusted, routes, all_time_stamps, clusters)


    # Extract the routes and time stamps
    new_route, new_time_stamps = utils.extract_single_route_n_times(G_e, len(V), x, t, new_to_old, edge_to_no, departure, arrival)

    routes[v_to_re_route] = new_route
    all_time_stamps[v_to_re_route] = new_time_stamps 
    
    
    # with open(dir + '/routes.json', 'w') as outfile:
    #     json.dump(routes, outfile)
    # with open(dir + '/all_time_stamps.json', 'w') as outfile:
    #     json.dump(all_time_stamps, outfile)
        
    return (v_to_re_route, G_adjusted, routes, all_time_stamps, clusters)
    