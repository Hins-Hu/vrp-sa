
"""
A modular of ultility functions 
"""

import numpy as np
import math
import random
import osmnx as ox
import networkx as nx
import os
import glob
import json
import secrets
import matplotlib.pyplot as plt


"""
Generate AV-enabled roads
"""

def generate_av_roads_random(G, rate):
    E = list(G.edges)
    av = random.sample(E, int(rate*len(E)))
    non_av = list(set(E) - set(av))
    
    for e in av:
        G.edges[e]['is_av'] = 1
    for e in non_av:
        G.edges[e]['is_av'] = 0
        
    return G

def generate_av_roads_primary(G):
    E = list(G.edges)
    for e in E:
        if G.edges[e]['highway'] == 'primary':
            G.edges[e]['is_av'] = 1
        else:
            G.edges[e]['is_av'] = 0
    
    return G
    
def generate_av_roads_primary_n_secondary(G):
    E = list(G.edges)
    for e in E:
        if G.edges[e]['highway'] == 'primary' or G.edges[e]['highway'] == 'secondary':
            G.edges[e]['is_av'] = 1
        else:
            G.edges[e]['is_av'] = 0
    
    return G



"""
Re-indexing nodes
"""
def re_index_nodes(V):
    
    old_to_new = {}
    new_to_old = {}
    index = 0
    for i in V:
        new_to_old[index] = i
        old_to_new[i] = index
        index += 1
    
    return (old_to_new, new_to_old)

"""
Re-indexing edges
"""
def re_index_edges(E):
    
    edge_to_no = {}
    no_to_edge = {}
    index = 0
    for e in E:
        no_to_edge[index] = e
        edge_to_no[e] = index
        index += 1

    return (edge_to_no, no_to_edge)


"""
Visualize a graph
"""

def visualize_graph(G, filepath, customer=None, depot=None):
    
    E = list(G.edges)
    V = list(G.nodes)

    index_auto_E = [i for i in range(len(E)) if G.edges[E[i]]['is_av'] == 1]
    node_color = ['w' for _ in range(len(V))]
    node_size = [15 for _ in range(len(V))]
    edge_color = ['#999999' for _ in range(len(E))]
    edge_linewidth = [1 for _ in range(len(E))]


    for idx in index_auto_E:
        edge_color[idx] = '#00FF00'
        edge_linewidth[idx] = 2
        
    if customer != None:
        
        old_to_new, _ = re_index_nodes(V)
        for i in customer:
            idx = old_to_new[i]
            node_color[idx] = 'r'
            node_size[idx] = 100
    
    if depot != None:
        
        old_to_new, _ = re_index_nodes(V)
        for i in depot:
            idx = old_to_new[i]
            node_color[idx] = 'b'
            node_size[idx] = 100


    ox.plot_graph(G, node_color = node_color, node_size = node_size, edge_color = edge_color, 
                    edge_linewidth = edge_linewidth, figsize = (30, 30), save = True, filepath = filepath)

"""
Visualize routes (for the generic VRP)
"""

def visualize_separate_routes(G, routes, num_av, dir):

    counter = 0
    for m in range(len(routes)):
        if len(routes[m]) <= 1:
            counter += 1
        if len(routes[m]) > 1:
            if m < num_av:
                color = '#FFC0CB'
            else:
                color = '#00FFFF'
            ox.plot.plot_graph_route(G, routes[m], route_color = color, route_linewidth=4,
                            route_alpha=0.8, orig_dest_size=100, figsize = (30, 30), save = True, filepath = dir + '/route-' + str(m - counter)+'.png')

    
"""
Generate a strongly-connected subgraph
"""
def subgraph_generator(G, num_nodes):
    V = list(G.nodes)
    E = list(G.edges)
    s = random.sample(V, 1)[0]

    visited = [s]
    frontier = [s]
    while len(visited) < num_nodes and len(frontier) > 0:
        i = frontier.pop(0)
        for j in G.neighbors(i):
            if j not in visited:
                visited.append(j)
                frontier.append(j)
    
    G = G.subgraph(visited)
    tmp_V = max(nx.strongly_connected_components(G), key = len)
    G = G.subgraph(list(tmp_V))
    
    return G


"""
Assign cost and travel time to edges
"""

def assign_cost_n_travel_time(G, discount_factor_av=0.5, inflated_factor_av = 1.2):
    
    
    # Average speed for urban vehicles in meter/sec
    speed = 1

    # Compute the travel time
    E = list(G.edges)
    for e in E:
        G.edges[e]['travel_time'] = G.edges[e]['length'] / speed

    # Compute the routing cost
    for e in E:
        if G.edges[e]['is_av'] == 1:
            G.edges[e]['av_cost'] = discount_factor_av * G.edges[e]['length']
            G.edges[e]['non_av_cost'] = G.edges[e]['length']
        else:
            G.edges[e]['av_cost'] = inflated_factor_av * G.edges[e]['length']
            G.edges[e]['non_av_cost'] = G.edges[e]['length']
            
            
"""
Clean a directory
"""
def clean_dir(dir):
    for item in glob.glob(dir+'/*'):
        os.remove(item)
        
        


"""
A time discretizer to produce a list of time intervals and two corresponding lists of start times and end times
"""

def time_discretizer(t_max, gran):
    
    t_max_sec = t_max * 3600
    num_interval = int((60 / gran) * t_max)
    time_intervals = np.arange(num_interval)
    start_time_list = np.arange(0, t_max_sec, gran * 60)
    end_time_list = start_time_list + gran * 60
    
    return (time_intervals, start_time_list, end_time_list)


"""
A time tracker for a given route
"""

def time_tracker(G, route):
    time_stamps = [(route[0], 0)]
    for i in range(1, len(route)):
        # time = time_stamp[i-1][1] + 2
        # print(r[i-1], r[i])
        time = time_stamps[i-1][1] + G[route[i-1]][route[i]]['travel_time']
        time_stamps.append((route[i], round(time)))
    
    return time_stamps



"""
A time tracker for a solution to a generic VRP
"""

def time_tracker_all_routes(G, routes):
    
    all_time_stamps = []
    for r in routes:
        if r == []:
            # all_time_stamps.append([])
            continue
        else:
            time_stamps = time_tracker(G, r)
            all_time_stamps.append(time_stamps)

    return all_time_stamps

"""
A control tracker
"""

def control_tracker(G, routes, fleet, num_av, t_max, gran):
    
    
    # Discrete the horizon 
    Q, a, b = time_discretizer(t_max, gran)

    control = np.zeros((len(fleet), len(Q)))
    for m in fleet[: num_av]:
        if routes[m] == []:
            continue

        time_stamp = time_tracker(G, routes[m])
        for q in Q:
            for i in range(len(time_stamp)-1):
                b1 = (a[q] <= time_stamp[i+1][1])
                b2 = (b[q] >= time_stamp[i][1])
                if b1 and b2:
                    s = time_stamp[i][0]
                    t = time_stamp[i+1][0]
                    if G[s][t][0]['is_av'] == 1:
                        control[m][q] = 0
                    else:
                        control[m][q] = 1
    
    return control


"""
Write down the experiment setting
"""

def write_input(seed, dir, num_node, num_customer, budget, num_av, num_non_av, capacity, t_max, gran, vehicle_cost, demand, num_layer=None):
    
    if num_layer == None:
        input = ['seed', 'num_node', 'num_customer', 'budget', 'num_av', 'num_non_av', 'capacity', 't_max', 'gran', 'vehicle_cost', 'demand']
        input_dict = dict(zip(input, [seed, num_node, num_customer, budget, num_av, num_non_av, capacity, t_max, gran, vehicle_cost, demand]))
    else:
        input = ['seed', 'num_node', 'num_customer', 'budget', 'num_av', 'num_non_av', 'capacity', 't_max', 'gran', 'demand', 'vehicle_cost', 'num_layer']
        input_dict = dict(zip(input, [seed, num_node, num_customer, budget, num_av, num_non_av, capacity, t_max, gran, demand, vehicle_cost, num_layer]))
    
    if not os.path.exists(dir):
        os.mkdir(dir)
    with open(dir + '/input.json', 'w') as outfile:
        json.dump(input_dict, outfile)
        
        
"""
Generate demand following lognormal distribution
"""
        
def demand_generator(num_customer, capacity):
    
    demand = np.random.lognormal(0.1, 1, size=num_customer)
    demand[demand >= capacity] = capacity
    
    return list(demand)


"""
Convert routes to a set of jobs
"""

def convert_route_to_job(G, routes):
    
    num_job = []
    idle_vehicle = []

    for key in routes:
        
        route = routes[key]
        
        if route == []:
            idle_vehicle.append(key)
            continue
        
        count = 0
        
        for ind in range(1, len(route)):
            
            if ind == len(route)-1:
                s_node = route[ind-1]
                m_node = route[ind]
                if G[s_node][m_node]['is_av'] == 0:
                    count += 1
            else:
                s_node = route[ind-1]
                m_node = route[ind]
                t_node = route[ind+1]
                
                if G[s_node][m_node]['is_av'] == 0:
                    count += 1
                    if G[m_node][t_node]['is_av'] == 0:
                        count -= 1
                
        num_job.append(count)
        
    return num_job


"""
Extract start times and end times from all time stamps 
"""

def extract_start_n_end_time(G, all_time_stamps):
    
    all_job_start = []
    all_job_end = []

    for key in all_time_stamps:
        
        time_stamps = all_time_stamps[key]
        
        if time_stamps == []:
            continue
        
        job_end_in_a_route = []
        job_start_in_a_route = []
        
        for i in range(len(time_stamps)):
            
            if i == 0:
                s_node, s_time = time_stamps[i]
                t_node = time_stamps[i+1][0]
                if G[s_node][t_node]['is_av'] == 0:
                    job_start_in_a_route.append(s_time)
                    
            elif i == len(time_stamps)-1:
                s_node = time_stamps[i-1][0]
                t_node, t_time = time_stamps[i]
                if G[s_node][t_node]['is_av'] == 0:
                    job_end_in_a_route.append(t_time)
                
            else:
                s_node, s_time = time_stamps[i-1]
                m_node, m_time = time_stamps[i]
                t_node, t_time = time_stamps[i+1]
                if G[s_node][m_node]['is_av'] == 0:
                    job_end_in_a_route.append(m_time)
                    if G[m_node][t_node]['is_av'] == 0:
                        job_end_in_a_route.pop()
                else:
                    job_start_in_a_route.append(m_time)
                    if G[m_node][t_node]['is_av'] == 1:
                        job_start_in_a_route.pop()
                
        all_job_end.append(job_end_in_a_route)
        all_job_start.append(job_start_in_a_route)
        
    return (all_job_start, all_job_end)




def contruct_time_expanded_graph(G, num_layer, depot, customer, demand):

    K = depot
    D = customer
    
    V = list(G.nodes)
    old_to_new, _ = re_index_nodes(V)
    

    # G_e = copy.deepcopy(G)
    G_e = nx.convert_node_labels_to_integers(G)
    D_e = [old_to_new[i] for i in D]
    for i in range(num_layer - 1):

        # Duplicate one more layer
        # layer = copy.deepcopy(G) 
        layer = nx.convert_node_labels_to_integers(G, first_label = (i+1) * len(V))
        G_e = nx.disjoint_union(G_e, layer)

        # Construct D_e
        for j in D:
            D_e.append(old_to_new[j] + (i+1) * len(V))

    # # Construct the demand dictionary
    demand_dict = dict(zip(D_e, np.tile(demand, num_layer)))

    # Construct the 2-D customer vertices
    D_e_2d = np.reshape(D_e, (num_layer, len(D)))

    # Add artificial edges of the depot
    departure = old_to_new[K[0]]
    K_e = []
    for i in range(num_layer - 1):
        idx_1 = departure + i * len(V)
        idx_2 = departure + (num_layer - 1) * len(V)
        K_e.append(idx_1)
        G_e.add_edge(idx_1, idx_2, length = 0, is_artificial = 1, is_av = 1, index = str(idx_1)+'-'+str(idx_2))
    arrival = departure + (num_layer - 1) * len(V)
    K_e.append(arrival)

    # Construct the set of complete undirected graphs for all customers
    for h in range(len(D)):
        D_h = D_e_2d[:, h]
        for u in D_h:
            for v in D_h:
                if u != v:
                    G_e.add_edge(u, v, length = 0, is_artificial = 1, is_av = 1, index = str(u)+'-'+str(v))


    return (G_e, departure, arrival, D_e, D_e_2d, demand_dict)
    # # Make new lists of nodes and edges in the extended graph
    # E_e = list(G_e.edges())
    # V_e = list(G_e.nodes)


    # # Mapping for edges
    # edge_to_no, no_to_edge = utils.re_index_edges(E_e)


    # # Sanity Check 
    # counter = 0
    # for j in range(num_layer - 1):

    #     tmp_k = old_to_new[K[0]]
    #     print(G_e[tmp_k + j * len(V)][tmp_k + (num_layer-1) * len(V)])
    #     print('++++++++++++++++++++++')
    #     for i in D:
    #         tmp_i = old_to_new[i]
    #         print(G_e[tmp_i + j * len(V)][tmp_i + (j+1) * len(V)])
    #         counter += 1
    #     print('======================')


    # # Assign costs and travel times
    # utils.assign_cost_n_travel_time(G_e, discount_factor_av = 0.5, inflated_factor_av = 1) 


    # # Sanity check
    # for h in range(len(D)):
    #     print("========"+str(h))
    #     D_h = D_e_2d[:, h]
    #     for u in D_h:
    #         for v in D_h:
    #             if u != v:
    #                 print(G_e[u][v])


    # print("The size of the extended graph is: ")
    # print(len(E_e))
    # print(len(V_e))



def extract_routes_n_times(G, num_node_base, x, t, new_to_old, edge_to_no, fleet, departure, arrival):
    
    routes = []
    all_time_stamps = []
    for m in fleet:
        
        departure_out_edge_no = [edge_to_no[e] for e in list(G.out_edges(departure))]
        departure_outflow = x[departure_out_edge_no, m].X
        
        if sum(departure_outflow) == 0:
            routes.append([])
            all_time_stamps.append([])
            continue
        
        route_m = [new_to_old[departure]]
        time_stamp_m = [(new_to_old[departure], round(t[departure, m].X))]
        
        i = departure
        while i != arrival:

            i_out_edges = list(G.out_edges(i))
            i_out_edges_no = [edge_to_no[e] for e in i_out_edges]
            i_outflow = x[i_out_edges_no, m].X
            next_edge = i_out_edges[np.where(i_outflow)[0][0]]
            next_node = next_edge[1]
            if next_node % num_node_base != i % num_node_base:
                route_m.append(new_to_old[next_node % num_node_base])
                time_stamp_m.append((new_to_old[next_node % num_node_base], round(t[next_node, m].X)))
            i = next_node
            # if m == 1:
                # tmp_route.append(i)
        routes.append(route_m)
        all_time_stamps.append(time_stamp_m)
        
    return (routes, all_time_stamps)


def extract_single_route_n_times(G, num_node_base, x, t, new_to_old, edge_to_no, departure, arrival):
    
        
    departure_out_edge_no = [edge_to_no[e] for e in list(G.out_edges(departure))]
    departure_outflow = x[departure_out_edge_no].X
    
    if sum(departure_outflow) == 0:
        return ([], [])
    
    route = [new_to_old[departure]]
    time_stamps = [(new_to_old[departure], round(t[departure].X))]
    
    i = departure
    while i != arrival:

        i_out_edges = list(G.out_edges(i))
        i_out_edges_no = [edge_to_no[e] for e in i_out_edges]
        i_outflow = x[i_out_edges_no].X
        next_edge = i_out_edges[np.where(i_outflow)[0][0]]
        next_node = next_edge[1]
        if next_node % num_node_base != i % num_node_base:
            route.append(new_to_old[next_node % num_node_base])
            time_stamps.append((new_to_old[next_node % num_node_base], round(t[next_node].X)))
        i = next_node
        # if m == 1:
            # tmp_route.append(i)

    return (route, time_stamps)    


def generate_random_colors_hex(length):
    
    color_list = []    
    for _ in range(length):
        s = "#"
        x = 0
        while x < 6:
            s += secrets.choice("0123456789ABCDEF")
            x += 1
        
        color_list.append(s)
    
    return color_list








def graph_edge_adjustment(G_s, G):
    
    E_s = list(G_s.edges)
    G_adjusted = G.copy()
    for e in E_s:
        i, j = e[:2]
        if not G_adjusted.has_edge(i, j):
            G_adjusted.add_edge(i, j, is_av = G_s[i][j]['is_av'], is_artificial = 0, length = G_s[i][j]['length'])

    return G_adjusted




# PART 2: ==========================

"""
Genenrate a graph with only grid points
"""

def av_grid_pts_generator(x_range, y_range, n):
    
    x_step = math.ceil((x_range[1] - x_range[0]) / n)
    y_step = math.ceil((y_range[1] - y_range[0]) / n)
    
    
    G = nx.DiGraph(n = n, x_range = x_range, y_range = y_range, x_step = x_step, y_step = y_step)
    
    x_coor = np.arange(x_range[0], x_range[0] + n*x_step + 1, x_step)
    y_coor = np.arange(y_range[0], y_range[0] + n*y_step + 1, y_step)
    
    for x in x_coor:
        for y in y_coor:
            G.add_node((x, y))
            
    return G


"""
Node_zoning helper function for grid graph generator
"""

def node_zoning(coordinates, G):
    
    zone_dict = {}
    x_min, _ = G.graph['x_range']
    y_min, _ = G.graph['y_range']
    x_step = G.graph['x_step']
    y_step = G.graph['y_step']
    n = G.graph['n']
    
    
    vertical = np.arange(x_min, x_min + (n + 1) * x_step, x_step)
    horizontal = np.arange(y_min, y_min + (n + 1) * y_step, y_step)
        
    for node in range(len(coordinates)):
        x, y = coordinates[node]
        
        
        if x in vertical or y in horizontal:
            
            if 'boundary' not in zone_dict.keys():
                zone_dict['boundary'] = [node]
            else:
                tmp_L = zone_dict['boundary']
                tmp_L += [node]
                zone_dict['boundary'] = tmp_L
            
        else:
            i = np.where(vertical > x)[0][0] - 1
            j = np.where(horizontal > y)[0][0] - 1
                        
            if (i, j) not in zone_dict.keys():
                zone_dict[(i, j)] = [node]
            else:
                tmp_L = zone_dict[(i, j)]
                tmp_L += [node]
                zone_dict[(i, j)] = tmp_L
    
    return zone_dict


"""
A helper function to generate a grid grapgh given a set of points
"""

def generate_grid_by_points(points):
    
    
    # Points: a list of coordintates
    
    
    pts = np.array(points)
    
    x_coor = np.unique(pts[:, 0])
    y_coor = np.unique(pts[:, 1])
    
    H_h = nx.DiGraph()
    for i in range(len(x_coor)-1):
        H_h.add_edge(x_coor[i], x_coor[i+1], length = x_coor[i+1] - x_coor[i], is_av = 0)
        H_h.add_edge(x_coor[i+1], x_coor[i], length = x_coor[i+1] - x_coor[i], is_av = 0)
    
    H_v = nx.DiGraph()
    for j in range(len(y_coor)-1):
        H_v.add_edge(y_coor[j], y_coor[j+1], length = y_coor[j+1] - y_coor[j], is_av = 0)
        H_v.add_edge(y_coor[j+1], y_coor[j], length = y_coor[j+1] - y_coor[j], is_av = 0)
        
    H = nx.cartesian_product(H_h, H_v)
    return H 


"""
A helper function to add customers and the depot to a grid-like graph
"""

def add_customer_n_depot_to_grid(G, coordinates, zone_dict):
    
    x_range = G.graph['x_range']
    y_range = G.graph['y_range']
    x_step = G.graph['x_step']
    y_step = G.graph['y_step']
    n = G.graph['n']
    
    # G = nx.DiGraph()
    
    # Add regional network into the grid
    for key in list(zone_dict.keys()):
        if key == 'boundary':
            continue
        else:
            subset = zone_dict[key]
            subset_coor = [coordinates[pt] for pt in subset]
            
            LL = [x_range[0] + key[0] * x_step, y_range[0] + key[1] * y_step]
            LR = [x_range[0] + (key[0] + 1) * x_step, y_range[0] + key[1] * y_step]
            UL = [x_range[0] + key[0] * x_step, y_range[0] + (key[1] + 1) * y_step]
            UR = [x_range[0] + (key[0] + 1) * x_step, y_range[0] + (key[1] + 1) * y_step]
            
            subset_coor += [LL, LR, UL, UR]
            
            H = generate_grid_by_points(subset_coor)
            # print('H:', nx.is_strongly_connected(H))
            G = nx.compose(H, G)
            # print('G:', nx.is_strongly_connected(G))
            


    # Add boundary customers (or depot) to the graph
    if 'boundary' in zone_dict.keys():
        for item in zone_dict['boundary']:
            x, y = coordinates[item]
            if not G.has_node((x, y)):
                G.add_node((x, y))



    # Regulate all edges aligned with the grid points     
    E = list(G.edges)
    V = list(G.nodes)
    
    for i in range(n+1):
        
        x_line = x_range[0] + i * x_step
        for e in E:
            if e[0][0] == x_line and e[1][0] == x_line and G.has_edge(e[0], e[1]):
                G.remove_edge(e[0], e[1])
        
        subset_i = []
        for node in V:
            if node[0] == x_line:      
                subset_i.append(node[1])  
        subset_i.sort()
        
        for j in range(len(subset_i) - 1):
            G.add_edge((x_line, subset_i[j]), (x_line, subset_i[j+1]), 
                       is_av = 1, length = subset_i[j+1] - subset_i[j])
            G.add_edge((x_line, subset_i[j+1]), (x_line, subset_i[j]), 
                       is_av = 1, length = subset_i[j+1] - subset_i[j])
        
    for j in range(n+1):
        
        y_line = y_range[0] + j * y_step
        for e in E:
            if e[0][1] == y_line and e[1][1] == y_line and G.has_edge(e[0], e[1]):
                G.remove_edge(e[0], e[1])
                
        subset_j = []
        for node in V:
            if node[1] == y_line:
                subset_j.append(node[0])
        subset_j.sort()
        
        for i in range(len(subset_j) - 1):
            G.add_edge((subset_j[i], y_line), (subset_j[i+1], y_line), 
                       is_av = 1, length = subset_j[i+1] - subset_j[i])
            G.add_edge((subset_j[i+1], y_line), (subset_j[i], y_line), 
                       is_av = 1, length = subset_j[i+1] - subset_j[i])
    
    return G



def mark_customer_n_depot(G, coordinates):
    
    for idx in range(len(coordinates)):
        x, y = coordinates[idx]
        G = nx.relabel_nodes(G, {(x, y): (x, y, 'c'+str(idx))})

    return G

def visualize_graph_nx(G, dir):

    node_color = []
    node_size = []
    for v in G.nodes:
        if type(v[-1]) != str:
            node_color.append('black')
            node_size.append(2)
        else:
            if v[-1] == 'c0':
                node_color.append('blue')
                node_size.append(25)
            else:
                node_color.append('red')
                node_size.append(25)
            # else:
            #     node_color.append('purple')
            #     node_size.append(25)
                

    edge_color = []
    width = []
    for e in G.edges:
        if G.edges[e]['is_av'] == 1:
            edge_color.append('green')
            width.append(3)
        else:
            edge_color.append('black')
            width.append(.5)
            
    pos = {}
    for node in G.nodes:
        pos[node] = [node[0], node[1]]
        
    fig = plt.figure()
    nx.draw(G, pos, arrowstyle = '-', node_color = node_color, edge_color = edge_color, 
            node_size = node_size, width = width, with_labels=False)
    fig.savefig(dir + '/network.png',  dpi=300)
    
    
def dist_n_path_matrices_4_experiments(G, I):
    
    
    K_D = [I.depot] + I.customers
    cost_non_av = np.zeros((len(K_D), len(K_D)))
    cost_av = np.zeros((len(K_D), len(K_D)))
    path_non_av = {}
    path_av = {}

    for i in range(len(K_D)):
        for j in range(len(K_D)):
            if i!=j:
                
                # if i == 0:
                #     s_node = (I.coordinates[i][0], I.coordinates[i][1], 'd')
                #     t_node = (I.coordinates[j][0], I.coordinates[j][1], 'c'+str(j))
                    
                # elif j == 0:
                #     s_node = (I.coordinates[i][0], I.coordinates[i][1], 'c'+str(i))
                #     t_node = (I.coordinates[j][0], I.coordinates[j][1], 'd')
                    
                s_node = (I.coordinates[i][0], I.coordinates[i][1], 'c'+str(i))
                t_node = (I.coordinates[j][0], I.coordinates[j][1], 'c'+str(j))
                    
                cost_non_av[i][j] = nx.shortest_path_length(G, s_node, t_node, weight='non_av_cost')
                cost_av[i][j] = nx.shortest_path_length(G, s_node, t_node, weight='av_cost')

                path_non_av[str(s_node)+'-'+str(t_node)] = nx.shortest_path(G, s_node, t_node, weight='non_av_cost')
                path_av[str(s_node)+'-'+str(t_node)] = nx.shortest_path(G, s_node, t_node, weight='av_cost')
    
    return cost_av, cost_non_av, path_av, path_non_av



def visualize_re_scheduling(dir):
    

    old_schedule_start = list(json.load(open(dir + '/old_schedule_start.json')).values())
    old_schedule_end = list(json.load(open(dir + '/old_schedule_end.json')).values())
    new_schedule_end = list(json.load(open(dir + '/new_schedule_end.json')).values())
    new_schedule_start = list(json.load(open(dir + '/new_schedule_start.json')).values())

    fig, ax = plt.subplots()

    old_tick_pos = [(v+1)*20 for v in range(len(old_schedule_start))]
    old_tick_name = [str(v) for v in range(len(old_schedule_start))]
    new_tick_pos = [(v+1)*20 - 6 for v in range(len(old_schedule_start))]
    new_tick_name = [str(v)+'-new' for v in range(len(old_schedule_start))]
    bar_width = 5


    # Old schedules
    for v in range(len(old_schedule_start)):
        old = [np.vstack((s, t)) for s, t in zip(old_schedule_start, old_schedule_end)][v]
        old = [x for x in zip(*old)]
        old = [(st[0], st[1] - st[0]) for st in old]
        ax.broken_barh(old, ((v+1)*20, bar_width), facecolors = 'tab:blue')

    # New schedules
    for v in range(len(new_schedule_start)):
        new = [np.vstack((s, t)) for s, t in zip(new_schedule_start, new_schedule_end)][v]
        new = [x for x in zip(*new)]
        new = [(st[0], st[1] - st[0]) for st in new]
        ax.broken_barh(new, ((v+1)*20-6, bar_width), facecolors = 'tab:red')

    ax.set_xlabel('Time intervals when vehicles need remote control')
    ax.set_yticks(old_tick_pos + new_tick_pos, labels = old_tick_name + new_tick_name)
    plt.savefig(dir + '/re-scheduling.png')
    
    
    
def visualize_re_routing(dir, G_adjusted, all_time_stamps):


    old_routes_start = list(json.load(open(dir + '/old_schedule_start.json')).values())
    old_routes_end = list(json.load(open(dir + '/old_schedule_end.json')).values())


    new_routes_start, new_routes_end = extract_start_n_end_time(G_adjusted, all_time_stamps)

    fig, ax = plt.subplots()

    old_tick_pos = [(v+1)*20 for v in range(len(old_routes_start))]
    old_tick_name = [str(v) for v in range(len(old_routes_start))]
    new_tick_pos = [(v+1)*20 - 6 for v in range(len(old_routes_start))]
    new_tick_name = [str(v)+'-new' for v in range(len(old_routes_start))]
    bar_width = 5


    # Old
    for v in range(len(old_routes_start)):
        old = [np.vstack((s, t)) for s, t in zip(old_routes_start, old_routes_end)][v]
        old = [x for x in zip(*old)]
        old = [(st[0], st[1] - st[0]) for st in old]
        ax.broken_barh(old, ((v+1)*20, bar_width), facecolors = 'tab:blue')


    # New    
    vehicles_mapping = list(all_time_stamps.keys())
    
    for v in range(len(new_routes_start)):
        new = [np.vstack((s, t)) for s, t in zip(new_routes_start, new_routes_end)][v]
        new = [x for x in zip(*new)]
        new = [(st[0], st[1] - st[0]) for st in new]
        
        v = vehicles_mapping[v]
        ax.broken_barh(new, ((v+1)*20-6, bar_width), facecolors = 'tab:red')

    ax.set_xlabel('Time intervals when vehicles need remote control')
    ax.set_yticks(old_tick_pos + new_tick_pos, labels = old_tick_name + new_tick_name)
    plt.savefig(dir + '/re-routing.png')




