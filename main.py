import cvrplib 
import utils
import random
import numpy as np
import algorithms as algo
import argparse
import csv
from termcolor import colored
import sys
import os


"""
Arguments
"""


parser = argparse.ArgumentParser(description ='The algorithmic pipeline')
 
parser.add_argument('instance', help = "The CVRP instance used") 
parser.add_argument('-p', dest = 'path', help = "The output path",
                    default = None, required = False)
parser.add_argument('-g', dest = "grid_size", help = "The # of neighbhoods divided by the grid in the underlying network",
                    default = 5, required = False, type = int, choices = range(1, 11))
parser.add_argument('-bf', dest = "budget_factor", help = "The budget factor", 
                    default = 2, required = False, type = int, choices = range(1, 4))
parser.add_argument('-b', dest = "budget", help = "The budget for the additional resources (i.e. remote control). Overwrite the function of budget factor.", 
                    default = None, required = False)
parser.add_argument('-t', dest = 't_max_factor', help = "The maximum delay factor for the makespan",
                    default = 1.5, required = False, type = float)
parser.add_argument('-g1', dest = 'gamma_1', help = "The minimum factor of travel time decrease",
                    default = 1.0, required = False, type = float) 
parser.add_argument('-g2', dest = 'gamma_2', help = "The maximum factor of travel time increase",
                    default = 1.0, required = False, type = float)
parser.add_argument('-d', dest = 'discount_factor_av', help = "The discount factor for AVs in AV-enabled roads",
                    default = 0.5, required = False, type = float)
parser.add_argument('-i', dest = 'inflated_factor_av', help = "The inflated factor for AVs in ridinary roads",
                    default = 1.2, required = False, type = float)
parser.add_argument('-nl', dest = 'num_layer', help = "The number of layers in the expanded graph",
                    default = 2, required = False, type = int)

args = parser.parse_args()




"""
Params
"""

# Link to the instance file
instance_path = 'Instances/' + args.instance + '.vrp'

# Create a dir

instance = args.instance
dir = args.path
if args.path == None:
    os.mkdir(instance)
    dir = instance
else:
    dir = args.path
    
# Set seeds
seed = 1
random.seed(seed)
np.random.seed(seed)

# The grid size
grid_size = args.grid_size

# The number of vehicles
# TODO: The naming system might be different in other instances
num_vehicle = int(args.instance.split('-')[-1][1:])

# The budget for the additional resources (i.e. remote control)
budget_factor = [1/3, 1/2, 2/3]
budget = max(1, round(num_vehicle * budget_factor[args.budget_factor - 1]))

if args.budget != None:
    budget = args.budget

# The maximum delay factor for the makespan
t_max_factor = args.t_max_factor

# The factors for flexible travel time assumption 
gamma = [args.gamma_1, args.gamma_2]

# The cost adjustment factors for AVs 
discount_factor_av = args.discount_factor_av
inflated_factor_av = args.inflated_factor_av

# The number of layers in the expanded graph
num_layer = args.num_layer





"""
Instance and Network
"""

# Load the CVRP instances
I = cvrplib.download(args.instance)

# Boundary of the instance
x_min, y_min = np.min(np.array(I.coordinates), axis=0)
x_max, y_max = np.max(np.array(I.coordinates), axis=0)

# Construct the grid-like network
G = algo.grid_like_graph_generator([x_min, x_max], [y_min, y_max], grid_size, I)

# Visualize the network and output it in the dir
utils.visualize_graph_nx(G, dir)

# Rename the customers and the depot
customer = [(I.coordinates[idx][0], I.coordinates[idx][1], 'c'+str(idx)) for idx in range(1, len(I.coordinates))]
depot = [(I.coordinates[0][0], I.coordinates[0][1], 'd')]

# Assign edge cost for ordinary roads and AV-enabled roads
utils.assign_cost_n_travel_time(G, discount_factor_av, inflated_factor_av)

# Compute the distance matrices 
print("Computing the distance matrix and the path matrix...")
cost_av, cost_non_av, path_av, path_non_av = utils.dist_n_path_matrices_4_experiments(G, I)


"""
Solve the CVRP
"""

print("Start to solve the CVRP instance with av_cost.")
sys.stdout.flush()

# TODO: Build another HCVRP solver to replace this one
hgs_LB = algo.UHGS_CVRPSolver(instance_path, dir, G, I, t = 5)
hgs_LB.add_explicit_dist_mtx(cost_av)
hgs_LB.solve('LB.sol')
LB, routes, all_time_stamps, clusters = hgs_LB.extract_results(path_av)

print("Start to solve the CVRP instance with non_av_cost.")
sys.stdout.flush()

hgs_UB = algo.UHGS_CVRPSolver(instance_path, dir, G, I, t = 5)
hgs_UB.add_explicit_dist_mtx(cost_non_av)
hgs_UB.solve('UB.sol')
UB, _, all_time_stamps_UB, _ = hgs_UB.extract_results(path_non_av)


# Extract the latest return time of the fleet
t_max = max([timestamp[-1] for schedule in all_time_stamps_UB.values() for timestamp in schedule])

#* DEBUG:
# t_max_tmp = max([timestamp[-1] for schedule in all_time_stamps.values() for timestamp in schedule])
# print("equal: ", all_time_stamps == all_time_stamps_UB)
# print("UB", all_time_stamps_UB)
# print("LB", all_time_stamps)
# print("UB T_max: " + str(t_max))
# print("LB T_max: " + str(t_max_tmp))



"""
Solve the re-scheduling FRP
"""

print("Start to solve the re-scheduling MILP.")
status = algo.frp_reschedule(G, routes, all_time_stamps, t_max * t_max_factor, budget, dir, gamma)

# Visualize the re-scheduling if it is feasible
if status == 1:
    utils.visualize_re_scheduling(dir)

    
    
"""
Solve the re-routing FRP
"""

print("Start to solve the re-routing FRP.")
re_routing_solver = algo.ReRoutingFRPSolver(G, I, routes, all_time_stamps, clusters, t_max * t_max_factor, budget, gamma, num_layer, dir)
cost, new_routes, new_all_time_stamps, D_bar, G_adjusted = re_routing_solver.solve()

# Visualize 
utils.visualize_re_routing(dir, G_adjusted, new_all_time_stamps)



"""
Replace unserved AV routes with HDV routes
"""

print("Dispatch HDVs to serve the unserved customers if necessary.")
sys.stdout.flush()

# Serve the unserved customers by HDVs 
if D_bar != []: 
    
    print('Unserved: ', D_bar)
    
    path = utils.write_vrp_instance(dir, I, D_bar, cost_non_av)
    hgs_unserved = algo.UHGS_CVRPSolver(path, dir, G, I, t = 5)
    hgs_unserved.tmp_path = path
    hgs_unserved.solve('HDV_replacement.sol')
    replaced_cost, _, _, _ = hgs_unserved.extract_results(path_non_av)    
    cost += replaced_cost
    

"""
Output the results
"""

print(colored("The LB, cost, and UB are: ", 'cyan'))
print(colored(LB, 'cyan'), colored(cost, 'cyan'), colored(UB, 'cyan'))

result_path = dir.split('/')[0] + '/result-'+ str(args.t_max_factor) + '-' + str(args.budget_factor) +'.csv'
with open (result_path, 'a', newline = '') as file:
    writer = csv.writer(file)
    writer.writerow([instance, LB, status, cost, UB])


#*: DEBUG
# tmp_cost = 0
# tmp_route = routes[2]
# for i, j in zip(tmp_route[:-1], tmp_route[1:]):
        # tmp_cost += G[i][j]['av_cost']
# print('tmp_cost: ' + str(tmp_cost))
# print('replaced_cost: ' + str(replaced_cost))


