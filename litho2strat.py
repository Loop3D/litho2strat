'''
  This program generates a set of plausible stratigraphies with uncertainties, for a given drillhole lithology log.
  It uses map data for distance and topology constraints, and several free parameters describing the solution complexity (level of deformation) constraints.

  Author: Vitaliy Ogarko, vogarko@gmail.com
  The University of Western Australia
'''

import numpy as np
import matplotlib.pylab as pl
import networkx as nx
import os
#import tracemalloc

from strata_solver import \
    generate_strat_routes, \
    generate_strata_table, \
    group_drillhole_litho_sequence, \
    StrataSolverParameters

from data_readers import \
    read_strat_data, \
    read_drillsample_data, \
    read_thickness_data, \
    read_topology_data, \
    read_ignore_list, \
    read_alternative_rock_names, \
    DrillSampleHeader

#==============================================================================
# Adding unit contacts topology (extracted from map data).
add_topology_constraints = True
# Ignore topology graph edge direction defining the unit age.
ignore_unit_age = True
#---------------------------------------------------------------------------
# The number of nearest units (for distance constraints).
number_nearest_units = 3
#---------------------------------------------------------------------------
# Minimum score for drillhole lithologies to use them.
min_drillhole_litho_score = 70
#---------------------------------------------------------------------------
# Group drillhole lithology sequence.
# Note: use this for max_num_unit_contacts_inside_litho > 0 to avoid the solution number to blow.
group_drillhole_lithos = False
#---------------------------------------------------------------------------
spar = StrataSolverParameters()

# 'Returning to the same unit' constraints.
spar.max_num_returns_per_unit = 2
#---------------------------------------------------------------------------
# The number of unit contacts inside the same litholgy sequence.
spar.max_num_unit_contacts_inside_litho = 0
#---------------------------------------------------------------------------
# Use the single closest unit for the top (first) lithology.
spar.single_top_unit = True

#---------------------------------------------------------------------------
# Adding thickness constraints. (Requires unit thickness data).
add_thickness_constraints = False

#==============================================================================
def write_routes_to_file(filename, drillsample_data, all_routes):
    '''
    Writing stratigraphic routes to file.
    '''
    f = open(filename, "w")
    num_rows = len(drillsample_data)
    for row in range(num_rows):
        depth = drillsample_data[row].depth_from
        f.write("%f " % depth)
        # Calculate the number of unique strata for this depth.
        unique_units = set([])
        for route in all_routes:
            unique_units.add(route.path[row])
        f.write("%d " % len(unique_units))

        for route in all_routes:
            f.write("%d " % route.path[row])
        f.write("\n")
    f.close()

#==============================================================================
def print_unique_routes(all_routes, num_print_paths):
    '''
    Print all unique routes (i.e., with unique strata sequence).
    '''
    unique_routes = set([])
    for route in all_routes:
        unique_routes.add(route.get_strata_sequence())

    print("Number of unique routes = ", len(unique_routes))
    if (num_print_paths > 0):
        num = 0
        for route in unique_routes:
            num += 1
            print(route)
            if (num >= num_print_paths):
                break

#=============================================================================
def plot_routes(drillsample_data, routes, strat_distr):
    '''
    Plot and display the routes.
    '''
    print("Plotting the routes...")

    # Using the "From" column.
    x_data = [d.depth_from for d in drillsample_data]

    for route in routes:
        pl.plot(x_data, route.path, '.-')

    pl.xlabel('Depth')
    pl.ylabel('Strata unit index')
    pl.show()

    #------------------------------------------
    # Plot route probabilities.
    num_rows = len(routes[0].path)

    for route in routes:
        route_proba = np.zeros(num_rows)
        for row in range(num_rows):
            unit_index = route.path[row]
            route_proba[row] = strat_distr[row, unit_index]
        pl.plot(x_data, route_proba, '.-')

    pl.xlabel('Depth')
    pl.ylabel('Probability')
    pl.show()

#=============================================================================
def get_strat_distr(all_routes, num_units):
    '''
    Returns the distribution of unit presence at every depth.
    '''
    num_rows = len(all_routes[0].path)
    strat_distr = np.zeros((num_rows, num_units))

    for route in all_routes:
        for row, unit_index in enumerate(route.path):
            strat_distr[row, unit_index] += 1

    # Normalize.
    strat_distr = strat_distr / float(len(all_routes))
    return strat_distr

#=============================================================================
def get_route_scores(all_routes, strat_distr):
    '''
    Returns the route scores (based on path probability).
    Needs strat_distr returned by get_strat_distr().
    '''
    num_rows = len(all_routes[0].path)
    route_scores = np.zeros(len(all_routes), dtype=float)

    for route_index, route in enumerate(all_routes):
        for row, unit_index in enumerate(route.path):
            route_scores[route_index] += strat_distr[row, unit_index]

    # Normalize.
    route_scores = route_scores / float(num_rows)
    return route_scores

#=============================================================================
def plot_unit_probabilities(all_routes, drillsample_data, unit_names):
    '''
    Generate a plot with probability of occurence for each unit.
    '''
    if (len(all_routes) == 0):
        return

    num_units = len(unit_names)

    # Building the distribution of unit presence at every depth.
    strat_distr = get_strat_distr(all_routes, num_units)

    #------------------------------------------
    # Plot distribution of the route scores (based on path probability).
    route_scores = get_route_scores(all_routes, strat_distr)

    pl.hist(route_scores, bins = 50)
    pl.xlabel('Route score')
    pl.ylabel('Frequency')
    pl.show()

    #------------------------------------------
    # Top scores.
    indexes_max = np.argsort(-route_scores) # A minus here to have largest to smallest score order.
    ntop = 10
    print('Top indexes: ', indexes_max[0:ntop])
    print('Top scores: ', route_scores[indexes_max[0:ntop]])

    index_max = indexes_max[0]

    #------------------------------------------
    # Print the most probable routes.
    plot_routes(drillsample_data, [all_routes[i] for i in indexes_max[0:ntop]], strat_distr)

    #------------------------------------------
    # Print if there are multiple best routes.
    multiple_best_routes = False
    if (len(indexes_max) > 1):
        if (route_scores[indexes_max[0]] == route_scores[indexes_max[1]]):
            multiple_best_routes = True

    print("Multiple best routes: ", multiple_best_routes)

    #------------------------------------------
    # Generating the plots.
    # Increasing the figure size.
    pl.rcParams["figure.figsize"] = (12.8, 9.6) # Default size = (6.4, 4.8)

    num_units_nonempty = 0
    for i in range(num_units):
        if (sum(strat_distr[:, i]) != 0):
            num_units_nonempty += 1
    fig, axs = pl.subplots(num_units_nonempty, sharey=True, squeeze=True)

    fig.suptitle('Probability of occurrence for every unit.', y=0.96)

    num_rows = len(all_routes[0].path)

    #-------------------------------------------------------------
    # Adding the "From" and "To" depths for visualisation.
    x_data = []
    for d in drillsample_data[0:num_rows]:
        x_data.append(d.depth_from)
        x_data.append(d.depth_to)

    # Duplicate each value, as the probability is the same between the "From" and "To" depths.
    strat_distr = np.repeat(strat_distr, 2, axis=0)
    #-------------------------------------------------------------

    # Skip empty units.
    nonempty_units = []
    for i in range(num_units):
        if (sum(strat_distr[:, i]) != 0):
            nonempty_units.append(i)

    j = 0
    for index, i in enumerate(nonempty_units):
        # Plot lines.
        axs[j].plot(x_data, strat_distr[:, i], zorder=1, c='blue')

        # Set red color for zero data.
        color = ['red' if p <= 0 else 'blue' for p in strat_distr[:, i]]

        # Plot dots.
        axs[j].scatter(x_data, strat_distr[:, i], s=5, c=color, zorder=2)

        axs[j].set_title(unit_names[i], size=9, y=0.97)
        axs[j].set_ylabel(str(j))

        if (index != len(nonempty_units) - 1):
            # Hide tick labels.
            axs[j].set_xticklabels([])

        # Add vertical lines.
        axs[j].xaxis.grid(True)

        j += 1
 
    #pl.tight_layout()
    pl.subplots_adjust(hspace = 0.5)
 
    pl.xlabel('Depth')
    pl.show()

#==============================================================================
def get_drillhole_lithos(drillsample_data):
    '''
    Returns the drillhole lithologies from drillsample data.
    '''
    all_lithos = set()
    for row in drillsample_data:
        lithos = row.lithos
        all_lithos.update(lithos)
    return all_lithos

#=============================================================================
def generate_missing_lithos():
    '''
    Generates a list of missing drillhole lithologies from unit data in all data files.
    '''

    # The Cover unit data file.
    cover_unit_filename = "data/real/cover_unit.txt"

    # The Ignore items list.
    ignore_list_filename = "data/real/ignore_list.txt"

    # Alternative rock names file.
    alternative_rock_names_file = "data/real/alternative_rock_names.txt"

    directory = "data/real/dist_files/litho_tables"
    #directory = "data/real/dh_files/litho_tables"

    # Drillsample data column names.
    drillsample_header = DrillSampleHeader('Fromdepth', 'Todepth', 'Lithologies', 'Scores')

    # Read the drillhole ignore items list.
    ignore_list = read_ignore_list(ignore_list_filename)

    # Read the alternative rock names.
    alternative_rock_names = read_alternative_rock_names(alternative_rock_names_file)

    # All missing lithologies.
    all_missing_lithos = set()

    counter = 0
    for file in os.listdir(directory):
        counter = counter + 1
        print("PROCESSING FILE NUMBER =", counter)

        filename = os.fsdecode(file)
        collarID = int(filename[6:-4])
        print(collarID)

        drillsample_filename = "data/real/dist_files/litho_tables_V2/litho_" + str(collarID) + ".csv"
        dist_table_filename = "data/real/dist_files/dist_tables/100_500k_map_near_" + str(collarID) + ".csv"

        #drillsample_filename = "data/real/dh_files/litho_tables/litho_" + str(collarID) + ".csv"
        #dist_table_filename = "data/real/dh_files/dist_tables/100k_map_near_" + str(collarID) + ".csv"

        # Drill sample data.
        drillsample_data = read_drillsample_data(drillsample_filename, drillsample_header, ignore_list, min_drillhole_litho_score)

        # Unit lithologies and distance data.
        strata_data = read_strat_data(dist_table_filename, alternative_rock_names)

        # Filter strat data.
        drillhole_lithos = get_drillhole_lithos(drillsample_data)
        strata_data.filter_strat_data_based_on_drillhole_lithos(drillhole_lithos)
        strata_data.filter_strat_data_based_on_distance(number_nearest_units)

        # Read the Cover unit lithologies.
        strata_data.add_cover_unit("Cover", cover_unit_filename)

        # Generating the table of possible strata paths.
        strata_table, missing_lithos = generate_strata_table(drillsample_data, strata_data, spar.single_top_unit)

        all_missing_lithos.update(missing_lithos)

    print("Missing lithologies list:")
    print("Number lithos missing =", len(all_missing_lithos))
    print(all_missing_lithos)

#=============================================================================
def main():
    print('Started litho2strat')

    #generate_missing_lithos()
    #exit()

    # Topology file.
    #topology_filename = "data/real/ASUD_strat.gml"
    topology_filename = "data/real/ASUD_strat2.gml"

    # The Cover unit data file.
    cover_unit_filename = "data/real/cover_unit.txt"

    # The Ignore items list file.
    ignore_list_filename = "data/real/ignore_list.txt"

    # Alternative rock names file.
    alternative_rock_names_file = "data/real/alternative_rock_names.txt"

    #----------------------------------------------------------------------------
    # # Mark's data.
    # collarID = 548917
    # drillsample_filename = "data/real/dh_files/litho_tables/litho_" + str(collarID) + ".csv"
    # dist_table_filename = "data/real/dh_files/dist_tables/100k_map_near_" + str(collarID) + ".csv"
    #
    # # Drillsample data column names.
    # drillsample_header = DrillSampleHeader('FromDepth', 'ToDepth', 'CET_Litho', 'Scores')

    #----------------------------------------------------------------------------
    # Mark's data with known solutions.
    #----------------------------------------------------------------------------

    # (!) Cannot pass topo constraints for row = 5 for mudstone-granite contact when both sanstone and mudstone are excluded from Cover.
    #collarID = 1209855

    # (!) Gravel issue for CollarId=1209857, at depth=244m
    #collarID = 1209857

    # Looks good using (number_nearest_units = 2, max_num_returns_per_unit = 2, add_topology_constraints = True, single_top_unit = True, max_num_unit_contacts_inside_litho = 0)
    #collarID = 353386

    # Looks good using (number_nearest_units = 2, max_num_returns_per_unit = 2, add_topology_constraints = True, single_top_unit = True, max_num_unit_contacts_inside_litho = 0)
    collarID = 2182301

    # Looks good using (number_nearest_units = 2, max_num_returns_per_unit = 2, add_topology_constraints = True, single_top_unit = True, max_num_unit_contacts_inside_litho = 0)
    #collarID = 2182076

    # (!) Has wacke at row = 28 which is not passing topology constraints! To discuss with Mark.
    #collarID = 810340

    # Confirmed results (using 1 closest unit & single top unit).
    #collarID = 2182336
    #collarID = 2182335
    #collarID = 2182340
    #collarID = 2182339
    #collarID = 2182338
    #collarID = 2182335
    #collarID = 2182334

    # (!) Strange gravel at 10m, which looks like real gravel, but there are rocks above...
    #collarID = 2470303
    # (!) Strange gravel at 4m, which looks like real gravel, but there are rocks above...
    #collarID = 2470304

    print('collarID =', collarID)

    drillsample_filename = "data/real/dist_files/litho_tables_V3/litho_" + str(collarID) + ".csv"
    dist_table_filename = "data/real/dist_files/dist_tables/100_500k_map_near_" + str(collarID) + ".csv"

    # Synthetic test.
    # Note: The tests #1 and #2 show very different probabilities for max_num_unit_contacts_inside_litho = 0 and 1, i.e., a constant and linear increasing transition.
    # IMPORTANT: For these tests, set max_num_returns_per_unit = 0.
    #drillsample_filename = "data/tests/litho_2.csv"
    #dist_table_filename = "data/tests/map_2.csv"

    # Drillsample data column names.
    drillsample_header = DrillSampleHeader('Fromdepth', 'Todepth', 'Lithologies', 'Scores')

    #--------------------------------------------------------------
    # Reading the input data.
    #--------------------------------------------------------------
    # Read the drillhole ignore items list.
    ignore_list = read_ignore_list(ignore_list_filename)

    # Read drill sample data.
    drillsample_data = read_drillsample_data(drillsample_filename, drillsample_header, ignore_list, min_drillhole_litho_score)

    if (group_drillhole_lithos):
        # Group the drillsample lithologies.
        drillsample_data = group_drillhole_litho_sequence(drillsample_data, spar.max_num_unit_contacts_inside_litho)

    # Read the alternative rock names.
    alternative_rock_names = read_alternative_rock_names(alternative_rock_names_file)

    # Read unit lithologies and distance data.
    strata_data = read_strat_data(dist_table_filename, alternative_rock_names)

    # Filter strat data.
    drillhole_lithos = get_drillhole_lithos(drillsample_data)
    strata_data.filter_strat_data_based_on_drillhole_lithos(drillhole_lithos)
    strata_data.filter_strat_data_based_on_distance(number_nearest_units)

    print("Filtered unit lithologies:", strata_data.get_unique_lithos())

    # Read the Cover unit lithologies.
    strata_data.add_cover_unit("Cover", cover_unit_filename)

    # Read thickness data.
    thickness_data = []
    if (add_thickness_constraints):
        thickness_data = read_thickness_data(thickness_filename)

    # Read topology data.
    graph = nx.Graph()
    if (add_topology_constraints):
        graph = read_topology_data(topology_filename, ignore_unit_age)
        # Sanity check: check that strata units exist in the graph.
        for unit_name in strata_data.unit2litho:
            if unit_name not in graph.nodes():
                print("WARNING: Not found graph unit: ", unit_name)

    #--------------------------------------------------------------
    # Generating the stratigraphies.
    #--------------------------------------------------------------
#    tracemalloc.start()

    all_routes, all_routes_number = generate_strat_routes(spar, strata_data, drillsample_data, thickness_data, graph)

    print("Total number of routes = ", len(all_routes))

#    current, peak = tracemalloc.get_traced_memory()
#    print("Current memory usage is {} MB; Peak was {} MB".format(current / 10**6, peak / 10**6))

    #--------------------------------------------------------------
    # Plot the results.
    #--------------------------------------------------------------
    # Plot the number of processed routes at each row.
    pl.xlabel('Row number')
    pl.ylabel('Number of routes')
    pl.plot(all_routes_number)
    pl.show()

    # Print all unique routes (i.e., unique strata sequence).
    print_unique_routes(all_routes, 10)

    # Write results to the file.
    #write_routes_to_file("strata.txt", drillsample_data, all_routes)

    # Plot unit probabilities.
    unit_names = strata_data.get_unit_names()
    plot_unit_probabilities(all_routes, drillsample_data, unit_names)

#=============================================================================
if __name__ == "__main__":
    main()

