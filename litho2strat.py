'''
  This program generates a set of plausible stratigraphies with uncertainties, for a given drillhole lithology log.
  It uses map data for distance and topology constraints, and several free parameters describing the solution complexity (level of deformation) constraints.

  Author: Vitaliy Ogarko, vogarko@gmail.com
  The University of Western Australia
'''

import os

from strata_solver import \
    generate_strat_routes, \
    generate_strata_table, \
    StrataSolverParameters

from data_readers import \
    read_strat_data, \
    read_drillsample_data, \
    read_thickness_data, \
    read_topology_data, \
    read_ignore_list, \
    read_alternative_rock_names, \
    DrillSampleDataHeader, \
    StrataDataHeader

from solution_utils import *

spar = StrataSolverParameters()

#==============================================================================
# Unit contact topology constraints (extracted from map data).
spar.add_topology_constraints = True
# 'Age alignment' constraints: the maximum number of times the age direction can flip.
spar.max_num_age_flips = 2
#---------------------------------------------------------------------------
# The number of nearest units (for distance constraints).
number_nearest_units = 3
#---------------------------------------------------------------------------
# Minimum score for drillhole lithologies to use them.
min_drillhole_litho_score = 80
#---------------------------------------------------------------------------
# Group drillhole lithology sequence.
# Note: use this for max_num_unit_contacts_inside_litho > 0 to avoid the solution number to blow.
group_drillhole_lithos = False
#---------------------------------------------------------------------------
# The cover ration threshold (relative length) for removing the cover.
cover_ratio_threshold = 0.65
#---------------------------------------------------------------------------

# 'Returning to the same unit' constraints.
spar.max_num_returns_per_unit = 1
#---------------------------------------------------------------------------
# The maximum number of unit contacts inside the same litholgy sequence.
spar.max_num_unit_contacts_inside_litho = 0
#---------------------------------------------------------------------------
# Use the single closest unit for the top (first) lithology.
spar.single_top_unit = True

#---------------------------------------------------------------------------
# Adding thickness constraints. (Requires unit thickness data).
add_thickness_constraints = False

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
    alternative_rock_names_filename = "data/real/alternative_rock_names.txt"

    directory = "data/real/dist_files/litho_tables"
    #directory = "data/real/dh_files/litho_tables"

    # Drillsample data column names.
    drillsample_header = DrillSampleDataHeader('Fromdepth', 'Todepth', 'Lithologies', 'Scores')

    # Strata data csv file column names.
    strata_data_header = StrataDataHeader('UNITNAME', 'lithos', 'distance', 'DESCRIPTN')

    # Read the drillhole ignore items list.
    ignore_list = read_ignore_list(ignore_list_filename)

    # Read the alternative rock names.
    alternative_rock_names = read_alternative_rock_names(alternative_rock_names_filename)

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
        drillsample_data = read_drillsample_data(drillsample_header, drillsample_filename, ignore_list, min_drillhole_litho_score)

        # Unit lithologies and distance data.
        strata_data = read_strat_data(strata_data_header, dist_table_filename, alternative_rock_names)

        # Filter strat data.
        drillhole_lithos = drillsample_data.get_drillhole_lithos()
        strata_data.filter_strat_data_based_on_drillhole_lithos(drillhole_lithos)
        strata_data.filter_strat_data_based_on_distance(number_nearest_units)

        # Generating the table of possible strata paths.
        strata_table, missing_lithos = generate_strata_table(drillsample_data, strata_data, spar.single_top_unit)

        all_missing_lithos.update(missing_lithos)

    # Reading the cover lithos list.
    with open(cover_unit_filename) as f:
        cover_lithos = f.read().splitlines()

    # Remove cover lithos from the missing lithos list.
    all_missing_lithos_list = list(all_missing_lithos)
    for litho in all_missing_lithos_list:
        if (litho in cover_lithos):
            all_missing_lithos.remove(litho)

    print("Missing lithologies list:")
    print("Number lithos missing =", len(all_missing_lithos))
    print(all_missing_lithos)

#=============================================================================
def main():
    print('Started litho2strat')

    #generate_missing_lithos()
    #exit()

    # Topology file.
    topology_filename = "data/real/ASUD_strat4.gml"

    # The Cover unit data file.
    cover_unit_filename = "data/real/cover_unit.txt"

    # The Ignore items list file.
    ignore_list_filename = "data/real/ignore_list.txt"

    # Alternative rock names file.
    alternative_rock_names_filename = "data/real/alternative_rock_names.txt"

    # Unit colours for drawing stratigraphy logs.
    unit_colors_filename = "data/real/500kibg_colours.csv"

    # Drillsample data csv file column names.
    drillsample_header = DrillSampleDataHeader('Fromdepth', 'Todepth', 'Lithologies', 'Scores')

    # Strata data csv file column names.
    strata_data_header = StrataDataHeader('UNITNAME', 'lithos', 'distance', 'DESCRIPTN')

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
    #collarID = 2182301

    # Looks good using (number_nearest_units = 2, max_num_returns_per_unit = 2, add_topology_constraints = True, single_top_unit = True, max_num_unit_contacts_inside_litho = 0)
    #collarID = 2182076

    # (!) Has wacke at row = 28 which is not passing topology constraints! To discuss with Mark.
    #collarID = 810340

    # Read topology data (the same for all collarIDs).
    graph = None
    if (spar.add_topology_constraints):
        graph = read_topology_data(topology_filename)

    #====================================================================================
    # Process the list of CollarIDs.
    #====================================================================================
    # First small cluster of 6.
    #collarIDs = [2182334, 2182335, 2182336, 2182338, 2182339, 2182340]

    # Second cluster of 14 (crossing the boundary of two units on the map).
    #collarIDs = [2182301, 2182306, 2182307, 2182308, 2182309, 2182310]
    #collarIDs.extend(list(range(2182312, 2182319 + 1)))

    # Third cluster of 18 (in the area with several units touching).
    collarIDs = [2182009, 2182013, 2182016, 2182029, 2182047, 2470197, 2470200, 2470301, 2470303, 2470304, 2470305]
    # Empty ones: 2182010, 2182035, 2470196
    # No solutions found: 2182017, 2182018
    # Solution number blows: 2182030
    # No data files: 2470193 

    # No correlation example.
    #collarIDs = [2470303, 2182308]

    collarIDs = [2470303]

    # Strong correlation.
    collarIDs = [2470303, 2182029]

    display_plots = True

    # Stores solutions for different drillholes.
    strat_solutions = []

    for collarID in collarIDs:

        print('collarID =', collarID)

        drillsample_filename = "data/real/dist_files/litho_tables_V3/litho_" + str(collarID) + ".csv"
        dist_table_filename = "data/real/dist_files/dist_tables/100_500k_map_near_" + str(collarID) + ".csv"

        # Synthetic test.
        # Note: The tests #1 and #2 show very different probabilities for max_num_unit_contacts_inside_litho = 0 and 1, i.e., a constant and linear increasing transition.
        # IMPORTANT: For these tests, set max_num_returns_per_unit = 0.
        #drillsample_filename = "data/tests/litho_2.csv"
        #dist_table_filename = "data/tests/map_2.csv"

        #--------------------------------------------------------------
        # Reading the input data.
        #--------------------------------------------------------------
        # Read the drillhole ignore items list.
        ignore_list = read_ignore_list(ignore_list_filename)

        # Read drill sample data.
        drillsample_data = read_drillsample_data(drillsample_header, drillsample_filename, ignore_list, min_drillhole_litho_score)

        # Remove the Cover.
        drillsample_data.remove_cover(cover_unit_filename, cover_ratio_threshold)

        if (group_drillhole_lithos):
            # Group the drillsample lithologies.
            drillsample_data.group_drillhole_litho_sequence(spar.max_num_unit_contacts_inside_litho)

        # Read the alternative rock names.
        alternative_rock_names = read_alternative_rock_names(alternative_rock_names_filename)

        # Read unit lithologies and distance data.
        strata_data = read_strat_data(strata_data_header, dist_table_filename, alternative_rock_names)

        # Filter strat data.
        drillhole_lithos = drillsample_data.get_drillhole_lithos()
        strata_data.filter_strat_data_based_on_drillhole_lithos(drillhole_lithos)
        strata_data.filter_strat_data_based_on_distance(number_nearest_units)

        filtered_lithos = strata_data.get_unique_lithos()
        print("The number of filtered unit lithologies:", len(filtered_lithos))
        print("Filtered unit lithologies:", sorted(filtered_lithos))

        # Read thickness data.
        thickness_data = []
        if (add_thickness_constraints):
            thickness_data = read_thickness_data(thickness_filename)

        #--------------------------------------------------------------
        # Generating the stratigraphies.
        #--------------------------------------------------------------
        strat_solution = generate_strat_routes(spar, strata_data, drillsample_data, thickness_data, graph)
        strat_solution.collarID = collarID

        print("Total number of routes = ", len(strat_solution.routes))

        strat_solutions.append(strat_solution)

        #--------------------------------------------------------------
        # Plot the results.
        #--------------------------------------------------------------

        # Print all unique routes (i.e., unique strata sequence).
        print_unique_routes(strat_solution.routes, 10)

        # Draw stratigraphy logs.
        draw_solution_logs(strat_solution, display_plots, 'strat', unit_colors_filename, True, None)

        # Draw probability logs.
        #draw_solution_logs(strat_solution, display_plots, 'proba', '', True, None)

        # Draw the age-rule logs.
        #draw_solution_logs(strat_solution, display_plots, 'age', '', True, graph)

        # Plot unit probabilities.
        #plot_unit_probabilities(strat_solution, display_plots)

        # Draw the topology graph of all solution routes.
        draw_solution_graph(strat_solution)

        # Plot histogram of solution scores.
        plot_route_scores(strat_solution.graph_route_scores)

        # Write the best routes to file.
        write_best_routes_to_file(strat_solution, 10)

    #-------------------------------------------------------------------------
    # Analyze solution correletion between different drillholes.
    #-------------------------------------------------------------------------

    # Loop over all solutipon pairs.
    for i in range(len(strat_solutions)):
        graph = strat_solutions[i].graph
        for j in range(len(strat_solutions)):
            if (i != j):
                # Calculate solution scores based on external graph.
                graph_route_scores = strat_solutions[j].calculate_graph_route_scores(graph)
                strat_solutions[j].external_graph_route_scores_list.append(graph_route_scores)

    for solution in strat_solutions:
        plot_solution_correlation(solution)

    #-------------------------------------------------------------------------
    # Calculate a new route score equal to the sum of graph scores ovel all drillholes.
    #-------------------------------------------------------------------------
    for solution in strat_solutions:
        solution.route_scores = solution.graph_route_scores
        for external_graph_route_scores in solution.external_graph_route_scores_list:
            solution.route_scores = solution.route_scores + external_graph_route_scores

    #-------------------------------------------------------------------------
    # Show the most correlated solution logs.
    for solution in strat_solutions:
        draw_solution_logs(solution, display_plots, 'strat', unit_colors_filename, False, None)

#=============================================================================
if __name__ == "__main__":
    main()

