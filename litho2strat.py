'''
  This program generates a set of plausible stratigraphies with uncertainties, for a given drillhole lithology log.
  It uses map data for distance and topology constraints, and several free parameters describing the solution complexity (level of deformation) constraints.

  Author: Vitaliy Ogarko, vogarko@gmail.com
  The University of Western Australia
'''

import matplotlib.pylab as pl
import os
#import tracemalloc

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

from solution_utils import \
    print_unique_routes, \
    plot_route_scores, \
    plot_top_routes, \
    plot_unit_probabilities, \
    write_best_routes_to_file

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
min_drillhole_litho_score = 80
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
    drillsample_header = DrillSampleDataHeader('Fromdepth', 'Todepth', 'Lithologies', 'Scores')

    # Strata data csv file column names.
    strata_data_header = StrataDataHeader('UNITNAME', 'lithos', 'distance', 'DESCRIPTN')

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
        drillsample_data = read_drillsample_data(drillsample_header, drillsample_filename, ignore_list, min_drillhole_litho_score)

        # Unit lithologies and distance data.
        strata_data = read_strat_data(strata_data_header, dist_table_filename, alternative_rock_names)

        # Filter strat data.
        drillhole_lithos = drillsample_data.get_drillhole_lithos()
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
    topology_filename = "data/real/ASUD_strat3.gml"

    # The Cover unit data file.
    cover_unit_filename = "data/real/cover_unit.txt"

    # The Ignore items list file.
    ignore_list_filename = "data/real/ignore_list.txt"

    # Alternative rock names file.
    alternative_rock_names_file = "data/real/alternative_rock_names.txt"

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
    collarID = 2182336
    #collarID = 2182335
    #collarID = 2182340
    #collarID = 2182339
    #collarID = 2182338
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

    # Drillsample data csv file column names.
    drillsample_header = DrillSampleDataHeader('Fromdepth', 'Todepth', 'Lithologies', 'Scores')

    # Strata data csv file column names.
    strata_data_header = StrataDataHeader('UNITNAME', 'lithos', 'distance', 'DESCRIPTN')

    #--------------------------------------------------------------
    # Reading the input data.
    #--------------------------------------------------------------
    # Read the drillhole ignore items list.
    ignore_list = read_ignore_list(ignore_list_filename)

    # Read drill sample data.
    drillsample_data = read_drillsample_data(drillsample_header, drillsample_filename, ignore_list, min_drillhole_litho_score)

    if (group_drillhole_lithos):
        # Group the drillsample lithologies.
        drillsample_data.group_drillhole_litho_sequence(spar.max_num_unit_contacts_inside_litho)

    # Read the alternative rock names.
    alternative_rock_names = read_alternative_rock_names(alternative_rock_names_file)

    # Read unit lithologies and distance data.
    strata_data = read_strat_data(strata_data_header, dist_table_filename, alternative_rock_names)

    # Filter strat data.
    drillhole_lithos = drillsample_data.get_drillhole_lithos()
    strata_data.filter_strat_data_based_on_drillhole_lithos(drillhole_lithos)
    strata_data.filter_strat_data_based_on_distance(number_nearest_units)

    filtered_lithos = strata_data.get_unique_lithos()
    print("The number of filtered unit lithologies:", len(filtered_lithos))
    print("Filtered unit lithologies:", sorted(filtered_lithos))

    # Read the Cover unit lithologies.
    strata_data.add_cover_unit("Cover", cover_unit_filename)

    # Read thickness data.
    thickness_data = []
    if (add_thickness_constraints):
        thickness_data = read_thickness_data(thickness_filename)

    # Read topology data.
    graph = None
    if (add_topology_constraints):
        graph = read_topology_data(topology_filename, ignore_unit_age)
        # Sanity check: check that strata units exist in the graph.
        unit_names = strata_data.get_unit_names()
        for unit_name in unit_names:
            if unit_name not in graph.nodes():
                print("WARNING: Not found graph unit: ", unit_name)

    #--------------------------------------------------------------
    # Generating the stratigraphies.
    #--------------------------------------------------------------
#    tracemalloc.start()

    strat_solution = generate_strat_routes(spar, strata_data, drillsample_data, thickness_data, graph)

    print("Total number of routes = ", len(strat_solution.routes))

#    current, peak = tracemalloc.get_traced_memory()
#    print("Current memory usage is {} MB; Peak was {} MB".format(current / 10**6, peak / 10**6))

    #--------------------------------------------------------------
    # Plot the results.
    #--------------------------------------------------------------
    # Plot the number of processed routes at each row.
    pl.xlabel('Row number')
    pl.ylabel('Number of routes')
    pl.plot(strat_solution.routes_number)
    pl.show()

    # Print all unique routes (i.e., unique strata sequence).
    print_unique_routes(strat_solution.routes, 10)

    # Plot route scores.
    plot_route_scores(strat_solution)

    # Plot top routes and their probabilities.
    plot_top_routes(strat_solution)

    # Plot unit probabilities.
    plot_unit_probabilities(strat_solution)

    # Write the best routes to file.
    filename = "best_routes_" + str(collarID) + ".txt"
    write_best_routes_to_file(strat_solution, filename, 10)

#=============================================================================
if __name__ == "__main__":
    main()

