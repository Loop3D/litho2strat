'''
  This program generates a set of plausible stratigraphies with uncertainties, for a given drillhole lithology log.
  It uses map data for distance and topology constraints, and several free parameters describing the solution complexity (level of deformation) constraints.

  Author: Vitaliy Ogarko, vogarko@gmail.com
  The University of Western Australia
'''

import os
import configparser

from strata_solver import \
    generate_strat_routes, \
    StrataSolverParameters
import data_readers as dr
from solution_utils import *
from solution_analysis import *

#=============================================================================
def main():
    print('Started litho2strat')

    #-----------------------------------------------------------------
    # Read input parameters.
    #-----------------------------------------------------------------
    config = configparser.ConfigParser()
    config.read('Parfile.txt')

    section = 'FilePaths'
    print(config.items(section))

    # Topology file.
    topology_filename = config.get(section, 'topology_filename')

    # The Cover unit data file.
    cover_unit_filename = config.get(section, 'cover_unit_filename')

    # The ignore items list file.
    ignore_list_filename = config.get(section, 'ignore_list_filename')

    # Alternative rock names file.
    alternative_rock_names_filename = config.get(section, 'alternative_rock_names_filename')

    # Unit colours for drawing stratigraphy logs.
    unit_colors_filename = config.get(section, 'unit_colors_filename')

    # Drillhole lithology data file. The $collarID$ in the file name will be replaced with the actual value.
    drillsample_filename_collarID = config.get(section, 'drillsample_filename')

    # Units near the collar with distances data file. The $collarID$ in the file name will be replaced with the actual value.
    dist_table_filename_collarID = config.get(section, 'dist_table_filename')

    #---------------------------------
    section = 'SolverParameters'
    print(config.items(section))

    spar = StrataSolverParameters()

    # Unit contact topology constraints (extracted from map data).
    spar.add_topology_constraints = config.getboolean(section, 'add_topology_constraints')

    # 'Age alignment' constraints: the maximum number of times the age direction can flip.
    spar.max_num_age_flips = config.getint(section, 'max_num_age_flips')

    # 'Returning to the same unit' constraints.
    spar.max_num_returns_per_unit = config.getint(section, 'max_num_returns_per_unit')

    # The maximum number of unit contacts inside the same litholgy sequence.
    spar.max_num_unit_contacts_inside_litho = config.getint(section, 'max_num_unit_contacts_inside_litho')

    # Use the single closest unit for the top (first) lithology.
    spar.single_top_unit = config.getboolean(section, 'single_top_unit')

    #---------------------------------
    section = 'DataPreprocessing'
    print(config.items(section))

    # The number of nearest units (for distance constraints).
    number_nearest_units = config.getint(section, 'number_nearest_units')

    # Minimum score for drillhole lithologies to use them.
    min_drillhole_litho_score = config.getint(section, 'min_drillhole_litho_score')

    # Group drillhole lithology sequence.
    # Note: use this for max_num_unit_contacts_inside_litho > 0 to avoid the solution number to blow.
    group_drillhole_lithos = config.getboolean(section, 'group_drillhole_lithos')

    # The cover ration threshold (relative length) for removing the cover.
    cover_ratio_threshold = config.getfloat(section, 'cover_ratio_threshold')

    #---------------------------------
    section = 'CollarIDs'
    print(config.items(section))

    collarIDs = config.get(section, 'collarIDs').split(",")
    collarIDs = [c.strip() for c in collarIDs]
    #-----------------------------------------------------------------------------------

    # Adding thickness constraints (requires unit thickness data which we do not have).
    add_thickness_constraints = False

    # Drillsample data csv file column names.
    drillsample_header = dr.DrillSampleDataHeader('Fromdepth', 'Todepth', 'Lithologies', 'Scores')

    # Strata data csv file column names.
    strata_data_header = dr.StrataDataHeader('UNITNAME', 'lithos', 'distance', 'DESCRIPTN')

    # Read topology data (the same for all collarIDs).
    map_graph = None
    if (spar.add_topology_constraints):
        map_graph = dr.read_topology_data(topology_filename)
        print("Number of units in the map graph =", map_graph.number_of_nodes())

    #====================================================================================
    # Process the list of CollarIDs.
    #====================================================================================
    display_plots = True

    # Stores solutions for different drillholes.
    strat_solutions = []

    for collarID in collarIDs:

        print('--------------------------------')
        print('collarID =', collarID)
        print('--------------------------------')

        drillsample_filename = drillsample_filename_collarID.replace("$collarID$", str(collarID))
        dist_table_filename = dist_table_filename_collarID.replace("$collarID$", str(collarID))

        #--------------------------------------------------------------
        # Reading the input data.
        #--------------------------------------------------------------
        # Read the drillhole ignore items list.
        ignore_list = dr.read_ignore_list(ignore_list_filename)

        # Read drill sample data.
        drillsample_data = dr.read_drillsample_data(drillsample_header, drillsample_filename, ignore_list, min_drillhole_litho_score)

        # Remove the Cover.
        drillsample_data.remove_cover(cover_unit_filename, cover_ratio_threshold)

        if (group_drillhole_lithos):
            # Group the drillsample lithologies.
            drillsample_data.group_drillhole_litho_sequence(spar.max_num_unit_contacts_inside_litho)

        # Read the alternative rock names.
        alternative_rock_names = dr.read_alternative_rock_names(alternative_rock_names_filename)

        # Read unit lithologies and distance data.
        strata_data = dr.read_strat_data(strata_data_header, dist_table_filename, alternative_rock_names)

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
            thickness_data = dr.read_thickness_data(thickness_filename)

        #--------------------------------------------------------------
        # Generating the stratigraphies.
        #--------------------------------------------------------------
        strat_solution = generate_strat_routes(spar, strata_data, drillsample_data, thickness_data, map_graph, alternative_rock_names)
        strat_solution.analyze_solution()
        strat_solution.collarID = collarID

        print("Total number of routes = ", len(strat_solution.routes))

        strat_solutions.append(strat_solution)

        #--------------------------------------------------------------
        # Plot the results.
        #--------------------------------------------------------------
        # Print all unique routes (i.e., unique strata sequence).
        print_unique_routes(strat_solution, 10)

        # Draw stratigraphy logs.
        draw_solution_logs(strat_solution, display_plots, 'strat', unit_colors_filename, True, None, [])

        # Draw probability logs.
        #draw_solution_logs(strat_solution, display_plots, 'proba', '', True, None, [])

        # Draw the age-rule logs.
        #draw_solution_logs(strat_solution, display_plots, 'age', '', True, map_graph, [])

        # Plot unit probabilities.
        #plot_unit_probabilities(strat_solution, display_plots)

        # Draw the topology graph of all solution routes.
        draw_solution_graph(strat_solution)

        # Plot histogram of solution scores.
        plot_route_scores(strat_solution.graph_route_scores)

        # Write the best routes to file.
        write_best_routes_to_file(strat_solution, 10)

    #--------------------------------------------------------------------
    # Solution correlation.
    #--------------------------------------------------------------------
    if (len(strat_solutions) > 1):
        compare_solution_graphs(strat_solutions[0].graph, strat_solutions[1].graph)

    correlate_solutions(strat_solutions)

    for solution in strat_solutions:
        plot_solution_correlation(solution)

    draw_correlated_solution_logs(strat_solutions, display_plots, unit_colors_filename, map_graph)

#=============================================================================
if __name__ == "__main__":
    main()

