'''
  This program generates a set of plausible stratigraphies with uncertainties, for a given drillhole lithology log.
  It uses map data for distance and topology constraints, and several free parameters describing the solution complexity (level of deformation) constraints.

  Author: Vitaliy Ogarko, vogarko@gmail.com
  The University of Western Australia
'''

from argparse import ArgumentParser
from dataclasses import fields

from strata_solver import \
    generate_strat_routes, \
    StrataSolverParameters
import data_readers as dr
from solution_utils import *
from solution_analysis import *
from input_params import *

#=============================================================================
def solve(par, drillsample_header, strata_data_header):
    '''
    Read data from files and run the solver.
    '''
    spar = StrataSolverParameters()

    # Copy solver parameters.
    for field in fields(StrataSolverParameters):
        setattr(spar, field.name, getattr(par, field.name))

    #====================================================================================
    # Read data common for all collarIDs.
    #====================================================================================
    # Read topology data (graph).
    map_graph = None
    if (par.add_topology_constraints):
        map_graph = dr.read_topology_data(par.topology_filename)
        print("Number of units in the map graph =", map_graph.number_of_nodes())

    # Read the drillhole ignore items list.
    ignore_list = dr.read_ignore_list(par.ignore_list_filename)

    # Read the alternative rock names.
    alternative_rock_names = dr.read_alternative_rock_names(par.alternative_rock_names_filename)

    # Read cover lithologies.
    cover_lithos = dr.read_cover_lithos(par.cover_unit_filename)

    # Read thickness data.
    thickness_data = []
    if (par.add_thickness_constraints):
        thickness_data = dr.read_thickness_data(thickness_filename)

    #====================================================================================
    # Process the list of CollarIDs.
    #====================================================================================
    display_plots = True

    # Stores solutions for different drillholes.
    strat_solutions = []

    for collarID in par.collarIDs:

        print('--------------------------------')
        print('collarID =', collarID)
        print('--------------------------------')

        drillsample_filename = par.drillsample_filename.replace("$collarID$", str(collarID))
        stratasample_filename = par.stratasample_filename.replace("$collarID$", str(collarID))
        dist_table_filename = par.dist_table_filename.replace("$collarID$", str(collarID))

        stratasample_header = dr.DrillSampleDataHeader('DEPTH_FROM_M', 'DEPTH_TO_M', 'STRAT_UNIT_NAME', '')

        #--------------------------------------------------------------
        # Read and preprocess the drillsample and unit data.
        #--------------------------------------------------------------
        # Read drill sample data.
        drillsample_data = dr.read_drillsample_data(drillsample_header, drillsample_filename, ignore_list, par.min_drillhole_litho_score)

        # Read the drillhole stratigraphy data.
        strata_sample_log = read_strata_sample_data(stratasample_header, stratasample_filename)
        strata_sample_log.collarID = collarID

        # Remove the Cover.
        if (len(cover_lithos) > 0):
            drillsample_data.remove_cover(cover_lithos, par.cover_ratio_threshold)

        if (par.group_drillhole_lithos):
            # Group the drillsample lithologies.
            drillsample_data.group_drillhole_litho_sequence(spar.max_num_unit_contacts_inside_litho)

        # Read unit lithologies and distance data.
        strata_data = dr.read_strat_data(strata_data_header, dist_table_filename, alternative_rock_names)

        # Filter strat data.
        drillhole_lithos = drillsample_data.get_drillhole_lithos()
        strata_data.filter_strat_data_based_on_drillhole_lithos(drillhole_lithos)
        strata_data.filter_strat_data_based_on_distance(par.number_nearest_units)

        filtered_lithos = strata_data.get_unique_lithos()
        print("The number of filtered unit lithologies:", len(filtered_lithos))
        print("Filtered unit lithologies:", sorted(filtered_lithos))

        #--------------------------------------------------------------
        # Generating the stratigraphies.
        #--------------------------------------------------------------
        strat_solution = generate_strat_routes(spar, strata_data, drillsample_data, thickness_data, map_graph, alternative_rock_names)
        strat_solution.analyze_solution(par.correlation_power)
        strat_solution.collarID = collarID

        print("Total number of routes = ", len(strat_solution.routes))

        strat_solutions.append(strat_solution)

        #--------------------------------------------------------------
        # Plot the results.
        #--------------------------------------------------------------
        # Print all unique routes (i.e., unique strata sequence).
        print_unique_routes(strat_solution, 10)

        # Draw strata sample log (known strata solution).
        draw_solution_logs(strata_sample_log, display_plots, 'strat', par.unit_colors_filename, True, None, [0])

        # Draw stratigraphy logs.
        draw_solution_logs(strat_solution, display_plots, 'strat', par.unit_colors_filename, True, None, [])

        # Draw probability logs.
        #draw_solution_logs(strat_solution, display_plots, 'proba', '', True, None, [])

        # Draw the age-rule logs.
        #draw_solution_logs(strat_solution, display_plots, 'age', '', True, map_graph, [])

        # Plot unit probabilities.
        plot_unit_probabilities(strat_solution, display_plots)

        # Draw the topology graph of all solution routes.
        draw_solution_graph(strat_solution)

        # Plot histogram of solution scores.
        plot_route_scores(strat_solution.graph_route_scores)

        # Write the best routes to file.
        write_best_routes_to_file(strat_solution, 10)

    #====================================================================================
    # Solution correlation.
    #====================================================================================
    if (len(strat_solutions) > 1):
        compare_solution_graphs(strat_solutions[0].graph, strat_solutions[1].graph)

    correlate_solutions(strat_solutions, par.correlation_power)

    for solution in strat_solutions:
        plot_solution_correlation(solution)

    # Draw the most correlated strata sequences.
    for solution in strat_solutions:
        draw_solution_logs(solution, display_plots, 'strat-seq', par.unit_colors_filename, False, None, [])

    #draw_correlated_solution_logs(strat_solutions, display_plots, par.unit_colors_filename, map_graph)

#=========================================================================================================
def read_strata_sample_data(stratasample_header, stratasample_filename):
    '''
    Reads the drillhole stratigraphy sample data.
    '''

    stratasample_data = dr.read_drillsample_data(stratasample_header, stratasample_filename, [], 0)
    depth_data = stratasample_data.get_depth_data()

    # Create route object.
    class Object(object):
        pass
    route = Object()
    route.path = []

    unit_names = []

    for row in stratasample_data.rows:
        unit_name = dr.align_unit_name(row.lithos[0])
        if unit_name not in unit_names:
            unit_names.append(unit_name)
        route.path.append(unit_names.index(unit_name))

    routes = [route]
    route_scores = [1.]

    strata_log = StrataLog(0, routes, route_scores, unit_names, depth_data)

    return strata_log

#=============================================================================
def main(parfile_path):
    print('Started litho2strat')

    # Read input parameters.
    par = read_input_parameters(parfile_path)

    # TODO: Define headers in the Parfile parameters.
    use_SA_data = True

    if use_SA_data:
        # Drillsample data csv file column names.
        drillsample_header = dr.DrillSampleDataHeader('DEPTH_FROM_M', 'DEPTH_TO_M', 'MAJOR_LITHOLOGY', '')

        # Strata data csv file column names.
        strata_data_header = dr.StrataDataHeader('strat', 'summary', 'distance', 'description')
    else:
        # Drillsample data csv file column names.
        drillsample_header = dr.DrillSampleDataHeader('Fromdepth', 'Todepth', 'Lithologies', 'Scores')

        # Strata data csv file column names.
        strata_data_header = dr.StrataDataHeader('UNITNAME', 'lithos', 'distance', 'DESCRIPTN')

    # Run the solver.
    solve(par, drillsample_header, strata_data_header)

#=============================================================================
if __name__ == "__main__":
    # Read command line arguments.
    parser = ArgumentParser()
    parser.add_argument("-p", "--parfile", dest="parfile_path",
                    help="path to the parameters file", default="parfiles/Parfile.txt")

    args = parser.parse_args()

    # The main program.
    main(args.parfile_path)

