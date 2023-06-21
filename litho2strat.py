'''
  This program generates a set of plausible stratigraphies with uncertainties, for a given drillhole lithology log.
  It uses map data for distance and topology constraints, and several free parameters describing the solution complexity (level of deformation) constraints.

  Author: Vitaliy Ogarko, vogarko@gmail.com
  The University of Western Australia
'''

import os
import configparser
from argparse import ArgumentParser
from dataclasses import dataclass, field, fields
from typing import List

from strata_solver import \
    generate_strat_routes, \
    StrataSolverParameters
import data_readers as dr
from solution_utils import *
from solution_analysis import *

#=============================================================================
@dataclass
class InputParameters:
    '''
    Contains all input parameters used in the Parfile.
    '''
    #-------------------------------
    # Section 'FilePaths'.
    #-------------------------------
    # Topology file.
    topology_filename: str = ""
    # The Cover unit data file.
    cover_unit_filename: str = ""
    # The ignore items list file.
    ignore_list_filename: str = ""
    # Alternative rock names file.
    alternative_rock_names_filename: str = ""
    # Unit colours for drawing stratigraphy logs.
    unit_colors_filename: str = ""
    # Drillhole lithology data file. The $collarID$ in the file name will be replaced with the actual value.
    drillsample_filename_collarID: str = ""
    # Units near the collar with distances data file. The $collarID$ in the file name will be replaced with the actual value.
    dist_table_filename_collarID: str = ""

    #-------------------------------
    # Section 'SolverParameters'.
    #-------------------------------
    # Unit contact topology constraints (extracted from map data).
    add_topology_constraints: bool = True
    # Jumps over units in topology graph to allow skipping units: i.e., one jump would allow contact A->C for the graph A->B->C.
    max_num_strata_jumps: int = 0
    # 'Age alignment' constraints: the maximum number of times the age direction can flip.
    max_num_age_flips: int = 2
    # 'Returning to the same unit' constraints.
    max_num_returns_per_unit: int = 1
    # The maximum number of unit contacts inside the same litholgy sequence.
    max_num_unit_contacts_inside_litho: int = 0
    # Use the single closest unit for the top (first) lithology.
    single_top_unit: bool = True

    #-------------------------------
    # Section 'Correlation'.
    #-------------------------------
    # Correlation score normalization power. Higher power leads to shorter strata sequence.
    correlation_power: float = 1.0

    #-------------------------------
    # Section 'DataPreprocessing'.
    #-------------------------------
    # The number of nearest units (for distance constraints).
    number_nearest_units: int = 3
    # Minimum score for drillhole lithologies to use them.
    min_drillhole_litho_score: int = 80
    # Group drillhole lithology sequence.
    # Note: use this for max_num_unit_contacts_inside_litho > 0 to avoid the solution number to blow.
    group_drillhole_lithos: bool = False
    # The cover ration threshold (relative length) for removing the cover.
    cover_ratio_threshold: float = 0.65

    #-------------------------------
    # Section 'CollarIDs'.
    #-------------------------------
    # List of collar IDs used.
    collarIDs: List[str] = field(default_factory=list)

    #-------------------------------
    # Other.
    #-------------------------------
    # Adding thickness constraints (requires unit thickness data which we do not have).
    add_thickness_constraints: bool = False

#=============================================================================
def read_input_parameters(parfile_path):
    '''
    Read input parameters from Parfile.
    '''
    config = configparser.ConfigParser()
    if (len(config.read(parfile_path)) == 0):
        raise ValueError("Failed to open/find a parameters file!")

    par = InputParameters()

    section = 'FilePaths'
    print(config.items(section))

    par.topology_filename = config.get(section, 'topology_filename', fallback = '')
    par.cover_unit_filename = config.get(section, 'cover_unit_filename', fallback = '')
    par.ignore_list_filename = config.get(section, 'ignore_list_filename', fallback = '')
    par.alternative_rock_names_filename = config.get(section, 'alternative_rock_names_filename', fallback = '')
    par.unit_colors_filename = config.get(section, 'unit_colors_filename', fallback = '')
    par.drillsample_filename_collarID = config.get(section, 'drillsample_filename')
    par.dist_table_filename_collarID = config.get(section, 'dist_table_filename')

    section = 'SolverParameters'
    print(config.items(section))

    par.add_topology_constraints = config.getboolean(section, 'add_topology_constraints', fallback = par.add_topology_constraints)
    par.max_num_strata_jumps = config.getint(section, 'max_num_strata_jumps', fallback = par.max_num_strata_jumps)
    par.max_num_age_flips = config.getint(section, 'max_num_age_flips', fallback = par.max_num_age_flips)
    par.max_num_returns_per_unit = config.getint(section, 'max_num_returns_per_unit', fallback = par.max_num_returns_per_unit)
    par.max_num_unit_contacts_inside_litho = config.getint(section, 'max_num_unit_contacts_inside_litho', fallback = par.max_num_unit_contacts_inside_litho)
    par.single_top_unit = config.getboolean(section, 'single_top_unit', fallback = par.single_top_unit)

    section = 'Correlation'
    print(config.items(section))

    par.correlation_power = config.getfloat(section, 'correlation_power', fallback = par.correlation_power)

    section = 'DataPreprocessing'
    print(config.items(section))

    par.number_nearest_units = config.getint(section, 'number_nearest_units', fallback = par.number_nearest_units)
    par.min_drillhole_litho_score = config.getint(section, 'min_drillhole_litho_score', fallback = par.min_drillhole_litho_score)
    par.group_drillhole_lithos = config.getboolean(section, 'group_drillhole_lithos', fallback = par.group_drillhole_lithos)
    par.cover_ratio_threshold = config.getfloat(section, 'cover_ratio_threshold', fallback = par.cover_ratio_threshold)

    section = 'CollarIDs'
    print(config.items(section))

    collarIDs = config.get(section, 'collarIDs').split(",")
    par.collarIDs = [c.strip() for c in collarIDs]

    # Hardcoded as we don't have the thickness data.
    par.add_thickness_constraints = False

    return par

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

        drillsample_filename = par.drillsample_filename_collarID.replace("$collarID$", str(collarID))
        dist_table_filename = par.dist_table_filename_collarID.replace("$collarID$", str(collarID))

        #--------------------------------------------------------------
        # Read and preprocess the drillsample and unit data.
        #--------------------------------------------------------------
        # Read drill sample data.
        drillsample_data = dr.read_drillsample_data(drillsample_header, drillsample_filename, ignore_list, par.min_drillhole_litho_score)

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

