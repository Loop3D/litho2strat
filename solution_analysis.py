'''
  This program generates a set of plausible stratigraphies with uncertainties, for a given drillhole lithology log.
  It uses map data for distance and topology constraints, and several free parameters describing the solution complexity (level of deformation) constraints.

  Author: Vitaliy Ogarko, vogarko@gmail.com
  The University of Western Australia
'''
from solution_utils import *

def correlate_solutions(strat_solutions):
    '''
    Correlete solutions from different drillholes.
    Build the correlated solution score.
    '''

    # Loop over all solution pairs.
    for i in range(len(strat_solutions)):
        solution_graph = strat_solutions[i].graph
        for j in range(len(strat_solutions)):
            if (i != j):
                # Calculate solution scores based on external solution graph.
                graph_route_scores = strat_solutions[j].calculate_graph_route_scores(solution_graph)
                strat_solutions[j].external_graph_route_scores_list.append(graph_route_scores)

    # Calculate a correlated route score (as the sum of graph scores ovel all drillholes).
    for solution in strat_solutions:
        solution.unique_route_scores = solution.graph_route_scores
        for external_graph_route_scores in solution.external_graph_route_scores_list:
            solution.unique_route_scores = solution.unique_route_scores + external_graph_route_scores

#====================================================================================================
def draw_correlated_solution_logs(strat_solutions, display_plots, unit_colors_filename, map_graph):
    '''
    Draw correlated solution logs.
    '''

    # Draw the most correlated strata sequences.
    for solution in strat_solutions:
        draw_solution_logs(solution, display_plots, 'strat-seq', unit_colors_filename, False, None, [])

    # Draw full routes (corresponding to the most correlated strata sequences).
    for solution in strat_solutions:
        # Sorted indexes by score.
        unique_route_indexes = np.argsort(-solution.unique_route_scores)
        # Retrieve corresponding full route indexes.
        route_indexes = []
        route_scores = []
        for unique_route_index in unique_route_indexes:
            # Select a corresponding full route.
            #route_index = solution.unique_routes[unique_route_index].route_indexes[0]

            # Choose the full route with the highest own route score.
            best_score = -1.
            for index in solution.unique_routes[unique_route_index].route_indexes:
                route_score = solution.route_scores[index]
                if (route_score > best_score):
                    route_index = index

            route_indexes.append(route_index)
            route_scores.append(solution.unique_route_scores[unique_route_index])

        # Draw full routes (corresponding to strata sequences).
        draw_solution_logs(solution, display_plots, 'strat', unit_colors_filename, False, None, route_indexes, route_scores)

        # Draw age alignment log (corresponding to strata sequences).
        draw_solution_logs(solution, display_plots, 'age', '', False, map_graph, route_indexes, route_scores)
