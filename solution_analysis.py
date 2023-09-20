'''
  This program generates a set of plausible stratigraphies with uncertainties, for a given drillhole lithology log.
  It uses map data for distance and topology constraints, and several free parameters describing the solution complexity (level of deformation) constraints.

  Author: Vitaliy Ogarko, vogarko@gmail.com
  The University of Western Australia
'''
import numpy as np
import networkx as nx

#=============================================================================
def build_solution_graph(solution):
    '''
    Builds the solution topology graph, with edges weighted by unit contact frequency (among all solution routes).
    '''
    G = nx.DiGraph()
    # Use unique routes here, as using all routes favours inclusion of units when it should not (SA collarID=265020 with diamictite added to Angepena Formation).
    routes = solution.unique_routes
    #routes = solution.routes

    for route in routes:
        route_edges = set()
        for i in range(len(route.path) - 1):
            unit_name1 = solution.unit_names[route.path[i]]
            unit_name2 = solution.unit_names[route.path[i + 1]]

            if (unit_name1 != unit_name2):
                e = (unit_name1, unit_name2)

                if (e not in route_edges):
                # Count only once each contact type on the route.
                # Note: they still will be included into the route score multiple times when we calculate the score based on all contacts on the route.
                    route_edges.add(e)
                    if (not G.has_edge(*e)):
                        # Adding a new graph edge.
                        G.add_edge(*e, weight=1)
                    else:
                        # Increase the edge weight.
                        G[e[0]][e[1]]['weight'] += 1

    # Normalize the edge weight.
    for u, v, d in G.edges(data=True):
        G[u][v]['weight'] = float(d['weight']) / float(len(routes))

    return G

#=============================================================================
def calculate_graph_route_scores(solution, graph, correlation_power):
    '''
    Calculate the route scores based on the solution topology graph.
    Note: the graph can be external, i.e., from another drillhole.
    '''
    num_routes = len(solution.unique_routes)
    graph_route_scores = np.zeros(num_routes, dtype=float)
    unit_names = solution.unit_names

    for route_index, route in enumerate(solution.unique_routes):
        unique_route = route.path
        score = 0.
        num_contacts = len(unique_route) - 1
        for i in range(num_contacts):
            # Graph edge.
            e = (unit_names[unique_route[i]], unit_names[unique_route[i + 1]])
            if (graph.has_edge(*e)):
            # When using the external graph (from another drillhole) it may not have this edge.
                # Adding the graph edge weight to the score.
                weight = graph[e[0]][e[1]]['weight']
                score = score + weight

        if (num_contacts > 0):
            # Normalize the score with the number of contacts.
            graph_route_scores[route_index] = float(score) / float(num_contacts)**correlation_power
        else:
            graph_route_scores[route_index] = 0.

    return graph_route_scores

#=============================================================================
def calculate_weighted_Jaccard_index(G1, G2):
    '''
    Calculate the weighted Jaccard index.
    '''
    # Union of all graph edges.
    edges = set(G1.edges()).union(G2.edges())

    mins = 0.
    maxs = 0.
    for edge in edges:
        weight1 = G1.get_edge_data(*edge, {}).get('weight', 0.)
        weight2 = G2.get_edge_data(*edge, {}).get('weight', 0.)
        mins += min(weight1, weight2)
        maxs += max(weight1, weight2)

    if (maxs != 0.):
        Jw = mins / maxs
    else:
        Jw = 0.

    return Jw

#=============================================================================
def calculate_correlation_matrix(solutions):
    '''
    Calculate correlatoin matrix between different drillholes (using weighted Jaccard index).
    '''
    num_solutions = len(solutions)
    Jw = np.zeros((num_solutions, num_solutions))
    for i in range(num_solutions):
        for j in range(num_solutions):
            Jw[i, j] = calculate_weighted_Jaccard_index(solutions[i].graph, solutions[j].graph)
    return Jw

#=============================================================================
def correlate_solutions(strat_solutions, correlation_power):
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
                graph_route_scores = calculate_graph_route_scores(strat_solutions[j], solution_graph, correlation_power)
                strat_solutions[j].external_graph_route_scores_list.append(graph_route_scores)

    # Calculate a correlated route score (as the sum of graph scores ovel all drillholes).
    for solution in strat_solutions:
        solution.unique_route_scores = np.copy(solution.graph_route_scores)
        for external_graph_route_scores in solution.external_graph_route_scores_list:
            solution.unique_route_scores += external_graph_route_scores

        solution.unique_route_scores /= float(len(strat_solutions))

#====================================================================================================
def get_correlated_solution_logs(solution):
    '''
    Retrieve the most correlated solution logs (corresponding to most correlated strata sequences).
    '''
    # Sorted indexes by the unique score.
    unique_route_indexes = np.argsort(-solution.unique_route_scores)

    # Retrieve the corresponding full routes.
    route_indexes = []
    route_scores = []
    for unique_route_index in unique_route_indexes:
        for route_index in solution.unique_routes[unique_route_index].route_indexes:
            route_indexes.append(route_index)
            route_scores.append(solution.unique_route_scores[unique_route_index])

    return route_indexes, route_scores
