'''
  This program generates a set of plausible stratigraphies with uncertainties, for a given drillhole lithology log.
  It uses map data for distance and topology constraints, and several free parameters describing the solution complexity (level of deformation) constraints.

  Author: Vitaliy Ogarko, vogarko@gmail.com
  The University of Western Australia
'''

import numpy as np
import networkx as nx

class StratSequence:
    def __init__(self, path):
        self.path = path

#========================================================================================================
class StrataSolution:
    '''
    Stores solutions of the Strata Solver, with their metadata and calculated scores.
    '''
    def __init__(self, routes, routes_number, unit_names, depth_data, unit2dist):
        # The strata solution paths.
        self.routes = routes
        # The number of solutions for every drillsample raw.
        self.routes_number = routes_number
        # The unit names.
        self.unit_names = unit_names
        # Depth data.
        self.depth_data = depth_data
        # Drillhole collar ID.
        self.collarID = 0

        num_units = len(self.unit_names)
        num_rows = len(depth_data.depth_from)

        # Building the distribution of unit presence at every depth.
        self.strat_distr = _calculate_strat_distr(routes, num_rows, num_units)

        # Calculate the route scores (based on path probability).
        self.route_scores = _calculate_route_scores(routes, self.strat_distr, self.depth_data)

        # Calcualte unique routes (e.g. two routes A-A-B-B-C and A-B-B-B-C become one route A-B-C).
        self.unique_routes = _calculate_unique_routes(routes)

        # Build the solution topology graph.
        self.graph = _build_solution_graph(self)

        # Calculate the route scores based on its own solution graph.
        self.graph_route_scores = self.calculate_graph_route_scores(self.graph)

        # Scores calculated using graphs from other drillholes.
        self.external_graph_route_scores_list = []

    #=====================================================================
    def unit_nonempty(self, unit_name):
        '''
        Checks if the unit has nonzero probability.
        '''
        unit_index = self.unit_names.index(unit_name)
        if (sum(self.strat_distr[:, unit_index]) != 0):
            return True
        else:
            return False

    #=============================================================================
    def calculate_graph_route_scores(self, graph):
        '''
        Calculate the route scores based on the solution topology graph.
        Note: it can use the external solution graph from other drillholes.
        '''
        num_routes = len(self.unique_routes)
        graph_route_scores = np.zeros(num_routes, dtype=float)
        unit_names = self.unit_names

        for route_index, route in enumerate(self.unique_routes):
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
                graph_route_scores[route_index] = float(score) / float(num_contacts)
            else:
                graph_route_scores[route_index] = 0.

        return graph_route_scores

    #=====================================================================
    def num_nonempty_units(self):
        '''
        Returns the number of non-empty units.
        '''
        counter = 0
        for unit_name in self.unit_names:
            if self.unit_nonempty(unit_name):
                counter += 1
        return counter

#=============================================================================
def _calculate_unique_routes(routes):
    unique_route_paths = set()
    for route in routes:
        unique_route_paths.add(route.get_strata_sequence())

    unique_routes = []
    for path in unique_route_paths:
        route = StratSequence(path)
        unique_routes.append(route)

    return unique_routes

#=============================================================================
def _calculate_strat_distr(all_routes, num_rows, num_units):
    '''
    Calculates the distribution of unit presence at every depth.
    '''
    strat_distr = np.zeros((num_rows, num_units))

    for route in all_routes:
        for row, unit_index in enumerate(route.path):
            strat_distr[row, unit_index] += 1

    if (len(all_routes) > 0):
        # Normalize.
        strat_distr = strat_distr / float(len(all_routes))

    return strat_distr

#=============================================================================
def _calculate_route_scores(all_routes, strat_distr, depth_data):
    '''
    Calculates the route scores (based on path probability).
    '''
    num_rows = strat_distr.shape[0]
    route_scores = np.zeros(len(all_routes), dtype=float)

    for route_index, route in enumerate(all_routes):
        total_length = 0.
        for row, unit_index in enumerate(route.path):
            length = depth_data.depth_to[row] - depth_data.depth_from[row]
            # Scale with the length as some data rows have different lenghts.
            route_scores[route_index] += strat_distr[row, unit_index] * length
            total_length += length
        # Normalize with the total drillhole length coverage.
        route_scores[route_index] /= total_length

    return route_scores

#=============================================================================
def _build_solution_graph(solution):
    '''
    Builds the solution topology graph, with edges weighted by unit contact frequency (among all solution routes).
    '''
    G = nx.DiGraph()
    #routes = solution.unique_routes
    routes = solution.routes

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
                        G[e[0]][e[1]]['weight'] = G[e[0]][e[1]]['weight'] + 1

    # Normalize the edge weight.
    for u, v, d in G.edges(data=True):
        G[u][v]['weight'] = float(d['weight']) / float(len(routes))

    return G
