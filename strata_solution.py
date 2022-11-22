'''
  This program generates a set of plausible stratigraphies with uncertainties, for a given drillhole lithology log.
  It uses map data for distance and topology constraints, and several free parameters describing the solution complexity (level of deformation) constraints.

  Author: Vitaliy Ogarko, vogarko@gmail.com
  The University of Western Australia
'''

import numpy as np
import networkx as nx

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

        # Calculate the route distance scores (based on RMS distance to units).
        #self.dist_scores = _calculate_distance_scores(routes, unit_names, unit2dist)

        # Build the solution topology graph.
        self.graph = _build_solution_graph(self)

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
def _calculate_distance_scores(all_routes, unit_names, unit2dist):
    '''
    Calculates the route scores (based on root mean square distance to units).
    '''
    dist_scores = np.zeros(len(all_routes), dtype=float)

    if (len(all_routes) == 0):
        return dist_scores

    num_rows = len(all_routes[0].path)

    for route_index, route in enumerate(all_routes):
        for row, unit_index in enumerate(route.path):
            unit_name = unit_names[unit_index]
            distance = unit2dist[unit_name]
            dist_scores[route_index] += distance**2

    # Normalize.
    if (num_rows > 0):
        dist_scores = np.sqrt(dist_scores / float(num_rows))

    return dist_scores

#=============================================================================
def _build_solution_graph(strat_solution):
    '''
    Builds the solution topology graph, with edges weighted by unit contact frequency (among all solution routes).
    '''
    G = nx.Graph()

    for route in strat_solution.routes:
        unique_route = route.get_strata_sequence()
        for i in range(len(unique_route) - 1):
            e = (unique_route[i], unique_route[i + 1])
            if (not G.has_edge(*e)):
                G.add_edge(*e, weight=1)
            else:
                # Increase the edge weight.
                G[e[0]][e[1]]['weight'] = G[e[0]][e[1]]['weight'] + 1

    # Define graph node names.
    nodes_mapping = dict()
    for unit_index, unit_name in enumerate(strat_solution.unit_names):
        nodes_mapping[unit_index] = unit_name

    # Label graph nodes with unit names.
    G = nx.relabel_nodes(G, nodes_mapping)

    return G

