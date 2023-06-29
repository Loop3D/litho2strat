'''
  This program generates a set of plausible stratigraphies with uncertainties, for a given drillhole lithology log.
  It uses map data for distance and topology constraints, and several free parameters describing the solution complexity (level of deformation) constraints.

  Author: Vitaliy Ogarko, vogarko@gmail.com
  The University of Western Australia
'''

import numpy as np
import networkx as nx

from solution_analysis import \
    build_solution_graph, \
    calculate_graph_route_scores

class StratSequence:
    def __init__(self, path, route_indexes):
        # Unique route path.
        self.path = path
        # Indexes of original full routes.
        self.route_indexes = route_indexes

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

        # Adding the cover unit.
        # Do this to be able to correlate solutions made of one unit only and thus have zero contacts.
        self.add_cover_unit()

        num_units = len(self.unit_names)
        num_rows = len(self.depth_data.depth_from)

        # Building the distribution of unit presence at every depth.
        self.strat_distr = _calculate_strat_distr(self.routes, num_rows, num_units)

        # Calculate the route scores (based on path probability).
        self.route_scores = _calculate_route_scores(self.routes, self.strat_distr, self.depth_data)

        # Calcualte unique routes (e.g. two routes A-A-B-B-C and A-B-B-B-C become one route A-B-C).
        self.unique_routes = _calculate_unique_routes(self.routes)

        # The solution topology graph.
        self.graph = None

        # The route scores based on its own solution graph.
        self.graph_route_scores = None

        # Scores calculated using graphs from other drillholes.
        self.external_graph_route_scores_list = []

        self.unique_route_scores = None

    #=====================================================================
    def add_cover_unit(self):
        '''
        Adding the cover contact.
        '''
        # Adjust unit names list.
        self.unit_names.append("cover")
        cover_index = self.unit_names.index("cover")

        # Adjust solution routes.
        for route in self.routes:
            route.path.insert(0, cover_index)

        # Adjust depth data.
        self.depth_data.depth_to.insert(0, self.depth_data.depth_from[0])
        self.depth_data.depth_from.insert(0, 0.)

    #=====================================================================
    def unit_nonempty(self, unit_name):
        '''
        Checks if the unit has nonzero probability.
        '''
        if (unit_name in self.unit_names):
            unit_index = self.unit_names.index(unit_name)
            if (sum(self.strat_distr[:, unit_index]) != 0):
                return True
            else:
                return False
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

    #=====================================================================
    def analyze_solution(self, correlation_power):
        '''
        Generates the solution graph and graph scores.
        '''
        self.graph = build_solution_graph(self)
        self.graph_route_scores = calculate_graph_route_scores(self, self.graph, correlation_power)

#=============================================================================
def _calculate_unique_routes(routes):
    # Building a dictionary with keys = unique route paths, and values = list of indexes of the original routes.
    unique2full = dict()
    for route_index, route in enumerate(routes):
        unique_path = route.get_strata_sequence()
        if (unique_path in unique2full):
            unique2full[unique_path].append(route_index)
        else:
            unique2full[unique_path] = [route_index]

    # Build a list of unique routes.
    unique_routes = []
    for key, value in unique2full.items():
        unique_route = StratSequence(key, value)
        unique_routes.append(unique_route)

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
