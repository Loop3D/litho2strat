'''
  This program generates a set of plausible stratigraphies with uncertainties, for a given drillhole lithology log.
  It uses map data for distance and topology constraints, and several free parameters describing the solution complexity (level of deformation) constraints.

  Author: Vitaliy Ogarko, vogarko@gmail.com
  The University of Western Australia
'''

import numpy as np

#========================================================================================================
class StrataSolution:
    '''
    The solutions of the Strata Solver, with their scores.
    '''
    def __init__(self, routes, routes_number, unit_names, depth_data):
        # The strata solution paths.
        self.routes = routes
        # The number of solutions for every drillsample raw.
        self.routes_number = routes_number
        # The unit names.
        self.unit_names = unit_names
        # Depth data.
        self.depth_data = depth_data

        num_units = len(self.unit_names)
        num_rows = len(depth_data.depth_from)

        # Building the distribution of unit presence at every depth.
        self.strat_distr = _get_strat_distr(routes, num_rows, num_units)

        # Calculate the route scores (based on path probability).
        self.route_scores = _get_route_scores(routes, self.strat_distr)

#=============================================================================
def _get_strat_distr(all_routes, num_rows, num_units):
    '''
    Returns the distribution of unit presence at every depth.
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
def _get_route_scores(all_routes, strat_distr):
    '''
    Returns the route scores (based on path probability).
    Needs strat_distr returned by get_strat_distr().
    '''
    num_rows = strat_distr.shape[0]
    route_scores = np.zeros(len(all_routes), dtype=float)

    for route_index, route in enumerate(all_routes):
        for row, unit_index in enumerate(route.path):
            route_scores[route_index] += strat_distr[row, unit_index]

    # Normalize.
    if (num_rows > 0):
        route_scores = route_scores / float(num_rows)

    return route_scores

