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
    def __init__(self, routes, routes_number, unit_names):
        # The strata solution paths.
        self.routes = routes
        # The number of solutions for every drillsample raw.
        self.routes_number = routes_number
        # The unit names.
        self.unit_names = unit_names

        num_units = len(self.unit_names)

        # Building the distribution of unit presence at every depth.
        self.strat_distr = get_strat_distr(routes, num_units)

        # Calculate the route scores (based on path probability).
        self.route_scores = get_route_scores(routes, self.strat_distr)

#=============================================================================
def get_strat_distr(all_routes, num_units):
    '''
    Returns the distribution of unit presence at every depth.
    '''
    num_rows = len(all_routes[0].path)
    strat_distr = np.zeros((num_rows, num_units))

    for route in all_routes:
        for row, unit_index in enumerate(route.path):
            strat_distr[row, unit_index] += 1

    # Normalize.
    strat_distr = strat_distr / float(len(all_routes))
    return strat_distr

#=============================================================================
def get_route_scores(all_routes, strat_distr):
    '''
    Returns the route scores (based on path probability).
    Needs strat_distr returned by get_strat_distr().
    '''
    num_rows = len(all_routes[0].path)
    route_scores = np.zeros(len(all_routes), dtype=float)

    for route_index, route in enumerate(all_routes):
        for row, unit_index in enumerate(route.path):
            route_scores[route_index] += strat_distr[row, unit_index]

    # Normalize.
    route_scores = route_scores / float(num_rows)
    return route_scores

