'''
  This program generates a set of plausible stratigraphies with uncertainties, for a given drillhole lithology log.
  It uses map data for distance and topology constraints, and several free parameters describing the solution complexity (level of deformation) constraints.

  Author: Vitaliy Ogarko, vogarko@gmail.com
  The University of Western Australia
'''

import numpy as np
from dataclasses import dataclass, field
from typing import List

from strata_solution import StrataSolution

#========================================================================================================
@dataclass
class StrataSolverParameters:
    '''
    The strata solver parameters holder.
    '''
    # 'Returning to the same unit' constraints.
    max_num_returns_per_unit: int = 0
    # The number of unit contacts inside the same litholgy sequence.
    max_num_unit_contacts_inside_litho: int = 0
    # Use the single closest unit for the top (first) lithology.
    single_top_unit: bool = False
    # Flag for using the unit age to determine connectivity.
    ignore_unit_age: bool = False

#========================================================================================================
@dataclass
class StrataTableElement:
    '''
    The element of a strata table.
    '''
    path_exists: bool = False
    lithos: List[str] = field(default_factory=list)

#=======================================================================================
def generate_strata_table(drillsample_data, strata_data, single_top_unit):
    '''
    Generates the stratigraphic table, and unit names list.
    '''
    num_rows = drillsample_data.get_num_rows()
    num_units = strata_data.get_num_units()
    unit_names = strata_data.get_unit_names()

    print("num_rows (before) = ", num_rows)
    print("num_units = ", num_units)

    strata_table = np.ndarray(shape=[num_rows, num_units], dtype=object)

    # Initialize with default objects.
    strata_table.flat = [StrataTableElement() for _ in strata_table.flat]

    new_row_index = 0

    # Missing lithologies.
    missing_lithos = set()

    for row in drillsample_data.rows[:]:
        any_litho_found = False

        # Loop over lithos in the current drillsample row.
        for litho in row.lithos[:]:
            litho_found = False

            if (new_row_index == 0 and single_top_unit):
            # Use only the closest unit for the top lithology.
                if (litho in strata_data.litho2dist):
                    # Sorted distance list for this lithology.
                    dist_list = strata_data.litho2dist[litho]

                    # Adding Cover if it is present for this litho.
                    add_cover = False
                    for item in dist_list:
                        if (item[1] == 'Cover'):
                            add_cover = True
                            closest_unit = 'Cover'
                            closest_unit_distance = item[0]

                    # Finding the closest unit that is not Cover.
                    for item in dist_list:
                        if (item[1] != 'Cover'):
                            closest_unit = item[1]
                            closest_unit_distance = item[0]
                            break

                    print('Closest top unit info (litho, unit, distance):', [litho, closest_unit, closest_unit_distance])

                    for unit_name in strata_data.unit2litho:
                        if (unit_name == closest_unit or (add_cover and unit_name == 'Cover')):
                            litho_found = True
                            unit_index = unit_names.index(unit_name)
                            strata_table[new_row_index, unit_index].path_exists = True
                            strata_table[new_row_index, unit_index].lithos.append(litho)

            else:
                for unit_name in strata_data.unit2litho:
                    if (litho in strata_data.unit2litho[unit_name]):
                        litho_found = True
                        unit_index = unit_names.index(unit_name)
                        strata_table[new_row_index, unit_index].path_exists = True
                        strata_table[new_row_index, unit_index].lithos.append(litho)

            if (not litho_found):
            # Drillhole lithology not found in units data.
                print("WARNING: Drillhole lithologies not found in units data: ", litho)

                # Remove this lithology and its score from drillsample data for the current row.
                del row.scores[row.lithos.index(litho)]
                row.lithos.remove(litho)

                # Update the missing lithos list.
                missing_lithos.add(litho)
            else:
                any_litho_found = True

        if (not any_litho_found):
            # Treat this as "no data".
            drillsample_data.rows.remove(row)
        else:
            new_row_index += 1

    num_rows = drillsample_data.get_num_rows()
    print("num_rows (after) = ", num_rows)

    if (num_rows == 0):
        print('No data left!!!')
        return np.ndarray(shape=[num_rows, num_units], dtype=object), set()

    # Remove rows due to missing lithologies.
    strata_table = strata_table[0:num_rows, :]

    return strata_table, missing_lithos

#==============================================================================
'''
A class for storing the stratigraphic route.
'''
class StrataRoute:
    __slots__ = 'to_remove', 'path', 'current_thickness', 'unit_visited', 'num_unit_contacts_inside_litho'

    def __init__(self, num_units):
        # Flag for removal.
        self.to_remove = False
        # Containts the number of times each unit was visited.
        self.unit_visited = np.zeros((num_units), dtype=int)
        # The number of unit contacts inside the same lithology sequence.
        self.num_unit_contacts_inside_litho = 0

    def add_first_position(self, strat, thickness_change):
        '''
        Adds the first position to the route.
        '''
        # Unit index for every drillhole data.
        self.path = [strat]
        # The thickness of the last strata unit.
        self.current_thickness = thickness_change
        # Mark this unit as 'visited'.
        self.unit_visited[strat] += 1

    def get_strata_sequence(self):
        '''
        Returns a sequence of strata units in the route excluding consequtive dublicates.
        '''
        return tuple([v for i, v in enumerate(self.path) if i == 0 or v != self.path[i - 1]])

#==============================================================================
def flatten(S):
    '''
    Flattens the multilevel list of lists.
    For example, it will convert [[[1,1],2,2],3,3] to [1,1,2,2,3,3].
    '''
    if S == []:
        return S
    if isinstance(S[0], list):
        return flatten(S[0]) + flatten(S[1:])
    return S[:1] + flatten(S[1:])

#==============================================================================
def apply_max_num_returns_constraint(route, strata_list, max_num_returns):
    '''
    Apply the "maximum number of returns to a unit" constraint:
        remove from the input unit list the units where the route cannot return anymore.
    '''
    # Apply the "max number of returns" constraint.
    for strat in strata_list[:]:
        if (route.unit_visited[strat] - 1 >= max_num_returns):
        # Reached the maximum numer of local returns (to this unit).
            strata_list.remove(strat)

#==============================================================================
def apply_topology_constraints(graph, unit_names, strat0, strata_list, ignore_unit_age):
    '''
    Apply unit topology (connectivity) constraints:
        remove from the input unit list the units not connected to a given one (strat0).
    '''
    # Cover can be in contact with any unit.
    if (unit_names[strat0] == 'Cover'):
        return

    if (unit_names[strat0] in graph.nodes()):
        for strat in strata_list[:]:
            e = (unit_names[strat0], unit_names[strat])
            if (not ignore_unit_age):
                if (not graph.has_edge(*e)):
                    # Units are not connected. Skip this unit.
                    strata_list.remove(strat)
            else:
            # NOTE: We test both edges here instead of converting graph to undirected, as we need to keep the age info in the graph.
                e2 = (unit_names[strat], unit_names[strat0])
                if (not graph.has_edge(*e) and not graph.has_edge(*e2)):
                    # Units are not connected. Skip this unit.
                    strata_list.remove(strat)

#=======================================================================================================
def generate_strat_routes(spar, strata_data, drillsample_data, thickness_data, graph):
    '''
    Generating stratigraphic routes.
    '''
    print('Starting strata solver with:', spar)

    if (len(drillsample_data.rows) == 0):
        # Empty drillsample data - return.
        return StrataSolution([], [], [], drillsample_data.get_depth_data())

    # Generating the table of possible strata paths.
    strata_table, missing_lithos = generate_strata_table(drillsample_data, strata_data, spar.single_top_unit)

    # Unit index to unit name mapping.
    unit_names = strata_data.get_unit_names()

    num_rows = strata_table.shape[0]
    num_units = strata_table.shape[1]

    # Local flags are based on the function input data.
    add_thickness_constraints = (len(thickness_data) > 0)
    add_topology_constraints = (graph != None)

    if (add_topology_constraints):
        # Sanity check: check that strata units exist in the graph.
        for unit_name in unit_names:
            if unit_name not in graph.nodes():
                print("WARNING: Not found graph unit: ", unit_name)

    # Extract strata thikcness to lists (faster data structures).
    min_strata_thickness = [d.mean - d.range for d in thickness_data]
    max_strata_thickness = [d.mean + d.range for d in thickness_data]

    all_routes = []
    all_routes_number = []

    # Set the initial routes.
    row = 0
    thickness_change = drillsample_data.get_thickness_change(row)

    for strat in range(num_units):
        if (strata_table[row, strat].path_exists):
            new_route = StrataRoute(num_units)
            new_route.add_first_position(strat, thickness_change)
            # Adding new route into the list.
            all_routes.append(new_route)

    print("Starting routes:")
    print([r.path for r in all_routes])

    row_max = num_rows
    print("row_max = ", row_max)

    litho_sequence_length = 1

    # Print the starting info.
    row = 0
    print("Processed row =", row, drillsample_data.rows[row].depth_from, drillsample_data.rows[row].lithos, len(all_routes))

    # Going through the strata table and generating the routes.
    for row in range(1, row_max):
        # Print the info.
        print("Processing row =", row, drillsample_data.rows[row].depth_from, drillsample_data.rows[row].lithos, end = "\r")

        # The drillhole lithos.
        # Note: we deliberately consider the full list of drillsample lithos instead of lithos- on the route.
        # Because considering the route lithos may lead to exponential growth of number of routes due to frequent unit change.
        current_lithos = drillsample_data.rows[row].lithos
        previous_lithos = drillsample_data.rows[row - 1].lithos

        thickness_change = drillsample_data.get_thickness_change(row)
        new_routes = []

        # Allowed strata units for a unit change.
        strata_allowed = [strat for strat in range(num_units) if (strata_table[row, strat].path_exists)]

        same_lithos = False
        if (set(current_lithos) == set(previous_lithos)):
            litho_sequence_length = litho_sequence_length + 1
            same_lithos = True
        else:
            litho_sequence_length = 1

        # Iterate over all routes.
        for route in all_routes:
            # The current strata index.
            strat0 = route.path[-1]
            current_thickness = route.current_thickness

            #--------------------------------------------------------------------
            # Check if we can go down in other stratas (and create new routes).
            #--------------------------------------------------------------------
            can_change = True

            # Add 'unit contacts inside the same litho' constraints.
            if (same_lithos):
                # Constrain the maximum number of unit contacts inside a litho sequence.
                if (route.num_unit_contacts_inside_litho >= spar.max_num_unit_contacts_inside_litho):
                    can_change = False
            else:
                # Reset.
                route.num_unit_contacts_inside_litho = 0

            # Apply unit thickness constraints.
            if (add_thickness_constraints):
                # Ignore thickness for the top unit
                if (len(route.path) > 1):
                    can_change = can_change and (current_thickness >= min_strata_thickness[strat0])

            if (can_change):
                # Strata units excluding the current one, and excluding the Cover (as we cannot change to Cover, but only can start from it).
                strata_list = [s for s in strata_allowed if (s != strat0 and unit_names[s] != "Cover")]

                # Applying the "maximum number of returns to a unit" constraint.
                apply_max_num_returns_constraint(route, strata_list, spar.max_num_returns_per_unit)

                # Apply unit topology constraints.
                if (add_topology_constraints):
                    apply_topology_constraints(graph, unit_names, strat0, strata_list, spar.ignore_unit_age)

                if (len(strata_list) != 0):
                    # Copy the route to create references to it below.
                    old_path = route.path.copy()

                    # Looking to which strata unit we can change.
                    for strat in strata_list:
                        # Making the new route.
                        new_route = StrataRoute(num_units)
                        # New path contains the reference to the old path, and the new route position.
                        # Note: we are not copying the full old path, but only store a reference to it to save memory.
                        new_route.path = [old_path, strat]
                        new_route.current_thickness = thickness_change
                        np.copyto(new_route.unit_visited, route.unit_visited)

                        # Count this unit as visited.
                        new_route.unit_visited[strat] += 1

                        # Processing the unit contact inside the same lithology sequence.
                        if (same_lithos):
                            new_route.num_unit_contacts_inside_litho = route.num_unit_contacts_inside_litho + 1
                        else:
                            new_route.num_unit_contacts_inside_litho = 0

                        # Adding new route into the list.
                        new_routes.append(new_route)

            #-----------------------------------------------------------------
            # Check if we can go down the same srata unit (if cannot -- remove the current route).
            #-----------------------------------------------------------------
            can_stay = True

            if (add_thickness_constraints):
                # Apply unit thickness constraints.
                can_stay = can_stay and (current_thickness < max_strata_thickness[strat0])

            path_exists = strata_table[row, strat0].path_exists

            can_stay = can_stay and path_exists

            # Processing the route.
            if (can_stay):
                # Adding new route position.
                route.path.append(strat0)
                route.current_thickness += thickness_change
            else:
            # Did not reach the end of a drillhole, and cannot go down the same unit.
                # Mark the route for removal.
                route.to_remove = True

        # Remove the routes marked for removal.
        all_routes = [route for route in all_routes if not route.to_remove]

        # Addig new routes.
        all_routes.extend(new_routes)

        # Update the number of routes.
        num_routes = len(all_routes)
        all_routes_number.append(num_routes)

        # Print the info.
        print("Processed row =", row, drillsample_data.rows[row].depth_from, drillsample_data.rows[row].lithos, num_routes)

        if (num_routes == 0):
            break

    #------------------------------------------------------------------
    # Adding the final number of routes.
    all_routes_number.append(len(all_routes))

    # Flatten the multilevel list of lists: convert [[[1,1],2,2],3,3] to [1,1,2,2,3,3].
    for route in all_routes:
        route.path = flatten(route.path)

    # Create the solution object.
    depth_data = drillsample_data.get_depth_data()
    solution = StrataSolution(all_routes, all_routes_number, unit_names, depth_data, strata_data.unit2dist)

    return solution
