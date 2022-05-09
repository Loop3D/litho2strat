'''
  This program generates a set of plausible stratigraphies with uncertainties, for a given drillhole lithology log.
  It uses map data for distance and topology constraints, and several free parameters describing the solution complexity (level of deformation) constraints.

  Author: Vitaliy Ogarko, vogarko@gmail.com
  The University of Western Australia
'''

import csv
import numpy as np
import matplotlib.pylab as pl
import networkx as nx
from dataclasses import dataclass, field
from typing import List
import tracemalloc
import os

# 'Returning to the same unit' constraints.
max_num_returns_per_unit = 1
#---------------------------------------------------------------------------
# Adding thickness constraints. (Requires unit thickness data).
add_thickness_constraints = False
#---------------------------------------------------------------------------
# Adding unit contacts topology (extracted from map data).
add_topology_constraints = False
# Ignore topology graph edge direction defining the unit age.
ignore_unit_age = True
#---------------------------------------------------------------------------
# The number of unit contacts inside the same litholgy sequence.
max_num_unit_contacts_inside_litho = 0
#---------------------------------------------------------------------------
# The number of nearest units (for distance constraints).
number_nearest_units = 2
#---------------------------------------------------------------------------
# Use the single closest unit for the top (first) lithology.
single_top_unit = True

#==============================================================================
# Missing lithos.
missing_lithos = set()

#==============================================================================
def fix_litho_name(litho):
    '''
    Convert the map lithology name to the 'CET lithology' name.
    '''
    # Convert things like metagranite to granite.
    if (litho[0:4] == 'meta'):
        litho = litho[4:]

    return litho

#==============================================================================
def get_unique_lithos_from_strat_data(strat_data):
    '''
    Returns a unique list of lithologies from strat_data.
    '''
    lithos = list()
    for unit_name in strat_data:
        for litho in strat_data[unit_name]:
            lithos.append(litho)
    unique_lithos = list(dict.fromkeys(lithos))

    return unique_lithos

#===========================================================================================
def add_unit_to_distance_map(unit_name, distance, lithos, litho2dist):
    '''
    Adds a unit with its lithologies to the distance map litho2dist.
    '''
    for litho in lithos:
        el = (distance, unit_name)
        if (litho in litho2dist):
            found_unit = False
            # Iterate over distance list.
            for tup in litho2dist[litho]:
                if (unit_name == tup[1]):
                # Found this unit in the list.
                    found_unit = True
                    if (distance < tup[0]):
                        # Update element to a smaller distance.
                        litho2dist[litho].remove(tup)
                        litho2dist[litho].append(el)
                    break
            if (not found_unit):
                litho2dist[litho].append(el)
        else:
            litho2dist[litho] = [el]
        # Sort the list by distance.
        litho2dist[litho].sort(key=lambda tup: tup[0])

#====================================================================================
def read_strat_data(dist_table_filename, drillhole_lithos, alternative_rock_names):
    '''
    Building lithologies list for every unit from csv file data.
    '''
    # Unit code column number.
    code_column = 1
    # Unit name column number.
    unitname_column = 2
    # Lithology list column number.
    lithos_column = 4
    # Distance column number.
    dist_column = 5
    # LITHNAME1 column.
    lithname1_column = 10
    # Description column.
    descript_column = 3

    # Converter from string to list.
    str2list = lambda x: x.strip("[]").replace("'", "").replace(" ", "").split(",")
    # To fix ? symbol in some names (gabbro?leucogabbro)
    str2list2 = lambda x: x.strip("[]").replace("'", "").replace(" ", "").replace("?", ",").split(",")

    # Reading the units table.
    strat_all = dict()

    # Lithology to distance and unitname mapping.
    litho2dist = dict()

    with open(dist_table_filename, 'r') as csvfile:
        # Reading the csv data.
        csvreader = csv.reader(csvfile, delimiter=',')
        # Skipping the header.
        next(csvreader)
        # Extracting the lighology list for every csv row (strata unit).
        for row in csvreader:
            # Extract the list of lithologies.
            lithos = str2list2(row[lithos_column])
            # Distance to the unit code.
            distance = float(row[dist_column])
            # Unit name.
            unit_name = row[unitname_column]

            # Convert the unitname to align it with format used in the topology graph.
            unit_name = unit_name.replace(" ", "_").replace(",", "_")

            #=========================================
            # Hard fixes for "mudstone" and "coal". 
            # TODO: Fix the original data.
            #=========================================
            lithname1 = row[lithname1_column]
            if ("mudstone" in lithname1):
                lithos.append("mudstone")

            description = row[descript_column]
            if ("coal" in description):
                lithos.append("coal")
                lithos.append("lignite")

            #=========================================
            # Fixing the litho names.
            #=========================================
            lithos = [fix_litho_name(l) for l in lithos]

            #=========================================
            # Adding the alternative litho names.
            #=========================================
            alt_lithos = []
            for litho in lithos:
                for alt_names in alternative_rock_names:
                    if litho in alt_names:
                        alt_lithos.extend(alt_names)

            lithos.extend(alt_lithos)

            # Remove duplicates.
            lithos = list(dict.fromkeys(lithos))

            #=========================================
            # Store the sorted distance map.
            #=========================================
            add_unit_to_distance_map(unit_name, distance, lithos, litho2dist)

            #=========================================
            # Adding the lithos to the dictionary (excluding the duplicates).
            if unit_name in strat_all:
                for litho in lithos:
                    if litho not in strat_all[unit_name]:
                        strat_all[unit_name].append(litho)
            else:
                strat_all[unit_name] = list(dict.fromkeys(lithos)) # Remove duplicates.

    print("The total number of units: " + str(len(strat_all)))

    #=====================================================================
    unique_lithos = get_unique_lithos_from_strat_data(strat_all)

    print(unique_lithos)
    print("The total number of lithologies: " + str(len(unique_lithos)))

    #=====================================================================
    # Filter units based on the drillhole lithologies: 
    # Remove lithologies, that are not present in the drillhole data.
    # Remove the units that do not contain the drillhole lithos.
    #=====================================================================
    strat_filtered = dict()
    for unit_name in strat_all:
        for litho in strat_all[unit_name]:
            # Only add lithologies that are present in drillhole data.
            if (litho in drillhole_lithos):
                if (unit_name in strat_filtered):
                    strat_filtered[unit_name].append(litho)
                else:
                    strat_filtered[unit_name] = [litho]

    print("The number of filtered units: " + str(len(strat_filtered)))

    #=====================================================================
    unique_lithos = get_unique_lithos_from_strat_data(strat_filtered)

    print(unique_lithos)
    print("The filtered number of lithologies: " + str(len(unique_lithos)))

    #=====================================================================
    # Filter units based on the distance from drillhole.
    #=====================================================================
    strat_dist = dict()
    for unit_name in strat_filtered:
        for litho in strat_filtered[unit_name]:
            # Sorted distance list for this lithology.
            dist_list = litho2dist[litho]

            # Consider only N closest codes.
            for el in dist_list[:number_nearest_units]:
                unit_name_nearest = el[1]
                if (unit_name == unit_name_nearest):
                    if (unit_name in strat_dist):
                        strat_dist[unit_name].append(litho)
                    else:
                        strat_dist[unit_name] = [litho]
                    break

    print("The number of filtered (by distance) units: " + str(len(strat_dist)))

    return strat_dist, litho2dist

#========================================================================================
def test_column_exist(column_name, fieldnames):
    '''
    Tests if the column name is present in the fieldnames list.
    '''
    if (column_name not in fieldnames):
        print("Error: The column name is not found in drillsample data header:", column_name)
        exit()

#========================================================================================================
@dataclass
class DrillSampleHeader:
    '''
    Contains the names of drillsample header data columns.
    '''
    depth_from: str
    depth_to: str
    lithos: str
    scores: str

#========================================================================================================
@dataclass
class DrillSampleDataRow:
    '''
    Contains the names of drillsample data row.
    '''
    depth_from: float = 0.
    depth_to: float = 0.
    lithos: List[str] = field(default_factory=list)
    scores: List[int] = field(default_factory=list)

#========================================================================================================
def read_drillsample_data(filename, header, ignore_list):
    '''
    Reading drill sample data from csv file.
    '''
    all_data = []
    all_lithos = set()

    with open(filename, 'r') as csvfile:
        reader = csv.DictReader(csvfile, delimiter=',')

        # The header column names.
        fieldnames = reader.fieldnames

        # Sanity check.
        test_column_exist(header.depth_from, fieldnames)
        test_column_exist(header.depth_to, fieldnames)
        test_column_exist(header.lithos, fieldnames)
        test_column_exist(header.scores, fieldnames)

        # Extracting the data for every csv row.
        for row in reader:
            # Extract the data columns.
            data = DrillSampleDataRow()
            data.depth_from = float(row[header.depth_from])
            data.depth_to = float(row[header.depth_to])

            if (row[header.lithos] == ''):
                # No data - skip the row.
                continue

            lithos = row[header.lithos].split(", ")
            scores = [int(s) for s in row[header.scores].split(", ")]

            # Sanity check.
            if (len(lithos) != len(scores)):
                print("Error: number of lithos differs from the number of scores for:", row)
                exit()

            # Filter the lithos and scores based on the Ignore list of lithos.
            for index, litho in enumerate(lithos):
                if (litho not in ignore_list):
                    data.lithos.append(litho)
                    data.scores.append(scores[index])
                    # Gather all unique lithos.
                    all_lithos.add(litho)

            all_data.append(data)

    print(all_lithos)
    print("The number of drillhole lithologies: " + str(len(all_lithos)))

    return all_data

#==============================================================================
def group_drillhole_litho_sequence(data):
    '''
    Group the litho sequence by name (inside the drillhole data), leaving at most N in each group,
    corresponding to number of contacts inside the litho sequence.
    '''
    data_grouped = []
    N = max_num_unit_contacts_inside_litho + 1
    current_litho = ""
    from_depth = float(data[0][0])
    to_depth = float(data[0][1])

    data_mod = [d for d in data]
    data_mod.append([0, 0, "last_litho_name_for_following_calculation"])

    for index, row in enumerate(data_mod):

        prev_litho = current_litho
        current_litho = row[2]

        if (current_litho != prev_litho and index > 0):
        # Change of litho name.
            total_thickness = to_depth - from_depth
            local_thickness = total_thickness / float(N)
            litho = prev_litho

            print("Grouping lithos for:", from_depth, to_depth, litho)

            # Generate local grouped lithos.
            for i in range(N):
                from_depth_local = from_depth + float(i) * local_thickness
                to_depth_local = from_depth_local + local_thickness

                row_grouped = list(range(3))
                row_grouped[0] = from_depth_local
                row_grouped[1] = to_depth_local
                row_grouped[2] = litho

                data_grouped.append(row_grouped)

            # Update the starting depth for the following grouped lithos.
            from_depth = float(row[0])
        to_depth = float(row[1])

    return data_grouped

#==============================================================================
def read_thickness_data(filename):
    '''
    Reading thickness data from csv file.
    '''
    data = []
    with open(filename, 'r') as csvfile:
        # Reading the csv data.
        csvreader = csv.reader(csvfile, delimiter=',')
        # Skipping the header.
        next(csvreader)
        # Extracting the data for every csv row.
        for row in csvreader:
            thickness = np.array([0, 0], dtype='f')
            thickness[0] = float(row[1]) # "thickness_mean".
            thickness[1] = float(row[2]) # "thickess_range".
            data.append(thickness)
    return data

#==============================================================================
def read_topology_data(topology_filename):
    '''
    Read topology data (gml format graph) from a file.
    '''
    print('Importing the graph data...')

    # Import graph from a file.
    Gf = nx.read_gml(topology_filename)

    # Modify the graph to have node names = unit names.
    for node in Gf.nodes():
        unit_name = Gf.nodes[node]['LabelGraphics']['text']
        mapping = {node:unit_name}
        Gf = nx.relabel_nodes(Gf, mapping)

    if (ignore_unit_age):
    # Ignore graph edge direction defining the unit age.
        Gf = Gf.to_undirected()

    print('Importing completed.')

    return Gf

#========================================================================================================
@dataclass
class StrataTableElement:
    '''
    The element of a strata table.
    '''
    path_exists: bool = False
    lithos: List[str] = field(default_factory=list)

#==============================================================================
def generate_strata_table(drillsample_data, strat_data, litho2dist):
    '''
    Generates the stratigraphic table, and unit names list.
    '''
    num_rows = len(drillsample_data)
    num_units = len(strat_data)
    unit_names = get_unit_names(strat_data)

    print("num_rows (before) = ", num_rows)
    print("num_units = ", num_units)

    strata_table = np.ndarray(shape=[num_rows, num_units], dtype=object)

    # Initialize with default objects.
    strata_table.flat = [StrataTableElement() for _ in strata_table.flat]

    new_row_index = 0

    for row in drillsample_data[:]:
        #litho = row.litho
        any_litho_found = False

        # Loop over lithos in the current drillsample row.
        for litho in row.lithos[:]:
            litho_found = False

            if (new_row_index == 0 and single_top_unit):
            # Use only the closest unit for the top lithology.
                if (litho in litho2dist):
                    # Sorted distance list for this lithology.
                    dist_list = litho2dist[litho]

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

                    if (closest_unit_distance > 0):
                        print("WARNING: The closest distance to the top unit > 0! Dist =", closest_unit_distance)

                    for unit_name in strat_data:
                        if (unit_name == closest_unit or (add_cover and unit_name == 'Cover')):
                            litho_found = True
                            unit_index = unit_names.index(unit_name)
                            strata_table[new_row_index, unit_index].path_exists = True
                            strata_table[new_row_index, unit_index].lithos.append(litho)
            else:
                for unit_name in strat_data:
                    if (litho in strat_data[unit_name]):
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
            drillsample_data.remove(row)
        else:
            new_row_index += 1

    num_rows = len(drillsample_data)
    print("num_rows (after) = ", num_rows)

    if (num_rows == 0):
        print('No data left!!!')
        return np.ndarray(shape=[num_rows, num_units], dtype=object)

    # Remove rows due to missing lithologies.
    strata_table = strata_table[0:num_rows, :]

    return strata_table

#==============================================================================
'''
A class for storing the stratigraphic route.
'''
class StrataRoute:
    __slots__ = 'to_remove', 'path', 'current_thickness', 'unit_visited', 'num_unit_contacts_inside_litho'

    def __init__(self):
        # Flag for removal.
        self.to_remove = False

    # Adding the first position to the route.
    # This method essentially initializes the route.
    def add_first_position(self, strat, thickness_change, num_units):
        # Unit index for every drillhole data.
        self.path = [strat]
        # The thickness of the last strata unit.
        self.current_thickness = thickness_change
        # Containts the number of times each unit was visited.
        self.unit_visited = np.zeros((num_units), dtype=int)
        # Mark this unit as 'visited'.
        self.unit_visited[strat] += 1
        # The number of unit contacts inside the same lithology sequence.
        self.num_unit_contacts_inside_litho = 0

    def __str__(self):
        return str(self.path)

    def __repr__(self):
        return str(self.path)

    def get_strata_sequence(self):
        '''
        Returns a sequence of strata units in the route excluding consequtive dublicates.
        '''
        return tuple([v for i, v in enumerate(self.path) if i == 0 or v != self.path[i - 1]])

#==============================================================================
def get_thickness_change(drillsample_data, row):
    '''
    Returns a thickness change for a given row in the drillhole sample.
    '''
    # "To" - "From"
    return drillsample_data[row].depth_to - drillsample_data[row].depth_from

#==============================================================================
def get_min_strata_thickness(thickness_data):
    # "thickness_mean" - "thickess_range"
    return [data[0] - data[1] for data in thickness_data]

#==============================================================================
def get_max_strata_thickness(thickness_data):
    # "thickness_mean" + "thickess_range"
    return [data[0] + data[1] for data in thickness_data]

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
def apply_max_num_returns_constraint(route, strata_list):
    '''
    Apply the "maximum number of returns to a unit" constraint:
        remove from the input unit list the units where the route cannot return anymore.
    '''
    # Apply the "max number of returns" constraint.
    for strat in strata_list[:]:
        if (route.unit_visited[strat] - 1 >= max_num_returns_per_unit):
        # Reached the maximum numer of local returns (to this unit).
            strata_list.remove(strat)

#==============================================================================
def apply_topology_constraints(graph, unit_names, strat0, strata_list):
    '''
    Apply unit topology (connectivity) constraints:
        remove from the input unit list the units not connected to a given one (strat0).
    '''
    if (unit_names[strat0] in graph.nodes()):
        for strat in strata_list[:]:
            if (not graph.has_edge(unit_names[strat0], unit_names[strat])):
                # Units are not connected. Skip this unit.
                strata_list.remove(strat)

#==============================================================================
def get_unit_names(strat_data):
    '''
    Defines the mapping between the unit index and unit name.
    '''
    unit_names = []
    for unit_name in strat_data:
        if (unit_name == 'Cover'):
            # Map the Cover's index to zero.
            unit_names.insert(0, unit_name)
        else:
            unit_names.append(unit_name)

    return unit_names

#==============================================================================
def get_drillhole_lithos(drillsample_data):
    '''
    Returns the drillhole lithologies from drillsample data.
    '''
    all_lithos = set()
    for row in drillsample_data:
        lithos = row.lithos
        all_lithos.update(lithos)
    return all_lithos

#==============================================================================
def generate_strat_routes(strat_data, litho2dist, drillsample_data, thickness_data, graph):
    '''
    Generating stratigraphic routes.
    '''
    # Generating the table of possible strata paths.
    strata_table = generate_strata_table(drillsample_data, strat_data, litho2dist)

    # Unit index to unit name mapping.
    unit_names = get_unit_names(strat_data)

    num_rows = strata_table.shape[0]
    num_units = strata_table.shape[1]

    # Extract strata thikcness to lists (faster data structures).
    min_strata_thickness = get_min_strata_thickness(thickness_data)
    max_strata_thickness = get_max_strata_thickness(thickness_data)

    all_routes = []
    all_routes_number = []

    # Set the initial routes.
    row = 0
    thickness_change = get_thickness_change(drillsample_data, row)

    for strat in range(num_units):
        if (strata_table[row, strat].path_exists):
            new_route = StrataRoute()
            new_route.add_first_position(strat, thickness_change, num_units)
            # Adding new route into the list.
            all_routes.append(new_route)

    print("Starting routes:")
    print(all_routes)

    row_max = num_rows
    print("row_max = ", row_max)

    litho_sequence_length = 1

    # Going through the strata table and generating the routes.
    for row in range(1, row_max):
        # The drillhole lithos.
        # Note: we deliberately consider the full list of drillsample lithos instead of lithos- on the route.
        # Because considering the route lithos may lead to exponential growth of number of routes due to frequent unit change.
        current_lithos = drillsample_data[row].lithos
        previous_lithos = drillsample_data[row - 1].lithos

        thickness_change = get_thickness_change(drillsample_data, row)
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
                if (route.num_unit_contacts_inside_litho >= max_num_unit_contacts_inside_litho):
                    can_change = False

                #-------------------------------------------------------------------------------------
                # Exluding superficial (duplicate) routes with unit contacts inside a litho sequence.
                #-------------------------------------------------------------------------------------
                # Example:
                # litho_sequence_length = 3
                # n = 1,2,3
                #     A-A-A
                # for n = 3, allow for unit change only if have 1 "unit contacts inside litho" already (between 1 and 2),
                # so when num_unit_contacts_inside_litho >= litho_sequence_length - 2
                # This way we skip (duplicate) cases with unit change at n = 3 for routes with zero "unit contacts inside litho".
                # All those cases should already have been built at n = 2.
                #-------------------------------------------------------------------------------------
                if (route.num_unit_contacts_inside_litho < litho_sequence_length - 2):
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
                apply_max_num_returns_constraint(route, strata_list)

                # Apply unit topology constraints.
                if (add_topology_constraints):
                    apply_topology_constraints(graph, unit_names, strat0, strata_list)

                if (len(strata_list) != 0):
                    # Copy the route to create references to it below.
                    old_path = route.path.copy()

                    # Looking to which strata unit we can change.
                    for strat in strata_list:
                        # Making the new route.
                        new_route = StrataRoute()
                        # New path contains the reference to the old path, and the new route position.
                        # Note: we are not copying the full old path, but only store a reference to it to save memory.
                        new_route.path = [old_path, strat]
                        new_route.current_thickness = thickness_change
                        new_route.unit_visited = np.array(route.unit_visited)

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
        print("ROW = ", row, drillsample_data[row].lithos, num_routes)

        if (num_routes == 0):
            break

    #------------------------------------------------------------------
    # Adding the final number of routes.
    all_routes_number.append(len(all_routes))

    # Flatten the multilevel list of lists: convert [[[1,1],2,2],3,3] to [1,1,2,2,3,3].
    for route in all_routes:
        route.path = flatten(route.path)

    return all_routes, all_routes_number

#==============================================================================
def write_routes_to_file(filename, drillsample_data, all_routes):
    '''
    Writing stratigraphic routes to file.
    '''
    f = open(filename, "w")
    num_rows = len(drillsample_data)
    for row in range(num_rows):
        depth = drillsample_data[row].depth_from
        f.write("%f " % depth)
        # Calculate the number of unique strata for this depth.
        unique_units = set([])
        for route in all_routes:
            unique_units.add(route.path[row])
        f.write("%d " % len(unique_units))

        for route in all_routes:
            f.write("%d " % route.path[row])
        f.write("\n")
    f.close()

#==============================================================================
def print_unique_routes(all_routes, num_print_paths):
    '''
    Print all unique routes (i.e., with unique strata sequence).
    '''
    unique_routes = set([])
    for route in all_routes:
        unique_routes.add(route.get_strata_sequence())

    print("Number of unique routes = ", len(unique_routes))
    if (num_print_paths > 0):
        num = 0
        for route in unique_routes:
            num += 1
            print(route)
            if (num >= num_print_paths):
                break

#=============================================================================
def plot_routes(drillsample_data, routes, strat_distr):
    '''
    Plot and display the routes.
    '''
    print("Plotting the routes...")

    # Using the "From" column.
    x_data = [d.depth_from for d in drillsample_data]

    for route in routes:
        pl.plot(x_data, route.path, '.-')

    pl.xlabel('Depth')
    pl.ylabel('Strata unit index')
    pl.show()

    #------------------------------------------
    # Plot route probabilities.
    num_rows = len(routes[0].path)

    for route in routes:
        route_proba = np.zeros(num_rows)
        for row in range(num_rows):
            unit_index = route.path[row]
            route_proba[row] = strat_distr[row, unit_index]
        pl.plot(x_data, route_proba, '.-')

    pl.xlabel('Depth')
    pl.ylabel('Probability')
    pl.show()

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

#=============================================================================
def plot_unit_probabilities(all_routes, drillsample_data, unit_names):
    '''
    Generate a plot with probability of occurence for each unit.
    '''
    if (len(all_routes) == 0):
        return

    num_units = len(unit_names)

    # Building the distribution of unit presence at every depth.
    strat_distr = get_strat_distr(all_routes, num_units)

    #------------------------------------------
    # Plot distribution of the route scores (based on path probability).
    route_scores = get_route_scores(all_routes, strat_distr)

    title_params = 'Max returns per unit = ' + str(max_num_returns_per_unit)

    pl.hist(route_scores, bins = 50)
    pl.title(title_params)
    pl.xlabel('Route score')
    pl.ylabel('Frequency')
    pl.show()

    #------------------------------------------
    # Top scores.
    indexes_max = np.argsort(-route_scores) # A minus here to have largest to smallest score order.
    ntop = 10
    print('Top indexes: ', indexes_max[0:ntop])
    print('Top scores: ', route_scores[indexes_max[0:ntop]])

    index_max = indexes_max[0]

    #------------------------------------------
    # Print the most probable routes.
    plot_routes(drillsample_data, [all_routes[i] for i in indexes_max[0:ntop]], strat_distr)

    #------------------------------------------
    # Print if there are multiple best routes.
    multiple_best_routes = False
    if (len(indexes_max) > 1):
        if (route_scores[indexes_max[0]] == route_scores[indexes_max[1]]):
            multiple_best_routes = True

    print("Multiple best routes: ", multiple_best_routes)

    #------------------------------------------
    # Generating the plots.
    # Increasing the figure size.
    pl.rcParams["figure.figsize"] = (12.8, 9.6) # Default size = (6.4, 4.8)

    num_units_nonempty = 0
    for i in range(num_units):
        if (sum(strat_distr[:, i]) != 0):
            num_units_nonempty += 1
    fig, axs = pl.subplots(num_units_nonempty, sharey=True, squeeze=False)

    fig.suptitle('Probability of occurrence for every unit.', y=0.96)

    num_rows = len(all_routes[0].path)

    #-------------------------------------------------------------
    # Adding the "From" and "To" depths for visualisation.
    x_data = []
    for d in drillsample_data[0:num_rows]:
        x_data.append(d.depth_from)
        x_data.append(d.depth_to)

    # Duplicate each value, as the probability is the same between the "From" and "To" depths.
    strat_distr = np.repeat(strat_distr, 2, axis=0)
    #-------------------------------------------------------------

    j = 0
    for i in range(num_units):
        if (sum(strat_distr[:, i]) == 0):
            # Skip empty units.
            continue

        # Plot lines.
        axs[j, 0].plot(x_data, strat_distr[:, i], zorder=1, c='blue')

        # Set red color for zero data.
        color = ['red' if p <= 0 else 'blue' for p in strat_distr[:, i]]

        # Plot dots.
        axs[j, 0].scatter(x_data, strat_distr[:, i], s=5, c=color, zorder=2)

        axs[j, 0].set_title(unit_names[i], size=9, y=0.91)
        axs[j, 0].set_ylabel(str(j))

        if (i != num_units - 1):
            # Hide tick labels.
            axs[j, 0].set_xticklabels([])

        # Add vertical lines.
        axs[j, 0].xaxis.grid(True)

        j += 1
 
    #pl.tight_layout()
 
    pl.xlabel('Depth')
    pl.show()

#=============================================================================
def add_cover_unit(unit_name, filename, strat_data, litho2dist):
    '''
    Adds the custom unit with its lithologies read from file.
    '''
    with open(filename) as f:
        lithos = f.read().splitlines()

    # Adding unit data to strat_data.
    if (unit_name in strat_data):
        strat_data[unit_name].append(lithos)
    else:
        strat_data[unit_name] = lithos

    # Adding unit data to the distance map litho2dist.
    distance = 0.
    add_unit_to_distance_map(unit_name, distance, lithos, litho2dist)

#=============================================================================
def read_ignore_list(filename):
    '''
    Read the drillhole items to ignore.
    '''
    with open(filename) as f:
        items = f.read().splitlines()

    return items

#=============================================================================
def read_alternative_rock_names(filename):
    '''
    Read the alternative rock names (synonyms).
    Returns a list with lists of alternative names.
    '''
    with open(filename) as f:
        items = f.read().splitlines()

    alternative_rock_names = []

    # Building a dictionary.
    for item in items:
        names_list = item.split(", ")
        alternative_rock_names.append(names_list)

    return alternative_rock_names

#=============================================================================
def generate_missing_lithos():
    '''
    Generates a list of missing drillhole lithologies from unit data in all data files.
    '''

    # The Cover unit data file.
    cover_unit_filename = "data/real/cover_unit.txt"

    # The Ignore items list.
    ignore_list_filename = "data/real/ignore_list.txt"

    # Alternative rock names file.
    alternative_rock_names_file = "data/real/alternative_rock_names.txt"

    directory = "data/real/dist_files/litho_tables"
    #directory = "data/real/dh_files/litho_tables"

    # Drillsample data column names.
    drillsample_header = DrillSampleHeader('Fromdepth', 'Todepth', 'Lithologies', 'Scores')

    # Read the drillhole ignore items list.
    ignore_list = read_ignore_list(ignore_list_filename)

    # Read the alternative rock names.
    alternative_rock_names = read_alternative_rock_names(alternative_rock_names_file)

    counter = 0
    for file in os.listdir(directory):
        counter = counter + 1
        print("PROCESSING FILE NUMBER =", counter)

        filename = os.fsdecode(file)
        collarID = int(filename[6:-4])
        print(collarID)

        drillsample_filename = "data/real/dist_files/litho_tables_V2/litho_" + str(collarID) + ".csv"
        dist_table_filename = "data/real/dist_files/dist_tables/100_500k_map_near_" + str(collarID) + ".csv"

        #drillsample_filename = "data/real/dh_files/litho_tables/litho_" + str(collarID) + ".csv"
        #dist_table_filename = "data/real/dh_files/dist_tables/100k_map_near_" + str(collarID) + ".csv"

        # Drill sample data.
        drillsample_data = read_drillsample_data(drillsample_filename, drillsample_header, ignore_list)

        # Unit lithologies and distance data.
        drillhole_lithos = get_drillhole_lithos(drillsample_data)
        strat_data, litho2dist = read_strat_data(dist_table_filename, drillhole_lithos, alternative_rock_names)

        # Read the Cover unit lithologies.
        add_cover_unit("Cover", cover_unit_filename, strat_data, litho2dist)

        # Generating the table of possible strata paths.
        strata_table = generate_strata_table(drillsample_data, strat_data, litho2dist)

    print("Missing lithologies list:")
    print("Number lithos missing =", len(missing_lithos))
    print(missing_lithos)

#=============================================================================
def main():
    print('Started litho2strat')

    #generate_missing_lithos()
    #exit()

    # Topology file.
    #topology_filename = "data/real/ASUD_strat.gml"
    topology_filename = "data/real/ASUD_strat2.gml"

    # The Cover unit data file.
    cover_unit_filename = "data/real/cover_unit.txt"

    # The Ignore items list file.
    ignore_list_filename = "data/real/ignore_list.txt"

    # Alternative rock names file.
    alternative_rock_names_file = "data/real/alternative_rock_names.txt"

    #----------------------------------------------------------------------------
    # # Mark's data.
    # collarID = 548917
    # drillsample_filename = "data/real/dh_files/litho_tables/litho_" + str(collarID) + ".csv"
    # dist_table_filename = "data/real/dh_files/dist_tables/100k_map_near_" + str(collarID) + ".csv"
    #
    # # Drillsample data column names.
    # drillsample_header = DrillSampleHeader('FromDepth', 'ToDepth', 'CET_Litho', 'Scores')

    #----------------------------------------------------------------------------
    # Mark's data with known solutions.
    #----------------------------------------------------------------------------
    # TODO: Discuss with Mark - we have here too long Cover...
    #collarID = 1209857
    #collarID = 353386
    #collarID = 2182301
    #collarID = 810340

    # Confirmed results (using 1 closest unit & single top unit).
    collarID = 2182336
    #collarID = 2182335
    #collarID = 2182340
    #collarID = 2182339
    #collarID = 2182338
    #collarID = 2182335
    #collarID = 2182334

    # TODO: Discuss with Mark - we have here conglomerate, which is rock, and then gravel, which we define as 'always Cover'. Thus we have the Cover below the rock here.
    #collarID = 1209855

    drillsample_filename = "data/real/dist_files/litho_tables_V2/litho_" + str(collarID) + ".csv"
    dist_table_filename = "data/real/dist_files/dist_tables/100_500k_map_near_" + str(collarID) + ".csv"

    # Synthetic test.
    #drillsample_filename = "data/real/test/litho_1.csv"
    #dist_table_filename = "data/real/test/map_1.csv"

    # Drillsample data column names.
    drillsample_header = DrillSampleHeader('Fromdepth', 'Todepth', 'Lithologies', 'Scores')

    #--------------------------------------------------------------
    # Reading the input data.
    #--------------------------------------------------------------
    # Read the drillhole ignore items list.
    ignore_list = read_ignore_list(ignore_list_filename)

    # Read drill sample data.
    drillsample_data = read_drillsample_data(drillsample_filename, drillsample_header, ignore_list)

    # Group the litho sequence inside the drillhole data.
    # Commented because it leads to many additional routes when splitting the existing single lithos.
    #drillsample_data = group_drillhole_litho_sequence(drillsample_data)

    # Read the alternative rock names.
    alternative_rock_names = read_alternative_rock_names(alternative_rock_names_file)

    # Read unit lithologies and distance data.
    drillhole_lithos = get_drillhole_lithos(drillsample_data)
    strat_data, litho2dist = read_strat_data(dist_table_filename, drillhole_lithos, alternative_rock_names)

    # Read the Cover unit lithologies.
    add_cover_unit("Cover", cover_unit_filename, strat_data, litho2dist)

    # Read thickness data.
    thickness_data = []
    if (add_thickness_constraints):
        thickness_data = read_thickness_data(thickness_filename)

    # Read topology data.
    graph = nx.Graph()
    if (add_topology_constraints):
        graph = read_topology_data(topology_filename)
        # Sanity check: check that strata units exist in the graph.
        for unit_name in strat_data:
            if unit_name not in graph.nodes():
                print("WARNING: Not found graph unit: ", unit_name)

    #--------------------------------------------------------------
    # Generating the stratigraphies.
    #--------------------------------------------------------------
#    tracemalloc.start()

    all_routes, all_routes_number = generate_strat_routes(strat_data, litho2dist, drillsample_data, thickness_data, graph)

    print("Total number of routes = ", len(all_routes))

#    current, peak = tracemalloc.get_traced_memory()
#    print("Current memory usage is {} MB; Peak was {} MB".format(current / 10**6, peak / 10**6))

    #--------------------------------------------------------------
    # Plot the results.
    #--------------------------------------------------------------
    # Plot the number of processed routes at each row.
    pl.xlabel('Row number')
    pl.ylabel('Number of routes')
    pl.plot(all_routes_number)
    pl.show()

    # Print all unique routes (i.e., unique strata sequence).
    print_unique_routes(all_routes, 10)

    # Write results to the file.
    #write_routes_to_file("strata.txt", drillsample_data, all_routes)

    # Plot unit probabilities.
    unit_names = get_unit_names(strat_data)
    plot_unit_probabilities(all_routes, drillsample_data, unit_names)

#=============================================================================
if __name__ == "__main__":
    main()

