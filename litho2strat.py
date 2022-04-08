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
import tracemalloc
import os

# 'Returning to the same unit' constraints.
max_num_returns_per_unit = 1
#---------------------------------------------------------------------------
# Adding thickness constraints. (Requires unit thickness data).
add_thickness_constraints = False
#---------------------------------------------------------------------------
# Adding unit contacts topology (extracted from map data).
add_topology_constraints = True
# Ignore topology graph edge direction defining the unit age.
ignore_unit_age = True
#---------------------------------------------------------------------------
# The number of unit contacts inside the same litholgy sequence.
max_num_unit_contacts_inside_litho = 0
#---------------------------------------------------------------------------
# The number of nearest units (for distance constraints).
number_nearest_units = 1
#---------------------------------------------------------------------------
# Use the single closest unit for the top (first) lithology.
single_top_unit = True

#==============================================================================
# List of all drillhole lithologies.
all_lithos = []
# Lithology to distance and unitname mapping.
litho2dist = dict()
# Unit names (filtered).
unit_names = []
# Missing lithos.
missing_lithos = set()

# Converter from string to list.
str2list = lambda x: x.strip("[]").replace("'", "").replace(" ", "").split(",")
str2list2 = lambda x: x.strip("[]").replace("'", "").replace(" ", "").replace("?", ",").split(",") # To fix ? symbol in some names (gabbro?leucogabbro)

#==============================================================================
def fix_litho_name(litho):
    '''
    Convert the map lithology name to the 'CET lithology' name.
    '''

    if (litho == 'iron-formation'):
        litho = 'banded_iron_formation'

    if (litho == 'ultramafic'):
        litho = 'ultramafic-rock'

    # Convert things like metagranite to granite.
    if (litho[0:4] == 'meta'):
        litho = litho[4:]

    return litho

#==============================================================================
def get_unique_lithos_from_strat_data(strat_data):
    lithos = list()
    for unit_name in strat_data:
        for litho in strat_data[unit_name]:
            lithos.append(litho)
    unique_lithos = list(dict.fromkeys(lithos))

    return unique_lithos

#==============================================================================
def read_strat_data(dist_table_filename):
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

    # Reading the units table.
    strat_all = dict()
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
            # Fix for "mudstone". 
            # TODO: Add mudstone to "lithos" column in the original data.
            #=========================================
            lithname1 = row[lithname1_column]
            if ("mudstone" in lithname1):
                lithos.append("mudstone")

            #=========================================
            # Fixing the litho names.
            #=========================================
            lithos2 = list()
            for litho in lithos:
                litho2 = fix_litho_name(litho)
                lithos2.append(litho2)

            lithos = [l for l in lithos2]

            #=========================================
            # Store the sorted distance map.
            #=========================================
            for litho in lithos:
                if (litho in all_lithos):
                # Lithology is present in drillhole data.
                    el = (distance, unit_name)
                    if (litho in litho2dist):
                        found_unit = False
                        # Iterate over distance list.
                        for tup in litho2dist[litho]:
                            if (unit_name == tup[1]):
                            # Found this unit in the list.
                                found_unit = True
                                if (distance < tup[0]):
                                    # Update element to smaller distance.
                                    litho2dist[litho].remove(tup)
                                    litho2dist[litho].append(el)
                                break
                        if (not found_unit):
                            litho2dist[litho].append(el)
                    else:
                        litho2dist[litho] = [el]
                    # Sort the list by distance.
                    litho2dist[litho].sort(key=lambda tup: tup[0])

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
            if (litho in all_lithos):
            # Lithology is present in drillhole data.
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

    return strat_dist

#==============================================================================
def read_drillsample_data(filename, column_from, column_to, column_litho):
    '''
    Reading drill sample data from csv file.
    '''
    data = []
    with open(filename, 'r') as csvfile:
        # Reading the csv data.
        csvreader = csv.reader(csvfile, delimiter=',')
        # Skipping the header.
        next(csvreader)
        # Extracting the data for every csv row.
        for row in csvreader:
            # Extract the data columns.
            row_formatted = list(range(3))
            row_formatted[0] = row[column_from] # FromDepth
            row_formatted[1] = row[column_to] # ToDepth
            row_formatted[2] = row[column_litho] # CET_Litho
            data.append(row_formatted)

            # Add lithology to the global list of all lithos.
            litho_name = row_formatted[2]
            if (litho_name not in all_lithos):
                all_lithos.append(litho_name)

    print(all_lithos)
    print("The number of drillhole lithologies: " + str(len(all_lithos)))

    return data

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

    return Gf

#==============================================================================
def generate_strata_table(drillsample_data, strat_data):
    '''
    Generates the stratigraphic table, and unit names list.
    '''
    num_rows = len(drillsample_data)
    num_units = len(strat_data)

    print("num_rows (before) = ", num_rows)
    print("num_units = ", num_units)

    strata_table = np.full((num_rows, num_units), False)
    new_row_index = 0

    for row in drillsample_data[:]:
        litho = row[2]
        litho_found = False

        if (new_row_index == 0 and single_top_unit):
        # Use only the closest unit for the top lithology.
            if (litho in litho2dist):
                # Sorted distance list for this lithology.
                dist_list = litho2dist[litho]
                closest_unit = dist_list[0][1]
                closest_unit_distance = dist_list[0][0]
    
                if (closest_unit_distance > 0):
                    print("WARNING: The closest distance to the top unit > 0! Dist =", closest_unit_distance)
    
                unit_index = 0
                for unit_name in strat_data:
                    if (unit_name == closest_unit):
                        litho_found = True
                        strata_table[new_row_index, unit_index] = True
                    unit_index += 1
        else:
            unit_index = 0
            for unit_name in strat_data:
                if (litho in strat_data[unit_name]):
                    litho_found = True
                    strata_table[new_row_index, unit_index] = True
                unit_index += 1

        if (not litho_found):
        # Drillhole lithology not found in units.
            print("WARNING: Drillhole lithology not found in units data: ", litho)
            missing_lithos.add(litho)

            # Treat this as "no data".
            drillsample_data.remove(row)
        else:
            new_row_index += 1

    num_rows = len(drillsample_data)
    print("num_rows (after) = ", num_rows)

    if (num_rows == 0):
        print('No data left!!!')
        return np.full((0, 0), False)

    # Remove rows due to missing lithologies.
    strata_table = strata_table[0:num_rows, :]

    # Map unit index to unit name.
    for unit_name in strat_data:
        unit_names.append(unit_name)

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
    return float(drillsample_data[row][1]) - float(drillsample_data[row][0])

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
def generate_strat_routes(strat_data, drillsample_data, thickness_data, graph):
    '''
    Generating stratigraphic routes.
    '''
    # Generating the table of possible strata paths.
    strata_table = generate_strata_table(drillsample_data, strat_data)

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
        if (strata_table[row, strat]):
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
        current_litho = drillsample_data[row][2]
        current_litho_index = all_lithos.index(current_litho)
        previous_litho = drillsample_data[row - 1][2]

        num_routes = len(all_routes)
        all_routes_number.append(num_routes)

        print("ROW = ", row - 1, previous_litho, num_routes)
        if (num_routes == 0):
            break

        thickness_change = get_thickness_change(drillsample_data, row)
        new_routes = []

        # Allowed strata units for a unit change.
        strata_allowed = [strat for strat in range(num_units) if (strata_table[row, strat])]

        if (current_litho == previous_litho):
            litho_sequence_length = litho_sequence_length + 1
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
            if (current_litho == previous_litho):
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
                # Strata units excluding the current one.
                strata_list = [s for s in strata_allowed if (s != strat0)]

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
                        if (current_litho == previous_litho):
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

            path_exists = strata_table[row, strat0]

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

        if (row == row_max - 1):
            # Print info for the last row.
            print("ROW = ", row, current_litho, len(all_routes))

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
        depth = float(drillsample_data[row][0])
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
    x_data = [float(d[0]) for d in drillsample_data]

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
def plot_unit_probabilities(all_routes, drillsample_data, num_units):
    '''
    Generate a plot with probability of occurence for each unit.
    '''
    if (len(all_routes) == 0):
        return

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

    # Using the "From" column.
    x_data = [float(d[0]) for d in drillsample_data[0:num_rows]]

    j = 0
    for i in range(num_units):
        if (sum(strat_distr[:, i]) == 0):
            # Skip empty units.
            continue

        # Plot lines.
        axs[j, 0].plot(x_data, strat_distr[:, i], zorder=1, c='blue')

        # Set red color for zero data.
        color = ['red' if d <= 0 else 'blue' for d in strat_distr[:, i]]

        # Set green color for the route with the highest score.
        for row in range(num_rows):
            if (all_routes[index_max].path[row] == i):
                color[row] = 'green'

        # Plot dots.
        axs[j, 0].scatter(x_data, strat_distr[:, i], s=5, c=color, zorder=2)

        axs[j, 0].set_title(unit_names[i], size=9, y=0.91)
        axs[j, 0].set_ylabel(str(j + 1))

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
def generate_missing_lithos():
    '''
    Generates a list of missing drillhole lithologies from unit data in all data files.
    '''

    directory = "data/real/dist_files/litho_tables"
    #directory = "data/real/dh_files/litho_tables"

    # Drillsample data columns.
    column_from = 3
    column_to = 4
    column_litho = 8

    counter = 0

    for file in os.listdir(directory):
        counter = counter + 1
        print("PROCESSING FILE NUMBER =", counter)

        filename = os.fsdecode(file)
        collarID = int(filename[6:-4])
        print(collarID)

        drillsample_filename = "data/real/dist_files/litho_tables/litho_" + str(collarID) + ".csv"
        dist_table_filename = "data/real/dist_files/dist_tables/100_500k_map_near_" + str(collarID) + ".csv"

        #drillsample_filename = "data/real/dh_files/litho_tables/litho_" + str(collarID) + ".csv"
        #dist_table_filename = "data/real/dh_files/dist_tables/100k_map_near_" + str(collarID) + ".csv"

        # Drill sample data.
        drillsample_data = read_drillsample_data(drillsample_filename, column_from, column_to, column_litho)

        # Unit lithologies data.
        strat_data = read_strat_data(dist_table_filename)

        # Generating the table of possible strata paths.
        strata_table = generate_strata_table(drillsample_data, strat_data)

    print("Missing lithologies list:")
    print("Number lithos missing =", len(missing_lithos))
    print(missing_lithos)

#=============================================================================
def main():
    print('Started litho2strat')

    #generate_missing_lithos()
    #exit()

    # Topology file.
    topology_filename = "data/real/ASUD_strat.gml"

    # # Mark's new data.
    # collarID = 548917
    # drillsample_filename = "data/real/dh_files/litho_tables/litho_" + str(collarID) + ".csv"
    # dist_table_filename = "data/real/dh_files/dist_tables/100k_map_near_" + str(collarID) + ".csv"
    #
    # # Drillsample data columns.
    # column_from = 4
    # column_to = 5
    # column_litho = 8

    # Mark's data with known solutions.
    #collarID = 1209857
    #collarID =  353386

    # Confirmed results (using 1 closest unit & single top unit).
    collarID = 2182336
    # collarID = 2182335
    # collarID = 2182340
    # collarID = 2182339
    # collarID = 2182338
    # collarID = 2182335
    # collarID = 2182334

    drillsample_filename = "data/real/dist_files/litho_tables/litho_" + str(collarID) + ".csv"
    dist_table_filename = "data/real/dist_files/dist_tables/100_500k_map_near_" + str(collarID) + ".csv"

    # Drillsample data columns.
    column_from = 3
    column_to = 4
    column_litho = 8

    #--------------------------------------------------------------
    # Reading the input data.
    #--------------------------------------------------------------
    # Drill sample data.
    drillsample_data = read_drillsample_data(drillsample_filename, column_from, column_to, column_litho)

    # Group the litho sequence inside the drillhole data.
    # Commented because it leads to many additional routes when splitting the existing single lithos.
    #drillsample_data = group_drillhole_litho_sequence(drillsample_data)

    # Unit lithologies data.
    strat_data = []
    strat_data = read_strat_data(dist_table_filename)

    # Thickness data.
    thickness_data = []
    if (add_thickness_constraints):
        thickness_data = read_thickness_data(thickness_filename)

    # Topology data.
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

    all_routes, all_routes_number = generate_strat_routes(strat_data, drillsample_data, thickness_data, graph)

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
    num_units = len(strat_data)
    plot_unit_probabilities(all_routes, drillsample_data, num_units)

#=============================================================================
if __name__ == "__main__":
    main()

