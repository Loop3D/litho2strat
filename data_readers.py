'''
  This program generates a set of plausible stratigraphies with uncertainties, for a given drillhole lithology log.
  It uses map data for distance and topology constraints, and several free parameters describing the solution complexity (level of deformation) constraints.

  Author: Vitaliy Ogarko, vogarko@gmail.com
  The University of Western Australia
'''

import csv
import numpy as np
import networkx as nx
from dataclasses import dataclass, field
from typing import List
import os

#---------------------------------------------------------------------------
# The number of nearest units (for distance constraints).
number_nearest_units = 3
#---------------------------------------------------------------------------
# Minimum score for drillhole lithologies to use them.
min_drillhole_litho_score = 70
#---------------------------------------------------------------------------
# Ignore topology graph edge direction defining the unit age.
ignore_unit_age = True

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

            # Filter the lithos and scores based on the ignore list and based on the minimum score.
            for index, litho in enumerate(lithos):
                if (litho not in ignore_list
                    and scores[index] >= min_drillhole_litho_score):
                    # Adding lithology and its score.
                    data.lithos.append(litho)
                    data.scores.append(scores[index])

                    # Gather all unique lithos.
                    all_lithos.add(litho)

            if (len(data.lithos) > 0):
                all_data.append(data)

    print(all_lithos)
    print("The number of drillhole lithologies: " + str(len(all_lithos)))

    return all_data

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

