'''
  This program generates a set of plausible stratigraphies with uncertainties, for a given drillhole lithology log.
  It uses map data for distance and topology constraints, and several free parameters describing the solution complexity (level of deformation) constraints.

  Author: Vitaliy Ogarko, vogarko@gmail.com
  The University of Western Australia
'''

import csv
import networkx as nx
from dataclasses import dataclass

from strata_data import StrataData
from drillsample_data import DrillSampleDataRow, DrillsampleData

#==============================================================================
def fix_litho_name(litho):
    '''
    Convert the map lithology name to the 'CET lithology' name.
    '''
    # Stripping the whitespaces.
    litho = litho.strip()

    # Convert things like metagranite to granite.
    if (litho[0:4] == 'meta'):
        litho = litho[4:]

    return litho

#====================================================================================
def read_strat_data(dist_table_filename, alternative_rock_names):
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

    # List of unique lithos in the map data (disregarding the alternative names).
    unique_lithos_original = set()

    # Create the object with strata data.
    strata_data = StrataData()

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

            #-----------------------------------------
            # Hard fixes for "mudstone" and "coal". 
            # TODO: Fix the original data.
            #-----------------------------------------
            lithname1 = row[lithname1_column]
            if ("mudstone" in lithname1):
                lithos.append("mudstone")

            description = row[descript_column]
            if ("coal" in description):
                lithos.append("coal")
                lithos.append("lignite")

            #-----------------------------------------
            # Fixing the litho names.
            #-----------------------------------------
            lithos = [fix_litho_name(l) for l in lithos]

            #-----------------------------------------
            unique_lithos_original.update(lithos)

            #-----------------------------------------
            # Adding the alternative litho names.
            #-----------------------------------------
            alt_lithos = []
            for litho in lithos:
                for alt_names in alternative_rock_names:
                    if litho in alt_names:
                        alt_lithos.extend(alt_names)

            lithos.extend(alt_lithos)

            # Remove duplicates.
            lithos = list(dict.fromkeys(lithos))

            #-----------------------------------------
            # Store unit in the sorted distance map.
            #-----------------------------------------
            strata_data.add_unit_to_distance_map(unit_name, distance, lithos)

            #-----------------------------------------
            # Adding the lithos to the dictionary (excluding the duplicates).
            #-----------------------------------------
            strata_data.add_unit_to_units_map(unit_name, lithos)

    print("The total number of units:", strata_data.get_num_units())

    #------------------------------------------------------------------------
    print("The total number of (original) lithologies:", len(unique_lithos_original))
    print(sorted(unique_lithos_original))

    return strata_data

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
def read_drillsample_data(filename, header, ignore_list, min_litho_score):
    '''
    Reading drill sample data from csv file.
    '''
    all_lithos = set()
    drillsample_data = DrillsampleData()

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
                # Stripping the whitespaces.
                litho = litho.strip()

                if (litho not in ignore_list
                    and scores[index] >= min_litho_score):
                    # Adding lithology and its score.
                    data.lithos.append(litho)
                    data.scores.append(scores[index])

                    # Gather all unique lithos.
                    all_lithos.add(litho)

            if (len(data.lithos) > 0):
                drillsample_data.rows.append(data)

    print(all_lithos)
    print("The number of drillhole lithologies: " + str(len(all_lithos)))

    return drillsample_data

#==============================================================================
@dataclass
class ThicknessDataElement:
    '''
    The thickness data element
    '''
    # Mean unit thickness.
    mean: float = 0.
    # The thickness variation range.
    range: float = 0.

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
            thickness = ThicknessDataElement()
            thickness.mean = float(row[1])
            thickness.range = float(row[2])
            data.append(thickness)
    return data

#==============================================================================
def read_topology_data(topology_filename, ignore_unit_age):
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

