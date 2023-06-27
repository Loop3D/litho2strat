'''
  This program generates a set of plausible stratigraphies with uncertainties, for a given drillhole lithology log.
  It uses map data for distance and topology constraints, and several free parameters describing the solution complexity (level of deformation) constraints.

  Author: Vitaliy Ogarko, vogarko@gmail.com
  The University of Western Australia
'''

import configparser
from dataclasses import dataclass

@dataclass
class InputParameters:
    '''
    Contains all input parameters used in the Parfile.
    '''
    #-------------------------------
    # Section 'FilePaths'.
    #-------------------------------
    # Topology file.
    topology_filename: str = ""
    # The Cover unit data file.
    cover_unit_filename: str = ""
    # The ignore items list file.
    ignore_list_filename: str = ""
    # Alternative rock names file.
    alternative_rock_names_filename: str = ""
    # Unit colours for drawing stratigraphy logs.
    unit_colors_filename: str = ""
    # Drillhole lithology data file. The $collarID$ in the file name will be replaced with the actual value.
    drillsample_filename: str = ""
    # Drillhole stratigraphy data file. The $collarID$ in the file name will be replaced with the actual value.
    stratasample_filename: str = ""
    # Units near the collar with distances data file. The $collarID$ in the file name will be replaced with the actual value.
    dist_table_filename: str = ""

    #-------------------------------
    # Section 'DataHeaders'.
    #-------------------------------
    # Drillsample lithology data csv file column names.
    drillsample_header: str = ""
    # Drillsample stratigraphy data csv file column names.
    stratasample_header: str = ""
    # Unit composition and distance data csv file column names.
    strata_data_header: str = ""

    #-------------------------------
    # Section 'SolverParameters'.
    #-------------------------------
    # Unit contact topology constraints (extracted from map data).
    add_topology_constraints: bool = True
    # Jumps over units in topology graph to allow skipping units: i.e., one jump would allow contact A->C for the graph A->B->C.
    max_num_strata_jumps: int = 0
    # 'Age alignment' constraints: the maximum number of times the age direction can flip.
    max_num_age_flips: int = 2
    # 'Returning to the same unit' constraints.
    max_num_returns_per_unit: int = 1
    # The maximum number of unit contacts inside the same litholgy sequence.
    max_num_unit_contacts_inside_litho: int = 0
    # Use the single closest unit for the top (first) lithology.
    single_top_unit: bool = True

    #-------------------------------
    # Section 'Correlation'.
    #-------------------------------
    # Correlation score normalization power. Higher power leads to shorter strata sequence.
    correlation_power: float = 1.0

    #-------------------------------
    # Section 'DataPreprocessing'.
    #-------------------------------
    # The number of nearest units (for distance constraints).
    number_nearest_units: int = 3
    # Minimum score for drillhole lithologies to use them.
    min_drillhole_litho_score: int = 80
    # Group drillhole lithology sequence.
    # Note: use this for max_num_unit_contacts_inside_litho > 0 to avoid the solution number to blow.
    group_drillhole_lithos: bool = False
    # The cover ration threshold (relative length) for removing the cover.
    cover_ratio_threshold: float = 0.65

    #-------------------------------
    # Section 'CollarIDs'.
    #-------------------------------
    # List of collar IDs used.
    collarIDs: str = ""

    #-------------------------------
    # Other.
    #-------------------------------
    # Adding thickness constraints (requires unit thickness data which we do not have).
    add_thickness_constraints: bool = False

#=============================================================================
def read_input_parameters(parfile_path):
    '''
    Read input parameters from Parfile.
    '''
    config = configparser.ConfigParser()
    if (len(config.read(parfile_path)) == 0):
        raise ValueError("Failed to open/find a parameters file!")

    par = InputParameters()

    section = 'FilePaths'
    print(config.items(section))

    par.topology_filename = config.get(section, 'topology_filename', fallback = '')
    par.cover_unit_filename = config.get(section, 'cover_unit_filename', fallback = '')
    par.ignore_list_filename = config.get(section, 'ignore_list_filename', fallback = '')
    par.alternative_rock_names_filename = config.get(section, 'alternative_rock_names_filename', fallback = '')
    par.unit_colors_filename = config.get(section, 'unit_colors_filename', fallback = '')
    par.drillsample_filename = config.get(section, 'drillsample_filename')
    par.stratasample_filename = config.get(section, 'stratasample_filename', fallback = '')
    par.dist_table_filename = config.get(section, 'dist_table_filename')

    section = 'DataHeaders'
    print(config.items(section))

    drillsample_header = config.get(section, 'drillsample_header')
    stratasample_header = config.get(section, 'stratasample_header', fallback = '')
    strata_data_header = config.get(section, 'strata_data_header')

    par.drillsample_header = drillsample_header.replace(" ","").split(",")
    par.stratasample_header = stratasample_header.replace(" ","").split(",")
    par.strata_data_header = strata_data_header.replace(" ","").split(",")

    section = 'SolverParameters'
    print(config.items(section))

    par.add_topology_constraints = config.getboolean(section, 'add_topology_constraints', fallback = par.add_topology_constraints)
    par.max_num_strata_jumps = config.getint(section, 'max_num_strata_jumps', fallback = par.max_num_strata_jumps)
    par.max_num_age_flips = config.getint(section, 'max_num_age_flips', fallback = par.max_num_age_flips)
    par.max_num_returns_per_unit = config.getint(section, 'max_num_returns_per_unit', fallback = par.max_num_returns_per_unit)
    par.max_num_unit_contacts_inside_litho = config.getint(section, 'max_num_unit_contacts_inside_litho', fallback = par.max_num_unit_contacts_inside_litho)
    par.single_top_unit = config.getboolean(section, 'single_top_unit', fallback = par.single_top_unit)

    section = 'Correlation'
    print(config.items(section))

    par.correlation_power = config.getfloat(section, 'correlation_power', fallback = par.correlation_power)

    section = 'DataPreprocessing'
    print(config.items(section))

    par.number_nearest_units = config.getint(section, 'number_nearest_units', fallback = par.number_nearest_units)
    par.min_drillhole_litho_score = config.getint(section, 'min_drillhole_litho_score', fallback = par.min_drillhole_litho_score)
    par.group_drillhole_lithos = config.getboolean(section, 'group_drillhole_lithos', fallback = par.group_drillhole_lithos)
    par.cover_ratio_threshold = config.getfloat(section, 'cover_ratio_threshold', fallback = par.cover_ratio_threshold)

    section = 'CollarIDs'
    print(config.items(section))

    collarIDs = config.get(section, 'collarIDs').split(",")
    par.collarIDs = [c.strip() for c in collarIDs]

    # Hardcoded as we don't have the thickness data.
    par.add_thickness_constraints = False

    return par
