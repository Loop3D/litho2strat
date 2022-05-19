'''
  This program generates a set of plausible stratigraphies with uncertainties, for a given drillhole lithology log.
  It uses map data for distance and topology constraints, and several free parameters describing the solution complexity (level of deformation) constraints.

  Author: Vitaliy Ogarko, vogarko@gmail.com
  The University of Western Australia
'''

from dataclasses import dataclass, field
from typing import Dict

#========================================================================================================
@dataclass
class StrataData:
    '''
    The strata data needed for the Strata Solver.
    '''
    # Maps the unit name to the list of its lithologies.
    unit2litho: Dict[str, list] = field(default_factory=dict)

    # Maps the lithology name to the sorted list of distances to units with corresponding unit names.
    litho2dist: Dict[str, list] = field(default_factory=dict)

    #========================================================================================================
    def get_unit_names(self):
        '''
        Defines the mapping between the unit index and unit name.
        '''
        unit_names = []
        for unit_name in self.unit2litho:
            if (unit_name == 'Cover'):
                # Map the Cover's index to zero.
                unit_names.insert(0, unit_name)
            else:
                unit_names.append(unit_name)

        return unit_names

    #========================================================================================================
    def filter_strat_data_based_on_drillhole_lithos(self, drillhole_lithos):
        '''
        Filter units based on the drillhole lithologies: 
            - Remove lithologies, that are not present in the drillhole data.
            - Remove the units that do not contain the drillhole lithos.
        '''
        strat_filtered = dict()
        for unit_name in self.unit2litho:
            for litho in self.unit2litho[unit_name]:
                # Only add lithologies that are present in drillhole data.
                if (litho in drillhole_lithos):
                    if (unit_name in strat_filtered):
                        strat_filtered[unit_name].append(litho)
                    else:
                        strat_filtered[unit_name] = [litho]

        print("The number of filtered units: " + str(len(strat_filtered)))

        self.unit2litho = strat_filtered

    #=======================================================================================
    def filter_strat_data_based_on_distance(self, number_nearest_units):
        '''
        Filter units based on the distance from drillhole.
        '''
        strat_dist = dict()
        for unit_name in self.unit2litho:
            for litho in self.unit2litho[unit_name]:
                # Sorted distance list for this lithology.
                dist_list = self.litho2dist[litho]

                # Consider only N closest unit codes.
                for el in dist_list[:number_nearest_units]:
                    unit_name_nearest = el[1]
                    if (unit_name == unit_name_nearest):
                        if (unit_name in strat_dist):
                            strat_dist[unit_name].append(litho)
                        else:
                            strat_dist[unit_name] = [litho]
                        break

        print("The number of filtered (by distance) units: " + str(len(strat_dist)))

        self.unit2litho = strat_dist

#========================================================================================================
