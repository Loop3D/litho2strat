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

    # Maps the lithology name to the sorted list of distances-to-units with corresponding unit names.
    litho2dist: Dict[str, list] = field(default_factory=dict)

    # Unit to shortest distance map.
    unit2dist: Dict[str, float] = field(default_factory=dict)

    #========================================================================================================
    def get_num_units(self):
        '''
        Returns the number of units.
        '''
        return len(self.unit2litho)

    #========================================================================================================
    def get_unit_names(self):
        '''
        Defines the mapping between the unit index and unit name.
        '''
        unit_names = []
        for unit_name in self.unit2litho:
            unit_names.append(unit_name)

        return unit_names

    #==============================================================================
    def get_unique_lithos(self):
        '''
        Returns a unique list of lithologies.
        '''
        lithos = list()
        for unit_name in self.unit2litho:
            for litho in self.unit2litho[unit_name]:
                lithos.append(litho)
        unique_lithos = list(dict.fromkeys(lithos))

        return unique_lithos

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

        self.unit2litho = strat_filtered

        print("The number of filtered (by drillhole lithos) units:", self.get_num_units())

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

        self.unit2litho = strat_dist

        print("The number of filtered (by distance) units:", self.get_num_units())

    #===========================================================================================
    def add_unit_to_units_map(self, unit_name, lithos):
        '''
        Adds a unit with its lithologies to the units map.
        '''
        if unit_name in self.unit2litho:
            for litho in lithos:
                if litho not in self.unit2litho[unit_name]:
                    self.unit2litho[unit_name].append(litho)
        else:
            # Adding a new unit, with unique list of lithologies.
            self.unit2litho[unit_name] = list(dict.fromkeys(lithos))

    #===========================================================================================
    def add_unit_to_distance_map(self, unit_name, distance, lithos):
        '''
        Adds a unit with its lithologies to the distance map.
        '''
        #---------------------------------------------------
        # Adding unit to unit2dist map.
        #---------------------------------------------------
        if (unit_name in self.unit2dist):
            if (distance < self.unit2dist[unit_name]):
                # Found a shorter distance - updating data.
                self.unit2dist[unit_name] = distance
        else:
            self.unit2dist[unit_name] = distance

        #---------------------------------------------------
        # Adding unit to litho2dist map.
        #---------------------------------------------------
        for litho in lithos:
            el = (distance, unit_name)
            if (litho in self.litho2dist):
                found_unit = False
                # Iterate over distance list.
                for tup in self.litho2dist[litho]:
                    if (unit_name == tup[1]):
                    # Found this unit in the list.
                        found_unit = True
                        if (distance < tup[0]):
                            # Update element to a smaller distance.
                            self.litho2dist[litho].remove(tup)
                            self.litho2dist[litho].append(el)
                        break
                if (not found_unit):
                    self.litho2dist[litho].append(el)
            else:
                self.litho2dist[litho] = [el]
            # Sort the list by distance.
            self.litho2dist[litho].sort(key=lambda tup: tup[0])

#========================================================================================================
