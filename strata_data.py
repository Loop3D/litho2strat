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

    #----------------------------------------------------------------------------------------------------
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
