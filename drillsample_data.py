'''
  This program generates a set of plausible stratigraphies with uncertainties, for a given drillhole lithology log.
  It uses map data for distance and topology constraints, and several free parameters describing the solution complexity (level of deformation) constraints.

  Author: Vitaliy Ogarko, vogarko@gmail.com
  The University of Western Australia
'''

from dataclasses import dataclass, field
from typing import List

#========================================================================================================
@dataclass
class DrillSampleDataRow:
    '''
    Contains the drillsample data attributes for one row.
    '''
    depth_from: float = 0.
    depth_to: float = 0.
    lithos: List[str] = field(default_factory=list)
    scores: List[int] = field(default_factory=list)

#========================================================================================================
@dataclass
class DrillsampleData:
    '''
    The drillsample data needed for the Strata Solver.
    '''
    # Maps the unit name to the list of its lithologies.
    rows: List[DrillSampleDataRow] = field(default_factory=list)

