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
class DepthData:
    '''
    Contains the drillsample depth data.
    '''
    depth_from: List[int] = field(default_factory=list)
    depth_to: List[int] = field(default_factory=list)

#========================================================================================================
@dataclass
class DrillsampleData:
    '''
    The drillsample data needed for the Strata Solver.
    '''
    # Maps the unit name to the list of its lithologies.
    rows: List[DrillSampleDataRow] = field(default_factory=list)

    # Cover part.
    rows_cover: List[DrillSampleDataRow] = field(default_factory=list)

    #==============================================================================
    def get_num_rows(self):
        '''
        Returns the number of rows.
        '''
        return len(self.rows)

    #==============================================================================
    def get_drillhole_lithos(self):
        '''
        Returns the drillhole lithologies.
        '''
        all_lithos = set()
        for row in self.rows:
            lithos = row.lithos
            all_lithos.update(lithos)
        return all_lithos

    #==============================================================================
    def get_depth_data(self):
        '''
        Returns the depth data.
        '''
        depth_data = DepthData()
        depth_data.depth_from = [d.depth_from for d in self.rows]
        depth_data.depth_to = [d.depth_to for d in self.rows]
        return depth_data

    #==============================================================================
    def get_thickness_change(self, row):
        '''
        Returns a thickness change for a given row in the drillhole sample.
        '''
        # "To" - "From"
        return self.rows[row].depth_to - self.rows[row].depth_from

    #==============================================================================
    def group_drillhole_litho_sequence(self, max_num_unit_contacts_inside_litho):
        '''
        Group the litho sequence by name (inside the drillhole data), leaving at most N in each group,
        corresponding to number of contacts inside the litho sequence.
        '''
        N = max_num_unit_contacts_inside_litho + 1
        current_litho = self.rows[0].lithos
        from_depth = self.rows[0].depth_from
        to_depth = self.rows[0].depth_to

        data_mod = [d for d in self.rows]
        data_mod.append(DrillSampleDataRow())

        num_same_lithos = 0
        data_grouped = []

        for index, row in enumerate(data_mod):

            prev_litho = current_litho
            current_litho = row.lithos

            if (set(current_litho) != set(prev_litho) and index > 0):
            # Detected a change of lithology name.
                if (num_same_lithos == 1):
                    # Don't split single data rows.
                    N_adj = 1
                else:
                    N_adj = N

                total_thickness = to_depth - from_depth
                local_thickness = total_thickness / float(N_adj)
                litho = prev_litho

                print("Grouping lithos for:", from_depth, to_depth, litho, num_same_lithos, N_adj)

                # Generate the local grouped lithos.
                for i in range(N_adj):
                    from_depth_local = from_depth + float(i) * local_thickness
                    to_depth_local = from_depth_local + local_thickness

                    row_grouped = DrillSampleDataRow()
                    row_grouped.depth_from = from_depth_local
                    row_grouped.depth_to = to_depth_local
                    row_grouped.lithos = litho

                    data_grouped.append(row_grouped)

                # Update the starting depth for the following grouped lithos.
                from_depth = row.depth_from
                # Reset the number of lithos.
                num_same_lithos = 1
            else:
                num_same_lithos = num_same_lithos + 1

            to_depth = row.depth_to

        self.rows = data_grouped

    #==============================================================================
    def remove_cover(self, cover_lithos_filename, threshold):
        '''
        Removes the cover.
        '''
        cover_index = self.identify_cover(cover_lithos_filename, threshold)

        if (cover_index >= 0):
            print("REMOVING THE COVER until depth =", self.rows[cover_index].depth_to)
 
            self.rows_cover = self.rows[0:cover_index + 1]
            self.rows = self.rows[cover_index + 1:]

            # Print the removed Cover lithologies.
            for index, row in reversed(list(enumerate(self.rows_cover))):
                print("Removed cover:", row.depth_to, row.lithos)

    #==============================================================================
    @staticmethod
    def __has_cover_litho(lithos, cover_lithos):
        for litho in lithos:
            if (litho in cover_lithos):
                return True
        return False

    #==============================================================================
    @staticmethod
    def __all_cover_lithos(lithos, cover_lithos):
        for litho in lithos:
            if (litho not in cover_lithos):
                return False
        return True

    #==============================================================================
    def identify_cover(self, cover_lithos_filename, threshold):
        '''
        Identifies the cover.
        '''
        print("Identifying the Cover...")

        # Read the Cover lithologies.
        with open(cover_lithos_filename) as f:
            cover_lithos = f.read().splitlines()

        # Iterate from the bottom to top, and search for the Cover.
        for index, row in reversed(list(enumerate(self.rows))):
            print(index, row.depth_to, row.lithos)

            found_cover_litho = self.__has_cover_litho(row.lithos, cover_lithos)

            if (found_cover_litho):
                cover_index = index
                total_length = 0.
                cover_length = 0.

                # Calcualte the total and cover lenghts.
                for i in reversed(range(cover_index + 1)):
                    row = self.rows[i]
                    thickness = row.depth_to - row.depth_from
                    total_length += thickness

                    if (self.__has_cover_litho(row.lithos, cover_lithos)):
                        cover_length += thickness

                cover_ratio = cover_length / total_length
                print("Cover ratio:", cover_ratio)

                if (cover_ratio >= threshold):
                    return cover_index
        return -1
