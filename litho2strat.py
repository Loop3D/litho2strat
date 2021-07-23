import numpy as np
import pandas as pd

#==============================================================================
'''
Helper functions for reading the input data.
'''
def read_strat_data(filename):
    # Converter from string to list.
    str2list = lambda x: x.strip("[]").replace("'","").split(", ")
    data = pd.read_csv(filename, converters={"lithologies": str2list})
    return data

def read_thickness_data(filename):
    data = pd.read_csv(filename)
    return data

def read_drillsample_data(filename):
    data = pd.read_csv(filename)
    return data

#==============================================================================
def generate_strata_table(drillsample_data, strat_data):
    '''
    Generates the stratigraphic table.
    '''
    nRows = drillsample_data.shape[0] # Corresponds to lithology in drillhole sample.
    nStrat = strat_data.shape[0]

    print("nRows = ", nRows)
    print("nStrat = ", nStrat)

    stratTable = np.full((nRows, nStrat), False)

    for row in range(nRows):
        litho = drillsample_data["lithology"][row]
        for strat in range(nStrat):
            if (litho in strat_data["lithologies"][strat]):
                stratTable[row][strat] = True

    return stratTable

#==============================================================================
'''
A class for storing the stratigraphic route.
'''
class StrataRoute:
    def __init__(self, orig=None):
        if orig is None:
        # Constructor.
            # Strata index for every drillhole depth (defined by list index).
            self.path = []
            # The thickness of the last strata unit.
            self.current_thickness = 0.
            # The number of strata units.
            self.num_units = 0
        else:
        # Copy constructor.
            self.path = orig.path.copy()
            self.current_thickness = orig.current_thickness
            self.num_units = orig.num_units

    def __str__(self):
        return str(self.path)

    def __repr__(self):
        return str(self.path)

    def get_last_position(self):
        return self.path[-1]

    def add_new_position(self, position):
        return self.path.append(position)

    def get_strata_sequence(self):
        return tuple(set(self.path))

#==============================================================================
def get_thickness_change(drillsample_data, row):
    '''
    Returns a thickness change for a given row in the drillhole sample..
    '''
    return drillsample_data["to"][row] - drillsample_data["from"][row]

#==============================================================================
def get_strata_thickness(thickness_data, strat_index):
    mean = thickness_data["thickness_mean"][strat_index]
    range = thickness_data["thickess_range"][strat_index]
    thickness = {
        "min": mean - range,
        "max": mean + range
    }
    return thickness

#==============================================================================
def can_stay_in_strata(current_thickness, thickness_data, strat_index):
    strata_thickness = get_strata_thickness(thickness_data, strat_index)
    if (current_thickness < strata_thickness["max"]):
        return True
    else:
        False

#==============================================================================
def can_change_strata(route, thickness_data, strat_index):
    if (route.num_units == 1):
        # Ignore thickness for the top unit.
        return True

    strata_thickness = get_strata_thickness(thickness_data, strat_index)

    if (route.current_thickness >= strata_thickness["min"]):
        return True
    else:
        False

#==============================================================================
def generate_strat_routes(stratTable, drillsample_data, thickness_data):
    '''
    Generating stratigraphic routes.
    '''
    nRows = stratTable.shape[0]
    nStrat = stratTable.shape[1]

    all_routes = []

    # Set the initial routes.
    row = 0
    thickness_change = get_thickness_change(drillsample_data, row)

    for strat in range(nStrat):
        if (stratTable[row][strat]):
            new_route = StrataRoute()
            new_route.path = [strat]
            new_route.current_thickness = thickness_change
            new_route.num_units = 1
            all_routes.append(new_route)

    print("Starting routes:")
    print(all_routes)

    rowMax = nRows
    print("rowMax = ", rowMax)

    # Going through the strata table and generating the routes.
    for row in range(1, rowMax):
        #print("ROW = ", row, drillsample_data["lithology"][row])

        # NOTE: We iterate over the COPY of the list (slice)!
        for route in all_routes[:]:
            route_old = StrataRoute(route)
            # The current strata index.
            strat0 = route.get_last_position()
            thickness_change = get_thickness_change(drillsample_data, row)

            #-----------------------------------------------------------------
            # First check if we can go down the same srata unit.
            #-----------------------------------------------------------------
            strat = strat0
            # Apply unit thickness constraints.
            can_stay = can_stay_in_strata(route.current_thickness, thickness_data, strat)

            if (can_stay and stratTable[row][strat]):
                # Adding new route position.
                route.add_new_position(strat)
                route.current_thickness += thickness_change
            else:
                if (row < rowMax - 1):
                # Did not reach the end of a drillhole, and cannot go down the same unit.
                    # REMOVE the route.
                    all_routes.remove(route)

            #-----------------------------------------------------------------
            # Check if we can go down in other stratas.
            #-----------------------------------------------------------------
            # Apply unit thickness constraints.
            can_change = can_change_strata(route_old, thickness_data, strat0)

            if (can_change): 
                for strat in range(strat0 + 1, nStrat):
                    if (stratTable[row][strat]):
                        # ADDING new route.
                        new_route = StrataRoute(route_old)
                        new_route.add_new_position(strat)
                        new_route.current_thickness = thickness_change
                        new_route.num_units += 1
                        all_routes.append(new_route)

    return all_routes

#==============================================================================
def write_routes_to_file(filename, drillsample_data, all_routes):
    '''
    Writing stratigraphic routes to file.
    '''
    f = open(filename, "w")
    nRows = drillsample_data.shape[0]
    for row in range(nRows):
        depth = drillsample_data["from"][row]
        f.write("%f " % depth)
        # Calculate the number of unique strata for this depth.
        unique_lithos = set([])
        for route in all_routes:
            unique_lithos.add(route.path[row])
        f.write("%d " % len(unique_lithos))

        for route in all_routes:
            unique_lithos
            f.write("%d " % route.path[row])
        f.write("\n")
    f.close()

#==============================================================================
def print_unique_routes(all_routes):
    '''
    Print all unique routes (i.e., with unique strata sequence).
    '''
    unique_routes = set([])
    for route in all_routes:
        unique_routes.add(route.get_strata_sequence())

    print("Number of unique routes = ", len(unique_routes))
    for route in unique_routes:
        print(route)
#=============================================================================

def main():
    print('Started litho2strat')

    # Paths to data.
    strat_filename = "data/simple/strat_1627022992.5507748.csv"
    thickness_filename = "data/simple/thickness_mean_1627022992.5507748.csv"
    drillsample_filename = "data/simple/drill_sample_1627022992.5507748.csv"

    strat_data = read_strat_data(strat_filename)
    thickness_data = read_thickness_data(thickness_filename)
    drillsample_data = read_drillsample_data(drillsample_filename)

    stratTable = generate_strata_table(drillsample_data, strat_data)
    all_routes = generate_strat_routes(stratTable, drillsample_data, thickness_data)

    print("Total number of routes = ", len(all_routes))

    # Print all unique routes (i.e., unique strata sequence).
    print_unique_routes(all_routes)

    # Write results to the file.
    write_routes_to_file("strata.txt", drillsample_data, all_routes)

if __name__ == "__main__":
    main()

