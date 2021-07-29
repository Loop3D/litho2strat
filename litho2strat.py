import numpy as np
import pandas as pd
import matplotlib.pylab as pl

max_num_foreign_litho = 4

#==============================================================================
'''
Helper functions for reading the input data.
'''
def read_strat_data(filename):
    # Converter from string to list.
    str2list = lambda x: x.strip("[]").replace("'", "").split(", ")
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
        litho_found = False
        for strat in range(nStrat):
            if (litho in strat_data["lithologies"][strat]):
                litho_found = True
                stratTable[row, strat] = True

        if (not litho_found):
        # Lithology not found in stratas.
            print("Not found litho: ", litho)
            # Allow this lithology to be present in any strata.
            for strat in range(nStrat):
                stratTable[row, strat] = True

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
            # The number of "foreign" lithologies.
            self.num_foreign_litho = 0
            # The list of "foreign" lithologies for current unit.
            self.current_foreign_lithos = []
            # Flag for removal.
            self.to_remove = False
        else:
        # Copy constructor.
            self.path = orig.path.copy()
            self.current_thickness = orig.current_thickness
            self.num_units = orig.num_units
            self.num_foreign_litho = orig.num_foreign_litho
            #self.current_foreign_lithos = []
            self.to_remove = orig.to_remove

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
    return drillsample_data["to"][row] - drillsample_data["from"][row]

#==============================================================================
def get_min_strata_thickness(thickness_data):
    nStrat = thickness_data.shape[0]
    return [thickness_data["thickness_mean"][i] - thickness_data["thickess_range"][i] for i in range(nStrat)]

#==============================================================================
def get_max_strata_thickness(thickness_data):
    nStrat = thickness_data.shape[0]
    return [thickness_data["thickness_mean"][i] + thickness_data["thickess_range"][i] for i in range(nStrat)]

#==============================================================================
def process_foreign_lithology(route, current_litho, row):
    '''
    Process the "foreign" lithology on the route.
    '''
    # Do not count lighology as 'new' if it was already present in this unit.
    if (current_litho not in route.current_foreign_lithos):
        # Count it as new "foreign" lithology.
        route.num_foreign_litho += 1
        route.current_foreign_lithos.append(current_litho)

#==============================================================================
def generate_strat_routes(stratTable, drillsample_data, thickness_data, 
                          add_thickness_constraints, add_foreign_lithology, can_return_to_unit):
    '''
    Generating stratigraphic routes.
    '''
    nRows = stratTable.shape[0]
    nStrat = stratTable.shape[1]

    # Extract strata thikcness to lists (faster data structures).
    min_strata_thickness = get_min_strata_thickness(thickness_data)
    max_strata_thickness = get_max_strata_thickness(thickness_data)

    all_routes = []

    # Set the initial routes.
    row = 0
    thickness_change = get_thickness_change(drillsample_data, row)

    for strat in range(nStrat):
        if (stratTable[row, strat]):
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
        current_litho = drillsample_data["lithology"][row]
        print("ROW = ", row, current_litho, len(all_routes))
        if (len(all_routes) == 0):
            break

        thickness_change = get_thickness_change(drillsample_data, row)
        new_routes = []

        # Iterate over all routes.
        for route in all_routes:
            # The current strata index.
            strat0 = route.path[-1]
            current_thickness = route.current_thickness
            can_add_foreign_litho = (route.num_foreign_litho < max_num_foreign_litho)

            #-----------------------------------------------------------------
            # Check if we can go down in other stratas (and create new routes).
            #-----------------------------------------------------------------
            can_change = True

            # Apply unit thickness constraints.
            if (add_thickness_constraints):
                # Ignore thickness for the top unit
                if (route.num_units > 1):
                    can_change = can_change and (current_thickness >= min_strata_thickness[strat0])

            if (can_change):
                strataList = []
                if can_return_to_unit:
                    # All strata units excluding the current one.
                    strataList = [strat for strat in range(nStrat) if strat != strat0]
                else:
                    # Only units older than the current one.
                    strataList = range(strat0 + 1, nStrat)

                # Looking to which strata unit we can change.
                for strat in strataList:
                    path_exists = stratTable[row, strat]

                    # Processing "foreign" litho constraints.
                    foreign_litho_added = False
                    if (add_foreign_lithology):
                        if (not path_exists and can_add_foreign_litho):
                            foreign_litho_added = True

                    if (path_exists or foreign_litho_added):
                        # Making the new route.
                        new_route = StrataRoute(route)
                        new_route.path.append(strat)
                        new_route.current_thickness = thickness_change
                        new_route.num_units += 1
                        new_route.current_foreign_lithos = []
                        # Processing the "foreign" lithology.
                        if (foreign_litho_added):
                            process_foreign_lithology(new_route, current_litho, row)
                        # Add new route into the list.
                        new_routes.append(new_route)

            #-----------------------------------------------------------------
            # Check if we can go down the same srata unit (if cannot -- remove the current route).
            #-----------------------------------------------------------------
            can_stay = True

            if (add_thickness_constraints):
                # Apply unit thickness constraints.
                can_stay = can_stay and (current_thickness < max_strata_thickness[strat0])

            path_exists = stratTable[row, strat0]

            # Processing "foreign" litho constraints.
            foreign_litho_added = False
            if (add_foreign_lithology):
                if (not path_exists and can_add_foreign_litho):
                    foreign_litho_added = True

            can_stay = can_stay and (path_exists or foreign_litho_added)

            if (can_stay):
                # Adding new route position.
                route.path.append(strat0)
                route.current_thickness += thickness_change
                # Processing the "foreign" lithology.
                if (foreign_litho_added):
                    process_foreign_lithology(route, current_litho, row)
            else:
            # Did not reach the end of a drillhole, and cannot go down the same unit.
                # Mark the route for removal.
                route.to_remove = True

        # Remove the routes marked for removal.
        all_routes = [route for route in all_routes if not route.to_remove]
        # Addig new routes.
        all_routes.extend(new_routes)

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
        unique_units = set([])
        for route in all_routes:
            unique_units.add(route.path[row])
        f.write("%d " % len(unique_units))

        for route in all_routes:
            f.write("%d " % route.path[row])
        f.write("\n")
    f.close()

#==============================================================================
def print_unique_routes(all_routes, print_paths):
    '''
    Print all unique routes (i.e., with unique strata sequence).
    '''
    unique_routes = set([])
    for route in all_routes:
        unique_routes.add(route.get_strata_sequence())

    print("Number of unique routes = ", len(unique_routes))
    if (print_paths):
        for route in unique_routes:
            print(route)

#=============================================================================
def plot_routes(depth_data, routes):
    '''
    Plot and display the routes.
    '''
    print("Plotting the routes...")
    for route in routes:
        pl.plot(depth_data, route.path, '.-')

    pl.xlabel('Depth')
    pl.ylabel('Strata unit index')
    pl.show()

#=============================================================================

def main():
    print('Started litho2strat')

    # Paths to data.
#     strat_filename = "data/simple/strat_1627022992.5507748.csv"
#     thickness_filename = "data/simple/thickness_mean_1627022992.5507748.csv"
#     drillsample_filename = "data/simple/drill_sample_1627022992.5507748.csv"

    strat_filename       = "data/foreign_litho/strat_1627025194.6300328.csv"
    thickness_filename   = "data/foreign_litho/thickness_mean_1627025194.6300328.csv"
    drillsample_filename = "data/foreign_litho/drill_sample_1627025194.6300328.csv"

    strat_data = read_strat_data(strat_filename)
    thickness_data = read_thickness_data(thickness_filename)
    drillsample_data = read_drillsample_data(drillsample_filename)

    stratTable = generate_strata_table(drillsample_data, strat_data)

    add_thickness_constraints = True
    add_foreign_lithology = True
    can_return_to_unit = True

    all_routes = generate_strat_routes(stratTable, drillsample_data, thickness_data, 
                                       add_thickness_constraints, add_foreign_lithology, can_return_to_unit)

    print("Total number of routes = ", len(all_routes))

    # Print all unique routes (i.e., unique strata sequence).
    #print_unique_routes(all_routes, False)

    # Write results to the file.
    #write_routes_to_file("strata.txt", drillsample_data, all_routes)

    # Plot the results.
    #plot_routes(drillsample_data["from"], all_routes)

    #--------------------------------------------------------------------
    # Plot route statistics.
    #'''
    nRows = drillsample_data.shape[0]
    nUnits = strat_data.shape[0]
    strat_distr = np.zeros((nRows, nUnits))
    # Building the distribution of unit presence at every depth.
    for route in all_routes:
        for row in range(nRows):
            unit_index = route.path[row]
            strat_distr[row, unit_index] += 1
    # Normalize.
    strat_distr = strat_distr / len(all_routes)

    # Generating the plots.
    # Increasing the figure size.
    pl.rcParams["figure.figsize"] = (12.8, 9.6) # Default size = (6.4, 4.8)
    x_data = drillsample_data["from"]
    fig, axs = pl.subplots(nUnits, sharey=True)
    fig.suptitle('Probability for each unit. Max foreign lithos = ' + str(max_num_foreign_litho))
    for i in range(nUnits):
        axs[i].plot(x_data, strat_distr[:, i])
        axs[i].set(ylabel=str(i))
        if (i != nUnits - 1):
            # Hide tick labels.
            axs[i].set_xticklabels([])
        axs[i].xaxis.grid(True)
 
    pl.xlabel('Depth')
    pl.show()
    #'''
    #--------------------------------------------------------------------

#=============================================================================
if __name__ == "__main__":
    main()

