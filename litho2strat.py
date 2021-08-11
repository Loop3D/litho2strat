import csv
import numpy as np
import matplotlib.pylab as pl
import tracemalloc

max_num_foreign_lithos = 4
max_num_foreign_lithos_per_unit = 1
max_num_returns = 0
max_num_returns_per_unit = 1


#==============================================================================
# Converter from string to list.
str2list = lambda x: x.strip("[]").replace("'", "").split(", ")

#==============================================================================
def read_strat_data(filename):
    '''
    Reading lithologies for every unit from csv file.
    '''
    strat_data = []
    with open(filename, 'r') as csvfile:
        # Reading the csv data.
        csvreader = csv.reader(csvfile, delimiter=',')
        # Skipping the header.
        next(csvreader)
        # Extracting the lighology list for every csv row (strata unit).
        for row in csvreader:
            # The lithologies are stored in the third column!
            lithos = str2list(row[2])
            strat_data.append(lithos)

    return strat_data

#==============================================================================
def read_strat_data2(all_strat_filename, unit_list_filename):
    '''
    Building lithologies list for every unit from csv files data.
    '''
    # Reading all units description.
    strat_all = dict()
    with open(all_strat_filename, 'r') as csvfile:
        # Reading the csv data.
        csvreader = csv.reader(csvfile, delimiter=',')
        # Skipping the header.
        next(csvreader)
        # Extracting the lighology list for every csv row (strata unit).
        for row in csvreader:
            # The lithologies are stored in the third column!
            lithos = str2list(row[1])
            strat_all[row[0]] = lithos

    # Building unit list for a partial set of units..
    strat_data = []
    with open(unit_list_filename, 'r') as csvfile:
        # Reading the csv data.
        csvreader = csv.reader(csvfile, delimiter=',')
        # Skipping the header.
        next(csvreader)
        # Go through the list of units and form its lithology list.
        for row in csvreader:
            unit_name = row[0]
            if (unit_name in strat_all):
                strat_data.append(strat_all[unit_name])
            else:
                print("Error! No unit description found in 'all strat' list for unit name: " + unit_name)
                exit()

    return strat_data

#==============================================================================
def read_drillsample_data(filename):
    '''
    Reading drill sample data from csv file.
    '''
    data = []
    lithos = set()
    with open(filename, 'r') as csvfile:
        # Reading the csv data.
        csvreader = csv.reader(csvfile, delimiter=',')
        # Skipping the header.
        next(csvreader)
        # Extracting the data for every csv row.
        for row in csvreader:
            data.append(row)
            lithos.add(row[2])
    print("The number of drillhole lithologies: " + str(len(lithos)))
    return data

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
def generate_strata_table(drillsample_data, strat_data):
    '''
    Generates the stratigraphic table.
    '''
    num_rows = len(drillsample_data)
    num_units = len(strat_data)

    print("num_rows (before) = ", num_rows)
    print("num_units = ", num_units)

    strata_table = np.full((num_rows, num_units), False)
    new_row_index = 0

    for row in drillsample_data[:]:
        litho = row[2]
        litho_found = False
        for strat in range(num_units):
            if (litho in strat_data[strat]):
                litho_found = True
                strata_table[new_row_index, strat] = True

        if (not litho_found):
        # Lithology not found in units.
            print("WARNING: Not found lithology: ", row)
            # Treat this as "no data".
            drillsample_data.remove(row)
        else:
            new_row_index += 1

    num_rows = len(drillsample_data)
    print("num_rows (after) = ", num_rows)

    # Remove rows due to missing lithologies.
    strata_table = strata_table[0:num_rows, :]

    return strata_table

#==============================================================================
'''
A class for storing the stratigraphic route.
'''
class StrataRoute:
    def __init__(self):
        # Flag for removal.
        self.to_remove = False

    # Adding the first position to the route.
    # This method essentially initializes the route.
    def add_first_position(self, strat, thickness_change, num_units):
        # Unit index for every drillhole data.
        self.path = [strat]
        # The thickness of the last strata unit.
        self.current_thickness = thickness_change
        # The number of strata units.
        self.num_units = 1
        # The number of "foreign" lithologies.
        self.num_foreign_lithos = 0
        # Containts the list of "foreign" lithologies for every unit.
        self.foreign_lithos = [[] for i in range(num_units)]
        # Containts the number of times each unit was visited.
        self.unit_visited = np.zeros((num_units), dtype=int)
        self.unit_visited[strat] += 1

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
    # "To" - "From"
    return float(drillsample_data[row][1]) - float(drillsample_data[row][0])

#==============================================================================
def get_min_strata_thickness(thickness_data):
    # "thickness_mean" - "thickess_range"
    return [data[0] - data[1] for data in thickness_data]

#==============================================================================
def get_max_strata_thickness(thickness_data):
    # "thickness_mean" + "thickess_range"
    return [data[0] + data[1] for data in thickness_data]

#==============================================================================
def can_add_foreign_lithology(path_exists, route, current_litho, strat):
    '''
    Checks if a 'foreign' lithology can be added to the route.
    '''
    if (path_exists):
        return False
    else:
        can_add = False
        # If this lithology was already added to this unit, 
        # or if the maximum number of foreign lithologies is not exceeded then set can_add to True.
        num_foreign_lithos_in_unit = len(route.foreign_lithos[strat])
        if ((route.num_foreign_lithos < max_num_foreign_lithos 
             and num_foreign_lithos_in_unit < max_num_foreign_lithos_per_unit)
            or current_litho in route.foreign_lithos[strat]):
            can_add = True

        return can_add

#==============================================================================
def process_foreign_lithology(route, current_litho, strat):
    '''
    Adding the "foreign" lithology to the route.
    '''
    # Do not count lighology as 'new' if it was already present in this unit.
    if (current_litho not in route.foreign_lithos[strat]):
        route.num_foreign_lithos += 1
        route.foreign_lithos[strat].append(current_litho)

#==============================================================================
def flatten(S):
    '''
    Flattens the multilevel list of lists.
    For example, it will convert [[[1,1],2,2],3,3] to [1,1,2,2,3,3].
    '''
    if S == []:
        return S
    if isinstance(S[0], list):
        return flatten(S[0]) + flatten(S[1:])
    return S[:1] + flatten(S[1:])

#==============================================================================
def apply_max_num_returns_constraint(route, strata_list):
    '''
    Applies the "maximum number of returns to a unit" constraint:
    removes from the input unit list the units where the route cannot return anymore.
    '''
    # Calculate the current global number of returns to the same unit.
    num_returns0 = 0
    for visited in route.unit_visited:
        if (visited > 0):
            num_returns0 += visited - 1

    # Apply the "max number of returns" constraint.
    for strat in strata_list[:]:
        if (route.unit_visited[strat] - 1 >= max_num_returns_per_unit):
        # Reached the maximum numer of local returns (to this unit).
            strata_list.remove(strat)
        else:
        # Calculate the global number of returns, if we visit this strat again.
            num_returns = num_returns0
            if (route.unit_visited[strat] > 0):
                num_returns += 1
            if (num_returns > max_num_returns):
                # Skip this unit.
                strata_list.remove(strat)

#==============================================================================
def generate_strat_routes(strat_data, drillsample_data, thickness_data,
                          add_thickness_constraints, keep_continuous_lithos):
    '''
    Generating stratigraphic routes.
    '''
    # Generating the table of possible strata paths.
    strata_table = generate_strata_table(drillsample_data, strat_data)

    num_rows = strata_table.shape[0]
    num_units = strata_table.shape[1]

    # Extract strata thikcness to lists (faster data structures).
    min_strata_thickness = get_min_strata_thickness(thickness_data)
    max_strata_thickness = get_max_strata_thickness(thickness_data)

    # Generate a list of all lithologies (to be able to map them to integers).
    all_lithos = []
    for row in drillsample_data:
        litho_name = row[2]
        if (litho_name not in all_lithos):
            all_lithos.append(litho_name)

    all_routes = []
    all_routes_number = []

    # Set the initial routes.
    row = 0
    thickness_change = get_thickness_change(drillsample_data, row)

    for strat in range(num_units):
        if (strata_table[row, strat]):
            new_route = StrataRoute()
            new_route.add_first_position(strat, thickness_change, num_units)
            # Adding new route into the list.
            all_routes.append(new_route)

    print("Starting routes:")
    print(all_routes)

    rowMax = num_rows
    print("rowMax = ", rowMax)

    # Going through the strata table and generating the routes.
    for row in range(1, rowMax):
        current_litho = drillsample_data[row][2]
        current_litho_index = all_lithos.index(current_litho)
        num_routes = len(all_routes)
        all_routes_number.append(num_routes)

        print("ROW = ", row, current_litho, num_routes)
        if (num_routes == 0):
            break

        thickness_change = get_thickness_change(drillsample_data, row)
        new_routes = []

        # Iterate over all routes.
        for route in all_routes:
            # The current strata index.
            strat0 = route.path[-1]
            current_thickness = route.current_thickness

            #-----------------------------------------------------------------
            # Check if we can go down in other stratas (and create new routes).
            #-----------------------------------------------------------------
            can_change = True

            # Add 'continous lithology' constraints.
            if (keep_continuous_lithos):
                previous_litho = drillsample_data[row - 1][2]
                if (current_litho == previous_litho):
                    can_change = False

            # Apply unit thickness constraints.
            if (add_thickness_constraints):
                # Ignore thickness for the top unit
                if (route.num_units > 1):
                    can_change = can_change and (current_thickness >= min_strata_thickness[strat0])

            if (can_change):
                # Strata units excluding the current one, and those visited the maximum number of times.
                strata_list = [strat for strat in range(num_units) if (strat != strat0)]
                # Applying the "maximum umber of returns to a unit" constraint.
                apply_max_num_returns_constraint(route, strata_list)

                if (len(strata_list) != 0):
                    # Copy the route to create references to it below.
                    old_path = route.path.copy()

                    # Looking to which strata unit we can change.
                    for strat in strata_list:
                        path_exists = strata_table[row, strat]

                        # Processing "foreign" litho constraints.
                        foreign_litho_added = can_add_foreign_lithology(path_exists, route, current_litho_index, strat)

                        if (path_exists or foreign_litho_added):
                            # Making the new route.
                            new_route = StrataRoute()
                            # New path contains the reference to the old path, and the new route position.
                            # Note: we are not copying the full old path, but only store a reference to it to save memory.
                            new_route.path = [old_path, strat]
                            new_route.current_thickness = thickness_change
                            new_route.num_units = route.num_units + 1
                            new_route.num_foreign_lithos = route.num_foreign_lithos
                            new_route.foreign_lithos = [x[:] for x in route.foreign_lithos] # Deep copy.
                            new_route.unit_visited = np.array(route.unit_visited)
                            # Count this unit as visited.
                            new_route.unit_visited[strat] += 1
                            # Processing the "foreign" lithology.
                            if (foreign_litho_added):
                                process_foreign_lithology(new_route, current_litho_index, strat)
                            # Adding new route into the list.
                            new_routes.append(new_route)

            #-----------------------------------------------------------------
            # Check if we can go down the same srata unit (if cannot -- remove the current route).
            #-----------------------------------------------------------------
            can_stay = True

            if (add_thickness_constraints):
                # Apply unit thickness constraints.
                can_stay = can_stay and (current_thickness < max_strata_thickness[strat0])

            path_exists = strata_table[row, strat0]

            # Processing "foreign" litho constraints.
            foreign_litho_added = can_add_foreign_lithology(path_exists, route, current_litho_index, strat0)

            can_stay = can_stay and (path_exists or foreign_litho_added)

            # Processing the route.
            if (can_stay):
                # Adding new route position.
                route.path.append(strat0)
                route.current_thickness += thickness_change
                # Processing the "foreign" lithology.
                if (foreign_litho_added):
                    process_foreign_lithology(route, current_litho_index, strat0)
            else:
            # Did not reach the end of a drillhole, and cannot go down the same unit.
                # Mark the route for removal.
                route.to_remove = True

        # Remove the routes marked for removal.
        all_routes = [route for route in all_routes if not route.to_remove]
        # Addig new routes.
        all_routes.extend(new_routes)

    # Adding the final number of routes.
    all_routes_number.append(len(all_routes))

    # Flatten the multilevel list of lists: convert [[[1,1],2,2],3,3] to [1,1,2,2,3,3].
    for route in all_routes:
        route.path = flatten(route.path)

    return all_routes, all_routes_number

#==============================================================================
def write_routes_to_file(filename, drillsample_data, all_routes):
    '''
    Writing stratigraphic routes to file.
    '''
    f = open(filename, "w")
    num_rows = len(drillsample_data)
    for row in range(num_rows):
        depth = float(drillsample_data[row][0])
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
def print_unique_routes(all_routes, num_print_paths):
    '''
    Print all unique routes (i.e., with unique strata sequence).
    '''
    unique_routes = set([])
    for route in all_routes:
        unique_routes.add(route.get_strata_sequence())

    print("Number of unique routes = ", len(unique_routes))
    if (num_print_paths > 0):
        num = 0
        for route in unique_routes:
            num += 1
            print(route)
            if (num >= num_print_paths):
                break

#=============================================================================
def plot_routes(drillsample_data, routes):
    '''
    Plot and display the routes.
    '''
    print("Plotting the routes...")

    # Using the "From" column.
    x_data = [float(d[0]) for d in drillsample_data]

    for route in routes:
        pl.plot(x_data, route.path, '.-')

    pl.xlabel('Depth')
    pl.ylabel('Strata unit index')
    pl.show()

#=============================================================================
def plot_unit_probabilities(all_routes, drillsample_data, num_units):
    '''
    Generate a plot with probability of occurence for each unit.
    '''
    if (len(all_routes) == 0):
        return

    num_rows = len(drillsample_data)
    strat_distr = np.zeros((num_rows, num_units))

    # Building the distribution of unit presence at every depth.
    for route in all_routes:
        for row in range(num_rows):
            unit_index = route.path[row]
            strat_distr[row, unit_index] += 1
    # Normalize.
    strat_distr = strat_distr / float(len(all_routes))

    #------------------------------------------
    # Plot distribution of the route scores (based on path probability).
    route_scores = np.zeros(len(all_routes), dtype=float)
    route_index = 0
    for route in all_routes:
        for row in range(num_rows):
            unit_index = route.path[row]
            route_scores[route_index] += strat_distr[row, unit_index]
        route_scores[route_index] = route_scores[route_index] / float(num_rows)
        route_index += 1

    title_params = 'Max foreign lithos = ' + str(max_num_foreign_lithos) + ', max returns = ' + str(max_num_returns)

    pl.hist(route_scores, bins = 50)
    pl.title(title_params)
    pl.xlabel('Route score')
    pl.ylabel('Frequency')
    pl.show()

    #------------------------------------------
    # Index of the route with highest score. (Note: can be several such routes)
    index_max = np.argmax(route_scores)
    print('Max index = ', index_max)
    print('Max score = ', route_scores[index_max])

    # Generating the plots.
    # Increasing the figure size.
    pl.rcParams["figure.figsize"] = (12.8, 9.6) # Default size = (6.4, 4.8)
    fig, axs = pl.subplots(num_units, sharey=True)
    fig.suptitle('Probability for each unit. ' + title_params)

    # Using the "From" column.
    x_data = [float(d[0]) for d in drillsample_data]

    for i in range(num_units):
        # Plot lines.
        axs[i].plot(x_data, strat_distr[:, i], zorder=1, c='lightblue')

        # Set red color for zero data.
        color = ['red' if d <= 0 else 'blue' for d in strat_distr[:, i]]

        # Set green color for the route with the highest score.
        for row in range(num_rows):
            if (all_routes[index_max].path[row] == i):
                color[row] = 'green'

        # Plot dots.
        axs[i].scatter(x_data, strat_distr[:, i], s=5, c=color, zorder=2)

        axs[i].set(ylabel=str(i))
        if (i != num_units - 1):
            # Hide tick labels.
            axs[i].set_xticklabels([])
        # Add vertical lines.
        axs[i].xaxis.grid(True)
 
    pl.xlabel('Depth')
    pl.show()

#=============================================================================
def main():
    print('Started litho2strat')

    # Paths to data.
    all_strat_filename = ""
    unit_list_filename = ""

#     strat_filename       = "data/simple/strat_1627022992.5507748.csv"
#     thickness_filename   = "data/simple/thickness_mean_1627022992.5507748.csv"
#     drillsample_filename = "data/simple/drill_sample_1627022992.5507748.csv"

    strat_filename       = "data/foreign_litho/strat_1627025194.6300328.csv"
    thickness_filename   = "data/foreign_litho/thickness_mean_1627025194.6300328.csv"
    drillsample_filename = "data/foreign_litho/drill_sample_1627025194.6300328.csv"

#     all_strat_filename   = "data/real/ALL_Strat descriptions.csv"
#     drillsample_filename = "data/real/MtGibson_drillhole.csv"
#     unit_list_filename   = "data/real/MtGibson_strat.csv"

#     drillsample_filename = "data/real/Hill_drillhole.csv"
#     unit_list_filename   = "data/real/Hill_strat.csv"

#     all_strat_filename = "data/test/ALL_strat.csv"
#     drillsample_filename = "data/test/drill.csv"
#     unit_list_filename = "data/test/ALL_strat.csv"

    #--------------------------------------------------------------
    new_file_format = False

    add_thickness_constraints = True

    # Flag for whether we should stay in the unit until the lithology name is changed.
    keep_continuous_lithos = False

    #--------------------------------------------------------------
    # Reading input data.
    #--------------------------------------------------------------
    # Drill sample data.
    drillsample_data = read_drillsample_data(drillsample_filename)

    # Unit lithologies data..
    strat_data = []
    if (new_file_format):
        strat_data = read_strat_data2(all_strat_filename, unit_list_filename)
    else:
        strat_data = read_strat_data(strat_filename)

    # Thickness data.
    thickness_data = []
    if (add_thickness_constraints):
        thickness_data = read_thickness_data(thickness_filename)

    #--------------------------------------------------------------
    # Generating the stratigraphy routes.
    tracemalloc.start()

    all_routes, all_routes_number = generate_strat_routes(strat_data, drillsample_data, thickness_data,
                                       add_thickness_constraints, keep_continuous_lithos)

    print("Total number of routes = ", len(all_routes))

    current, peak = tracemalloc.get_traced_memory()
    print("Current memory usage is {} MB; Peak was {} MB".format(current / 10**6, peak / 10**6))

    #--------------------------------------------------------------
    # Plot the number of processed routes at each row.
    pl.xlabel('Row number')
    pl.ylabel('Number of routes')
    pl.plot(all_routes_number)
    pl.show()

    # Print all unique routes (i.e., unique strata sequence).
    print_unique_routes(all_routes, 10)

    # Write results to the file.
    #write_routes_to_file("strata.txt", drillsample_data, all_routes)

    # Plot the routes.
    #plot_routes(drillsample_data, all_routes[0:1000])

    # Plot unit probabilities.
    num_units = len(strat_data)
    plot_unit_probabilities(all_routes, drillsample_data, num_units)

#=============================================================================
if __name__ == "__main__":
    main()

