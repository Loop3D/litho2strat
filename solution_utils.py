'''
  This program generates a set of plausible stratigraphies with uncertainties, for a given drillhole lithology log.
  It uses map data for distance and topology constraints, and several free parameters describing the solution complexity (level of deformation) constraints.

  Author: Vitaliy Ogarko, vogarko@gmail.com
  The University of Western Australia
'''

import numpy as np
import matplotlib.pylab as pl
from matplotlib.patches import Rectangle
import os

output_folder = "output"

#==============================================================================
def draw_strata_logs(strat_solution, display_plot):
    '''
    Drawing the strata logs.
    '''
    print("Drawing the strata logs...")

    # Define qualitative palette.
    colors = [pl.cm.tab20(i) for i in range(20)]

    # Map nonempty units to continous index.
    # We use that instead of original index, as we do not have many qualitative colors in the colormap.
    counter = 0
    unit_index_nonempty = dict()
    for index, unit_name in enumerate(strat_solution.unit_names):
        if strat_solution.unit_nonempty(unit_name):
            unit_index_nonempty[index] = counter
            counter += 1

    print("Num nonempty units:", counter)
    if (counter > 20):
        print("Too many units! Adjust the color map.")
        return

    # Calculate the figure size.
    num_routes = min(len(strat_solution.routes), 200)
    x_max = strat_solution.depth_data.depth_to[-1]
    y_max = float(num_routes) + 0.5

    # Define figure dimensions.
    fig = pl.figure()
    pl.xlim(0, x_max)
    pl.ylim(0.5, y_max)

    currentAxis = pl.gca()

    for i in range(num_routes):
        for row, unit_index in enumerate(strat_solution.routes[i].path):
            x1 = strat_solution.depth_data.depth_from[row]
            x2 = strat_solution.depth_data.depth_to[row]
            y1 = 0.5 + float(i)
            y2 = 0.5 + float(i + 1)
            dx = x2 - x1
            dy = y2 - y1

            # Adding rectangle.
            color_index = unit_index_nonempty[unit_index]
            currentAxis.add_patch(Rectangle((x1, y1), dx, dy, facecolor=colors[color_index]))

    pl.xlabel('Depth')
    pl.ylabel('Stratigraphy')

    # Save image.
    filename = output_folder + "/strata_logs_" + str(strat_solution.collarID) + ".png"
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    pl.savefig(filename)

    if display_plot:
        pl.show()

    pl.close(pl.gcf())

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
def plot_routes(routes, x_data):
    '''
    Plot and display the routes.
    '''
    for route in routes:
        pl.plot(x_data, route.path, '.-')

    pl.xlabel('Depth')
    pl.ylabel('Strata unit index')
    pl.show()

#=============================================================================
def plot_route_probability(routes, x_data, strat_distr):
    '''
    Plot and display the route probability.
    '''
    num_rows = len(routes[0].path)

    for route in routes:
        route_proba = np.zeros(num_rows)
        for row in range(num_rows):
            unit_index = route.path[row]
            route_proba[row] = strat_distr[row, unit_index]
        pl.plot(x_data, route_proba, '.-')

    pl.xlabel('Depth')
    pl.ylabel('Probability')
    pl.show()

#=============================================================================
def plot_route_scores(strat_solution):
    '''
    Plot distribution of the route scores (based on path probability).
    '''
    pl.hist(strat_solution.route_scores, bins = 50)
    pl.xlabel('Route score')
    pl.ylabel('Frequency')
    pl.show()

#==============================================================================
def write_best_routes_to_file(strat_solution, ntop):
    '''
    Write the best ntop routes to file.
    '''
    if (len(strat_solution.route_scores) == 0):
        return

    filename = output_folder + "/best_routes_" + str(strat_solution.collarID) + ".txt"
    os.makedirs(os.path.dirname(filename), exist_ok=True)

    # Extract the indexes of the best routes.
    route_indexes = np.argsort(-strat_solution.route_scores)
    num_routes = len(route_indexes)

    if (ntop > num_routes):
        # Adjust ntop if we have less routes.
        ntop = num_routes

    with open(filename, "w") as file:
        # Write the number of units.
        num_units = strat_solution.num_nonempty_units()
        file.write("%d\n" % num_units)

        # Write unit names.
        for index, unit_name in enumerate(strat_solution.unit_names):
            if (strat_solution.unit_nonempty(unit_name)):
                file.write("%d,%s\n" % (index, unit_name))

        num_rows = len(strat_solution.depth_data.depth_from)

        file.write("%d,%d,%d\n" % (num_rows, ntop, num_routes))

        # Write stratigraphy.
        for row in range(num_rows):
            # Extract depths for this row.
            depth_from = strat_solution.depth_data.depth_from[row]
            depth_to = strat_solution.depth_data.depth_to[row]

            # Write depth data.
            file.write("%f,%f" % (depth_from, depth_to))

            # Write unit indexes.
            for route_index in route_indexes[0:ntop]:
                unit_index = strat_solution.routes[route_index].path[row]
                file.write(",%d" % unit_index)

            # Write probabilities.
            for route_index in route_indexes[0:ntop]:
                unit_index = strat_solution.routes[route_index].path[row]
                route_proba = strat_solution.strat_distr[row, unit_index]
                file.write(",%.3f" % route_proba)
            file.write("\n")

#=============================================================================
def plot_top_routes(strat_solution):
    '''
    Plot top routes and their probability.
    '''
    if (len(strat_solution.routes) == 0):
        return

    num_units = len(strat_solution.unit_names)
    route_scores = strat_solution.route_scores

    #------------------------------------------
    # Top scores.
    indexes_max = np.argsort(-route_scores) # A minus here to have largest to smallest score order.
    ntop = 10
    print('Top indexes: ', indexes_max[0:ntop])
    print('Top scores: ', route_scores[indexes_max[0:ntop]])

    index_max = indexes_max[0]

    #------------------------------------------
    top_routes = [strat_solution.routes[i] for i in indexes_max[0:ntop]]
    x_data = strat_solution.depth_data.depth_from

    # Print the most probable routes.
    plot_routes(top_routes, x_data)

    # Print the probability of the most probable routes.
    plot_route_probability(top_routes, x_data, strat_solution.strat_distr)

#=============================================================================
def plot_unit_probabilities(strat_solution, display_plot):
    '''
    Generate a plot with probability of occurence for each unit.
    '''
    if (len(strat_solution.routes) == 0):
        return

    # Increasing the figure size.
    pl.rcParams["figure.figsize"] = (12.8, 9.6) # Default size = (6.4, 4.8)

    num_units = len(strat_solution.unit_names)
    num_units_nonempty = strat_solution.num_nonempty_units()

    fig, axs = pl.subplots(nrows=num_units_nonempty, ncols=1, sharey=True, squeeze=False)

    fig.suptitle('Probability of occurrence for every unit.', y=0.96)

    num_rows = len(strat_solution.routes[0].path)

    #-------------------------------------------------------------
    # Adding the "From" and "To" depths for visualisation.
    x_data = []
    for i in range(num_rows):
        x_data.append(strat_solution.depth_data.depth_from[i])
        x_data.append(strat_solution.depth_data.depth_to[i])

    # Duplicate each value, as the probability is the same between the "From" and "To" depths.
    strat_distr = np.repeat(strat_solution.strat_distr, 2, axis=0)
    #-------------------------------------------------------------

    # Skip empty units.
    nonempty_units = []
    for index, unit_name in enumerate(strat_solution.unit_names):
        if (strat_solution.unit_nonempty(unit_name)):
            nonempty_units.append(index)

    j = 0
    for index, i in enumerate(nonempty_units):
        # Plot lines.
        axs[j, 0].plot(x_data, strat_distr[:, i], zorder=1, c='blue')

        # Set red color for zero data.
        color = ['red' if p <= 0 else 'blue' for p in strat_distr[:, i]]

        # Plot dots.
        axs[j, 0].scatter(x_data, strat_distr[:, i], s=5, c=color, zorder=2)

        axs[j, 0].set_title(strat_solution.unit_names[i], size=9, y=0.97)
        axs[j, 0].set_ylabel(str(j))

        if (index != len(nonempty_units) - 1):
            # Hide tick labels.
            axs[j, 0].set_xticklabels([])

        # Add vertical lines.
        axs[j, 0].xaxis.grid(True)

        j += 1
 
    #pl.tight_layout()
    pl.subplots_adjust(hspace = 0.5)
    pl.xlabel('Depth')

    # Save image.
    filename = output_folder + "/unit_proba_" + str(strat_solution.collarID) + ".png"
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    pl.savefig(filename)

    if display_plot:
        pl.show()

    pl.close(pl.gcf())

#=============================================================================

