'''
  This program generates a set of plausible stratigraphies with uncertainties, for a given drillhole lithology log.
  It uses map data for distance and topology constraints, and several free parameters describing the solution complexity (level of deformation) constraints.

  Author: Vitaliy Ogarko, vogarko@gmail.com
  The University of Western Australia
'''

import numpy as np
import matplotlib.pylab as pl

#==============================================================================
def write_routes_to_file(filename, drillsample_data, all_routes):
    '''
    Writing stratigraphic routes to file.
    '''
    f = open(filename, "w")
    num_rows = drillsample_data.get_num_rows()
    for row in range(num_rows):
        depth = drillsample_data.rows[row].depth_from
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
def plot_unit_probabilities(strat_solution):
    '''
    Generate a plot with probability of occurence for each unit.
    '''
    if (len(strat_solution.routes) == 0):
        return

    # Increasing the figure size.
    pl.rcParams["figure.figsize"] = (12.8, 9.6) # Default size = (6.4, 4.8)

    num_units = len(strat_solution.unit_names)

    # Count the number of non-empty units.
    num_units_nonempty = 0
    for i in range(num_units):
        if (sum(strat_solution.strat_distr[:, i]) != 0):
            num_units_nonempty += 1

    fig, axs = pl.subplots(num_units_nonempty, sharey=True, squeeze=True)

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
    for i in range(num_units):
        if (sum(strat_distr[:, i]) != 0):
            nonempty_units.append(i)

    j = 0
    for index, i in enumerate(nonempty_units):
        # Plot lines.
        axs[j].plot(x_data, strat_distr[:, i], zorder=1, c='blue')

        # Set red color for zero data.
        color = ['red' if p <= 0 else 'blue' for p in strat_distr[:, i]]

        # Plot dots.
        axs[j].scatter(x_data, strat_distr[:, i], s=5, c=color, zorder=2)

        axs[j].set_title(strat_solution.unit_names[i], size=9, y=0.97)
        axs[j].set_ylabel(str(j))

        if (index != len(nonempty_units) - 1):
            # Hide tick labels.
            axs[j].set_xticklabels([])

        # Add vertical lines.
        axs[j].xaxis.grid(True)

        j += 1
 
    #pl.tight_layout()
    pl.subplots_adjust(hspace = 0.5)
 
    pl.xlabel('Depth')
    pl.show()

#=============================================================================

