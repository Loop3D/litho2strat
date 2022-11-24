'''
  This program generates a set of plausible stratigraphies with uncertainties, for a given drillhole lithology log.
  It uses map data for distance and topology constraints, and several free parameters describing the solution complexity (level of deformation) constraints.

  Author: Vitaliy Ogarko, vogarko@gmail.com
  The University of Western Australia
'''

import numpy as np
import matplotlib.pyplot as pl
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.colors import ListedColormap
import matplotlib.colorbar as cbar
import networkx as nx
import csv
import os

output_folder = "output"

#==============================================================================
def draw_solution_graph(strat_solution):
    '''
    Drawing solution topology graph, with edges weighted by unit contact frequency (among all solutions).
    '''
    G = strat_solution.graph

    pos = nx.spring_layout(G)

    edges = nx.get_edge_attributes(G, 'weight')
    nodelist = G.nodes()

    # Scale the edges weight with the minimum weight.
    min_weight = float(min(edges.values()))
    edges = {key: value / min_weight for key, value in edges.items()}

    pl.figure(figsize=(12, 8))

    # Draw the graph.
    nx.draw_networkx_nodes(G, pos,
                           nodelist=nodelist,
                           node_size=400,
                           node_color='black',
                           alpha=0.5)

    nx.draw_networkx_edges(G, pos,
                           edgelist=edges.keys(),
                           width=list(edges.values()),
                           edge_color='lightblue',
                           alpha=0.9)

    nx.draw_networkx_labels(G, pos=pos,
                            labels=dict(zip(nodelist,nodelist)),
                            font_color='black')

    # Set margins for the axes so that nodes aren't clipped.
    ax = pl.gca()
    ax.margins(0.20)
    pl.axis("off")
    pl.show()

#==============================================================================
def draw_solution_logs(strat_solution, display_plot, type, unit_colors_filename,
                       sample_scores_uniformly, graph):
    '''
    Drawing solution logs.
    '''
    num_routes = len(strat_solution.routes)
    if (num_routes == 0):
        return

    if (type == 'age' and graph == None):
        # No graph supplied.
        return

    # Determine the number of routes to display (cannot show too many routes due to pixel size limitations).
    max_routes_displayed = 100
    if (num_routes > max_routes_displayed):
        num_routes_displayed = max_routes_displayed
    else:
        num_routes_displayed = num_routes

    print("Drawing solution logs, type =", type)

    pl.rcParams["figure.figsize"] = (12.8, 9.6) # Default size = (6.4, 4.8)

    # Qualitative palette.
    colors = [pl.cm.tab20(i) for i in range(20)]
    # Gradient palette.
    cmap = pl.get_cmap('viridis')

    if (type == 'strat'):
        if (unit_colors_filename != ""):
        # Define colors using the colormap provided in the file.
            unit_colors = dict()
            with open(unit_colors_filename, 'r') as csvfile:
                reader = csv.DictReader(csvfile, delimiter=',')
                # Extracting the data for every csv row.
                for row in reader:
                    unit_name = row['UNITNAME']
                    # Convert the unitname to align it with format used in the topology graph.
                    unit_name = unit_name.replace(" ", "_").replace(",", "_").replace("-", "_").lower()
                    # Add colour to dictionary.
                    unit_colors[unit_name] = row['colour']
        else:
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
    x_max = strat_solution.depth_data.depth_to[-1]
    y_max = float(num_routes_displayed) + 0.5

    # Define figure dimensions.
    pl.xlim(0, x_max)
    pl.ylim(0.5, y_max)

    currentAxis = pl.gca()
    currentAxis.set_title(str(strat_solution.collarID))

    # Top scores (a minus here to have the largest-to-smallest score order).
    indexes_max = np.argsort(-strat_solution.route_scores)

    patches = []
    color_list = []

    for i in range(num_routes_displayed):
        if (sample_scores_uniformly):
            # Sample the index uniformly (to show high and low score routes).
            ind = int(float(i) / float(num_routes_displayed) * float(num_routes))
        else:
            # Showing the top score routes.
            ind = i

        # Select the sorted by score solution index.
        route_index = indexes_max[ind]

        #print("score = ", strat_solution.route_scores[route_index])

        path_size = len(strat_solution.routes[route_index].path)
        for row in range(path_size):
            unit_index = strat_solution.routes[route_index].path[row]

            x1 = strat_solution.depth_data.depth_from[row]
            x2 = strat_solution.depth_data.depth_to[row]
            y1 = 0.5 + float(i)
            y2 = 0.5 + float(i + 1)
            dx = x2 - x1
            dy = y2 - y1

            # Define the rectangle color.
            if (type == 'strat'):
            # Draw stratigraphy log.
                if (unit_colors_filename != ""):
                    unit_name = strat_solution.unit_names[unit_index]
                    if (unit_name in unit_colors):
                        color = unit_colors[unit_name]
                    else:
                        print("WARNING: No color in the color map found for unit name =", unit_name)
                        color = "#000000"
                else:
                    color_index = unit_index_nonempty[unit_index]
                    color = colors[color_index]

            elif (type == 'proba'):
            # Draw route probabilities.
                route_proba = strat_solution.strat_distr[row, unit_index]
                color = cmap(route_proba)

            elif (type == 'age'):
            # Draw age alignment.
                # Find the next unit in the log.
                unit_index2 = unit_index
                for j in range(row + 1, path_size):
                    unit_index2 = strat_solution.routes[route_index].path[j]
                    # Returns the first unit change.
                    if (unit_index2 != unit_index):
                        break

                last_unit = False
                if (unit_index2 == unit_index):
                    last_unit = True

                unit_name = strat_solution.unit_names[unit_index]
                unit_name2 = strat_solution.unit_names[unit_index2]

                # Graph edge.
                e = (unit_name, unit_name2)

                if (last_unit):
                # Mark the last unit with black color.
                    color = "#000000"
                elif (graph.has_edge(*e)):
                # Unit contact is aligned with the age.
                    color = '#00FF00'
                else:
                # Not aligned - draw with red color.
                    color = '#FF0000'

            # Adding rectangle.
            patches.append(Rectangle((x1, y1), dx, dy))
            color_list.append(color)

    # Define patches collection with colormap.
    patches_cmap = ListedColormap(color_list)
    patches_collection = PatchCollection(patches, cmap=patches_cmap)
    patches_collection.set_array(np.arange(len(patches)))

    # Add rectangle collection to the figure.
    currentAxis.add_collection(patches_collection)

    if (type == 'strat'):
        ylabel = 'Stratigraphy'
        file_prefix = "strata_logs_"
        add_colorbar = False
    elif (type == 'proba'):
        ylabel = 'Probability'
        file_prefix = "proba_logs_"
        add_colorbar = True
    elif (type == 'age'):
        ylabel = 'Age alignment'
        file_prefix = "age_logs_"
        add_colorbar = False

    pl.xlabel('Depth')
    pl.ylabel(ylabel)

    if (add_colorbar):
        # Show the colorbar.
        cax, _ = cbar.make_axes(currentAxis) 
        cb2 = cbar.ColorbarBase(cax, cmap=cmap)

    # Save image.
    filename = output_folder + "/" + file_prefix + str(strat_solution.collarID) + ".png"
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
def plot_route_scores(scores):
    '''
    Plot distribution of the route scores.
    '''
    pl.hist(scores, bins = 50)
    pl.xlabel('Route score')
    pl.ylabel('Frequency')
    pl.show()

#=============================================================================
def plot_solution_correlation(solution):
    '''
    Plot the correlation of solution scores from different drillholes.
    '''
    if (len(solution.external_graph_route_scores_list) == 0):
        return

    x = solution.graph_route_scores
    y = solution.external_graph_route_scores_list[0]

    # Pearson correlation coefficient.
    rho = np.corrcoef(x, y)[0][1]
    print("Correlation coeff rho = ", rho)

    pl.scatter(x, y, s=5)

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

    title = "Probability of occurrence for every unit. CollarID = " + str(strat_solution.collarID)
    fig.suptitle(title, y=0.96)

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
        y_data = strat_distr[:, i]

        # Plot lines.
        axs[j, 0].plot(x_data, y_data, zorder=1, c='blue')

        # Plot missing data segments.
        for k in range(1, num_rows):
            pair_x = [x_data[2 * k - 1], x_data[2 * k]]
            pair_y = [y_data[2 * k - 1], y_data[2 * k]]
            if (pair_x[1] > pair_x[0]):
                # Plot line segment.
                axs[j, 0].plot(pair_x, pair_y, zorder=2, c='red')

        # Set red color for zero data.
        #color = ['red' if p <= 0 else 'blue' for p in strat_distr[:, i]]
        color = 'blue'

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

