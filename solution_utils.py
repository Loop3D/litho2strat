'''
  This program generates a set of plausible stratigraphies with uncertainties, for a given drillhole lithology log.
  It uses map data for distance and topology constraints, and several free parameters describing the solution complexity (level of deformation) constraints.

  Author: Vitaliy Ogarko, vogarko@gmail.com
  The University of Western Australia
'''

import numpy as np
import matplotlib.pyplot as pl
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.colors import ListedColormap
import matplotlib.colorbar as cbar
import networkx as nx
from dataclasses import dataclass
import csv
import os

output_folder = "output"

#==============================================================================
def draw_topology_graph(G, title):
    '''
    Drawing topology graph, with edges weighted by unit contact frequency.
    '''
    pos = nx.circular_layout(G)

    edges = nx.get_edge_attributes(G, 'weight')
    nodelist = G.nodes()

    pl.figure(figsize=(12, 8))

    # Draw the graph.
    nx.draw_networkx_nodes(G, pos,
                           nodelist=nodelist,
                           node_size=300,
                           node_color='black',
                           alpha=0.5)

    nx.draw_networkx_edges(G, pos,
                           edgelist=edges.keys(),
                           width=10.0,
                           edge_color=list(edges.values()),
                           edge_cmap=pl.cm.Blues,
                           edge_vmin=0.0,
                           edge_vmax=1.0,
                           alpha=0.8,
                           arrowsize=20)

    nx.draw_networkx_labels(G, pos=pos,
                            labels=dict(zip(nodelist,nodelist)),
                            font_color='black')

    # Set margins for the axes so that nodes aren't clipped.
    ax = pl.gca()
    ax.margins(0.20)

    ax.set_title(title)

    pl.axis("off")
    pl.show()

    # Write to graph to GML format.
    export_graph_to_GML(G, title)

#==============================================================================
def export_graph_to_GML(G, title):
    '''
    Write graph to file in GML format. For nicer visualisation in yEd.
    '''
    edges = nx.get_edge_attributes(G, 'weight')

    # Determine if all edge weights are the same.
    if len(set(edges.values())) == 1:
        all_weights_equal = True
    else:
        all_weights_equal = False

    # Define graph edge width.
    for u, v in G.edges:
        if (not all_weights_equal):
            G[u][v]["graphics"]["width"] = G[u][v]["weight"] * 7
        else:
            G[u][v]["graphics"]["width"] = 3

        G[u][v]["graphics"]["arrow"] = "last"

    # Adjust node style.
    for v in G.nodes():
        # Increase the node width.
        G.nodes[v]["graphics"] = {'w': 200., "fill": "#FFCC00"}
        # Adjust the node label.
        nodeLabel = v.replace("_", " ").title()
        G.nodes[v]["LabelGraphics"] = {'text': nodeLabel}

    # Form a file name.
    filename = output_folder + "/GML/" + str(title).replace(" ", "_") + ".gml"
    os.makedirs(os.path.dirname(filename), exist_ok=True)

    print("Write graph to:", filename)

    # Export the graph.
    nx.write_gml(G, filename)

#==============================================================================
def _get_unit_colors(strat_solution, unit_colors_filename):
    '''
    Define unit colors. Use a colormap from a file or a qualitative colormap.
    '''
    unit_colors = dict()
    if (unit_colors_filename != ""):
    # Use custom colors from a file.
        with open(unit_colors_filename, 'r') as csvfile:
            reader = csv.DictReader(csvfile, delimiter=',')
            # Extracting the data for every csv row.
            for row in reader:
                unit_name = row['UNITNAME']
                # Convert the unitname to align it with format used in the topology graph.
                unit_name = unit_name.replace(" ", "_").replace(",", "_").replace("-", "_").lower()
                if strat_solution.unit_nonempty(unit_name):
                    # Add colour to dictionary.
                    unit_colors[unit_name] = row['colour']
    else:
    # Use qualitative palette.
        colors = [pl.cm.tab20(i) for i in range(20)]

        # Define qualitative color map for nonempty units.
        color_index = 0
        for index, unit_name in enumerate(strat_solution.unit_names):
            if strat_solution.unit_nonempty(unit_name):
                unit_colors[unit_name] = colors[color_index]
                color_index += 1

        if (color_index > 20):
            print("Too many units! Adjust the color map.")
            exit(0)

    return unit_colors

#=========================================================================================
def _get_rectangular_color(type, unit_colors, strat_solution, route, row, graph, cmap):
    '''
    Retrieve the rectangular color.
    '''
    unit_index = route.path[row]

    if (type == 'strat' or type == 'strat-seq'):
    # Draw stratigraphy log.
        unit_name = strat_solution.unit_names[unit_index]
        if (unit_name in unit_colors):
            color = unit_colors[unit_name]
        else:
            print("WARNING: No color in the color map found for unit name =", unit_name)
            color = "#000000"

    elif (type == 'proba'):
    # Draw route probabilities.
        route_proba = strat_solution.strat_distr[row, unit_index]
        color = cmap(route_proba)

    elif (type == 'age'):
    # Draw age alignment.
        # Find the next unit in the log.
        unit_index2 = unit_index
        for j in range(row + 1, len(route.path)):
            unit_index2 = route.path[j]
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
        elif (graph.has_edge(*e) or unit_name == "cover"):
        # Unit contact is aligned with the age, or a cover.
            color = '#00FF00'
        else:
        # Not aligned - draw with red color.
            color = '#FF0000'

    return color

#==============================================================================
def draw_solution_logs(strat_solution, display_plot, type, unit_colors_filename,
                       sample_scores_uniformly, graph, custom_route_indexes, custom_route_scores=[], title=''):
    '''
    Drawing solution logs.
    '''
    if (type == 'age' and graph == None):
        # No graph supplied.
        return

    if (type == 'strat-seq'):
        routes = strat_solution.unique_routes
        route_scores = strat_solution.unique_route_scores
    else:
        routes = strat_solution.routes
        route_scores = strat_solution.route_scores

    # Select routes to display.
    if (len(custom_route_indexes) == 0):
        # Top scores (a minus here to have the largest-to-smallest score order).
        route_indexes = np.argsort(-route_scores)
    else:
        route_indexes = custom_route_indexes

    num_routes = len(route_indexes)
    if (num_routes == 0):
        return

    # Determine the number of routes to display (cannot show too many routes due to pixel size limitations).
    max_routes_displayed = 20
    if (num_routes > max_routes_displayed):
        num_routes_displayed = max_routes_displayed
    else:
        num_routes_displayed = num_routes

    print("Drawing solution logs, type =", type)

    if (num_routes_displayed == 1):
        fig_height = 4.8
    else:
        fig_height = 9.6
    pl.rcParams["figure.figsize"] = (12.8, fig_height) # Default size = (6.4, 4.8)

    # Gradient palette.
    cmap = pl.get_cmap('viridis')

    # Retrieve the unit colormap.
    unit_colors = None
    if (type == 'strat' or type == 'strat-seq'):
        unit_colors = _get_unit_colors(strat_solution, unit_colors_filename)

    # Calculate the figure size.
    if (type == 'strat-seq'):
        max_num_units = max(len(routes[i].path) for i in route_indexes[:num_routes_displayed])
        x_max = float(max_num_units) + 0.5
    else:
        x_max = strat_solution.depth_data.depth_to[-1]
    y_max = float(num_routes_displayed) + 0.5

    # Define figure dimensions.
    pl.xlim(0, x_max)
    pl.ylim(0.5, y_max)

    currentAxis = pl.gca()

    if (num_routes_displayed == 1):
    # Adjust figure aspect ratio for a single route.
        aspect = np.diff(currentAxis.get_xlim()) / np.diff(currentAxis.get_ylim())
        pl.gca().set_aspect(aspect / 20.)

    patches = []
    color_list = []

    yticklabels = []

    for i in range(num_routes_displayed):
        if (sample_scores_uniformly):
            # Sample the index uniformly (to show high and low score routes).
            ind = int(float(i) / float(num_routes_displayed) * float(num_routes))
        else:
            # Showing the top score routes.
            ind = i

        route_index = route_indexes[ind]
        route = routes[route_index]

        # Display route score on the y-axis.
        if (len(custom_route_scores) == 0):
            route_score = route_scores[route_index]
        else:
            route_score = custom_route_scores[ind]

        yticklabels.append(str(route_score)[0:5])

        for row in range(len(route.path)):
            if (type == 'strat-seq'):
                x1 = 0.5 + float(row)
                x2 = 0.5 + float(row + 1)
            else:
                x1 = strat_solution.depth_data.depth_from[row]
                x2 = strat_solution.depth_data.depth_to[row]

            y1 = 0.5 + float(i)
            y2 = 0.5 + float(i + 1)
            dx = x2 - x1
            dy = y2 - y1

            # Define the rectangle color.
            color = _get_rectangular_color(type, unit_colors, strat_solution, route, row, graph, cmap)

            # Adding rectangle.
            patches.append(Rectangle((x1, y1), dx, dy))
            color_list.append(color)

            if (type == 'strat-seq'):
                display_unit_name = True
            else:
                display_unit_name = False
            if (display_unit_name):
                pl.text(x1, y1, strat_solution.unit_names[route.path[row]][0:12])

    # Define patches collection with colormap.
    patches_cmap = ListedColormap(color_list)
    patches_collection = PatchCollection(patches, cmap=patches_cmap)
    patches_collection.set_array(np.arange(len(patches)))

    # Add rectangle collection to the figure.
    currentAxis.add_collection(patches_collection)

    # Add y-ticks.
    currentAxis.set_yticks([k + 1 for k in range(num_routes_displayed)])
    currentAxis.set_yticklabels(yticklabels)

    # Define labels and other things.
    skip_zero_cover = False
    add_colorbar = False
    add_legend = False
    if (type == 'strat'):
        default_title = "Most probable stratigraphies for collarID=" + str(strat_solution.collarID)
        xlabel = 'Depth'
        ylabel = 'Score'
        file_prefix = "strata_logs_"
        add_legend = True
        skip_zero_cover = True
    elif (type == 'strat-seq'):
        default_title = "Most probable strata sequences for collarID=" + str(strat_solution.collarID)
        xlabel = 'Unit number'
        ylabel = 'Score'
        file_prefix = "strata_seq_"
        add_legend = True
    elif (type == 'proba'):
        default_title = "Unit probability for collarID=" + str(strat_solution.collarID)
        xlabel = 'Depth'
        ylabel = 'Score'
        file_prefix = "proba_logs_"
        add_colorbar = True
    elif (type == 'age'):
        default_title = "Age alignment for collarID=" + str(strat_solution.collarID)
        xlabel = 'Depth'
        ylabel = 'Score'
        file_prefix = "age_logs_"

    if (title != ''):
        currentAxis.set_title(title)
    else:
        currentAxis.set_title(default_title)

    pl.xlabel(xlabel)
    pl.ylabel(ylabel)

    if (add_colorbar):
        # Show the colorbar.
        cax, _ = cbar.make_axes(currentAxis) 
        cb2 = cbar.ColorbarBase(cax, cmap=cmap)

    #-----------------------------------------------------------------
    # Adding legend.
    #-----------------------------------------------------------------
    if (add_legend):
        # Flag if the cover width is zero.
        zero_cover = strat_solution.depth_data.depth_to[0] == 0.0

        legend_patches = []
        for unit_name, unit_color in unit_colors.items():
            if (skip_zero_cover and zero_cover and unit_name == 'cover'):
                # Skip the Cover unit when it has zero width.
                continue
            patch = mpatches.Patch(color=unit_color, label=unit_name)
            legend_patches.append(patch)

        # Shrink current axis's height by 10% on the bottom.
        box = currentAxis.get_position()
        currentAxis.set_position([box.x0, box.y0 + box.height * 0.1,
                                  box.width, box.height * 0.9])

        if (num_routes_displayed == 1):
            by = -1.0
        else:
            by = -0.1
        pl.legend(handles=legend_patches, loc='upper center', bbox_to_anchor=(0.5, by),
                  fancybox=True, shadow=False, ncol=4)
    #-----------------------------------------------------------------

    # Save image.
    filename = output_folder + "/" + file_prefix + str(strat_solution.collarID) + ".png"
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    pl.savefig(filename)

    if display_plot:
        pl.show()

    pl.close(pl.gcf())

#==============================================================================
def print_unique_routes(solution, num_print_paths):
    '''
    Print all unique routes (i.e., with unique strata sequence).
    '''
    unique_routes = solution.unique_routes

    print("Number of unique routes = ", len(unique_routes))
    if (num_print_paths > 0):
        for index, route in enumerate(unique_routes):
            print(route.path)
            if (index >= num_print_paths - 1):
                break

#=============================================================================
def plot_route_scores(scores):
    '''
    Plot distribution of the route scores.
    '''
    # Set figure size.
    pl.rcParams["figure.figsize"] = (6.4, 4.8)

    pl.hist(scores, bins = 50)
    pl.xlabel('Score')
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

    if (np.std(x) != 0. and np.std(y) != 0.):
        # Pearson correlation coefficient.
        rho = np.corrcoef(x, y)[0][1]
        print("Correlation coeff rho =", rho)
    else:
        print("Correlation coeff rho is undefined!")

    # Set figure size.
    pl.rcParams["figure.figsize"] = (6.4, 6.4)

    pl.scatter(x, y, s=5)

    pl.show()

#=============================================================================
def plot_correlation_matrix(corr_matrix):
    '''
    Plot the correlation matrix.
    '''
    pl.matshow(corr_matrix, vmin=0., vmax=1.)
    pl.colorbar()
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

