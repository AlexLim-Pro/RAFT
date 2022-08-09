#!/usr/bin/env python3
#******************************************************************************
#rrr_swo_flo.py
#******************************************************************************

# Purpose:
# Given a river data coordinate csv file, a river connectivity csv file,
# a river network shapefile, and a river model output netCDF file, this
# program creates an interactive model to better understand and communicate the
# propagation of water through space and time within rivers.
# Authors:
# Alex Christopher Lim, Cedric H. David, 2022-2022


##### Import Python modules #####
import csv
import itertools
import sys
import time
from datetime import datetime

import matplotlib as mpl
import matplotlib.backend_bases
import matplotlib.font_manager as font_manager
import matplotlib.pyplot as plt
import matplotlib.style as mplstyle
import matplotlib.widgets
import numpy as np
import progressbar
import shapefile as shp
import statsmodels.api as sm
from scipy.io import netcdf_file

from rrr_msc_fnc import DistFncs


##### Styling the UI #####
### GUI Styling ###
mplstyle.use("fast")
mplstyle.use("seaborn-darkgrid")
mpl.rcParams["path.simplify"] = True
mpl.rcParams["path.simplify_threshold"] = 1.0
mpl.rcParams["text.color"] = "#FFFFFF"
mpl.rcParams["axes.labelcolor"] = "#FFFFFF"
mpl.rcParams["xtick.color"] = "#FFFFFF"
mpl.rcParams["ytick.color"] = "#FFFFFF"
mpl.rcParams["axes.facecolor"] = "#333333"
mpl.rcParams["patch.facecolor"] = "#F0852D"
mpl.rcParams["figure.facecolor"] = "#333333"
mpl.rcParams["boxplot.flierprops.marker"] = "o"
font_path = "./assets/fonts/Inter.ttf"
font_manager.fontManager.addfont(font_path)
prop = font_manager.FontProperties(fname=font_path)
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = prop.get_name()

widget_bg = "#F0852D"
widget_bg_active = "#F8B723"
quit_bg = "#FF0000"
quit_bg_active = "#CC0202"
peak_duration_color = "#82EEFD"
default_point_color = "#FFFFFF"
cycol = itertools.cycle([
    "#CF0DFF",
    "#4E0CE8",
    "#0033FF",
    "#0CA2E8",
    "#00FFD0",
])

### UI Styling ###
progressbarWidgets = [
    "    [",
    progressbar.Timer(format="Elapsed Time: %(elapsed)s"),
    "] ",
    progressbar.GranularBar(), " (",
    progressbar.ETA(), ") ",
    progressbar.Percentage(),
    " "
]
"""Standard widgets for use in progress bars"""


##### Miscellaneous Globals #####
coords_f_path = "../San_Gaud_data/coords_San_Guad.csv"
connectivity_f_path = "../San_Gaud_data/rapid_connect_San_Guad.csv"
sf_path = "../San_Gaud_data/NHDFlowline_San_Guad/NHDFlowline_San_Guad.shp"
Qout_f_path = "../San_Gaud_data/Qout_San_Guad_exp00.nc"

show_errors = True
mpl_backend_bases = matplotlib.backend_bases
fig, ax = plt.subplots(label="River Network")
coords = list()
coords_dict = dict()
id_ind = dict()
ind_id = dict()
connectivity = dict()
connectivity_rev = dict()
p = dict()
rivers = dict()
x_vals = dict()
x_vals_list = list()
y_vals = dict()
y_vals_list = list()
river_lengths = dict()
"""The length of each river in degrees"""
discharge_graph_rivers = dict()
river_colors = dict()
num_downstream = 5
currently_selected_river = ""
downstream_rivers_list = list()
truncated_downstream = False
select_multiple = False
selected_rivers_list = list()
offsets = list()
open_in_new_window = False
skip_rerender = False
max_river_selections = 10
enlarged_rivers = list()
hidden_alpha = 0.1
point_scaling = 5
num_reaches_str = "5"
num_reaches = int(num_reaches_str)
reach_dist_str = "5"
reach_dist = float(reach_dist_str)
reach_dist_units = "km"
have_default_view = False
have_cleared_view = False
shown_rivers_list = list()
default_xlims = list()
default_ylims = list()
threshold_level = "1"
sleep_time = 0.01
showing_river_discharges = False


##### Declaration of variables (given as command line arguments) #####
# (1) - rrr_coords_file
# (2) - rrr_connect_file
# (3) - rrr_shp_file
# (4) - rrr_Qout_file


##### Get command line arguments #####
IS_arg = len(sys.argv)
if IS_arg > 5:
    print("Error: A maximum of 4 arguments can be used.")
    raise SystemExit(22)

if IS_arg > 1:
    coords_f_path = sys.argv[1]
    if IS_arg > 2:
        connectivity_f_path = sys.argv[2]
        if IS_arg > 3:
            sf_path = sys.argv[3]
            if IS_arg > 4:
                Qout_f_path = sys.argv[4]

Qout_f = netcdf_file(Qout_f_path, "r")

Qout_data = Qout_f.variables["Qout"][:]*1
"""average river water discharge downstream of each river reach"""

Qout_data_normed = (Qout_data.T / abs(Qout_data).max(axis=1)).T
"""normalized average river water discharge downstream of each river reach"""

# y_vals_array = Qout_f.variables["lat"][:]*1
# """latitude of river data points"""
#
# x_vals_array = Qout_f.variables["lon"][:]*1
# """longitude of river data points"""

Qout_data_ids = Qout_f.variables["rivid"][:]*1
"""unique identifier for each river reach"""

Qout_data_max = max(abs(Qout_data[0]))

num_rows = 0
with open(coords_f_path, newline="\n") as f:
    num_rows += sum(1 for row in csv.reader(f, delimiter=","))

i = 0
with open(coords_f_path, newline="\n") as f:
    print("Loading data points:")
    time.sleep(sleep_time)
    bar = progressbar.ProgressBar(
        maxval=num_rows,
        widgets=progressbarWidgets).start()
    for row in csv.reader(f, delimiter=","):
        x_vals[float(row[-2].replace(" ", ""))] = row[0]
        y_vals[float(row[-1].replace(" ", ""))] = row[0]
        x_vals_list.append(float(row[-2].replace(" ", "")))
        y_vals_list.append(float(row[-1].replace(" ", "")))
        coords.append([float(row[-2]), float(row[-1])])
        coords_dict[row[0]] = [float(row[-2]), float(row[-1])]
        idx = np.where(Qout_data_ids == int(row[0]))
        size = Qout_data_normed[0][idx] * point_scaling
        color = default_point_color
        p[row[0]] = plt.scatter(float(row[-2]), float(row[-1]),
                                s=size, picker=5, c=color)
        id_ind[row[0]] = i
        ind_id[i] = row[0]
        i += 1
        bar.update(i)
        sys.stdout.flush()
time.sleep(sleep_time)
bar.finish()
time.sleep(sleep_time)
print("Finished loading data points")


x_vals_array = np.array(x_vals_list)
y_vals_array = np.array(y_vals_list)

num_rows = 0
with open(connectivity_f_path, newline="\n") as f:
    num_rows += sum(1 for row in csv.reader(f, delimiter=","))

with open(connectivity_f_path, newline="\n") as f:
    print("Loading connectivity information:")
    time.sleep(sleep_time)
    bar = progressbar.ProgressBar(maxval=num_rows,
                                  widgets=progressbarWidgets).start()
    j = 0
    for row in csv.reader(f, delimiter=","):
        data = list()
        data_rev = list()
        for d in row[1:]:
            data.append(int(d))
        connectivity[row[0].replace(" ", "")] = data
        k = row[0].replace(" ", "")
        if k in connectivity_rev:
            connectivity_rev[k] = connectivity_rev[k] + data[1:]
        else:
            connectivity_rev[k] = data[1:]
        j += 1
        bar.update(j)
        sys.stdout.flush()
time.sleep(sleep_time)
bar.finish()
time.sleep(sleep_time)
print("Finished loading connectivity information")


sf = shp.Reader(sf_path)
num_rows = len(sf.shapeRecords())
print("Loading river morphology:")
time.sleep(sleep_time)
bar = progressbar.ProgressBar(
        maxval=num_rows,
        widgets=progressbarWidgets).start()
j = 0
for shape in sf.shapeRecords():
    xy = [k for k in shape.shape.points[:]]
    id = ""
    x, y = zip(*[(k[0], k[1]) for k in xy])
    nearest_vals = list()
    for xj in x:
        temp_array = np.abs(x_vals_array - float(xj))
        close = np.isclose(temp_array, np.zeros(temp_array.shape),
                           atol=1e-13)
        if np.any(close):
            idx = np.abs(temp_array).argmin()
            id = x_vals[x_vals_list[idx]]
            break
    if id == "":
        for yj in y:
            temp_array = y_vals_array - float(yj)
            close = np.isclose(temp_array, np.zeros(temp_array.shape),
                               atol=1e-13)
            if np.any(close):
                idx = np.abs(temp_array).argmin()
                id = y_vals[y_vals_list[idx]]
                break
    if id == "":
        continue
    elif id in rivers:
        continue
    i = 0
    for k in xy[1:]:
        d = DistFncs.great_circle_dist(k[0], xy[i][0], k[1], xy[i][1])
        if id in river_lengths:
            river_lengths[id] += d
        else:
            river_lengths[id] = d
        i += 1
    color = next(cycol)
    rivers[id] = plt.plot(x, y, linewidth=0.5, alpha=1, c=color)
    river_colors[str(id)] = color
    j += 1
    bar.update(j)
    sys.stdout.flush()
time.sleep(sleep_time)
bar.finish()
time.sleep(sleep_time)
print("Finished loading river morphology")


##### River Network Window #####
plt.title("Please click on one river reach")
plt.xlabel("Longitude")
plt.ylabel("Latitude")

### Reset Option ###
ax_reset = plt.axes([0.875, 1 - 0.05 - 0.075 * 2 - 0.01, 0.1, 0.075])
b_reset = plt.Button(
    ax_reset,
    label="Reset",
    color=widget_bg,
    hovercolor=widget_bg_active,
)

### Quit Option ###
ax_quit = plt.axes([0.875, 1 - 0.05 - 0.075, 0.1, 0.075])
b_quit = plt.Button(
    ax_quit,
    label="Quit",
    color=quit_bg,
    hovercolor=quit_bg_active,
)


print("Loading river discharge data:")
time.sleep(sleep_time)
bar = progressbar.ProgressBar(
        maxval=Qout_data.shape[1],
        widgets=progressbarWidgets).start()
fig_graph, ax_graph = plt.subplots(label="fig_graph")
for s in range(Qout_data.shape[1]):
    x = list(range(Qout_data.shape[0]))
    y = Qout_data[:, s]
    if str(Qout_data_ids[s]) not in river_colors:
        continue
    id = str(Qout_data_ids[s])
    color = river_colors[id]
    discharge_graph_rivers[id] = plt.plot(x, y,
                                          linewidth=0.5,
                                          alpha=0.5,
                                          c=color)
    bar.update(s)
    sys.stdout.flush()
time.sleep(sleep_time)
bar.finish()
time.sleep(sleep_time)
print("Finished loading river discharge data")


plt.title("River Discharge Over Time")
plt.xlabel("Time (3 hours)")
plt.ylabel("Average River Discharge (m^3/s)")


def on_pick(event):
    """
    Handles selection of points.

    - Callback function

    :param event: Automatically generated by the event manager.
    :type event: mpl_backend_bases.PickEvent
    """
    global currently_selected_river
    global selected_rivers_list
    global offsets
    offsets_temp = event.artist.get_offsets()[event.ind][0]
    if np.any(offsets):
        if np.all(np.isclose(offsets_temp, offsets, rtol=1e-2)):
            return
    offsets = offsets_temp
    i = 0
    for c in coords:
        if c[0] == offsets[0] and c[1] == offsets[1]:
            break
        i += 1
    id = str(ind_id[i])
    if select_multiple:
        print("Selecting river", id)
        selected_rivers_list.append(id)
        return
    currently_selected_river = id
    paths = list()
    done_paths = list()
    start_time = datetime.now()
    print("Showing downstream for", id)
    paths.append(connectivity[id][0])
    clear_plots()
    p[str(id)].set_alpha(1)
    selected_rivers_list.append(str(id))
    if str(id) in rivers:
        j = 0
        for s in rivers[str(id)]:
            rivers[str(id)][j].set_alpha(1)
            j += 1
    done_paths.append(str(id))
    while paths:
        curr_id = str(paths.pop())
        if str(curr_id) not in p:
            continue
        paths.append(str(connectivity[str(curr_id)][0]))
        downstream_rivers_list.append(str(connectivity[str(curr_id)][0]))
        p[str(curr_id)].set_alpha(1)
        selected_rivers_list.append(str(curr_id))
        done_paths.append(str(curr_id))
    paths, done_paths = done_paths, paths
    while paths:
        try:
            curr_id = str(paths.pop())
            if str(curr_id) not in rivers:
                continue
            j = 0
            for s in rivers[str(curr_id)]:
                try:
                    rivers[str(curr_id)][j].set_alpha(1)
                except:
                    pass
                j += 1
            done_paths.append(str(curr_id))
        except:
            pass
    paths, done_paths = done_paths, paths
    while paths:
        curr_id = str(paths.pop())
        if str(curr_id) not in discharge_graph_rivers:
            continue
        j = 0
        for s in discharge_graph_rivers[str(curr_id)]:
            if str(curr_id) in discharge_graph_rivers:
                discharge_graph_rivers[str(curr_id)][j].set_alpha(1)
            j += 1
        done_paths.append(str(curr_id))
    print("Took",
          (datetime.now() - start_time).total_seconds(),
          "seconds to select")
    create_config_window()
    redraw_canvases()
    print("Took",
          (datetime.now() - start_time).total_seconds(),
          "seconds to select and draw")


def apply_multiple(*args, **kwargs):
    """
    Shows the selected rivers
    
    :param args: Unused parameter to allow function to work as a callback
    :param kwargs: Unused parameter to allow function to work as a callback
    """
    update_idle()
    print("Drawing graphs")
    i = 0
    for r in selected_rivers_list:
        fig_temp = plt.figure()
        plt.title("River " + str(r) + " Discharge Over Time")
        plt.xlabel("Time (3 hours)")
        plt.ylabel(r"Average River Discharge (m^3/s)")
        if str(r) in discharge_graph_rivers:
            j = 0
            for s in discharge_graph_rivers[str(r)]:
                r_temp = discharge_graph_rivers[str(r)][j]
                plt.plot(*r_temp.get_data(),
                         c=r_temp.get_color(),
                         linewidth=r_temp.get_linewidth(),
                         alpha=r_temp.get_alpha())
                j += 1
        if i > max_river_selections:
            break
        i += 1


def show_together(*args, **kwargs):
    """
    Shows the selected rivers in the same window

    :param args: Unused parameter to allow function to work as a callback
    :param kwargs: Unused parameter to allow function to work as a callback
    """
    print("Drawing graphs")
    i = 0
    for r in selected_rivers_list:
        if str(r) in discharge_graph_rivers:
            fig_temp = plt.figure()
            plt.title("River " + str(r) + " Discharge Over Time")
            plt.xlabel("Time (3 hours)")
            plt.ylabel(r"Average River Discharge ($m^3/s$)")
            j = 0
            for s in discharge_graph_rivers[str(r)]:
                r_temp = discharge_graph_rivers[str(r)][j]
                plt.plot(*r_temp.get_data(),
                         c=r_temp.get_color(),
                         linewidth=r_temp.get_linewidth(),
                         alpha=r_temp.get_alpha(),
                         label=str(r))
                j += 1
            plt.legend(loc="upper right")
        if i > max_river_selections:
            break
        i += 1


def set_threshold_level(*args, **kwargs):
    """
    Updates the number of river reaches to show

    :param args: Unused parameter to allow function to work as a callback
    :param kwargs: Unused parameter to allow function to work as a callback
    """
    global threshold_level
    threshold_level = b_threshold_level.text
    if float(threshold_level) > 100:
        print("Error: Cannot set threshold level above 100%.")
        threshold_level = "100"
    elif float(threshold_level) < 0:
        print("Error: Cannot set threshold level to below 0%.")
        threshold_level = "0"
    print("Setting the river threshold level to " + threshold_level + "%.")
    update_idle()


def update_num_reaches(*args, **kwargs):
    """
    Updates the number of river reaches to show

    :param args: Unused parameter to allow function to work as a callback
    :param kwargs: Unused parameter to allow function to work as a callback
    """
    global num_reaches_str
    global num_reaches
    num_reaches_str = b_num_reaches.text
    try:
        num_reaches = int(num_reaches_str)
    except ValueError:
        num_reaches = 0
    print("Showing", num_reaches, "river reaches")
    update_idle()


def update_reach_dist(*args, **kwargs):
    """
    Updates the distance between river reaches

    :param args: Unused parameter to allow function to work as a callback
    :param kwargs: Unused parameter to allow function to work as a callback
    """
    global reach_dist_str
    global reach_dist
    reach_dist_str = b_reach_dist.text
    try:
        reach_dist = float(reach_dist_str)
    except ValueError:
        reach_dist = 5
    if reach_dist == 0:
        print("Error: Cannot have a distance between reaches of 0.\n\t"
              "Using default value of 5", reach_dist_units, ".")
        reach_dist = 5
    print("Maintaining a distance of",
          reach_dist, reach_dist_units,
          "between river reaches")
    update_idle()


def show_propagation(*args, **kwargs):
    """
    Shows the propagation over time.

    :param args: Unused parameter to allow function to work as a callback
    :param kwargs: Unused parameter to allow function to work as a callback
    """
    update_idle()
    fig_name = "River Propagation Time"
    if fig_name in plt.get_figlabels():
        return
    rivers_list = downstream_rivers_list
    if not downstream_rivers_list:
        return
    print("Drawing the event duration:")
    time.sleep(sleep_time)
    bar = progressbar.ProgressBar(
        maxval=num_reaches,
        widgets=progressbarWidgets).start()
    i = 0
    fig_prop = plt.figure(label=fig_name)
    fig_prop.canvas.manager.set_window_title(fig_name)
    plt.title(fig_name)
    plt.xlabel("Time (3 hours)")
    plt.ylabel("Distance (" + reach_dist_units + ")")
    net_distance = 0
    shown = False
    k = 0
    x_graphs = []
    y_graphs = []
    c_graphs = []
    last_index = 0
    last_x = 0
    for r in rivers_list:
        if r == 0:
            continue
        print("Graphing River", str(r))
        if str(r) in rivers_list:
            if not shown:
                j = 0
                try:
                    for s in discharge_graph_rivers[str(r)]:
                        index = list(Qout_data_ids).index(int(r))
                        color = river_colors[str(r)]
                        if k == 0:
                            plt.scatter(0, 0, s=1, alpha=1, c=color)
                            x_graphs.append(last_x)
                            y_graphs.append(net_distance)
                            c_graphs.append(color)
                            last_index = index
                            k += 1
                            break
                        x = Qout_data[:, last_index]
                        y = Qout_data[:, index]
                        corr = sm.tsa.stattools.ccf(x, y,
                                                    adjusted=False, fft=True)
                        last_x += list(corr).index(min(corr, key=abs))
                        plt.scatter(last_x, net_distance,
                                    s=1, alpha=1, c=color)
                        plt.plot([x_graphs[-1], last_x],
                                 [y_graphs[-1], net_distance],
                                 linewidth=0.5, alpha=1, c=color)
                        x_graphs.append(last_x)
                        y_graphs.append(net_distance)
                        c_graphs.append(color)
                        last_index = index
                        k += 1
                        break
                        j += 1
                    shown = True
                except:
                    pass
            if str(r) in river_lengths:
                net_distance += river_lengths[str(r)]
                if net_distance >= k * get_reach_dist_deg():
                    shown = False
        if k >= num_reaches:
            break
        i += 1
        bar.update(k)
        sys.stdout.flush()
    time.sleep(sleep_time)
    bar.finish()
    time.sleep(sleep_time)
    print("Finished drawing the event duration")
    ax_temp = plt.gca()
    ax_temp.set_ylim(ax_temp.get_ylim()[::-1])
    plt.show()


def show_event_duration(*args, **kwargs):
    """
    Shows the propagation over time

    :param args: Unused parameter to allow function to work as a callback
    :param kwargs: Unused parameter to allow function to work as a callback
    """
    update_idle()
    fig_name = "Peak Duration"
    if fig_name in plt.get_figlabels():
        return
    rivers_list = downstream_rivers_list
    if not downstream_rivers_list:
        return
    i = 0
    print("Drawing the event duration:")
    time.sleep(sleep_time)
    bar = progressbar.ProgressBar(
        maxval=len(rivers_list),
        widgets=progressbarWidgets).start()
    fig_temp, axs = plt.subplots(nrows=num_reaches,
                                 sharex="all", sharey="all",
                                 label=fig_name)
    fig_temp.suptitle(fig_name)
    fig_temp.supxlabel("Time (3 hours)")
    fig_temp.supylabel(r"Magnitude (m^3/s)")
    net_distance = 0
    shown = False
    threshold_min = ax_graph.get_ylim()[-1] * float(threshold_level) / 100
    k = 0
    for r in rivers_list:
        if r == 0:
            continue
        if str(r) in rivers_list:
            if not shown:
                above_threshold_x = list()
                above_threshold_y = list()
                is_last_above_threshold = False
                j = 0
                if str(r) in discharge_graph_rivers:
                    for s in discharge_graph_rivers[str(r)]:
                        r_temp = discharge_graph_rivers[str(r)][j]
                        x_temp, y_temp = r_temp.get_data()
                        for y in y_temp:
                            if y >= threshold_min:
                                x = x_temp[list(y_temp).index(y)]
                                if is_last_above_threshold:
                                    above_threshold_x[-1].append(x)
                                    above_threshold_y[-1].append(y)
                                else:
                                    above_threshold_x.append([x])
                                    above_threshold_y.append([y])
                                is_last_above_threshold = True
                            else:
                                is_last_above_threshold = False
                        if len(list(above_threshold_y)) > 0:
                            for xi, yi in zip(above_threshold_x,
                                              above_threshold_y):
                                axs[k].fill_between(
                                    xi, [max(yi) for i in yi],
                                    [threshold_min for i in yi],
                                    color=peak_duration_color)
                        j += 1
                    shown = True
            if str(r) in river_lengths:
                net_distance += river_lengths[str(r)]
            if net_distance >= k * get_reach_dist_deg():
                k += 1
                shown = False
        if k >= num_reaches:
            break
        i += 1
        bar.update(i)
        sys.stdout.flush()
    time.sleep(sleep_time)
    bar.finish()
    time.sleep(sleep_time)
    print("Finished drawing the event duration")
    plt.show()


def get_reach_dist_deg():
    """
    Gets the current river reach distance in degrees

    :return: The current river reach distance in degrees
    :rtype: float
    """
    if reach_dist_units == "Kilometers":
        return DistFncs.km2deg(reach_dist)
    elif reach_dist_units == "Nautical Miles":
        return DistFncs.M2deg(reach_dist)
    elif reach_dist_units == "Degrees":
        return reach_dist
    elif reach_dist_units == "Radians":
        return np.rad2deg(reach_dist)
    elif reach_dist_units == "Miles":
        return DistFncs.mi2deg(reach_dist)
    return DistFncs.km2deg(5)


def show_discharge_over_time(*args, **kwargs):
    """
    Shows the discharge over time based on the number of reaches and the
    distance between reaches

    :param args: Unused parameter to allow function to work as a callback
    :param kwargs: Unused parameter to allow function to work as a callback
    """
    discharge_title = "River Discharges Over Time"
    if discharge_title in plt.get_figlabels():
        return
    rivers_list = downstream_rivers_list
    if not downstream_rivers_list:
        return
    print("Drawing the river discharges:")
    time.sleep(sleep_time)
    bar = progressbar.ProgressBar(
        maxval=num_reaches,
        widgets=progressbarWidgets).start()
    i = 0
    fig_temp, axs = plt.subplots(nrows=num_reaches,
                                 sharex="all", sharey="all")
    fig_temp.suptitle("River Discharges Over Time")
    fig_temp.supxlabel("Time (3 hours)")
    fig_temp.supylabel(r"Average River Discharge (m^3/s)")
    fig_temp.canvas.manager.set_window_title(discharge_title)
    net_distance = 0
    shown = False
    k = 0
    for r in rivers_list:
        if r == 0:
            continue
        if str(r) in rivers_list:
            if not shown:
                above_threshold_x = list()
                above_threshold_y = list()
                is_last_above_threshold = False
                j = 0
                if str(r) in discharge_graph_rivers:
                    for s in discharge_graph_rivers[str(r)]:
                        r_temp = discharge_graph_rivers[str(r)][j]
                        x_temp, y_temp = r_temp.get_data()
                        axs[k].plot(*r_temp.get_data(),
                                    c=r_temp.get_color(),
                                    linewidth=r_temp.get_linewidth(),
                                    alpha=1,
                                    label=str(r))
                        axs[k].legend(loc="upper right")
                        threshold_min = axs[k].get_ylim()[-1] * float(
                            threshold_level) / 100
                        for y in y_temp:
                            if y >= threshold_min:
                                x = x_temp[list(y_temp).index(y)]
                                if is_last_above_threshold:
                                    above_threshold_x[-1].append(x)
                                    above_threshold_y[-1].append(y)
                                else:
                                    above_threshold_x.append([x])
                                    above_threshold_y.append([y])
                                is_last_above_threshold = True
                            else:
                                is_last_above_threshold = False
                        if len(list(above_threshold_y)) > 0:
                            for xi, yi in zip(above_threshold_x,
                                              above_threshold_y):
                                axs[k].fill_between(
                                    xi, yi,
                                    [threshold_min for i in yi],
                                    color=default_point_color)
                        j += 1
                    shown = True
            if str(r) in river_lengths:
                net_distance += river_lengths[str(r)]
            if net_distance >= k * get_reach_dist_deg():
                k += 1
                shown = False
        if k >= num_reaches:
            break
        i += 1
        bar.update(k)
        sys.stdout.flush()
    time.sleep(sleep_time)
    bar.finish()
    time.sleep(sleep_time)
    print("Finished drawing the river discharges")
    plt.show()


def get_total_distance(*args, **kwargs):
    """
    Gets the total distance of the downstream river reaches

    :param args: Unused parameter to allow function to work as a callback
    :param kwargs: Unused parameter to allow function to work as a callback
    :return: The total distance of the downstream river reaches
    :rtype: float
    """
    if not downstream_rivers_list:
        return
    net_distance = 0
    for r in downstream_rivers_list:
        if r == 0:
            continue
        if str(r) in downstream_rivers_list:
            if str(r) in river_lengths:
                net_distance += river_lengths[str(r)]
    return net_distance


def redraw_canvases(*args, **kwargs):
    """
    Redraws all canvases

    :param args: Unused parameter to allow function to work as a callback
    :param kwargs: Unused parameter to allow function to work as a callback
    """
    if not skip_rerender:
        fig.canvas.draw_idle()
        fig_config.canvas.draw_idle()


def update_idle(*args, **kwargs):
    """
    Updates idle elements for canvases

    :param args: Unused parameter to allow function to work as a callback
    :param kwargs: Unused parameter to allow function to work as a callback
    """
    fig_config.canvas.draw_idle()


def fix_config_size(*args, **kwargs):
    """
    Prevents the config canvas from being resized

    :param args: Unused parameter to allow function to work as a callback
    :param kwargs: Unused parameter to allow function to work as a callback
    """
    fig_config.set_size_inches(5.5, 2.5)


def clear_plots(*args, **kwargs):
    """
    Clears plots dynamically linked to the main interactive map

    :param args: Unused parameter to allow function to work as a callback
    :param kwargs: Unused parameter to allow function to work as a callback
    """
    global shown_rivers_list
    for path in p:
        if path in p:
            try:
                p[path].set_alpha(hidden_alpha)
            except:
                pass
    for path in rivers:
        if path in rivers:
            j = 0
            while True:
                try:
                    rivers[path][j].set_alpha(hidden_alpha)
                except IndexError:
                    break
                except:
                    continue
                j += 1
    for path in discharge_graph_rivers:
        if path in discharge_graph_rivers:
            j = 0
            while True:
                try:
                    discharge_graph_rivers[path][j].set_alpha(0)
                except IndexError:
                    break
                except:
                    continue
                j += 1
    shown_rivers_list = list()


def reset(*args, **kwargs):
    """
    Resets the view.

    :param args: Unused parameter to allow function to work as a callback
    :param kwargs: Unused parameter to allow function to work as a callback
    """
    global truncated_downstream
    global selected_rivers_list
    global currently_selected_river
    global downstream_rivers_list
    update_idle()
    print("Resetting view")
    start_time = datetime.now()
    for path in p:
        p[path].set_alpha(1)
    for path in rivers:
        j = 0
        for s in rivers[path]:
            rivers[path][j].set_alpha(1)
            j += 1
    for path in discharge_graph_rivers:
        j = 0
        while True:
            try:
                discharge_graph_rivers[path][j].set_alpha(1)
            except IndexError:
                break
            except:
                continue
            j += 1
    ax.set_xlim(default_xlims)
    ax.set_ylim(default_ylims)
    truncated_downstream = False
    selected_rivers_list = list()
    currently_selected_river = ""
    downstream_rivers_list = list()
    for i in plt.get_fignums():
        if i != 1:
            plt.close(i)
    redraw_canvases()
    print("Finished resetting view")
    print("Took", (datetime.now() - start_time).total_seconds(), "to reset")


def close_all(*args, **kwargs):
    """
    Closes all figures and ends the program

    :param args: Unused parameter to allow function to work as a callback
    :param kwargs: Unused parameter to allow function to work as a callback
    """
    print("Quitting SWOF.\nSee you next time!")
    plt.close("all")


def create_config_window():
    global fig_config
    global ax_downstream_dist
    global ax_num_downstream
    global b_num_reaches
    global ax_reach_dist
    global b_reach_dist
    global ax_threshold_level
    global b_threshold_level
    global ax_reach_dist
    global b_show_discharge
    global ax_propagation
    global b_propagation
    global ax_event_dur
    global b_event_dur
    ##### Config Window #####
    fig_config = plt.figure(figsize=(5.5, 2.5),
                            label="Control Room")
    fig_config.canvas.mpl_connect("pick_event", update_idle)
    fig_config.canvas.mpl_connect("button_press_event", update_idle)
    fig_config.canvas.mpl_connect("button_release_event", update_idle)
    fig_config.canvas.mpl_connect("key_press_event", update_idle)
    fig_config.canvas.mpl_connect("key_release_event", update_idle)
    fig_config.canvas.mpl_connect("motion_notify_event", update_idle)
    fig_config.canvas.mpl_connect("resize_event", fix_config_size)

    ### Total Distance of Downstream Selected River Reach ###
    ax_downstream_dist = plt.axes([0.025, 1 - 0.085, 0.9, 0.075])
    plt.axis("off")
    plt.text(
        0, 0,
        r"Total distance downstream of selected river reach (km):"
        + r"                          "
        + str(round(DistFncs.deg2km(get_total_distance()), 2))
    )

    ### Number of River Reaches ###
    ax_num_downstream = plt.axes([0.875, 1 - 0.1 - 0.095, 0.09, 0.075])
    b_num_reaches = matplotlib.widgets.TextBox(
        ax_num_downstream,
        label="Enter number of downstream river reaches to be analyzed:",
        initial=num_reaches_str,
        textalignment="center",
        label_pad=1.64,
        color=widget_bg,
        hovercolor=widget_bg_active,
    )
    b_num_reaches.on_text_change(update_idle)
    b_num_reaches.on_submit(update_num_reaches)

    ### Distance Between Reaches ###
    ax_reach_dist = plt.axes([0.875, 1 - 0.1 - 0.095 * 2, 0.09, 0.075])
    b_reach_dist = matplotlib.widgets.TextBox(
        ax_reach_dist,
        label="Enter distance separating river reaches to be analyzed (km):",
        initial=reach_dist_str,
        textalignment="center",
        label_pad=1.43,
        color=widget_bg,
        hovercolor=widget_bg_active,
    )
    b_reach_dist.on_text_change(update_idle)
    b_reach_dist.on_submit(update_reach_dist)

    ### Threshold Level ###
    ax_threshold_level = plt.axes([0.875, 1 - 0.1 - 0.095 * 3, 0.09, 0.075])
    b_threshold_level = matplotlib.widgets.TextBox(
        ax_threshold_level,
        label=r"Enter percentile for event detection in each river reach (%):",
        initial=threshold_level,
        textalignment="center",
        label_pad=1.59,
        color=widget_bg,
        hovercolor=widget_bg_active,
    )
    b_threshold_level.on_text_change(update_idle())
    b_threshold_level.on_submit(set_threshold_level)

    ### Distance Based Discharge Over Time ###
    ax_reach_dist = plt.axes([0.025, 0.025, 0.31, 0.4])
    b_show_discharge = plt.Button(
        ax_reach_dist,
        label="Show\nDischarge Timeseries",
        color=widget_bg,
        hovercolor=widget_bg_active,
    )
    b_show_discharge.on_clicked(show_discharge_over_time)

    ### Propagation ###
    ax_propagation = plt.axes([0.025 + 0.32, 0.025, 0.31, 0.4])
    b_propagation = plt.Button(
        ax_propagation,
        label="Show\nFlow Wave Propagation",
        color=widget_bg,
        hovercolor=widget_bg_active,
    )
    b_propagation.on_clicked(show_propagation)

    ### Event Duration ###
    ax_event_dur = plt.axes([0.025 + 0.32 * 2, 0.025, 0.31, 0.4])
    b_event_dur = plt.Button(
        ax_event_dur,
        label="Show\nEvent Duration",
        color=widget_bg,
        hovercolor=widget_bg_active,
    )
    b_event_dur.on_clicked(show_event_duration)
    fig_config.canvas.manager.set_window_title("Control Room")
    plt.show(block=False)


##### Config Window #####
fig_config = plt.figure(figsize=(5.5, 2.5))
fig_config.canvas.mpl_connect("pick_event", update_idle)
fig_config.canvas.mpl_connect("button_press_event", update_idle)
fig_config.canvas.mpl_connect("button_release_event", update_idle)
fig_config.canvas.mpl_connect("key_press_event", update_idle)
fig_config.canvas.mpl_connect("key_release_event", update_idle)
fig_config.canvas.mpl_connect("motion_notify_event", update_idle)
fig_config.canvas.mpl_connect("resize_event", fix_config_size)

### Total Distance of Downstream Selected River Reach ###
ax_downstream_dist = plt.axes([0.025, 1 - 0.085, 0.9, 0.075])
plt.axis("off")
plt.text(
    0, 0,
    r"Total distance downstream of selected river reach (km):"
    + r"                             "
    + str(300)
)

### Number of River Reaches ###
ax_num_downstream = plt.axes([0.875, 1 - 0.1 - 0.095, 0.09, 0.075])
b_num_reaches = matplotlib.widgets.TextBox(
    ax_num_downstream,
    label="Enter number of downstream river reaches to be analyzed:",
    initial=num_reaches_str,
    textalignment="center",
    label_pad=1.64,
    color=widget_bg,
    hovercolor=widget_bg_active,
)
b_num_reaches.on_text_change(update_idle)
b_num_reaches.on_submit(update_num_reaches)

### Distance Between Reaches ###
ax_reach_dist = plt.axes([0.875, 1 - 0.1 - 0.095 * 2, 0.09, 0.075])
b_reach_dist = matplotlib.widgets.TextBox(
    ax_reach_dist,
    label="Enter distance separating river reaches to be analyzed (km):",
    initial=reach_dist_str,
    textalignment="center",
    label_pad=1.43,
    color=widget_bg,
    hovercolor=widget_bg_active,
)
b_reach_dist.on_text_change(update_idle)
b_reach_dist.on_submit(update_reach_dist)

### Threshold Level ###
ax_threshold_level = plt.axes([0.875, 1 - 0.1 - 0.095 * 3, 0.09, 0.075])
b_threshold_level = matplotlib.widgets.TextBox(
    ax_threshold_level,
    label=r"Enter percentile for event detection in each river reach (%):",
    initial=threshold_level,
    textalignment="center",
    label_pad=1.59,
    color=widget_bg,
    hovercolor=widget_bg_active,
)
b_threshold_level.on_text_change(update_idle())
b_threshold_level.on_submit(set_threshold_level)

### Distance Based Discharge Over Time ###
ax_reach_dist = plt.axes([0.025, 0.025, 0.31, 0.4])
b_show_discharge = plt.Button(
    ax_reach_dist,
    label="Show\nDischarge Timeseries",
    color=widget_bg,
    hovercolor=widget_bg_active,
)
b_show_discharge.on_clicked(show_discharge_over_time)

### Propagation ###
ax_propagation = plt.axes([0.025 + 0.32, 0.025, 0.31, 0.4])
b_propagation = plt.Button(
    ax_propagation,
    label="Show\nFlow Wave Propagation",
    color=widget_bg,
    hovercolor=widget_bg_active,
)
b_propagation.on_clicked(show_propagation)

### Event Duration ###
ax_event_dur = plt.axes([0.025 + 0.32 * 2, 0.025, 0.31, 0.4])
b_event_dur = plt.Button(
    ax_event_dur,
    label="Show\nEvent Duration",
    color=widget_bg,
    hovercolor=widget_bg_active,
)
b_event_dur.on_clicked(show_event_duration)


##### River Network Window Final Processing #####
### Allowing for Resetting ###
b_reset.on_clicked(reset)

### Allowing for Quiting the Program ###
b_quit.on_clicked(close_all)

### Allowing for River Reaches to be Clicked on ###
fig.canvas.mpl_connect("pick_event", on_pick)


##### Getting default map values #####
default_xlims = ax.get_xlim()
default_ylims = ax.get_ylim()


##### Finished Preprocessing #####
print("Finished preprocessing")
fig.canvas.draw()
fig_config.canvas.draw()
fig.canvas.manager.set_window_title("River Network")
fig_config.canvas.manager.set_window_title("Control Room")
plt.close(fig_graph)
plt.close(fig_config)
plt.show(block=True)
