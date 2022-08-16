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
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
import matplotlib.style as mplstyle
import matplotlib.widgets
import numpy as np
import progressbar
# import shapefile as shp
from descartes import PolygonPatch
from shapely.geometry import MultiLineString, LineString
import fiona
import statsmodels.api as sm
from scipy.io import netcdf_file

from rrr_msc_fnc import DistFncs

from svgpathtools import svg2paths
from svgpath2mpl import parse_path


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

grid_alpha = 0.1
widget_bg = "#F0852D"
widget_bg_active = "#F8B723"
quit_bg = "#CC0202"
quit_bg_active = "#FF0000"
b_c = "#0492C2"
b_c_a = "#48AAAD"
g_c = "#018A0F"
g_c_a = "#36A300"
peak_duration_color = "#82EEFD"
default_point_color = "#FFFFFF"
swapping_colors = [
    "#CF0DFF",
    "#4E0CE8",
    "#0033FF",
    "#0CA2E8",
    "#00FFD0",
]
cycol = itertools.cycle(swapping_colors)

waypoint_path, attributes = svg2paths("./assets/svgs/waypoint.svg")
waypoint = parse_path(attributes[0]["d"])
waypoint.vertices -= waypoint.vertices.mean(axis=0)
waypoint = waypoint.transformed(mpl.transforms.Affine2D().rotate_deg(180))
waypoint = waypoint.transformed(mpl.transforms.Affine2D().scale(-1, 1))

# waypoint_path, attributes = svg2paths("./assets/svgs/waypoint.svg")
# waypoints = [
#     parse_path(attributes[i]["d"]) for i in range(len(attributes) - 1)
# ]
# for i in range(len(waypoints)):
#     waypoints[i].vertices -= waypoints[i].vertices.mean(axis=0)
#     waypoints[i] = waypoints[i].transformed(mpl.transforms.Affine2D().rotate_deg(180))
#     waypoints[i] = waypoints[i].transformed(mpl.transforms.Affine2D().scale(-1, 1))

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
connectivity_f_path = "../San_Gaud_data/rapid_connect_San_Guad.csv"
sf_path = "../San_Gaud_data/NHDFlowline_San_Guad/NHDFlowline_San_Guad.shp"
Qout_f_path = "../San_Gaud_data/Qout_San_Guad_exp00.nc"

show_errors = True
mpl_backend_bases = matplotlib.backend_bases
fig, ax = plt.subplots(label="River Network")
ax.grid(alpha=0.1)
coords = list()
coords_dict = dict()
id_ind = dict()
ind_id = dict()
connectivity = dict()
p = dict()
rivers = dict()
x_vals = dict()
x_vals_list = list()
y_vals = dict()
y_vals_list = list()
discharge_graph_rivers = dict()
river_colors = dict()
num_downstream = 5
downstream_rivers_list = list()
offsets = list()
open_in_new_window = False
skip_rerender = False
max_river_selections = 10
enlarged_rivers = list()
hidden_alpha = 0.1
point_scaling = 2
num_reaches_str = "3"
num_reaches = int(num_reaches_str)
reach_dist_str = "117"
reach_dist = float(reach_dist_str)
default_reach_dist = reach_dist
reach_dist_units = "km"
have_default_view = False
have_cleared_view = False
shown_rivers_list = list()
default_xlims = list()
default_ylims = list()
threshold_level = "90"
sleep_time = 0.01
showing_river_discharges = False
default_linewidth = 1


##### Declaration of variables (given as command line arguments) #####
# (1) - rrr_connect_file
# (2) - rrr_shp_file
# (3) - rrr_Qout_file


##### Get command line arguments #####
IS_arg = len(sys.argv)
if IS_arg > 4:
    print("Error: A maximum of 3 arguments can be used.")
    raise SystemExit(22)

if IS_arg > 1:
    connectivity_f_path = sys.argv[2]
    if IS_arg > 2:
        sf_path = sys.argv[3]
        if IS_arg > 3:
            Qout_f_path = sys.argv[4]

Qout_f = netcdf_file(Qout_f_path, "r")

Qout_data = Qout_f.variables["Qout"][:]*1
"""average river water discharge downstream of each river reach"""

Qout_data_normed = (Qout_data.T / abs(Qout_data).max(axis=1)).T
"""normalized average river water discharge downstream of each river reach"""

y_vals_array = Qout_f.variables["lat"][:]*1
"""latitude of river data points"""

x_vals_array = Qout_f.variables["lon"][:]*1
"""longitude of river data points"""

Qout_data_ids = Qout_f.variables["rivid"][:]*1
"""unique identifier for each river reach"""

Qout_data_max = max(abs(Qout_data[0]))


plt.scatter(
    x_vals_array,
    y_vals_array,
    s=point_scaling,
    picker=5,
    color=mpl.colors.to_rgba(default_point_color, 1)
)


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
        data = [int(d) for d in row[1:]]
        connectivity[int(row[0].replace(" ", ""))] = data
        j += 1
        bar.update(j)
        sys.stdout.flush()
time.sleep(sleep_time)
bar.finish()
time.sleep(sleep_time)
print("Finished loading connectivity information")

with fiona.open(sf_path) as f:
    river_lengths = [i["properties"]["LENGTHKM"] for i in f]
    shp_comids = [i["properties"]["COMID"] for i in f]
    mp = MultiLineString([LineString(i["geometry"]["coordinates"]) for i in f])
    coord_arr = [[[x for x, y, z in i["geometry"]["coordinates"]],
                  [y for x, y, z in i["geometry"]["coordinates"]]] for i in f]
    [plt.plot(xi, yi,
              color=b_c,
              linewidth=default_linewidth) for xi, yi in coord_arr]


##### River Network Window #####
plt.title("Please click on one river reach")
plt.xlabel("Longitude")
plt.ylabel("Latitude")

### Quit Option ###
ax_quit = plt.axes([0.875, 1 - 0.05 - 0.075, 0.1, 0.075])
b_quit = plt.Button(
    ax_quit,
    label="Quit",
    color=quit_bg,
    hovercolor=quit_bg_active,
)
b_quit.label.set_fontsize(14)

### Reset Option ###
ax_reset = plt.axes([0.875, 1 - 0.05 - 0.075 * 2 - 0.01, 0.1, 0.075])
b_reset = plt.Button(
    ax_reset,
    label="Reset",
    color=widget_bg,
    hovercolor=widget_bg_active,
)
b_reset.label.set_fontsize(14)

### Option to Save All ###
ax_save = plt.axes([0.875, 1 - 0.05 - 0.075 * 3 - 0.02, 0.1, 0.075])
b_save = plt.Button(
    ax_save,
    label="Save\nFigures",
    color=g_c,
    hovercolor=g_c_a,
)


def get_rivids_along_downstream(num_rivids, dist):
    """
    Gets the river ids downstream at dist distance intervals

    :param num_rivids: The maximum number of river ids to return
    :type num_rivids: int
    :param dist: The distance between river ids
    :type dist: float
    :return: The river ids downstream and the distances thereof
    :rtype: tuple[list[int], list[float]]
    """
    global downstream_rivers_list
    rivids = list()
    dists = list()
    net_distance = 0
    for rivid in downstream_rivers_list:
        if len(rivids) >= num_rivids:
            break
        net_distance += river_lengths[shp_comids.index(rivid)]
        if net_distance > len(rivids) * dist:
            rivids.append(rivid)
            dists.append(net_distance)
    return rivids, dists


def choose_downstream(rivid):
    """
    Chooses the new downstream as being rivid.

    :param rivid: The id of the head of the river downstream
    :type rivid: int or float or str
    :return: The list of downstream rivers
    :rtype: list[int]
    """
    global downstream_rivers_list
    rivid = int(rivid)
    downstream_rivers_list.append(rivid)
    while rivid in connectivity:
        rivid = connectivity[rivid][0]
        downstream_rivers_list.append(rivid)
    return downstream_rivers_list


def on_pick(event):
    """
    Handles selection of points.

    - Callback function

    :param event: Automatically generated by the event manager.
    :type event: mpl_backend_bases.PickEvent
    """
    global offsets
    global downstream_rivers_list
    offsets_temp = event.artist.get_offsets()[event.ind][0]
    if np.any(offsets):
        if np.all(np.isclose(offsets_temp, offsets, rtol=1e-3)):
            return
    offsets = offsets_temp
    i = 0
    for j in range(len(x_vals_array)):
        if x_vals_array[i] == offsets[0] and y_vals_array[i] == offsets[1]:
            break
        i += 1
    id = Qout_data_ids[i]
    start_time = datetime.now()
    print("Showing downstream for", id)
    choose_downstream(id)
    rivids, dists = get_rivids_along_downstream(num_reaches, reach_dist)
    ax.cla()
    ax.set_title("Please click on one river reach")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    c_p = [mpl.colors.to_rgba(default_point_color, 1)
           if i in downstream_rivers_list
           else mpl.colors.to_rgba(default_point_color, hidden_alpha)
           for i in Qout_data_ids]

    # m_p = [waypoint
    #        if i in rivids
    #        else "o"
    #        for i in Qout_data_ids]
    ax.scatter(
        x_vals_array,
        y_vals_array,
        s=point_scaling,
        picker=5,
        color=c_p,
    )
    c_l = [mpl.colors.to_rgba(b_c, 1)
           if i in downstream_rivers_list
           else mpl.colors.to_rgba(b_c, hidden_alpha)
           for i in shp_comids]
    [ax.plot(xy[0], xy[1], color=c, linewidth=default_linewidth)
     for xy, c in zip(coord_arr, c_l)]

    x_p = [xi for xi, idx in zip(x_vals_array, Qout_data_ids) if idx in rivids]
    y_p = [yi for yi, idx in zip(y_vals_array, Qout_data_ids) if idx in rivids]
    ax.scatter(
        x_p,
        y_p,
        s=point_scaling * 1000,
        picker=5,
        color=mpl.colors.to_rgba("#FF0000", 1),
        marker=waypoint,
    )
    # [[ax.scatter(
    #     xi,
    #     yi,
    #     s=point_scaling * 1000,
    #     picker=5,
    #     color=mpl.colors.to_rgba("#FF0000", 1),
    #     marker=w,
    # )
    #  if idx in rivids else 0
    #  for xi, yi, idx in zip(x_vals_array, y_vals_array, Qout_data_ids)]
    #     for w in waypoints
    #  ]
    print("Took",
          (datetime.now() - start_time).total_seconds(),
          "seconds to select")
    create_config_window()
    redraw_canvases()
    print("Took",
          (datetime.now() - start_time).total_seconds(),
          "seconds to select and draw")


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
        reach_dist = default_reach_dist
    if reach_dist == 0:
        print("Error: Cannot have a distance between reaches of 0.\n\t"
              "Using default value of",
              default_reach_dist, str(reach_dist_units) + ".")
        reach_dist = default_reach_dist
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
    if not downstream_rivers_list:
        return
    rivids, dists = get_rivids_along_downstream(num_reaches, reach_dist)
    print("Drawing the propagation timeseries:")
    time.sleep(sleep_time)
    bar = progressbar.ProgressBar(
        maxval=len(rivids),
        widgets=progressbarWidgets).start()
    i = 0
    fig_prop = plt.figure(label=fig_name)
    fig_prop.canvas.manager.set_window_title(fig_name)
    plt.title(fig_name)
    plt.xlabel("Time (3 hours)")
    plt.ylabel("Distance (" + reach_dist_units + ")")

    first_s = None
    last_dist = None
    prev_x = None
    for s, i, idx, d in zip(
            [list(Qout_data_ids).index(rivid) for rivid in rivids],
            list(range(len(rivids))),
            rivids,
            dists):
        if i == 0:
            first_s = s
            last_dist = d
            prev_x = 0
            continue
        x = Qout_data[:, first_s]
        y = Qout_data[:, s]
        corr = sm.tsa.stattools.ccf(x, y,
                                    adjusted=True, fft=True)
        last_x = list(corr).index(max(corr))
        plt.grid(alpha=grid_alpha)
        plt.scatter(
            prev_x + last_x, last_dist + d,
            s=point_scaling,
            alpha=1,
            # c=b_c,
        )
        plt.plot(
            [prev_x, prev_x + last_x],
            [last_dist, last_dist + d],
            linewidth=default_linewidth,
            alpha=1,
            # color=b_c,
            label=str(idx)
        )
        prev_x += last_x
        last_dist += d
        bar.update(i)
        sys.stdout.flush()
    time.sleep(sleep_time)
    bar.finish()
    time.sleep(sleep_time)
    print("Finished drawing the propagation timeseries")
    ax_temp = plt.gca()
    ax_temp.legend(loc="upper right")
    ax_temp.set_ylim(ax_temp.get_ylim()[::-1])
    fig_prop.canvas.draw_idle()
    plt.show(block=False)


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
    rivids, dists = get_rivids_along_downstream(num_reaches, reach_dist)
    if not downstream_rivers_list:
        return
    i = 0
    print("Drawing the event duration:")
    time.sleep(sleep_time)
    bar = progressbar.ProgressBar(
        maxval=len(rivids),
        widgets=progressbarWidgets).start()
    fig_temp = plt.figure(label=fig_name)
    plt.title(fig_name)
    plt.xlabel("Time (3 hours)")
    plt.ylabel(r"Magnitude (m^3/s)")
    plt.grid(alpha=grid_alpha)
    fig_temp.canvas.manager.set_window_title(fig_name)
    color_swap = itertools.cycle(swapping_colors)
    for s, i, idx in zip(
            [list(Qout_data_ids).index(rivid) for rivid in rivids],
            list(range(len(rivids))),
            rivids):
        color = next(color_swap)
        y_percentile = np.percentile(Qout_data[:, s], float(threshold_level))
        x_p = [xi for xi, yi in zip(
            list(range(len(list(Qout_data[:, s])))),
            list(Qout_data[:, s])) if yi >= y_percentile]
        y_p = [yi for xi, yi in zip(
            list(range(len(list(Qout_data[:, s])))),
            list(Qout_data[:, s])) if yi >= y_percentile]
        [plt.fill_between(xi, [yi - 20], [yi],
                          alpha=1,
                          color=color,
                          label=str(idx)) if xi == x_p[0]
         else plt.fill_between(xi, [yi - 20], [yi],
                               alpha=1,
                               color=color)
         for xi, yi in zip(x_p, y_p)]
        bar.update(i)
        sys.stdout.flush()
    time.sleep(sleep_time)
    bar.finish()
    time.sleep(sleep_time)
    print("Finished drawing the river discharges")
    plt.legend(loc="upper right")
    fig_temp.canvas.draw_idle()
    plt.show(block=False)


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
    return DistFncs.km2deg(80)


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
    if not downstream_rivers_list:
        return
    rivids, dists = get_rivids_along_downstream(num_reaches, reach_dist)
    print("Drawing the river discharges:")
    time.sleep(sleep_time)
    bar = progressbar.ProgressBar(
        maxval=len(rivids),
        widgets=progressbarWidgets).start()
    fig_temp, axs = plt.subplots(nrows=len(rivids),
                                 sharex="all", sharey="all",
                                 label=discharge_title)
    fig_temp.suptitle("River Discharges Over Time")
    fig_temp.supxlabel("Time (3 hours)")
    fig_temp.supylabel(r"Average River Discharge (m^3/s)")
    fig_temp.canvas.manager.set_window_title(discharge_title)
    for s, i, idx in zip(
            [list(Qout_data_ids).index(rivid) for rivid in rivids],
            list(range(len(rivids))),
            rivids):
        axs[i].grid(alpha=grid_alpha)
        axs[i].plot(
            list(range(Qout_data.shape[0])),
            list(Qout_data[:, s]),
            linewidth=default_linewidth,
            alpha=1,
            color=b_c,
            label=str(idx)
        )
        y_percentile = np.percentile(Qout_data[:, s], float(threshold_level))
        axs[i].fill_between(x=list(range(len(list(Qout_data[:, s])))),
                            y1=[y_percentile] * len(list(Qout_data[:, s])),
                            y2=list(Qout_data[:, s]),
                            where=[yi >= y_percentile
                                   for yi in list(Qout_data[:, s])],
                            alpha=1,
                            color="#FFFFFF")
        axs[i].legend(loc="upper right")
        bar.update(i)
        sys.stdout.flush()
    time.sleep(sleep_time)
    bar.finish()
    time.sleep(sleep_time)
    print("Finished drawing the river discharges")
    fig_temp.canvas.draw_idle()
    plt.show(block=False)


def get_total_distance(*args, **kwargs):
    """
    Gets the total distance of the downstream river reaches

    :param args: Unused parameter to allow function to work as a callback
    :param kwargs: Unused parameter to allow function to work as a callback
    :return: The total distance of the downstream river reaches
    :rtype: float
    """
    if not downstream_rivers_list:
        return 0
    return sum([rl if rivid in downstream_rivers_list else 0
                for rl, rivid in zip(river_lengths, shp_comids)])


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


def close_all(*args, **kwargs):
    """
    Closes all figures and ends the program

    :param args: Unused parameter to allow function to work as a callback
    :param kwargs: Unused parameter to allow function to work as a callback
    """
    print("Quitting RAFT.\nSee you next time!")
    plt.close("all")


def reset(*args, **kwargs):
    """
    Resets the view.

    :param args: Unused parameter to allow function to work as a callback
    :param kwargs: Unused parameter to allow function to work as a callback
    """
    global downstream_rivers_list
    global offsets
    update_idle()
    print("Resetting view")
    start_time = datetime.now()
    ax.cla()
    ax.set_title("Please click on one river reach")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    [ax.plot(xi, yi, color=b_c, linewidth=default_linewidth)
     for xi, yi in coord_arr]
    ax.scatter(
        x_vals_array,
        y_vals_array,
        s=point_scaling,
        picker=5,
        color=mpl.colors.to_rgba(default_point_color, 1)
    )
    downstream_rivers_list = list()
    offsets = list()
    for i in plt.get_fignums():
        if i != 1:
            plt.close(i)
    fig.canvas.draw_idle()
    print("Finished resetting view")
    print("Took", (datetime.now() - start_time).total_seconds(),
          "seconds to reset")


def save_all(*args, **kwargs):
    """
    Saves all figures as svgs, pdf, and pngs.

    :param args: Unused parameter to allow function to work as a callback
    :param kwargs: Unused parameter to allow function to work as a callback
    """
    print("Saving all figures as svgs, pdf, and pngs")
    time.sleep(sleep_time)
    fig_labels = plt.get_figlabels()
    bar = progressbar.ProgressBar(
        maxval=len(fig_labels),
        widgets=progressbarWidgets).start()
    i = 0
    for l in fig_labels:
        plt.figure(l).savefig("./saved_outputs/svgs/" + str(l) + ".svg",
                              format="svg")
        plt.figure(l).savefig("./saved_outputs/pdfs/" + str(l) + ".pdf",
                              format="pdf")
        plt.figure(l).savefig("./saved_outputs/pngs/" + str(l) + ".png",
                              format="png")
        i += 1
        bar.update(i)
        sys.stdout.flush()
    time.sleep(sleep_time)
    bar.finish()
    time.sleep(sleep_time)
    print("Finished saving all figures as svgs, pdf, and pngs")


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
        x=0, y=0,
        s=r"Total distance downstream of selected river reach (km):"
          # + r"                          "
          # + str(round(get_total_distance(), 2))
    )
    ax_downstream_dist_num = plt.axes([0.875 + 0.09 / 2, 1 - 0.1, 0.09, 0.075])
    plt.axis("off")
    plt.text(
        x=0, y=0,
        s=str(round(get_total_distance(), 2)),
        ha="center",
        va="bottom",
        transform=ax_downstream_dist_num.transAxes,
    )


    ### Number of River Reaches ###
    ax_num_downstream = plt.axes([0.875, 1 - 0.1 - 0.095, 0.09, 0.075])
    b_num_reaches = matplotlib.widgets.TextBox(
        ax_num_downstream,
        label="Enter number of downstream river reaches to be analyzed:",
        initial=num_reaches_str,
        textalignment="center",
        label_pad=1.656,
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
        label_pad=1.45,
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
        label_pad=1.612,
        color=widget_bg,
        hovercolor=widget_bg_active,
    )
    b_threshold_level.on_text_change(update_idle())
    b_threshold_level.on_submit(set_threshold_level)

    ### Distance Based Discharge Over Time ###
    ax_show_discharge = plt.axes([0.025, 0.025, 0.31, 0.4])
    b_show_discharge = plt.Button(
        ax_show_discharge,
        label="Show\nDischarge Timeseries",
        color=b_c,
        hovercolor=b_c_a,
    )
    b_show_discharge.on_clicked(show_discharge_over_time)

    ### Propagation ###
    ax_propagation = plt.axes([0.025 + 0.32, 0.025, 0.31, 0.4])
    b_propagation = plt.Button(
        ax_propagation,
        label="Show\nFlow Wave Propagation",
        color=b_c,
        hovercolor=b_c_a,
    )
    b_propagation.on_clicked(show_propagation)

    ### Event Duration ###
    ax_event_dur = plt.axes([0.025 + 0.32 * 2, 0.025, 0.31, 0.4])
    b_event_dur = plt.Button(
        ax_event_dur,
        label="Show\nEvent Duration",
        color=b_c,
        hovercolor=b_c_a,
    )
    b_event_dur.on_clicked(show_event_duration)
    fig_config.canvas.manager.set_window_title("Control Room")
    plt.show(block=False)


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

fig.canvas.mpl_connect("motion_notify_event", redraw_canvases)

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
ax_show_discharge = plt.axes([0.025, 0.025, 0.31, 0.4])
b_show_discharge = plt.Button(
    ax_show_discharge,
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
### Allowing for Quiting the Program ###
b_quit.on_clicked(close_all)

### Allowing for Resetting ###
b_reset.on_clicked(reset)

### Allowing for Saving All Figures ###
b_save.on_clicked(save_all)

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
plt.close(fig_config)
plt.show(block=True)
