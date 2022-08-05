#******************************************************************************
#coord_viewer.py
#******************************************************************************

# Purpose:
# Given a river model output netCDF file, and a river network shapefile, this
# program creates an interactive model to better understand and communicate the
# propagation of water through space and time within rivers.
# Authors:
# Alex Christopher Lim, Cedric H. David, 2022-2022


##### Import Python modules #####
import csv
import itertools
import pickle
from datetime import datetime
from cycler import cycler

import matplotlib as mpl
import matplotlib.backend_bases
import matplotlib.colors
import matplotlib.pyplot as plt
import matplotlib.style as mplstyle
import matplotlib.widgets
import numpy as np
import shapefile as shp
from PIL import Image
from scipy.io import netcdf_file

from misc_functions import DistFuncs

mplstyle.use('fast')
mpl.rcParams['path.simplify'] = True


##### Miscellaneous Globals #####
backend_bases = matplotlib.backend_bases
fig, ax = plt.subplots()
coords_f_path = "../San_Gaud_data/coords_San_Guad.csv"
coords = list()
coords_dict = dict()
id_ind = dict()
ind_id = dict()
connectivity_f_path = "../San_Gaud_data/rapid_connect_San_Guad.csv"
sf_path = "../San_Gaud_data/NHDFlowline_San_Guad/NHDFlowline_San_Guad.shp"
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
Qout_f_path = "../San_Gaud_data/Qout_San_Guad_exp00.nc"
default_point_color = "b"
discharge_graph_rivers = dict()
river_colors = dict()
num_downstream = 5
currently_selected_river = ""
downstream_rivers_list = list()
truncated_downstream = False
select_multiple = False
selected_rivers_list = list()
open_in_new_window = False
max_river_selections = 10
enlarged_rivers = list()
hidden_alpha = 0.1
num_reaches_str = "2"
num_reaches = int(num_reaches_str)
reach_dist_str = "5"
reach_dist = float(reach_dist_str)
reach_dist_units = "km"
have_default_view = False
have_cleared_view = False
shown_rivers_list = list()

cycol = itertools.cycle("bgrcmk")


Qout_f = netcdf_file(Qout_f_path, "r")

Qout_data = Qout_f.variables["Qout"][:]*1
"""average river water discharge downstream of each river reach"""

Qout_data_normed = (Qout_data.T / abs(Qout_data).max(axis=1)).T
print("Max value:", abs(Qout_data).max())
"""normalized average river water discharge downstream of each river reach"""

Qout_data_ids = Qout_f.variables["rivid"][:]*1
"""unique identifier for each river reach"""

Qout_data_max = max(abs(Qout_data[0]))

i = 0
with open(coords_f_path, newline="\n") as f:
    for row in csv.reader(f, delimiter=","):
        print("Loading data point", i)
        draw_point = True
        x_vals[float(row[-2].replace(" ", ""))] = row[0]
        y_vals[float(row[-1].replace(" ", ""))] = row[0]
        x_vals_list.append(float(row[-2].replace(" ", "")))
        y_vals_list.append(float(row[-1].replace(" ", "")))
        coords.append([float(row[-2]), float(row[-1])])
        coords_dict[row[0]] = [float(row[-2]), float(row[-1])]
        idx = np.where(Qout_data_ids == int(row[0]))
        size = Qout_data_normed[0][idx]
        color = default_point_color
        if draw_point:
            p[row[0]] = plt.scatter(float(row[-2]), float(row[-1]),
                                    s=size, picker=5, c=color,
                                    alpha=hidden_alpha)
        id_ind[row[0]] = i
        ind_id[i] = row[0]
        i += 1

x_vals_array = np.array(x_vals_list)
y_vals_array = np.array(y_vals_list)

with open(connectivity_f_path, newline="\n") as f:
    j = 0
    for row in csv.reader(f, delimiter=","):
        print("Loading connectivity information for data point", j)
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


print("Reading shapefile", sf_path)
sf = shp.Reader(sf_path)
print("Finished reading shapefile", sf_path)
j = 0
for shape in sf.shapeRecords():
    print("Drawing shape", j)
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
        d = DistFuncs.great_circle_dist(k[0], xy[i][0], k[1], xy[i][1])
        if id in river_lengths:
            river_lengths[id] += d
        else:
            river_lengths[id] = d
        i += 1
    color = next(cycol)
    rivers[id] = plt.plot(x, y, linewidth=0.5, alpha=hidden_alpha, c=color)
    river_colors[str(id)] = color
    j += 1


plt.title("Raw River Coordinate Data")
plt.xlabel("Longitude")
plt.ylabel("Latitude")

Qout_graph_x_vals = list()
Qout_graph_y_vals = list()
fig_graph, ax_graph = plt.subplots()
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
plt.rcParams.update({
    "text.usetex": True,
})
plt.title("River Discharge Over Time")
plt.xlabel("Time (3 hours)")
plt.ylabel(r"Average River Discharge ($m^3/s$)")


def on_pick(event):
    """
    Handles selection of points.

    - Callback function

    :param event: Automatically generated by the event manager.
    :type event: matplotlib.backend_bases.PickEvent
    """
    global currently_selected_river
    global selected_rivers_list
    offsets = event.artist.get_offsets()[event.ind][0]
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
    redraw_canvases()


def apply_multiple(*args, **kwargs):
    """
    Shows the selected rivers
    
    :param args: Unused parameter to allow function to work as a callback
    :param kwargs: Unused parameter to allow function to work as a callback
    """
    print("Drawing graphs")
    i = 0
    for r in selected_rivers_list:
        fig_temp = plt.figure()
        plt.title("River " + str(r) + " Discharge Over Time")
        plt.xlabel("Time (3 hours)")
        plt.ylabel(r"Average River Discharge ($m^3/s$)")
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
    plt.show(block=False)


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
            plt.legend()
        if i > max_river_selections:
            break
        i += 1
    plt.show(block=False)


def isolate_evenly_spaced(*args, **kwargs):
    """
    Shows the evenly spaced downstream rivers in one window

    :param args: Unused parameter to allow function to work as a callback
    :param kwargs: Unused parameter to allow function to work as a callback
    """
    if open_in_new_window:
        isolate_evenly_spaced_new_windows(rivers_list=discharge_graph_rivers)
    else:
        isolate_evenly_spaced_one_window(rivers_list=discharge_graph_rivers)


def isolate_evenly_spaced_new_windows(rivers_list=None, *args, **kwargs):
    """
    Shows the evenly spaced downstream rivers in multiple windows

    :param rivers_list: The names of the downstream rivers, defaults to the
        average river water discharge downstream of each river reach
    :param args: Unused parameter to allow function to work as a callback
    :param kwargs: Unused parameter to allow function to work as a callback
    """
    if not rivers_list:
        rivers_list = discharge_graph_rivers
    print("Drawing graphs")
    i = 0
    m = int(len(rivers_list) / max_river_selections)
    k = 0
    for r in rivers_list:
        if i % m != 0:
            i += 1
            continue
        print("Graphing River", str(r))
        if str(r) in rivers_list:
            fig_temp = plt.figure()
            plt.title("River " + str(r) + " Discharges Over Time")
            plt.xlabel("Time (3 hours)")
            plt.ylabel(r"Average River Discharge ($m^3/s$)")
            j = 0
            for s in rivers_list[str(r)]:
                r_temp = rivers_list[str(r)][j]
                plt.plot(*r_temp.get_data(),
                         c=r_temp.get_color(),
                         linewidth=r_temp.get_linewidth(),
                         alpha=r_temp.get_alpha(),
                         label=str(r))
                plt.legend()
                j += 1
            k += 1
        if k >= max_river_selections:
            break
        i += 1
    plt.show(block=False)


def isolate_evenly_spaced_one_window(rivers_list=None, *args, **kwargs):
    """
    Shows the evenly spaced downstream rivers in one window

    :param rivers_list: The names of the downstream rivers, defaults to the
        average river water discharge downstream of each river reach
    :param args: Unused parameter to allow function to work as a callback
    :param kwargs: Unused parameter to allow function to work as a callback
    """
    if not rivers_list:
        rivers_list = discharge_graph_rivers
    print("Drawing graphs")
    i = 0
    fig_temp, axs = plt.subplots(nrows=max_river_selections,
                                 sharex="all", sharey="all")
    fig_temp.suptitle("River Discharges Over Time")
    fig_temp.supxlabel("Time (3 hours)")
    fig_temp.supylabel(r"Average River Discharge ($m^3/s$)")
    m = int(len(rivers_list) / max_river_selections)
    k = 0
    for r in rivers_list:
        if i % m != 0:
            i += 1
            continue
        print("Graphing River", str(r))
        if str(r) in rivers_list:
            j = 0
            for s in rivers_list[str(r)]:
                r_temp = rivers_list[str(r)][j]
                axs[k].plot(*r_temp.get_data(),
                            c=r_temp.get_color(),
                            linewidth=r_temp.get_linewidth(),
                            alpha=r_temp.get_alpha(),
                            label=str(r))
                axs[k].legend()
                j += 1
            k += 1
        if k >= max_river_selections:
            break
        i += 1
    plt.show(block=False)


def all_one_window_saved(rivers_list=None, *args, **kwargs):
    """
    Shows all downstream rivers in one window

    :param rivers_list: The names of the downstream rivers, defaults to the
        average river water discharge downstream of each river reach
    :param args: Unused parameter to allow function to work as a callback
    :param kwargs: Unused parameter to allow function to work as a callback
    """
    rivers_list = downstream_rivers_list
    if not downstream_rivers_list:
        return
    print("Drawing graphs")
    i = 0
    fig_temp, axs = plt.subplots(nrows=num_downstream,
                                 sharex="all", sharey="all")
    fig_temp.suptitle("River Discharges Over Time")
    fig_temp.supxlabel("Time (3 hours)")
    fig_temp.supylabel(r"Average River Discharge ($m^3/s$)")
    k = 0
    for r in rivers_list:
        print("Graphing River", str(r))
        if str(r) in rivers_list:
            j = 0
            for s in discharge_graph_rivers[str(r)]:
                r_temp = discharge_graph_rivers[str(r)][j]
                axs[k].plot(*r_temp.get_data(),
                            c=r_temp.get_color(),
                            linewidth=r_temp.get_linewidth(),
                            alpha=r_temp.get_alpha(),
                            label=str(r))
                axs[k].legend()
                j += 1
            k += 1
        if k >= num_downstream:
            break
        i += 1
    redraw_canvases()
    plt.show(block=False)


def checkboxes_callback(*args, **kwargs):
    """
    Handles changes to checkboxes

    :param args: Unused parameter to allow function to work as a callback
    :param kwargs: Unused parameter to allow function to work as a callback
    """
    global select_multiple
    global selected_rivers_list
    global open_in_new_window
    select_multiple_temp, open_in_new_window_temp = b_checkboxes.get_status()
    if select_multiple_temp != select_multiple:
        select_multiple = select_multiple_temp
        if select_multiple:
            print("Selecting Multiple")
        else:
            print("Selecting Entire Downstream")
            selected_rivers_list = list()
    elif open_in_new_window_temp != open_in_new_window:
        open_in_new_window = open_in_new_window_temp
        if open_in_new_window:
            print("Selecting to Open in New Windows")
        else:
            print("Selecting to Open in Same Window")


def update_num_downstream(*args, **kwargs):
    """
    Updates the number of downstream rivers to show

    :param args: Unused parameter to allow function to work as a callback
    :param kwargs: Unused parameter to allow function to work as a callback
    """
    global num_downstream
    global truncated_downstream
    if not truncated_downstream:
        num_downstream = len(downstream_rivers_list)
        truncated_downstream = True
    num_downstream_temp = num_downstream_slider.val
    print("Showing", num_downstream_temp, "downstream rivers")
    while num_downstream < num_downstream_temp:
        id = downstream_rivers_list[num_downstream]
        if str(id) in p:
            p[str(id)].set_alpha(1)
        if str(id) in rivers:
            j = 0
            for s in rivers[str(id)]:
                rivers[str(id)][j].set_alpha(1)
                j += 1
        if str(id) in discharge_graph_rivers:
            j = 0
            for s in discharge_graph_rivers[str(id)]:
                discharge_graph_rivers[str(id)][j].set_alpha(1)
                j += 1
        num_downstream += 1
    while num_downstream > num_downstream_temp:
        id = downstream_rivers_list[num_downstream - 1]
        if str(id) in p:
            p[str(id)].set_alpha(hidden_alpha * 2)
        if str(id) in rivers:
            j = 0
            for s in rivers[str(id)]:
                rivers[str(id)][j].set_alpha(hidden_alpha * 2)
                j += 1
        if str(id) in discharge_graph_rivers:
            j = 0
            for s in discharge_graph_rivers[str(id)]:
                discharge_graph_rivers[str(id)][j].set_alpha(0)
                j += 1
        num_downstream -= 1
    redraw_canvases()


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


# TODO: add real data to the implementation
def show_propagation(*args, **kwargs):
    """
    Shows the propagation over time

    :param args: Unused parameter to allow function to work as a callback
    :param kwargs: Unused parameter to allow function to work as a callback
    """
    fig_temp = plt.figure()
    Image.open("propagation_time_model.png").show()


# TODO: add real data to the implementation
def show_event_duration(*args, **kwargs):
    """
    Shows the propagation over time

    :param args: Unused parameter to allow function to work as a callback
    :param kwargs: Unused parameter to allow function to work as a callback
    """
    fig_temp = plt.figure()
    Image.open("event_duration_model.png").show()


def update_units(*args, **kwargs):
    """
    Updates what units are being used for distance

    :param args: Unused parameter to allow function to work as a callback
    :param kwargs: Unused parameter to allow function to work as a callback
    """
    global reach_dist_units
    new_units = b_radio_buttons.value_selected
    if new_units == "Kilometers":
        reach_dist_units = "Kilometers"
    elif new_units == "Nautical Miles":
        reach_dist_units = "Nautical Miles"
    elif new_units == "Degrees":
        reach_dist_units = "Degrees"
    elif new_units == "Radians":
        reach_dist_units = "Radians"
    elif new_units == "Miles":
        reach_dist_units = "Miles"
    else:
        print("Error:", new_units, "has not been implemented as a unit option")
    print("Now using", reach_dist_units)
    print("Maintaining a distance of",
          reach_dist, reach_dist_units,
          "between river reaches")


def get_reach_dist_deg():
    """
    Gets the current river reach distance in degrees

    :return: The current river reach distance in degrees
    :rtype: float
    """
    if reach_dist_units == "Kilometers":
        return DistFuncs.km2deg(reach_dist)
    elif reach_dist_units == "Nautical Miles":
        return DistFuncs.M2deg(reach_dist)
    elif reach_dist_units == "Degrees":
        return reach_dist
    elif reach_dist_units == "Radians":
        return np.rad2deg(reach_dist)
    elif reach_dist_units == "Miles":
        return DistFuncs.mi2deg(reach_dist)
    return DistFuncs.km2deg(5)

def show_discharge_over_time(*args, **kwargs):
    """
    Shows the discharge over time based on the number of reaches and the
    distance between reaches

    :param args: Unused parameter to allow function to work as a callback
    :param kwargs: Unused parameter to allow function to work as a callback
    """
    rivers_list = downstream_rivers_list
    if not downstream_rivers_list:
        return
    print("Drawing graphs")
    i = 0
    fig_temp, axs = plt.subplots(nrows=num_reaches,
                                 sharex="all", sharey="all")
    fig_temp.suptitle("River Discharges Over Time")
    fig_temp.supxlabel("Time (3 hours)")
    fig_temp.supylabel(r"Average River Discharge ($m^3/s$)")
    net_distance = 0
    shown = False
    k = 0
    for r in rivers_list:
        if r == 0:
            continue
        print("Graphing River", str(r))
        if str(r) in rivers_list:
            if not shown:
                j = 0
                for s in discharge_graph_rivers[str(r)]:
                    r_temp = discharge_graph_rivers[str(r)][j]
                    axs[k].plot(*r_temp.get_data(),
                                c=r_temp.get_color(),
                                linewidth=r_temp.get_linewidth(),
                                alpha=1,
                                label=str(r))
                    axs[k].legend()
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
    redraw_canvases()
    plt.show(block=False)


# TODO: implement time scrubbing
def fast_forward(*args, **kwargs):
    """
    Moves the main interactive map forward in time

    :param args: Unused parameter to allow function to work as a callback
    :param kwargs: Unused parameter to allow function to work as a callback
    """
    pass


def redraw_canvases(*args, **kwargs):
    """
    Redraws all canvases

    :param args: Unused parameter to allow function to work as a callback
    :param kwargs: Unused parameter to allow function to work as a callback
    """
    fig.canvas.draw()
    fig_graph.canvas.draw()
    fig.canvas.flush_events()
    fig_graph.canvas.flush_events()
    plt.show(block=False)


def clear_plots(*args, **kwargs):
    """
    Clears plots dynamically linked to the main interactive map

    :param args: Unused parameter to allow function to work as a callback
    :param kwargs: Unused parameter to allow function to work as a callback
    """
    global shown_rivers_list
    # for path in p:
    for path in shown_rivers_list:
        if path in p:
            try:
                p[path].set_alpha(hidden_alpha)
            except:
                pass
    # for path in rivers:
    for path in shown_rivers_list:
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
    # for path in discharge_graph_rivers:
    for path in shown_rivers_list:
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
    print("Resetting view")
    start_time = datetime.now()
    # for path in p:
    #     p[path].set_alpha(1)
    # for path in rivers:
    #     j = 0
    #     for s in rivers[path]:
    #         rivers[path][j].set_alpha(1)
    #         j += 1
    # for path in discharge_graph_rivers:
    #     j = 0
    #     while True:
    #         try:
    #             discharge_graph_rivers[path][j].set_alpha(1)
    #         except IndexError:
    #             break
    #         except:
    #             continue
    #         j += 1
    clear_plots()
    truncated_downstream = False
    selected_rivers_list = list()
    currently_selected_river = ""
    downstream_rivers_list = list()
    print("Finished resetting view")
    print("Took", (datetime.now() - start_time).total_seconds(), "to reset")


##### Config Window #####
fig.canvas.mpl_connect("pick_event", on_pick)

fig_config = plt.figure()

### Reset Option ###
ax_reset = plt.axes([0.81, 0.05, 0.1, 0.075])
b_reset = plt.Button(ax_reset, "Reset")
b_reset.on_clicked(reset)

### Show Selected Rivers ###
ax_apply_multiple = plt.axes([0.65, 0.05, 0.15, 0.075])
b_apply_multiple = plt.Button(ax_apply_multiple, "Show Selected")
b_apply_multiple.on_clicked(apply_multiple)

### Show All Downstream Rivers ###
ax_all = plt.axes([0.65, 0.135, 0.15, 0.075])
b_all = plt.Button(ax_all, "All")
b_all.on_clicked(all_one_window_saved)

### Show Downstream Rivers at Evenly Spaced Intervals ###
ax_evenly_spaced = plt.axes([0.65, 0.22, 0.15, 0.075])
b_evenly_spaced = plt.Button(ax_evenly_spaced, "Evenly Spaced")
b_evenly_spaced.on_clicked(isolate_evenly_spaced)

### Number of River Reaches ###
ax_num_reaches = plt.axes([0.65, 0.305, 0.15, 0.075])
b_num_reaches = matplotlib.widgets.TextBox(
    ax_num_reaches,
    label=r"Enter \# of River Reaches",
    initial=num_reaches_str,
)
b_num_reaches.on_submit(update_num_reaches)

### Distance Between Reaches ###
ax_reach_dist = plt.axes([0.65, 0.39, 0.15, 0.075])
b_reach_dist = matplotlib.widgets.TextBox(
    ax_reach_dist,
    label="Enter a Distance Between Reaches",
    initial=reach_dist_str,
)
b_reach_dist.on_submit(update_reach_dist)

### Propagation ###
ax_propagation = plt.axes([0.65, 0.475, 0.15, 0.075])
b_propagation = plt.Button(ax_propagation, "Propagation")
b_propagation.on_clicked(show_propagation)

### Event Duration ###
ax_reach_dist = plt.axes([0.65, 0.56, 0.15, 0.075])
b_reach_dist = plt.Button(ax_reach_dist, "Event Duration")
b_reach_dist.on_clicked(show_event_duration)

### Distance Based Discharge Over Time ###
ax_reach_dist = plt.axes([0.65, 0.645, 0.15, 0.075])
b_reach_dist = plt.Button(ax_reach_dist, "Discharge")
b_reach_dist.on_clicked(show_discharge_over_time)

### Distance Units Options ###
ax_radio_buttons = plt.axes([0.34, 0.645, 0.3, 0.16])
b_radio_buttons = matplotlib.widgets.RadioButtons(
    ax_radio_buttons,
    labels=["Kilometers",
            "Nautical Miles",
            "Degrees",
            "Radians",
            "Miles"],
    active=0)
b_radio_buttons.on_clicked(update_units)

### Checkbox Options ###
ax_checkboxes = plt.axes([0.34, 0.05, 0.3, 0.16])
b_checkboxes = matplotlib.widgets.CheckButtons(
    ax_checkboxes, labels=["Select Multiple",
                           "Open in New Window"])
b_checkboxes.on_clicked(checkboxes_callback)

### Change the Number of Downstream Rivers ###
ax_num_downstream = plt.axes([0.2, 0.25, 0.0225, 0.63])
num_downstream_slider = matplotlib.widgets.Slider(
    ax=ax_num_downstream,
    label="Number of Downstream Rivers",
    valmin=0,
    valmax=max_river_selections,
    valinit=num_downstream,
    valstep=1,
    orientation="vertical",
)
num_downstream_slider.on_changed(update_num_downstream)


##### Finished Preprocessing #####
print("Finished preprocessing")
plt.show(block=True)
