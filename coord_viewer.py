import csv
import itertools

import matplotlib.pyplot as plt
import matplotlib.widgets
import numpy as np
import shapefile as shp
from scipy.io import netcdf_file

kfac_scaling = False
discharge_scaling = True
river_length_scaling = False


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
river_lengths_list = list()
river_lengths_max = 0
upstream = False
background_path = "../Animation_script/San_Guad_Imagery.tif"
kfac_path = "../San_Gaud_data/kfac_San_Guad_1km_hour.csv"
kfacs = list()
max_kfacs = 0
Qout_f_path = "../San_Gaud_data/Qout_San_Guad_exp00.nc"
default_point_color = "b"
discharge_graph_rivers = dict()
river_colors = dict()
num_downstream = 5
currently_selected_river = ""
downstream_rivers_list = list()
truncated_downstream = False
select_multiple = False
select_multiple_label = "Select Multiple"
selected_rivers_list = list()

cycol = itertools.cycle("bgrcmk")


if kfac_scaling:
    i = 0
    with open(kfac_path, newline="\n") as f:
        for row in csv.reader(f, delimiter=","):
            print("Loading Muskingum K value", i)
            kfacs.append(float(row[0]))
            i += 1

    max_kfacs = max(kfacs)

Qout_f = netcdf_file(Qout_f_path, "r")

temp = Qout_f.variables["Qout"]
Qout_data = temp[:]*1
"""average river water discharge downstream of each river reach"""

temp = Qout_f.variables["Qout_err"]
Qout_data_err = temp[:]*1

temp = Qout_f.variables["rivid"]
"""average river water discharge uncertainty downstream of each river reach"""
Qout_data_ids = temp[:]*1
"""unique identifier for each river reach"""

Qout_data_max = max(abs(Qout_data[0]))
Qout_data_err_max = max(abs(Qout_data_err[0]))

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
        if kfac_scaling:
            size = kfacs[i] / max_kfacs
            color = default_point_color
        elif discharge_scaling:
            idx = np.where(Qout_data_ids == int(row[0]))
            size = Qout_data[0][idx] / Qout_data_max
            color = default_point_color
        elif river_length_scaling:
            draw_point = False
        else:
            size = 1
            color = default_point_color
        if draw_point:
            p[row[0]] = plt.scatter(float(row[-2]), float(row[-1]),
                                    s=size, picker=5, c=color)
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
        d = np.sqrt((k[0] - xy[i][0]) ** 2 + (k[1] - xy[i][1]) ** 2)
        river_lengths[id] = d
        river_lengths_list.append(d)
        i += 1
    color = next(cycol)
    rivers[id] = plt.plot(x, y, linewidth=0.5, alpha=1, c=color)
    river_colors[str(id)] = color
    j += 1


river_lengths_max = max(river_lengths_list)
if river_length_scaling:
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
            if row[0] in river_lengths:
                size = river_lengths[row[0]] / river_lengths_max
                color = default_point_color
                p[row[0]] = plt.scatter(float(row[-2]), float(row[-1]),
                                        s=size, picker=5, c=color)
            id_ind[row[0]] = i
            ind_id[i] = row[0]
            i += 1

plt.title("Raw River Coordinate Data")
plt.xlabel("Longitude")
plt.ylabel("Latitude")

Qout_graph_x_vals = list()
Qout_graph_y_vals = list()
fig_graph = plt.figure()
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
fig.canvas.draw()


def on_pick(event):
    """
    Handles selection of points
    """
    global currently_selected_river
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
    if upstream:  # This is no longer supported
        print("Showing upstream for", id)
        for d in connectivity[id][1:]:
            paths.append(str(d))
    else:
        print("Showing downstream for", id)
        paths.append(connectivity[id][0])
    clear_plots()
    p[str(id)].set_visible(True)
    j = 0
    if str(id) in rivers:
        for s in rivers[str(id)]:
            rivers[str(id)][j].set_visible(True)
            j += 1
    done_paths.append(str(id))
    while paths:
        curr_id = str(paths.pop())
        if str(curr_id) not in p:
            continue
        if upstream:  # This is no longer supported
            for r in connectivity[str(curr_id)][1:]:
                paths.append(str(r))
                downstream_rivers_list.append(str(r))
        else:
            paths.append(str(connectivity[str(curr_id)][0]))
            downstream_rivers_list.append(str(connectivity[str(curr_id)][0]))
        p[str(curr_id)].set_visible(True)
        done_paths.append(str(curr_id))

    if upstream:  # This is no longer supported
        print("Showing upstream for", id)
        for d in connectivity[id][1:]:
            paths.append(str(d))
    else:
        print("Showing downstream for", id)
        paths.append(connectivity[id][0])
    paths = done_paths
    done_paths = list()
    while paths:
        try:
            curr_id = str(paths.pop())
            if str(curr_id) not in rivers:
                continue
            j = 0
            for s in rivers[str(curr_id)]:
                try:
                    rivers[str(curr_id)][j].set_visible(True)
                except:
                    pass
                j += 1
            done_paths.append(str(curr_id))
        except:
            pass

    if upstream:  # This is no longer supported
        print("Showing upstream for", id)
        for d in connectivity[id][1:]:
            paths.append(str(d))
    else:
        print("Showing downstream for", id)
        paths.append(connectivity[id][0])
    paths = done_paths
    done_paths = list()
    while paths:
        curr_id = str(paths.pop())
        if str(curr_id) not in discharge_graph_rivers:
            continue
        j = 0
        for s in discharge_graph_rivers[str(curr_id)]:
            if str(curr_id) in discharge_graph_rivers:
                discharge_graph_rivers[str(curr_id)][j].set_visible(True)
            j += 1
        done_paths.append(str(curr_id))
    redraw_canvases()


def set_downstream(event):
    """
    Sets the option to downstream
    """
    global upstream
    upstream = False
    print("Set to show downstream")


def set_upstream(event):
    """
    Sets the option to upstream
    """
    global upstream
    upstream = True
    print("Set to show upstream")


def apply_multiple(event):
    """
    Applies the selected rivers
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
        if i > 10:
            break
        i += 1
    plt.show(block=False)


def set_select_multiple(event):
    """
    Toggles the select multiple rivers setting
    """
    global select_multiple
    global selected_rivers_list
    select_multiple = not select_multiple
    if select_multiple:
        print("Selecting Multiple")
    else:
        print("Selecting Entire Downstream")
        selected_rivers_list = list()


def update_num_downstream(event):
    """
    Updates the number of downstream rivers to show
    """
    global num_downstream
    global truncated_downstream
    if not truncated_downstream:
        num_downstream = len(downstream_rivers_list)
        truncated_downstream = True
    num_downstream_temp = num_downstream_slider.val
    print(num_downstream_temp)
    while num_downstream < num_downstream_temp:
        id = downstream_rivers_list[num_downstream]
        if str(id) in p:
            p[str(id)].set_visible(True)
        if str(id) in rivers:
            j = 0
            for s in rivers[str(id)]:
                rivers[str(id)][j].set_visible(True)
                j += 1
        if str(id) in discharge_graph_rivers:
            j = 0
            for s in discharge_graph_rivers[str(id)]:
                discharge_graph_rivers[str(id)][j].set_visible(True)
                j += 1
        num_downstream += 1
    while num_downstream > num_downstream_temp:
        id = downstream_rivers_list[num_downstream - 1]
        if str(id) in p:
            p[str(id)].set_visible(False)
        if str(id) in rivers:
            j = 0
            for s in rivers[str(id)]:
                rivers[str(id)][j].set_visible(False)
                j += 1
        if str(id) in discharge_graph_rivers:
            j = 0
            for s in discharge_graph_rivers[str(id)]:
                discharge_graph_rivers[str(id)][j].set_visible(False)
                j += 1
        num_downstream -= 1
    redraw_canvases()


def redraw_canvases():
    """
    Redraws all canvases
    """
    fig.canvas.draw()
    fig_graph.canvas.draw()
    fig_config.canvas.draw()


def clear_plots():
    for path in p:
        try:
            p[path].set_visible(False)
        except:
            pass
    for path in rivers:
        j = 0
        while True:
            try:
                rivers[path][j].set_visible(False)
            except IndexError:
                break
            except:
                continue
            j += 1
    for path in discharge_graph_rivers:
        j = 0
        while True:
            try:
                discharge_graph_rivers[path][j].set_visible(False)
            except IndexError:
                break
            except:
                continue
            j += 1


def reset(event):
    """
    Resets the view
    """
    print("Resetting view")
    for path in p:
        p[path].set_visible(True)
    for path in rivers:
        j = 0
        for s in rivers[path]:
            rivers[path][j].set_visible(True)
            j += 1
    for path in discharge_graph_rivers:
        j = 0
        while True:
            try:
                discharge_graph_rivers[path][j].set_visible(True)
            except IndexError:
                break
            except:
                continue
            j += 1
    fig.canvas.draw()
    fig_graph.canvas.draw()
    print("Finished resetting view")


fig.canvas.mpl_connect("pick_event", on_pick)

fig_config = plt.figure()
ax_reset = plt.axes([0.81, 0.05, 0.1, 0.075])
b_reset = plt.Button(ax_reset, "Reset")
b_reset.on_clicked(reset)
ax_apply_multiple = plt.axes([0.65, 0.05, 0.15, 0.075])
b_apply_multiple = plt.Button(ax_apply_multiple, "Apply")
b_apply_multiple.on_clicked(apply_multiple)
ax_select_multiple = plt.axes([0.44, 0.05, 0.2, 0.075])
b_select_multiple = matplotlib.widgets.CheckButtons(
    ax_select_multiple, labels=[select_multiple_label])
b_select_multiple.on_clicked(set_select_multiple)
ax_num_downstream = plt.axes([0.2, 0.25, 0.0225, 0.63])
num_downstream_slider = matplotlib.widgets.Slider(
    ax=ax_num_downstream,
    label="Number of Downstream Rivers",
    valmin=0,
    valmax=10,
    valinit=num_downstream,
    valstep=1,
    orientation="vertical",
)
num_downstream_slider.on_changed(update_num_downstream)

print("Finished preprocessing")


plt.show(block=True)
