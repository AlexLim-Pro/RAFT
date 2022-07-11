import csv
import shapefile as shp
import matplotlib.pyplot as plt
from math import isclose
import numpy as np


fig, ax = plt.subplots()
coords_f_path = "../San_Gaud_data/coords_San_Guad.csv"
coords = list()
id_ind = dict()
ind_id = dict()
connectivity_f_path = "../San_Gaud_data/rapid_connect_San_Guad.csv"
# sf_path = "../San_Gaud_data/catchment_Reg12/catchment_Reg12.shp"
sf_path = "../San_Gaud_data/NHDFlowline_San_Guad/NHDFlowline_San_Guad.shp"
connectivity = dict()
connectivity_rev = dict()
p = dict()
rivers = dict()
x_vals = dict()
x_vals_list = list()
y_vals = dict()
y_vals_list = list()
upstream = True
background_path = "../Animation_script/San_Guad_Imagery.tif"

i = 0
with open(coords_f_path, newline="\n") as f:
    for row in csv.reader(f, delimiter=","):
        print("Loading data point", i)
        x_vals[float(row[-2].replace(" ", ""))] = row[0]
        y_vals[float(row[-1].replace(" ", ""))] = row[0]
        x_vals_list.append(float(row[-2].replace(" ", "")))
        y_vals_list.append(float(row[-1].replace(" ", "")))
        coords.append([float(row[-2]), float(row[-1])])
        p[row[0]] = plt.scatter(float(row[-2]), float(row[-1]),
                                s=0.5, picker=5, c="b")
        id_ind[row[0]] = i
        ind_id[i] = row[0]
        i += 1

x_vals_array = np.array(x_vals_list)
y_vals_array = np.array(y_vals_list)


def records(filename):
    # generator
    reader = shp.Reader(filename)
    for sr in reader.shapeRecords():
        geom = sr.shape.__geo_interface__
        yield geom


with open(connectivity_f_path, newline="\n") as f:
    j = 0
    for row in csv.reader(f, delimiter=","):
        print("Loading connectivity information for data point", j)
        data = list()
        data_rev = list()
        for d in row[1:]:
            data.append(int(d))
        connectivity[row[0].replace(" ", "")] = data
        if row[0].replace(" ", "") in connectivity_rev:
            connectivity_rev[row[0].replace(" ", "")] = connectivity_rev[row[0].replace(" ", "")] + data[1:]  # [row[0].replace(" ", "")]
        else:
            connectivity_rev[row[0].replace(" ", "")] = data[1:]
        j += 1

print("Reading shapefile", sf_path)
sf = shp.Reader(sf_path)
print("Finished reading shapefile", sf_path)
j = 0
for shape in sf.shapeRecords():
    print("Drawing shape", j)
    xy = [k for k in shape.shape.points[:]]
    id = ""
    # for c in xy:
    #     if c in coords:
    #         id = str(ind_id[coords.index(c)])
    x, y = zip(*[(k[0], k[1]) for k in xy])
    nearest_vals = list()
    for xj in x:
        temp_array = np.abs(x_vals_array - float(xj))
        close = np.isclose(temp_array, np.zeros(temp_array.shape), atol=1e-13)
        if np.any(close):
            # idx = np.where(close)
            idx = np.abs(temp_array).argmin()
            id = x_vals[x_vals_list[idx]]
            print("idx:", idx)
            print("id:", id)
            print("xj:", xj)
            print("x_vals_list[idx]:", x_vals_list[idx])
            break
    if id == "":
        for yj in y:
            temp_array = y_vals_array - float(yj)
            close = np.isclose(temp_array, np.zeros(temp_array.shape), atol=1e-13)
            if np.any(close):
                # idx = np.where(close)
                idx = np.abs(temp_array).argmin()
                id = y_vals[y_vals_list[idx]]
                print("idx:", idx)
                print("id:", id)
                print("yj:", yj)
                print("y_vals_list[idx]:", y_vals_list[idx])
                break
    # for xi in x_vals:
    #     for xj in x:
    #         if isclose(float(xi), float(xj), abs_tol=1e-13):
    #             id = x_vals[xi]
    #             break
    # if id == "":
    #     for yi in y_vals:
    #         for yj in y:
    #             if isclose(float(yi), float(yj), abs_tol=1e-13):
    #                 id = y_vals[yi]
    #                 break
    if id == "":
        continue
    elif id in rivers:
        continue
    rivers[id] = plt.plot(x, y, linewidth=0.5, alpha=1)
    j += 1


def on_pick(event):
    offsets = event.artist.get_offsets()[event.ind][0]
    i = 0
    for c in coords:
        if c[0] == offsets[0] and c[1] == offsets[1]:
            break
        i += 1
    id = str(ind_id[i])
    paths = list()
    done_paths = list()
    if upstream:
        print("Showing upstream for", id)
        for d in connectivity[id][1:]:
            paths.append(str(d))
    else:
        print("Showing downstream for", id)
        paths.append(connectivity[id][0])
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
        # for s in rivers[path]:
        #     rivers[path][j].set_visible(False)
            j += 1
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
        if upstream:
            for r in connectivity[str(curr_id)][1:]:
                paths.append(str(r))
        else:
            paths.append(str(connectivity[str(curr_id)][0]))
        p[str(curr_id)].set_visible(True)
        done_paths.append(str(curr_id))

    if upstream:
        print("Showing upstream for", id)
        for d in connectivity[id][1:]:
            paths.append(str(d))
    else:
        print("Showing downstream for", id)
        paths.append(connectivity[id][0])
    done_paths = list()
    while paths:
        curr_id = str(paths.pop())
        if str(curr_id) not in rivers:
            continue
        if upstream:
            for r in connectivity[str(curr_id)][1:]:
                paths.append(str(r))
        else:
            paths.append(str(connectivity[str(curr_id)][0]))
        j = 0
        for s in rivers[str(curr_id)]:
            rivers[str(curr_id)][j].set_visible(True)
            j += 1
        done_paths.append(str(curr_id))

    fig.canvas.draw()


def set_downstream(event):
    global upstream
    upstream = False
    print("Set to show downstream")


def set_upstream(event):
    global upstream
    upstream = True
    print("Set to show upstream")


def reset(event):
    print("Resetting view")
    for path in p:
        p[path].set_visible(True)
    for path in rivers:
        j = 0
        for s in rivers[path]:
            rivers[path][j].set_visible(True)
            j += 1
    fig.canvas.draw()
    print("Finished resetting view")


fig.canvas.mpl_connect("pick_event", on_pick)

plt.title("Raw River Coordinate Data")
plt.xlabel("Longitude")
plt.ylabel("Latitude")

plt.figure()
ax_reset = plt.axes([0.81, 0.05, 0.1, 0.075])
b_reset = plt.Button(ax_reset, 'Reset')
b_reset.on_clicked(reset)
ax_downstream = plt.axes([0.49, 0.05, 0.15, 0.075])
b_downstream = plt.Button(ax_downstream, 'Downstream')
b_downstream.on_clicked(set_downstream)
ax_upstream = plt.axes([0.65, 0.05, 0.15, 0.075])
b_upstream = plt.Button(ax_upstream, 'Upstream')
b_upstream.on_clicked(set_upstream)

print("Finished preprocessing")


plt.show()
