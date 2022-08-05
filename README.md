# Surface Water and Observable Flow

---

## To Better Understand and Communicate the Propagation of Water Through Space and Time Within Rivers

_Computer Science Lead: Alex Christopher Lim_

_Design Leads: Isabel Li and Racquel Fygenson_

_Mentors: Scott Davidoff, Santiago Lombeyda, Hillary Mushkin, and Maggie Hendrie_

_PI: Cédric David_

---

### Abstract

Despite over 332 million cubic miles of water on Earth, over 97 percent is saline.  Of that three percent of freshwater, 69 percent is locked up in ice and glaciers and 30 percent is underground.  Fresh surface water makes up the remaining one percent; however, it is this 1/150th of one percent of Earth’s water that makes up the large majority of our freshwater supply.  To better understand and observe Earth’s surface water, NASA has established and is developing satellites such as the future SWOT and NISAR missions.  This project–SWOF (Surface Water and Observable Flow)–seeks to use current and future data to better understand and communicate the propagation of water through space and time within rivers.

---

### Runtime Instructions

This code has been tested to function on Python 3.7, 3.8, and 3.9.  Other Python 3 versions should work as well; however, they are not recommended.

Please install the required packages in [requirements.txt](requirements.txt).

#### [rrr_swo_flo.py](rrr_swo_flo.py)

Given a river data coordinate csv file, a river connectivity csv file, a river network shapefile, and a river model output netCDF file, this program creates an interactive model to better understand and communicate the propagation of water through space and time within rivers.

By default, these files are listed as being inside a folder called "San_Guad_data" in this file's parent directory.

If you would like to use different files or file paths, you may specify thereof using the command line.

**Example:**

`python3 rrr_swo_flow.py ../San_Gaud_data/coords_San_Guad.csv ../San_Gaud_data/rapid_connect_San_Guad.csv ../San_Gaud_data/NHDFlowline_San_Guad/NHDFlowline_San_Guad.shp ../San_Gaud_data/Qout_San_Guad_exp00.nc`

#### [rrr_msc_fnc.py](rrr_msc_fnc.py)

Miscellaneous useful functions.