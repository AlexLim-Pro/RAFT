# Rivers and Flow Traversal

---

## To Better Understand and Communicate the Propagation of Water Through Space and Time Within Rivers

_**Computer Science Lead:**
Alex Christopher Lim (University of California Santa Barbara in the College of Creative Studies)_

_**Design Leads:**
Isabel (zongzhan) Li (ArtCenter College of Design)
and Racquel Fygenson (Northeastern University)_

_**Mentors:**
Santiago Lombeyda (California Institute of Technology),
Maggie Hendrie (ArtCenter College of Design),
Hillary Mushkin (California Institute of Technology),
and Scott Davidoff (Jet Propulsion Laboratory, California Institute of Technology)_

_**PI:**
Cédric David (Jet Propulsion Laboratory, California Institute of Technology)_

_**Computer Science Collaborators:**
Matt-Heun Hong (University of North Carolina at Chapel Hill)
and Kazi Jawad (Carnegie Mellon University)_

### Acknowledgements

The research was carried out[, in part,] by the Jet
Propulsion Laboratory, California Institute of Technology,
under a contract with the National Aeronautics and Space
Administration (80NM0018D0004).

This work originated from the 2022 JPL/Caltech/ArtCenter
Data to Discovery Program.

---

### Abstract

Despite over 332 million cubic miles of water on Earth, over 97 percent is saline.  Of that three percent of freshwater, 69 percent is locked up in ice and glaciers and 30 percent is underground.  Fresh surface water makes up the remaining one percent; however, it is this 1/150th of one percent of Earth’s water that makes up the large majority of our freshwater supply.  To better understand and observe Earth’s surface water, NASA has established and is developing satellites such as the future SWOT and NISAR missions.  This project–RAFT (Rivers and Flow Traversal)–seeks to use current and future data to better understand and communicate the propagation of water through space and time within rivers.

---

### Runtime Instructions

This code has been tested to function on Python 3.7, 3.8, and 3.9.  Other Python 3 versions and even Python 2 versions should work as well; however, they are not recommended.

Please install the required packages in [requirements.txt](requirements.txt).

Sample data can be found at [https://zenodo.org/record/3688691](https://zenodo.org/record/3688691) or [https://github.com/AlexLim-Pro/RAFT_Sample_Data.git](https://github.com/AlexLim-Pro/RAFT_Sample_Data.git).

Important assets are included in [assets](assets).

When running RAFT, please click on one river reach to begin data analysis.

#### [rrr_raf_tal.py](rrr_raf_tal.py)

Given a river data coordinate csv file, a river connectivity csv file, a river network shapefile, and a river model output netCDF file, this program creates an interactive model to better understand and communicate the propagation of water through space and time within rivers.

By default, these files are listed as being inside a folder called "San_Guad_data" in this file's parent directory.

If you would like to use different files or file paths, you may specify thereof using the command line.

**Example:**

`python3 rrr_raf_tal.py ../San_Guad_data/rapid_connect_San_Guad.csv ../San_Guad_data/NHDFlowline_San_Guad/NHDFlowline_San_Guad.shp ../San_Guad_data/Qout_San_Guad_exp00.nc`
