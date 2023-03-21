# GETM-setup_2DV-fjord

**Setup files of a 2D-vertical model of the 79°N Glacier fjord (Greenland)
implemented in GETM, the General Estuarine Transport Model**


Go to the latest release:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7755515.svg)](https://doi.org/10.5281/zenodo.7755515)

*Please cite the DOI of the version you use.
The DOI badge above always redirects to the latest version.*


## Reference

Markus Reinert, Marvin Lorenz, Knut Klingbeil, Bjarne Büchmann, and
Hans Burchard (submitted to JAMES, 2023)


## Instructions to run the model

### Required software

* GETM source code from
  [DOI: 10.5281/zenodo.7741925](https://doi.org/10.5281/zenodo.7741925),
  see [the GETM website](https://getm.eu/) for compilation instructions

* [Flexible output manager](https://github.com/BoldingBruggeman/flexout)
  for a simplifed handling of the model output

* [Editscenario](https://github.com/BoldingBruggeman/editscenario)
  to create the namelist files

* [Ncmerge](https://sourceforge.net/p/getm-utils/wiki/ncmerge/)
  to merge the output files from parallel GETM runs into a single file

* [Python](https://www.python.org/) in version 3 with
  [Xarray](https://xarray.dev/) and [SciPy](https://scipy.org/)
  to create the input files for the model

### Preparation

1. Download this repository and go to the downloaded directory.
2. Create a sub-directory `bin` and copy the GETM executable to `bin/getm`.
3. Look at the script [run.sh](run.sh), check the two settings at the
   beginning of the file, and correct them if needed:
    1. The path `GETMDIR` must point to the code-directory of GETM.
    2. The number of CPU cores in `nCPU` must be available on your computer.

### Running the simulation

1. Modify the parameters in [fjord_322.xml](fjord_322.xml) for the
   experiment you want to run,
   or keep the file as it is for the default scenario.
2. Run the script [run.sh](run.sh), for example with `nohup ./run.sh &`
   to start the model in the background and to keep the model running
   even if you disconnect from the server.
3. Look at the output files in the folder `store`.


## Description

This is a 150 km-long, 20 km-wide 2D-vertical estuary with a maximum
depth of 900 m and a minimum bottom depth of 300 m.  The model has
301×3 grid cells, of which j=1, j=3, and i=1 are land points.
The model is initialized with a horizontally homogeneous temperature
and salinity stratification.  It has an open boundary on the right and
a glacial discharge modeled as river input on the left.
The left part of the domain is covered by glacial ice.  The lower edge
of the ice tongue decreases monotonically from a depth of 600 m at the
grounding line (x = 0) to a depth of 75 m at the calving front
(x = 75 km), then it decreases linearly with a 1%-slope to sea level
(z = 0).
The system reaches an almost steady state after about 6 months.
