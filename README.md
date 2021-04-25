# Anna Du's Ocean Particle Simulation System

## Introduction

This is my science project to simulate microplastics movement in the ocean environment. The work is based on [NASA's ECCO Dataset](https://data.nas.nasa.gov/ecco/).

## Install

### Prerequisite

1. Anaconda Python 3

1. [Nasa ECCO Python Package](https://ecco-v4-python-tutorial.readthedocs.io/Installing_Python_and_Python_Packages.html). `ECCOv4-py` should be installed under home directory or any directory in system environment variable PYTHONPATH.

1. [NASA ECCO Dataset](https://ecco-v4-python-tutorial.readthedocs.io/Downloading_the_ECCO_v4_state_estimate.html). This simulator uses Version 4 Release 4. The default location is in `~/eccodata/`

## Simulate

### Input file

The simulator takes an input file that specifies inital particle location in `.csv` format. Each row represents a particle with format of `tile,x,y,k` where k can be omitted (by default 0), where tile is 0-12, x and y are 1-89 and k is 0-49. See [ECCO Illustration](https://ecco-group.org/images/ecco_tiles.png)
[![ECCO Illustration](https://ecco-group.org/images/ecco_tiles.png)](https://ecco-group.org/analysis-tools.htm)

### Simulate 1027 particles and plot the global map

`python3 ocean_particle_simulator.py input-1027p.csv --plot-all-lonlat`

### Simulate 50 particles in tile 10 and plot only tile 10

`python3 ocean_particle_simulator.py input-50p --plot-1tile 10`

### Simulate 1027 particles and adjust initial k speed to 0.2, and fudge factor 20%, without plotting

`python3 ocean_particle_simulator.py input-1027p.csv --kvel 0.2 --fudge-pct 20 --only-compute`

### Use a previously computed result to plot a tile 2 video

`python ocean_particle_simulator.py results_input-1027p_f30kv0.15 --only-plot --plot-1tile 2`
