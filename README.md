### Loop Extrusion Simulations
Before starting, install the envorinment with `conda env create -f env.yml`. \
Clone this repository to your local machine with `git clone`
#### 1-D Trajectory
1. Navigate to and open `1D_polychrom_simulation.py`
2. Change the following variables if necessary:
    - `RUN_NAME` (line 21): A descriptive name for the simulation run
    - `N1_pol` (line 23): The size of the polymer in monomers
    - `front_buffer` (line 25): The number of monomers to add to beginning of polymer as a buffer zone
    - `end_buffer` (line 26): Same as front_buffer, but for end of polymer
    - `LIFETIME` (line 42): The number of monomers a LEF can extrude on average. Probability of unloading is `1/LIFETIME` (probably won't need to change this)
    - Blocking region arrays (lines 55-60): An array of integers defining where blockers should be placed. The lines look like this: `np.arange(181+front_buffer, 185+front_buffer)`, where the blocking region runs from 181 to 185.
    - Blocking region dictionary `blockingRegions` (line 65): Maps a name for each blocking region to the array containing the monomers it spans. If we want the blocking region to be bidirectional, we add `"_EBF1"` to the end of the name.
    - `cap` and `rel` (lines 68-91): Capture and release probabilities $[0..1]$. Each blocking region has one of each, and are denoted by the `if` statements.
    - Cohesin _loading_ regions (lines 137-142): Only for loading scenario; similar to the blocking regions, we set lower and upper bounds in the arrays.
    - `loading_regions` (line 148) can be set to wither an array of all possible loading regions (again, _only for loading scenario_), or set to just the array with coordinates spanning the entire array (_blocking scenario only_)
    - `loading_region_freqs` (line 151): One integer per loading region denoting how many LEFs should be loaded there. In _blocking_, this is a single number.
  
  3. Save changes and quit
  4. Execute `python3 1D_polychrom_simulation`. This should take about 10-20 seconds.

#### 3-D trajectory
1. The outputs from the 1-D trajectory are now in `1D_trajectory/trajectory`. **You do not need to move this.**
2. Navigate from `1D_trajectory/` to `3D_simulation/`
3. Execute `python3 3D_polychrom_simulation`. This will take ~45 minutes to complete on an NVIDIA T400 4GB GPU.
4. After this completes, you will have the directory `3D_simulation/sim_outs/`

#### Making contact matrix
1. Ensure you are still in `3D_trajectory/`
2. Execute `./trajectory_to_txt.py 10000 ./sim_outs` (`10000` = number of conformations)
3. This will generate `confs_txt/`
4. Execute `./make_contactMap.py ./confs_txt`
5. This will generate `matrix.txt`, the contact matrix
6. Plot `matrix.txt` in R, and remember to account for the buffer zones when plotting / analyzing.
