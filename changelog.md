# Changelog for JSBACH4+

## jsbach-4.20p5

_Release Date: 20 June 2018_

- Timers
  - Implement constructor function for `t_jsb_process_task` abstract type and add timers for each process task
  - Namelist switches are now `debug_level` and `timer_level`, both can take values between 0 and 3 and are independent of each other
  - Move main JSBACH timer from ICON/ECHAM to JSBACH
- Use proper MPI communicator in `p_bcast` calls (created problem on CRAY with asynchronous output)
- Carbon module
  - Add the summation of `faPAR` from `faPAR_cl`
  - Changed calculation of fuel from canopy to potential vegetation area and adapted fire_litter_threshold
  - Correction of fire bug due to fuel and typo correction of subroutine name
- Some re-factoring in `update_radiation_par`
- Deletion of `l_compat401` switch, changes for accompanying tuning of model and pre-processing script for input data
  - Delete `l_compat401` since standard experiments will now use new snow model and freezing of soil water
  - Change canopy albedo in lctlib
  - Set default for switch `use_alb_canopy` to false in order to actually use canopy albedo values from lctlib (instead of from ini file)
  - Change albedo parameters for glaciers in `mo_rad_constants.f90`
  - Use remapdis for ICON grids above R2B5 instead of remapycon in pre-processing script `jsbach4_ini_files_from_gauss.sh`
  - Add function to modify soil and root depth as part of tuning the model in pre-processing script

    Soil water holding capacity is modified by changing root depth and soil depth.
    They are decreased by 30% in the range 30S - 30N to warm the tropical ontinents
    and increased by 0.2m for the other land grid points to cool the extra-tropical
    continents.
  - Use revision for input data; put data for revision r0002 into new directory in pool
- Further smaller fixes and code cleanup


## jsbach-4.20p4

_Release Date: 17 April 2018_

- Reading of carbon pools from ini file
  - New function `jsbach_init_after_restart` which now reads the carbon pools from the ini file **after** the restart files have been read, but **before** the output is initialized. In ICON, the output of the first step now contains the correct pools that have just been read in. Note that the run script has to take care of re-setting `read_cpools` to `.FALSE.` for the jobs following the first one!
  - Some minor changes (comment for `read_cpools` config parameter and changes in script `cpool_file_from_restart_file.sh` for generation of cpool file from restart file)
- Re-structuring of aggregator

  Leads to a few percent of speed gain in aggregation
- Change in dsl4jsb.py pre-preprocessing script

  If the content of the source file hasn't changed compared to the target file, reset the modification time of the source file to the one of the target file, so that make doesn't think the source file is newer and therefore tries to process the file (but `dsl4jsb.py` doesn't create a new target since the content hasn't changed ... this can lead to an infinite loop in make)
- Some changes to `jsbach4_ini_files_from_gauss.sh` script for generation of input files
  - Set URL of svn repo to new address `svn.mpimet.mpg.de`
  - Allow generation of files for coupled setup for non-r2b4 resolutions (might still need some work)

## jsbach-4.20p3

_Release Date: 10 April 2018_

- Give back new saturated specific humidity over lakes (water and ice) to the atmosphere

  With the corresponding fix in the ICON atmosphere which puts `qsat_lwtr` and `qsat_lice` into the `bb_btm` array, this closes a water leakage bug.
- Some bug fixes in the carbon module
- Some other bug fixes
- Some changes to `jsbach4_ini_files_from_gauss.sh` script for generation of input files
  - updated module versions (cdo etc)
  - output files for debugging are in T255 (icon)
  - updated files for coupled setup
  - changed default cdo remap scheme from "con" to "ycon"

## jsbach-4.20p2

_Release Date: 29 March 2018_

- Bug fixes and changes for carbon/fuel/disturbance modules
- Reading of carbon pools from external file at start or restart
- Variable name changes
- Other bug fixes
- Implement (again) setup for terraplanet simulations
- `dsl4jsb.py`: only process source file if it has changed

  If a source file is processed, the md5 hash value of the file is now written into the header of the target file. The next time `dsl4jsb.py` is called on the source file, a new target file is only generated if the hash value of the source file and the hash value stored in the target file are different.
- Re-structure creation and handling of vertical axes in order to work with new features in ICON infrastructure
  - Works now with asynchronous I/O.
  - In preparation for use with CDI/PIO.
  - Initialization of JSBACH had to be separated into several subroutines because these have to be called at different stages in ICON.

## jsbach-4.20p1

_Release Date: 13 February 2018_

## jsbach-4.20

_Release Date: 16 December 2017_
