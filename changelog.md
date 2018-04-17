#Changelog for JSBACH4+

##jsbach-4.20p4

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

##jsbach-4.20p3

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

##jsbach-4.20p2

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

##jsbach-4.20p1

_Release Date: 13 February 2018_

##jsbach-4.20
_Release Date: 16 December 2017_
