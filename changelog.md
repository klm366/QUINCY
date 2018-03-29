#Changelog for JSBACH4+

##jsbach-4.20p2

_Release Date: 29 March 2018_

- Bug fixes and changes for carbon/fuel/disturbance modules
- Reading of carbon pools from external file at start or restart
- Variable name changes
- Other bug fixes
- Implement (again) setup for terraplanet simulations
- dsl4jsb.py: only process source file if it has changed

  If a source file is processed, the md5 hash value of the file is now written into the header of the target file. The next time dsl4jsb.py is called on the source file, a new target file is only generated if the hash value of the source file and the hash value stored in the target file are different.

- Re-structure creation and handling of vertical axes in order to work with new features in ICON infrastructure

  - Works now with asynchronous I/O.
  - In preparation for use with CDI/PIO.
  - Initialization of JSBACH had to be separated into several subroutines because these have to be called at different stages in ICON.

jsbach-4.20p1
-------------
_Release Date: 13 Feb 2018_

jsbach-4.20
-----------
_Release Date: 16 December 2017_
