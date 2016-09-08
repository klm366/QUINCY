.. include:: defs.hrst


Runtime
=======


Fetching Data
-------------

**Specialized experiment type handles ...**

* provide output streams from

  * remote server via passwordless ssh 
  * archive by means of pftp, storing it to :math:`\mathtt{downloads}`-directory

* select variables and periods needed for diagnostic

  * using e.g. |cdo|
  * write it to :math:`\mathtt{processed}`-directory using key/infix from configuration
	    
.. figure::  ./graphic/workingdir_vga0186.jpg
   :width: 50 %	   
   :align: center

   **Figure:** working directory expanded to focus attention on fetched and selected data, example

      
