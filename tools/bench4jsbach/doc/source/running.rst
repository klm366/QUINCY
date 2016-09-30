.. include:: defs.hrst

	  
Running
=======


Getting Started
---------------

**Steps**


* check out branch from git repository
* install python, pycmbs, ...
* set up passwordless ssh
* update path variable or change directory
* run bench4jsbach.py


**directory structure**


:bin:     source files
:configs: specification and templates for configuration files
:doc:     Sphinx documentation source files and Makefile
:prod:    scripts and configurations for diagnostic fields
:scripts: essential parts of |b4j|
:tests:   different testing environments
:work:    currently hardcoded working directory


.. figure::  ./graphic/b4j-basedir.jpg
   :width: 35 %
   :align: center

   **Figure:** \<mypath\>/jsbach/tools/bench4jsbach


Running
-------

**Workflow**

:startup:   choose a keyname for run and execute |b4j| with it. This
	    initial run creates
	    
	    * working environment
	    * configuration file
:adapt:     generated standard configuration to your requirements
	    
	    * specify candidate and reference experiments
	    * use all standart diagnostic included or select some tests
	    * (developers: supply a new diagnostic field)
:run:       |b4j| providing the configuration file
	    
	    * fetches required output files of the experiments
	    * interpolates the data according to requirements
	    * runs preprocessing, diagnostics, plotting in loops over experiments
	    * writes report
:study:     report, if required check other |b4j| output, rethink, enjoy :)
:keep:      |b4j| report, delete working directory




Experiments
-----------

**Candidates and references**

* define them by a subsection key

    * you are done for references from wikipage |wiki|
    * otherwise provide type, path, and host/archive inside the dkrz/zmaw domain

* multiple experiments possible
* experiments could serve both in one run

 
.. raw:: html

    <font color=#D8D8D8>Todo: Find a method to grey the now commented lines ('@') in the listing.</font>


::
  
    [experiments]
        [[tra0163]]
	    @host=mistral@
            @path=/mnt/lustre01/work/mh0469/m212070/Cosmos/mpiesm-1.2.00/experiments@
            type=MPIESM1	    
            candidate=True
            reference=False
            report=...analyse XXX albedo scheme...
	    
**Listing:** experiment definition in main configuration file  :math:`\mathtt{main.cfg}`
    

References
----------

You can find reference experiments on the wikipage |wiki|.

.. figure::  ./graphic/b4j-wiki-ref.jpg
   :width: 90 %
   :align: center

   **Figure:** wikipage |wiki|, detail
   
