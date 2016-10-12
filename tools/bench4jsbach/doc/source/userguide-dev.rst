.. include:: defs.hrst


User guide for developing |diag-f| s
====================================

In the following we show how a diagnostic for |exp| comparison from
external sources can be prepared to serve as |diag-f| that may be triggered by
a |b4j| run. A |diag-f| consists of 
scripts for the processing chain (preprocessing, analysing, plotting), declarations
also for data retrieval and interpolation, descriptive textblocks for generating reports
and a configuration file to bring everything together. The scripts may be
written in any programming language that is supported by your platform.
A script has to be layed out for one candidate and one reference
so that it is scalable, i.e. it can be called in loops over multiple |exp| s. It
will work on local data which the tool has selected and interpolated
from the |exp| output beforehand as
it had been indicated by the user and the |exp| type. The
interface of each script has to be declared in a respective section of the
|diag-f| s' configuration file. The
information to provide include input and output filenames and possibly hints
on format, grid,
spatial resolution, tiling and global fragmentation and time resolution and
period. 


We use the example :math:`\mathtt{myb4jrun}` from the
last section :ref:`userguide-com` again.
In the run we want to include some
diagnostic and benchmarking for a clear cut out example scientific issue called
:math:`\mathtt{myfield}`.

Declaration in main configuration
---------------------------------
      
The central entry point for a
|diag-f| into the run is the :math:`\mathtt{fields}` :ref:`diag-f-section` in the main
configuration file, in this example :math:`\mathtt{myb4jrun.cfg}`. We declare a subsection
called :math:`\mathtt{myfield}` and fill in the path to a new field configuration file
as value for the key :math:`\mathtt{config\_filename}`, here
:math:`\mathtt{mypath/myfield.cfg}`.

.. literalinclude::  ./graphic/mydiagfield.cfg
	   
**Listing:** Declaration of :math:`\mathtt{myfield}` in main configuration file :math:`\mathtt{myb4jrun.cfg}`

Configuration file :math:`\mathtt{myfield.cfg}`
-----------------------------------------------

Now we want to develop this field configuration file. It has to contain all the settings and
communication parameters that are needed for a |diag-f| from the |exp| output data to the
benchmark results in a final report.

Top level config
++++++++++++++++

For verification it will be checked if the value for the key
:math:`\mathtt{field\_name}` is consistent with the name given for the
:math:`\mathtt{field}` in the main configuration.
      
+------------------------------------------------+-----------------------------------+
| key                                            | value                             |
+-----------------+------------------------------+-----------------+-----------------+
| name            | description                  | datatype        | default         |
+=================+==============================+=================+=================+
| field_name      | |long-string-4|              | string          | \-              |
+-----------------+------------------------------+-----------------+-----------------+
| suffix          | |long-string-5|              | string          | \-              |
+-----------------+------------------------------+-----------------+-----------------+
| frequency       | |long-string-6|              | |long-string-7| | 'monthly'       |
+-----------------+------------------------------+-----------------+-----------------+
| documentation   | |long-string-8|              | string          | None            |
+-----------------+------------------------------+-----------------+-----------------+

.. |long-string-4| replace:: field id conform to id/name of diagnostic field in main configuration
.. |long-string-5| replace:: output format suffix for processed data of the run
.. |long-string-6| replace:: output frequence to use in the run
.. |long-string-7| replace:: option('yearly','monthly','6-hourly')
.. |long-string-8| replace:: self documentation

In our example we fill in the respective values.

::

   field_name=myfield
   suffix=nc
   frequence=monthly
   documentation="""
   NAME 
	   myfield_diagnostic

   DESCRIPTION
	   This configuration ...

   DETAILS
	   For all experiments ...
	   """
	   
**Listing:** Example top level configuration

Processing chain
++++++++++++++++

The processes involved in the |diag-f| are
interconnected by I/O dependencies. For the time being they have a binding
order: 1. :math:`\mathtt{fetch}`,
2. :math:`\mathtt{interpolate}`, 3. :math:`\mathtt{preprocess}`,
4. :math:`\mathtt{diagnose}`, 5. :math:`\mathtt{plot}`. |b4j| works internally
with an indexing system in which files are addressed by unique infixes,
f.e. :math:`\mathtt{FILE-INFIX}`.  For each
file that shall be written by a process the appropriate section has to hold a subsection :math:`\mathtt{[[FILE-INFIX]]}`
wich contains the detailed configuration. This infix
will also be used as substring for the respective filename. If the file
has to be read by a subsequent process the infix is used in the respective
configuration as value for its
:math:`\mathtt{input}` item. For illustration  
the listing at the bottom of the page :ref:`io-dependencies-html` shows a detail from the configuration file for the
|diag-f| :math:`\mathtt{myfield}`. In HTML clicking on the graphic opens another picture,
showing how the I/O dependencies are taken into account.


Fetch section
+++++++++++++

For each variable that the |diag-f| requires from the |exp| output the |b4j| run will
write a file containing the selected data field. You have to define a subsection f.e.
:math:`\mathtt{FILE-INFIX}` for it.


Single subsection in section  :math:`\mathtt{fetch}`:

+------------------------------------------------+-----------------------------------+
| key                                            | value                             |
+-----------------+------------------------------+-----------------+-----------------+
| name            | description                  | datatype        | default         |
+=================+==============================+=================+=================+
| longname        | complete name                | string          | None            |
+-----------------+------------------------------+-----------------+-----------------+
| unit            | physical unit                | string          | \-              |
+-----------------+------------------------------+-----------------+-----------------+

Underneath there has to be one subsubsection for each :math:`\mathtt{type}` given
for any item in the |exp| section of the main configuration. It must be labeled with
the name of the :math:`\mathtt{experiment\ type}`.

Subsubsection for the :math:`\mathtt{MPIESM1}` :math:`\mathtt{experiment\ type}`:


+------------------------------------------------+-----------------------------------+
| key                                            | value                             |
+-----------------+------------------------------+-----------------+-----------------+
| name            | description                  | datatype        | default         |
+=================+==============================+=================+=================+
| name            | |long-string-9a|             | string          | \-              |
+-----------------+------------------------------+-----------------+-----------------+
| code            | |long-string-9b|             | option          | \-              |
+-----------------+------------------------------+-----------------+-----------------+
| submodel        | |long-string-9c|             | string          | None            |
+-----------------+------------------------------+-----------------+-----------------+
| stream          | |long-string-10|             | string          | None            |
+-----------------+------------------------------+-----------------+-----------------+

.. |long-string-9a| replace:: name of variable to select
.. |long-string-9b| replace:: code number of variable to select 
.. |long-string-9c| replace:: submodel identifier in the original context
.. |long-string-10| replace:: stream identifier in the original context

			      
For clarification we show  one variable subsection labeled
‘:math:`\mathtt{D\_selcode}`’  of the :math:`\mathtt{fetch}` section from
our example :math:`\mathtt{mypath/myfield.cfg}`.

::

   [fetch]
       [[D_selcode]]
	   longname=Surface Downwelling Solar Radiation [W m-2]
	   unit=Wm**(-2)
	   [[[MPIESM1]]]
	       name=swdown_acc
	       code=21
	       submodel=jsbach
	       stream=jsbach
	       
**Listing:** Example :math:`\mathtt{fetch}` section configuration


Preprocessing and diagnostic sections
+++++++++++++++++++++++++++++++++++++

To do preprocessing, diagnostic and plotting |b4j| calls executable scripts. In order to
run a new diagnostic The name of the script is sufficient if it
is located in :math:`\mathtt{bench4jsbach/prod/scripts}`.

Top level in this sections:

+------------------------------------------------+-----------------------------------+
| key                                            | value                             |
+-----------------+------------------------------+-----------------+-----------------+
| name            | description                  | datatype        | default         |
+=================+==============================+=================+=================+
| exec            | |long-string-11|             | string          | None            |
+-----------------+------------------------------+-----------------+-----------------+
| input           | |long-string-12|             | stringlist      | \-              |
+-----------------+------------------------------+-----------------+-----------------+

.. |long-string-11| replace:: executable script to call in b4j run
.. |long-string-12| replace:: input files needed by script


In the example below the executable script :math:`\mathtt{myfield\_prep.py}`
shall run for :math:`\mathtt{preprocessing}`.  It has to read a data file 
::
   
   [preprocessing]
       exec=mypath/myfield_prep.py
       input=D_selcode,
       [[D_seasmean]]
	       
**Listing:** Example :math:`\mathtt{preprocessing}` section configuration


.. _io-dependencies-html:

I/O dependencies
++++++++++++++++

.. figure:: ./graphic/diagfield-conf.jpg
   :target: {{ workdir }}/graphic/diagfield-conf-overlay.jpg

**Listing:** myfield.cfg, example
