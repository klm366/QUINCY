.. include:: defs.hrst


User guide for developing |diag-f| s
====================================

In the following you learn how to develop a |diag-f| and to integrate it
in a |b4j| run. 


A diagnostic for experiment comparison from external sources can be
included as |diag-f| into a |b4j| run. Therefore it needs
scripts for the processing chain (preprocessing, analysing, plotting), a
configuration file and descriptive textblocks for the report. A script may be
written in any programming language that is supported by your platform.
It has to be layed out for one candidate and one reference
so that it can be called in loops over multiple experiments. The scripts work
on local data which have been selected and interpolated from the experiment output as
indicated by the user and the experiment type. They
interface of each script has to be declared in a respective section of the
|diag-f| s' configuration file. The
information to provide include input and output filenames and possibly hints
on format, grid,
spatial resolution, tiling and global fragmentation and time resolution and
period. 


Configuration
-------------
The
configuration for a |diag-f| defines a processing chain that is
interconnected by I/O dependencies. For the time being the processes
have a binding order: 1. :math:`\mathtt{fetch}`,
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
the listing below shows a detail from the configuration file for the
|diag-f| :math:`\mathtt{Albedo}`. Clicking on the graphic opens another picture,
showing how the I/O dependencies are taken into account.

.. figure:: ./graphic/diagfield-conf.jpg
   :target: {{ workdir }}/graphic/diagfield-conf-overlay.jpg

**Listing:** albedo.cfg, example

As can be seen in detail in the :ref:`userguide-c` a  |diag-f| has to be introduced
by a labeled subsection in the :ref:`diag-f` of the main
configuration file providing in particular the filename respectively
the whole local path for the specialized configuration.


Fetch section
+++++++++++++

Each variable that has to be read from the experiment's dataset is represented by a 
subsection labeled with a run internal unique name. The data fetched and
selected from the experiment output will be written to the subdirectory
:math:`\mathtt{interpolated\_selected}` in :math:`\mathtt{work}`, f.e.
:math:`\mathtt{work/myb4jrun/interpolated\_selected}` with a filename containing
this name as infix. The single variable subsection level has subsubsections for
each experiment type holding specialized configuration.

Single variable subsection:

+------------------------------------------------+-----------------------------------+
| key                                            | value                             |
+-----------------+------------------------------+-----------------+-----------------+
| name            | description                  | datatype        | default         |
+=================+==============================+=================+=================+
| longname        | long descriptive name        | string          | None            |
+-----------------+------------------------------+-----------------+-----------------+
| unit            | physical unit                | string          | \-              |
+-----------------+------------------------------+-----------------+-----------------+


:math:`\mathtt{MPIESM1}` subsubsection:

+------------------------------------------------+-----------------------------------+
| key                                            | value                             |
+-----------------+------------------------------+-----------------+-----------------+
| name            | description                  | datatype        | default         |
+=================+==============================+=================+=================+
| name            | variable name in datafile    | string          | \-              |
+-----------------+------------------------------+-----------------+-----------------+
| code            | variable code number in file | string          | \-              |
+-----------------+------------------------------+-----------------+-----------------+
| submodel        | hosting submodel in coupled  | string          | None            |
|                 | model                        |                 |                 |
+-----------------+------------------------------+-----------------+-----------------+
| stream          | stream, hosting datafile     | string          | None            |
+-----------------+------------------------------+-----------------+-----------------+


Processing sections, :math:`\mathtt{preprocessing, diagnose, plotting}`:
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

+------------------------------------------------+-----------------------------------+
| key                                            | value                             |
+-----------------+------------------------------+-----------------+-----------------+
| name            | description                  | datatype        | default         |
+=================+==============================+=================+=================+
| suffix          | format for output files      | option ( nc )   | nc              |
+-----------------+------------------------------+-----------------+-----------------+
| exec            | executable script            | string          | \-              |
+-----------------+------------------------------+-----------------+-----------------+
| input           | data needed from precursor   | string          | \-              |
+-----------------+------------------------------+-----------------+-----------------+


Specialized configuration for the processing section :math:`\mathtt{plotting}`:

+------------------------------------------------+-----------------------------------+
| key                                            | value                             |
+-----------------+------------------------------+-----------------+-----------------+
| name            | description                  | datatype        | default         |
+=================+==============================+=================+=================+
| suffix          | format for output files      | option ( png )  | png             |
+-----------------+------------------------------+-----------------+-----------------+
| nref            | members in comparison        | |long_string-2| | single          |
+-----------------+------------------------------+-----------------+-----------------+
| plot_title      | title for the plot           | string          | None            |
+-----------------+------------------------------+-----------------+-----------------+

.. |long_string-2| replace:: option (single, multiple)


Processing scripts
------------------

For integrated |diag-f| s the processing scripts can be found under
:math:`\mathtt{prod/scripts}`.  If you want to use your own |diag-f| you can
provide the paths for the executable scripts under the item :math:`\mathtt{exec}` in
the respective process' subsection of the configuration file.



As this is work in progress it has to be reckoned with major changes while
more |diag-f| s and other experiment types are integrated. 


