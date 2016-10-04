.. include:: defs.hrst

Introduction
============

Main idea
---------

With |b4j| we started developing a tool for the scientists of
the land department in the |mpi| to make the
comparison of |jsbach| model simulation output data [ReJSBA]_ with one another and with data
products straight forward and user friendly. In the first instance the tool
shall serve for the benchmarking of land models, it has to be able to quantify
the deviation of a simulation from reference experiments in the most crucial and/or
influential parameters out-of-the-box. This will help developing models
as well as managing a large number of them. It is desirable to cover most
research areas of the land department and to keep the tool well adapted
to the local technical conditions at the |mpi|.

In order to have a living, dynamical tool, the possibility for individual
extension is provided to the user from the beginning. This will help the
benchmarking process to
retain its validity but also help the scientists to handle future questions
bringing new values into focus. It may foster the communication between
the groups as individual developments can easily be made usable for the
whole department and thus serve as a basis for a fruitfuly discussion on the most
significant diagnostic for the generally applicable benchmarking in the
department.

Concepts
--------

A guiding idea for the tool is the separation of issues that strongly depend on the
origin of the climate data from issues directly related with their
comparison.

For the first we use :math:`\mathtt{experiment}`
as a generic term, independent from the data origin being an earth system
model simulation or an
observation. Issues in this context are f.e. locating and extracting
necessary raw data on/from  different hosts, filesystems and datafiles but also
data interpolation. It is intended to
solve these tasks
'under the hood' so that a scientist can
focus on the diagnostic and benchmarking methods.

The term |diag-f| clusters
processing with a defined set of variables and a particular scientific question.
It should contain at least one diagnostic and one benchmarking method.  In this
context the term
benchmarking refers to a method that condenses the deviation of two
experiments to one number, i.e. yields a global and time average value.


Context
-------

The following diagram shows |b4j| in context with other actors ( users, software, ... ). 

.. figure:: ./graphic/context-diagram.jpg
   :width: 40 %
   :align: center

   **Figure:** Context diagram of |b4j|

========================    ===============
**Benchmarking user**       * scientist managing many models
                            * model developer making  still-in-line-check  
**Diagnostic developer**    dedicated scientist integrating a |diag-f|
**Database**                database for reference |exp| s
**Redmine wiki**            communication platform for |b4j|
========================    ===============




:Benchmarking user:             * scientist managing many models
                                * model developer making  still-in-line-check  
:Diagnostic developer:          dedicated scientist integrating a |diag-f|
:Database:    \ \               database for reference |exp| s
:Redmine wiki:   \              communication platform for |b4j|

Status
------

This documentation is based on the vertical prototype from a not yet tagged
version ( `0.0.1` ) of |b4j|. A basic 
operating concept, the program structure and the workflow is implemented. The
tool handles |exp| s of type  MPIESM1 (Max Planck
Institute for Meteorology Earth System Model, [GiESM1]_) and a |diag-f| analysing
albedo.
