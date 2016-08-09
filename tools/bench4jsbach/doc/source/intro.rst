.. include:: defs.hrst

Introduction
============

The tool |b4j| is being developed in order to make the comparison of 
|jsbach| land model runs with one another and with data products straight
forward and user friendly. For this purpose it provides a framework to
separate the associated tasks

1. fetch, select and interpolate data
2. perform analyses, benchmarks and plotting
3. document the results in a report

It is capable to cope with multiple experiments, candidates as well as
references. Each candidate experiment is compared with all references,
pairwise or as a group depending on the diagnostic.

|b4j| is flexible and extensible to grow inline with
requirements. It is written in object oriented Python and uses
specialized classes for the tasks enumerated above.


Whole department knowhow for your model evaluation
--------------------------------------------------

**From the guiding vision ...**

* enable the groups to contribute their diagnostic
* make the most significant diagnostic usable for the whole department
* solve technical issues once, in particular concerning data retrieval and preparation
* develop a library for diagnostic, benchmark and plotting
  
**... to the usable potential**

* support your model development
* benchmark model output


Status
------

The basic functionality of |b4j| is available. The 

**Experiment Type:**

* MPIESM1

**Diagnostic Field:**

* Albedo
* Hydrology

are placed at your disposal.
