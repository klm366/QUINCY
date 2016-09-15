.. include:: defs.hrst

Configuration
=============

Configuration API Reference
---------------------------

The configuration files for a |b4j| run are validated with
corresponding specification files to check integrity and types of items and
values.

**Specification for the main configuration**

.. literalinclude::  ../../configs/b4j_config-NEXT.spec
	   
**Listing:** Specification file :math:`\mathtt{b4j\_config.spec}`



**Specification for the diagnostic field configurations**

.. literalinclude::  ../../configs/b4j_field_config-NEXT.spec
	   
**Listing:** Specification file :math:`\mathtt{b4j\_field\_config.spec}`

Example configuration for the report section of the diagnostic field "albedo"
-----------------------------------------------------------------------------

::

   [report]
       title="diagnostic and benchmarks"
       abstract="""The analysis and the benchmarking ..."""
       [[albedo_wds]]
           title="Weighted difference of seasonal means"
	   name=delta
	   type=matrix
	   technical_teXed="""The difference between candid ..."""
	   abstract="""The four figures provide an overview ..."""

**Listing:** Example configuration file :math:`\mathtt{albedo.cfg}`, abridged


.. figure::  ./graphic/b4j_vga0187_201605261205.jpg
   :width: 80 %
   :align: center

   **Figure:** Resulting albedo section in the report 

   
**Remarks on report configuration**
   
* report related items are gathered at the end of configuration files to improve
  readability
* mapping via key correspondence with diagnostic sections

  
**Advantage of static report texts**

* documentation on board
* evaluation criteria for benchmarks supplied directly

