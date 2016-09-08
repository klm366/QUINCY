.. include:: defs.hrst

Configuration
=============

Configuration API Reference
---------------------------

The configuration files for a |b4j| run are validated with
corresponding specification files to check integrity and types of items and
values.

**main configuration and its specification**

.. literalinclude::  ../../configs/b4j_config-NEXT.spec
	   
**Listing:** :specification file for main configuration, :math:`\mathtt{b4j\_config.spec}`



**diagnostic field configuration and its specification**

.. literalinclude::  ../../configs/b4j_field_config-NEXT.spec
	   
**Listing:** :specification file for diagnostic field configuration, :math:`\mathtt{b4j\_field\_config.spec}`

Example configuration for the report section of the diagnostic field `albedo´
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

**Listing:** :math:`\mathtt{albedo.cfg}`, abridged


.. figure::  ./graphic/b4j_vga0187_201605261205.jpg
   :width: 80 %
   :align: center

   **Figure:** resulting albedo section in the report 

**Advantage of static report texts**

* documentation on board
* evaluation criteria for benchmarks supplied directly

