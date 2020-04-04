:tocdepth: 3

.. _protocols:

Protocols
----------------

This section describes a number of basic protocols for using Icepack in other models.

.. _calling:

Access
~~~~~~~~~~~~~~~~~~~

Icepack provides several public interfaces.  These are defined in **columnphysics/icepack\_intfc.F90**.  
Icepack interfaces all contain the icepack\_ prefix.
Icepack interfaces follow a general design where data is passed in on a gridpoint by gridpoint
basis, that data is updated and returned to the driver, and the data is not stored within Icepack.  
Additional information about the interfaces can be found in :ref:`sequence_and_interface`.

Icepack interfaces can have long argument lists.  These are documented in :ref:`docintfc`.  In
some cases, arguments are required for optional features (i.e. biogeochemistry) even when that
feature is turned off in Icepack.  The Icepack
development team continues to work towards having more optional arguments.  If an argument is 
required for the interface but not needed, the driver will still have to pass a (dummy) variable 
thru the interface to meet the interface specification.

.. _initialization:

Initialization
~~~~~~~~~~~~~~~~~~~~~~~~~~

The subroutine **icepack_configure** should be called before any other icepack interfaces are called.
This subroutine initializes the abort flag and a few other important defaults.  We recommend that
call be implemented as::

      call icepack_configure()  ! initialize icepack
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call my_abort_method()

The 2nd and 3rd line above are described further in :ref:`aborts`.


.. _aborts:

Error Messages and Aborts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Icepack does not generally handle I/O (file units), the parallel computing environment (MPI, etc),
or model aborts.  Icepack generates and buffers error messages that can be accessed by the
driver.  In addition, if Icepack fails, it will set an abort flag that can be queried by the driver.
To best use those features, it's recommended that after every icepack interface call, the user
add the following::

      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call my_abort_method()

**icepack_warnings_flush** is a public interface in icepack that writes any warning or error messages
generated in icepack to the driver file unit number defined by nu_diag.  
The function **icepack_warnings_aborted** queries the internal icepack abort flag and
returns true if icepack generated an abort error.  
my_abort_method represents a method in the driver that will abort the model cleanly.

In addition to writing Icepack messages thru the icepack_warnings_flush interface,
there are also several methods in icepack that write general information to a file.  
The various **icepack_write_** interfaces accept a unit number provided by the driver
and then document internal Icepack values.

.. _setinternal:

Setting Internal Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

While Icepack does not generally store the model state, there are several Icepack interfaces
that allow the driver to set various scientific and technical parameters internally in Icepack
for later use.  Those interfaces generally have a set (init), get (query) and write method that
allows the driver to set, get, or write the values defined internally.  Some parameters
are required to be set by the driver and others take on defaults.  The following table
defines the available interfaces that fit into this category.

.. table:: *Init, Query, and Write Interfaces* 

   +----------------+---------------------------------+----------------------------------+----------------------------------+-------------------------------------------------+
   | type           | init                            |                      query       |                      write       |           notes                                 |
   +================+=================================+==================================+==================================+=================================================+
   | orbital        | icepack\_init\_orbit            | icepack\_query\_orbit            |                                  | orbital settings                                |
   +----------------+---------------------------------+----------------------------------+----------------------------------+-------------------------------------------------+
   | parameters     | icepack\_init\_parameters       | icepack\_query\_parameters       | icepack\_write\_parameters       | scientific parameters                           |
   +----------------+---------------------------------+----------------------------------+----------------------------------+-------------------------------------------------+
   | tracer flags   | icepack\_init\_tracer\_flags    | icepack\_query\_tracer\_flags    | icepack\_write\_tracer\_flags    | tracer flags                                    |
   +----------------+---------------------------------+----------------------------------+----------------------------------+-------------------------------------------------+
   | tracer sizes   |                                 | icepack\_query\_tracer\_sizes    | icepack\_write\_tracer\_sizes    | tracer counts and tracer maximum sizes          |
   +----------------+---------------------------------+----------------------------------+----------------------------------+-------------------------------------------------+
   | tracer indices | icepack\_init\_tracer\_indices  | icepack\_query\_tracer\_indices  | icepack\_write\_tracer\_indices  | tracer indexing in a broader tracer array       |
   +----------------+---------------------------------+----------------------------------+----------------------------------+-------------------------------------------------+

Many of these interfaces are related to tracers and in particular, tracer indexing in broader arrays.  This is further explained in :ref:`tracerindex`.

.. _tracerindex:

Tracer Indexing
~~~~~~~~~~~~~~~~~~~

Tracers are really just variables associated with the model state.  Some of the tracers are
prognostic, vary each timestep, and are updated in Icepack.  Other tracers are just used by
Icepack to evolve other tracers.

One of the most complicated aspects of the Icepack usage are managing tracers.  Some tracers (i.e.
Tsfc, qice, qsno) are required, while other tracers (i.e. FY or bgc tracers) are optional and used 
only when certain features are triggered.  As a general rule, Icepack is aware of only a specific set
of tracers and each tracer takes on multiple properties including counts, dependencies (:ref:`pondtr`), 
and indexing in a broader tracer array.  The following table summarize the various types of 
tracers understood by Icepack and lists some of their properties.  See also :ref:`tab-bio-tracer`.

.. table:: *Tracer Types and Properties* 

   +------------+----------+---------------+---------+---------+-----------------------------------------------------------------------------------+
   | name       | status   | optional flag | number  | count   | notes                                                                             |
   +============+==========+===============+=========+=========+===================================================================================+
   | Tsfc       | required |               | 1       | 1       | ice/snow temperature                                                              |
   +------------+----------+---------------+---------+---------+-----------------------------------------------------------------------------------+
   | qice       | required |               | 1       | nilyr   | ice enthalpy                                                                      |
   +------------+----------+---------------+---------+---------+-----------------------------------------------------------------------------------+
   | qsno       | required |               | 1       | nslyr   | snow enthalpy                                                                     |
   +------------+----------+---------------+---------+---------+-----------------------------------------------------------------------------------+
   | sice       | required |               | 1       | nilyr   | ice bulk salinity                                                                 |
   +------------+----------+---------------+---------+---------+-----------------------------------------------------------------------------------+
   | iage       | optional | tr_iage       | 1       | 1       | ice age                                                                           |
   +------------+----------+---------------+---------+---------+-----------------------------------------------------------------------------------+
   | FY         | optional | tr_FY         | 1       | 1       | first year ice                                                                    |
   +------------+----------+---------------+---------+---------+-----------------------------------------------------------------------------------+
   | alvl       | optional | tr_lvl        | 1       | 1       | level ice area fraction                                                           |
   +------------+----------+---------------+---------+---------+-----------------------------------------------------------------------------------+
   | vlvl       | optional | tr_lvl        | 1       | 1       | level ice area volume                                                             |
   +------------+----------+---------------+---------+---------+-----------------------------------------------------------------------------------+
   | apnd       | optional | tr_pond       | 1       | 1       | melt pond area fraction                                                           |
   +------------+----------+---------------+---------+---------+-----------------------------------------------------------------------------------+
   | hpnd       | optional | tr_pond       | 1       | 1       | melt pond depth                                                                   |
   +------------+----------+---------------+---------+---------+-----------------------------------------------------------------------------------+
   | ipnd       | optional | tr_pond       | 1       | 1       | melt pond refrozen thickness                                                      |
   +------------+----------+---------------+---------+---------+-----------------------------------------------------------------------------------+
   | fsd        | optional | tr_fsd        | 1       | nfsd    | floe size distribution                                                            |
   +------------+----------+---------------+---------+---------+-----------------------------------------------------------------------------------+
   | iso        | optional | tr_iso        | n_iso   | 2       | water isotopes (snow, sea ice)                                                    |
   +------------+----------+---------------+---------+---------+-----------------------------------------------------------------------------------+
   | aero       | optional | tr_aero       | n_aero* | 4       | aerosols (snow SSL, snow below SSL, sea ice SSL, sea ice below SSL in that order) |
   +------------+----------+---------------+---------+---------+-----------------------------------------------------------------------------------+
   | fbri       | optional | tr_brine      | 1       | 1       |                                                                                   |
   +------------+----------+---------------+---------+---------+-----------------------------------------------------------------------------------+
   | bgc_S      | optional |               | 1       | nblyr   | bulk salinity in fraction ice                                                     |
   +------------+----------+---------------+---------+---------+-----------------------------------------------------------------------------------+
   | bgc_N      | optional | tr_bgc_N      | n_algae | nblyr+3 | nutrients                                                                         |
   +------------+----------+---------------+---------+---------+-----------------------------------------------------------------------------------+
   | bgc_Nit    | optional |               | 1       | nblyr+3 | diatoms, phaeocystis, pico/small                                                  |
   +------------+----------+---------------+---------+---------+-----------------------------------------------------------------------------------+
   | bgc_DOC    | optional | tr_bgc_DOC    | n_doc   | nblyr+3 | dissolved organic carbon                                                          |
   +------------+----------+---------------+---------+---------+-----------------------------------------------------------------------------------+
   | bgc_DIC    | optional |               | n_dic   | nblyr+3 | dissolved inorganic carbon                                                        |
   +------------+----------+---------------+---------+---------+-----------------------------------------------------------------------------------+
   | bgc_chl    | optional |               | n_algae | nblyr+3 | algal chlorophyll                                                                 |
   +------------+----------+---------------+---------+---------+-----------------------------------------------------------------------------------+
   | bgc_Am     | optional | tr_bgc_Am     | 1       | nblyr+3 | ammonia                                                                           |
   +------------+----------+---------------+---------+---------+-----------------------------------------------------------------------------------+
   | bgc_Sil    | optional | tr_bgc_Sil    | 1       | nblyr+3 | silicon                                                                           |
   +------------+----------+---------------+---------+---------+-----------------------------------------------------------------------------------+
   | bgc_DMSPp  | optional | tr_bgc_DMS    | 1       | nblyr+3 |                                                                                   |
   +------------+----------+---------------+---------+---------+-----------------------------------------------------------------------------------+
   | bgc_DMSPd  | optional | tr_bgc_DMS    | 1       | nblyr+3 |                                                                                   |
   +------------+----------+---------------+---------+---------+-----------------------------------------------------------------------------------+
   | bgc_DMS    | optional | tr_bgc_DMS    | 1       | nblyr+3 |                                                                                   |
   +------------+----------+---------------+---------+---------+-----------------------------------------------------------------------------------+
   | bgc_PON    | optional | tr_bgc_PON    | 1       | nblyr+3 | zooplankton and detritus                                                          |
   +------------+----------+---------------+---------+---------+-----------------------------------------------------------------------------------+
   | bgc_DON    | optional | tr_bgc_DON    | n_don   | nblyr+3 | dissolved organic nitrogen                                                        |
   +------------+----------+---------------+---------+---------+-----------------------------------------------------------------------------------+
   | bgc_Fed    | optional | tr_bgc_Fe     | n_fed   | nblyr+3 | dissolved iron                                                                    |
   +------------+----------+---------------+---------+---------+-----------------------------------------------------------------------------------+
   | bgc_Fep    | optional | tr_bgc_Fe     | n_fep   | nblyr+3 | particulate iron                                                                  |
   +------------+----------+---------------+---------+---------+-----------------------------------------------------------------------------------+
   | bgc_hum    | optional | tr_bgc_hum    | 1       | nblyr+3 | humic material                                                                    |
   +------------+----------+---------------+---------+---------+-----------------------------------------------------------------------------------+
   | zaero      | optional | tr_zaero      | n_zaero | nblyr+3 | bgc aerosols like black carbon                                                    |
   +------------+----------+---------------+---------+---------+-----------------------------------------------------------------------------------+
   | zbgc_frac  | optional |               | 1       | nbtrcr  | fraction of tracer in mobile phase                                                |
   +------------+----------+---------------+---------+---------+-----------------------------------------------------------------------------------+

* NOTE the aero tracer indexing is a little more complicated depending which aero option is chosen.

The nt\_ start index in a full tracer array is the start index associated with tracer
relative to the number*count.  The nlt\_ start index in a bgc array is the start index 
associated with the tracer relative to the number only and it generally contains only
bgc tracers.

Generally, tracers are passed into the Icepack interfaces by type where each type is a separate
argument.  There are some cases where an array of tracers is required and this is where the
tracer indexing is particularly important.  Below is a list of the various tracer indexing used

  - nt\_ references the tracer start index in a broader tracer array
  - nlt\_ references a bgc specific tracer start index for a different bgc array with different indexing from the nt\_ indexing
  - trcrn_depend/strata/etc defines dependency properties for tracers associated with the full array reference by nt\_ indexing
  - bio_index and bio_index_o is something else

In **icepack_aggregate**, the arguments
*trcr_depend*, *trcr_base*, *n_trcr_strata*, and *nt_strata* are passed into the interface, and they
provide information on dependencies between tracers.  This information needs to be initialized in
the driving code.  In the bgc implementation, there are arrays *bio_index* and *bio_index_o* which
also need to be initialized in the driving code and passed to Icepack.

