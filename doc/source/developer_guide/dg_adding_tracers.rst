:tocdepth: 3 

.. _addtrcr:

Adding tracers
====================

We require that any changes made to the code be implemented in such a way that they can
be "turned off" through namelist flags.  In most cases, code run with such changes should 
be bit-for-bit identical with the unmodified code.  Occasionally, non-bit-for-bit changes
are necessary, e.g. associated with an unavoidable change in the order of operations. In
these cases, changes should be made in stages to isolate the non-bit-for-bit changes, 
so that those that should be bit-for-bit can be tested separately.

Tracers added to Icepack will also require extensive modifications to the host
sea ice model, including initialization on the horizontal grid, namelist flags 
and restart capabilities.  Modifications to the Icepack driver should reflect
the modifications needed in the host model but are not expected to match completely.
We recommend that the logical namelist variable
``tr_[tracer]`` be used for all calls involving the new tracer outside of
**ice\_[tracer].F90**, in case other users do not want to use that
tracer.

A number of optional tracers are available in the code, including ice
age, first-year ice area, melt pond area and volume, brine height,
aerosols, and level ice area and volume (from which ridged ice
quantities are derived). Salinity, enthalpies, age, aerosols, level-ice
volume, brine height and most melt pond quantities are volume-weighted
tracers, while first-year area, pond area, level-ice area and all of the
biogeochemistry tracers in this release are area-weighted tracers. In
the absence of sources and sinks, the total mass of a volume-weighted
tracer such as aerosol (kg) is conserved under transport in horizontal
and thickness space (the mass in a given grid cell will change), whereas
the aerosol concentration (kg/m) is unchanged following the motion, and
in particular, the concentration is unchanged when there is surface or
basal melting. The proper units for a volume-weighted mass tracer in the
tracer array are kg/m.

In several places in the code, tracer computations must be performed on
the conserved "tracer volume" rather than the tracer itself; for
example, the conserved quantity is :math:`h_{pnd}a_{pnd}a_{lvl}a_{i}`,
not :math:`h_{pnd}`. Conserved quantities are thus computed according to
the tracer dependencies, and code must be included to account for new
dependencies (e.g., :math:`a_{lvl}` and :math:`a_{pnd}` in
**ice\_itd.F90** and **ice\_mechred.F90**).

To add a tracer, follow these steps using one of the existing tracers as
a pattern.

#. **icedrv\_domain\_size.F90**: increase ``max_ntrcr`` (can also add option
   to **icepack.settings** and **icepack.build**)

#. **icedrv\_state.F90**: declare ``nt_[tracer]`` and ``tr_[tracer]``

#. **icepack\_[tracer].F90**: create initialization, physics routines

#. **ice\_drv\_init.F90**: (some of this may be done in **ice\_[tracer].F90**
   instead)

   -  add new module and ``tr_[tracer]`` to list of used modules and
      variables

   -  add logical namelist variable ``tr_[tracer]``

   -  initialize namelist variable

   -  print namelist variable to diagnostic output file

   -  increment number of tracers in use based on namelist input (``ntrcr``)

   -  define tracer types (``trcr_depend`` = 0 for ice area tracers, 1 for
      ice volume, 2 for snow volume, 2+``nt_``[tracer] for dependence on
      other tracers)

#. **icepack\_itd.F90**, **icepack\_mechred.F90**: Account for new dependencies
   if needed.

#. **icedrv\_InitMod.F90**: initialize tracer (includes reading restart
   file)

#. **icedrv\_RunMod.F90**, **icedrv\_step\_mod.F90**:

   -  call routine to write tracer restart data

   -  call physics routines in **icepack\_[tracer].F90** (often called from
      **icedrv\_step\_mod.F90**)

#. **icedrv\_restart.F90**: define restart variables

#. **icepack\_in**: add namelist variables to *tracer\_nml* and
   *icefields\_nml*

#. If strict conservation is necessary, add diagnostics as noted for
   topo ponds in Section :ref:`ponds`.
