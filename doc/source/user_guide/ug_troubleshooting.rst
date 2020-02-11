:tocdepth: 3

.. _troubleshooting:

Troubleshooting 
================

Check the FAQ: https://github.com/CICE-Consortium/Icepack/wiki

.. _setup:

Initial setup
-------------

If there are problems, you can manually edit 
the env, Macros, and **icepack.run** files in the case directory until things are 
working properly.  Then you can copy the env and Macros files back to 
**configuration/scripts/machines**.  

- Changes made directly in the run directory, e.g. to the namelist file, will be overwritten
  if scripts in the case directory are run again later.

- If changes are needed in the **icepack.run.setup.csh** script, it must be manually modified.

.. _restarttrouble:

Restarts
--------

- Manual restart tests require the path to the restart file be included in ``ice_in`` in the 
  namelist file.

- Ensure that ``kcatbound`` is the same as that used to create the restart file.  
  Other configuration parameters, such as ``NICELYR``, must also be consistent between runs.

.. _testtrouble:

..
      this is commented out now
    Underflows
    -----------
    - Tests using a debug flag that traps underflows will fail unless a "flush-to-zero" flag 
  is set in the Macros file.  This is due to very small exponential values in the delta-Eddington
      radiation scheme.

.. _debugging:

Debugging hints
---------------

Icepack has a warning package (**/columnphysics/icepack_warnings.F90**) where icepack 
stores information not set in write routines. Details about the package can be found 
in :ref:`aborts`. This package can be useful to detect an abort  

A printing utility is available in the driver that can be helpful when debugging the
code. Not all of these will work everywhere in the code, due to possible
conflicts in module dependencies.

*debug\_icepack* (**configuration/driver/ice\_diagnostics.F90**)
    A wrapper for *print\_state* that is easily called from numerous
    points during initialization and the timestepping loop

*print\_state* (**configuration/driver/ice\_diagnostics.F90**)
    Print the ice state and forcing fields for a given grid cell.

.. _bugs:

Known bugs and other issues
---------------------------

-   With the old CCSM radiative scheme (``shortwave`` = ‘default’ or
    ‘ccsm3’), a sizable fraction (more than 10%) of the total shortwave
    radiation is absorbed at the surface but should be penetrating into
    the ice interior instead. This is due to use of the aggregated,
    effective albedo rather than the bare ice albedo 
    when ``snowpatch`` < 1.

-   The linear remapping algorithm for thickness is not monotonic for tracers.

-   The form drag parameterization assumes a fixed ridge shape where both the 
    macroscopic ridge porosity and the angle of repose are specified parameters.  
    One high-resolution coupled model that uses the CICE column physics package
    has been unable to make this parameterization work in its current form.
    Development of a new, variational approach for ridging is underway 
    that will generate ridge shapes differently from
    the current parameterization, and is expected to alleviate the reported
    problem (Roberts, A., et al, in prep. 2018). 

Interpretation of albedos
-------------------------

The snow-and-ice albedo, ``albsni``, and diagnostic albedos ``albice``, ``albsno``,
and ``albpnd`` are merged over categories but not scaled (divided) by the
total ice area. (This is a change from CICE v4.1 for ``albsni``.) The latter
three history variables represent completely bare or completely snow- or
melt-pond-covered ice; that is, they do not take into account the snow
or melt pond fraction (``albsni`` does, as does the code itself during
thermodyamic computations). This is to facilitate comparison with
typical values in measurements or other albedo parameterizations. The
melt pond albedo ``albpnd`` is only computed for the Delta-Eddington
shortwave case.

With the Delta-Eddington parameterization, the albedo depends on the
cosine of the zenith angle (:math:`\cos\varphi`, ``coszen``) and is one if
the sun is below the horizon (:math:`\cos\varphi < 0`). Thus, the albedos
will be one in the dark, polar winter hemisphere. However, the
time-averaged albedo fields will be high if a diurnal solar cycle is
used, because values of one would be included in the average for half of
each 24-hour period. To rectify this, a separate counter should be used for the
averaging that is incremented only when :math:`\cos\varphi > 0`. However, this is
still a work in progress.

Interpretation of general results
---------------------------------

Icepack releases are "functional releases" in the sense that the code runs, 
does not crash, passes various tests, and requires further work to establish 
its scientific validity.  In general, users are not encouraged to use any of the
CICE Consortium's model configurations to obtain "scientific" results.  The
test configurations are useful for model development, but sea ice models must
be evaluated from a physical standpoint in a couple system because simplified
configurations do not necessarily represent what is actually happening in the
fully coupled system that includes interactive ocean and atmosphere components.


Proliferating subprocess parameterizations
------------------------------------------

With the addition of several alternative parameterizations for sea ice
processes, a number of subprocesses now appear in multiple parts of the
code with differing descriptions. For instance, sea ice porosity and
permeability, along with associated flushing and flooding, are
calculated separately for mushy thermodynamics, topo and level-ice melt
ponds, and for the brine height tracer, each employing its own
equations. Likewise, the Bitz99 and mushy thermodynamics compute freeboard
and snow–ice formation differently, and the topo and level-ice melt pond
schemes both allow fresh ice to grow atop melt ponds, using slightly
different formulations for Stefan freezing. These various process
parameterizations will be compared and their subprocess descriptions
possibly unified in the future.
