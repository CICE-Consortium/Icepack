:tocdepth: 3 

.. _addtrcr:

Adding tracers
====================

We require that any changes made to the code be implemented in such a way that they can
be "turned off" through namelist flags.  In most cases, code run with such changes should 
be bit-for-bit identical with the unmodified code.  Occasionally, non-bit-for-bit changes
are necessary or unavoidable due to a change in the order of operations. In
these cases, changes should be made in stages to isolate the non-bit-for-bit changes, 
so that those that should be bit-for-bit can be tested separately.

Tracers added to Icepack will also require extensive modifications to the host
sea ice model, including initialization on the horizontal grid, namelist flags 
and restart capabilities.  Modifications to the Icepack driver should reflect
the modifications needed in the host model but are not expected to match completely.
We recommend that a logical namelist variable
``tr_[tracer]`` be added and used for all calls involving the new tracer to provide
clean backwards compatibility.

A number of optional tracers are available in the code, including ice
age, first-year ice area, melt pond area and volume, brine height,
aerosols, and level ice area and volume (from which ridged ice
quantities are derived). Salinity, enthalpies, age, aerosols, level-ice
volume, brine height and most melt pond quantities are volume-weighted
tracers, while first-year area, pond area, and level-ice area are area-weighted 
tracers. Biogeochemistry tracers in the skeletal layer are area-weighted,
and vertical biogeochemistry tracers are volume-weighted.  In
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
the tracer dependencies (weights), which are tracked using the arrays
``trcr_depend`` (indicates dependency on area, ice volume or snow volume),
``trcr_base`` (a dependency mask), ``n_trcr_strata`` (the number of
underlying tracer layers), and ``nt_strata`` (indices of underlying layers). 
See subroutine *icepack_compute_tracers* in **icepack_tracers.F90**.

To add a tracer, follow these steps using one of the existing tracers, for example 
age or isotopes, as an example.  Lets call the new tracer xyz.  Note that many
of the changes are defined in the Icepack driver or scripts.  For CICE or other models
using Icepack, adding a new tracer may be done differently.  However, changes to the
Icepack columnphyics should be similar.  As you make changes, we recommend that you
build and run with the tracer off to ensure the code modifications are working properly.
Changes to setup, namelist, and build files either need to be tested by generating
a new case or by modifying the same files in an existing case before running.
Code changes can be tested by rebuilding and rerunning an existing case.

The appendix has a tutorial activity, :ref:`tutorialfluff`, where a tracer is added to 
Icepack.  A set of code differences are provided there as an example of how this could
be done in practice.

#. Define a new tracer.  First, setup the tracer cpp in **configuration/scripts/icepack.settings**
   and **configuration/scripts/icepack.build**

   - In **icepack.settings**, add::

        setenv NXYZ 1  # number of xyz tracers

   - In **icepack.build**, add::

        -DNXYZ=${NXYZ}

     to the setenv ICE_CPPDEFS line

#. Define the new tracer dimension in **icedrv_domain_size.F90**.  

   - Add a new dimension::

        n_xyz = NXYZ, & ! number of xyz tracers

   - Update the size of ``max_ntrcr`` by adding::

        + n_xyz  & ! number of xyz tracers

#. Add the new tracer to **icepack_tracers.F90**: 

   - Define new variables::

        tr_xyz=.false.
        n_xyz=0
        nt_xyz=0

     By default, these are turned off.

   - Add the new variables to the tracer init, query, and write subroutine arguments
     (tracer_flags, tracer_sizes, and tracer_indices).  The driver of Icepack will turn
     this tracer on as needed.

#. Update the driver code to initialize the tracer flag and size and allocate space
   for the tracer.  You will also
   define dependency, whether the variable is dependent on ice area (aice), ice
   volume (vice), snow volume (vsno), or something else.  Edit **icedrv_init.F90**.

   - Add ``n_xyz`` to use icedrv_domain_size in subroutine input_data

   - Define ``tr_xyz`` and ``nt_xyz`` variables in subroutine input_data

   - Add a logical namelist variable ``tr_xyz` to tracer_nml

   - Add a line to initialize ``tr_xyz=.false.``

   - Add a check that ``tr_xyz`` and ``n_xyz`` are consistent

   - Add code to increment the number of tracers in use, ``ntrcr``, based on namelist input

   - Add some code to print xyz tracer info to the output file

   - Add calls to subroutines icepack_init_tracer_sizes, icepack_init_tracer_flags, and 
     icepack_init_tracer_indices to initialize the xyz tracer information in Icepack.

   - Define the tracer dependency in subroutine init_state.  Area tracers can be thought of in terms of averages across 
     the ice area (the tracer value itself) or averages over the grid cell area (tracer * aice).  
     Volume tracers can be considered in a similar way, but it can also be helpful to think of 
     the tracer as the “density” of the variable, which when multiplied by the ice or snow volume 
     gives the total content of the variable in that volume of the ice or snow.  The “total content” 
     is often a conserved quantity that can be checked during the calculation.

#. Add the tracer to the default namelist **configuration/scripts/icepack_in**

   - Add the namelist flag ``tr_xyz`` to *tracer_nml*.
     Best practice is to set the default namelist values so that the 
     new capability is turned off.  You can create an option file with your preferred
     configuration in **configuration/scripts/options**.

#. If your tracer depends on ocean or atmosphere forcing

   - Initialize the sources and sinks in **icedrv_flux.F90**

   - Add a subroutine to generate these sources or sinks in **icedrv_forcing.F90** 
     or **icedrv_forcing_bgc.F90**.

   - Add the new tracer forcing calls in **icedrv_InitMod.F90** and **icedrv_RunMod.F90**

   - Always use the flag ``tr_xyz`` to determine whether to call these routines.

#. Create a new physics file/module for your tracer, **icepack_xyz.F90**.
   This subroutine handles a number of the processes within the 
   sea ice and snow which affect the distribution of the tracer. These processes could include 
   congelation growth, snow melt, sea ice melt, evaporation / sublimation, or snow-ice formation. 
   Note there is a more sophisticated vertical distribution of tracer under the zbgc package. 
   If your tracer is impacted by frazil ice growth or lateral melt, this is discussed later.
   It’s often helpful to copy and modify existing modules such as icepack_age.F90 or icepack_isotope.F90.

#. Add the physics calls to Icepack or the driver.  

   - Depending on the physics implementation, the
     new tracer physics calls might be done in **icepack_therm_vertical**, **icedrv_step.F90**, and/or
     elsewhere.  See use of subroutines ``update_aerosol``, ``update_isotope``, or ``increment_age``.

   - Pass tracer array into Icepack via public interfaces as needed
   
   - Always use the flag ``tr_xyz`` to determine whether to call these routines.

#. Add the new tracer to the restart files.  Edit **icedrv_restart.F90**,

   -  define restart variables

   -  call routines to read and write tracer restart data

#. If strict conservation is necessary, add conservation diagnostics using the 
   topographical ponds as an example, :ref:`ponds`

#. Update documentation, including **icepack_index.rst** and **ug_case_settings.rst**

#. Test and validate.  Verify backwards compatibility.

