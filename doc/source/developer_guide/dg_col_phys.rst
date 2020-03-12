:tocdepth: 3

.. _dev_colphys:

Icepack Column Physics
========================

File List
------------------------------------

The column physics source code contains the following files

| **columnphysics/**   the column physics code
|    **icepack_aerosol.F90**       handles most work associated with the aerosol tracers
|    **icepack_age.F90**           handles most work associated with the age tracer
|    **icepack_algae.F90**         biogeochemistry
|    **icepack_atmo.F90**          stability-based parameterization for calculation of turbulent ice–atmosphere fluxes
|    **icepack_brine.F90**         evolves the brine height tracer
|    **icepack_firstyear.F90**     handles most work associated with the first-year ice area tracer
|    **icepack_flux.F90**          fluxes needed/produced by the model
|    **icepack_fsd.F90**           supports floe size distribution
|    **icepack_intfc.F90**         interface routines for linking Icepack with a host sea ice model
|    **icepack_itd.F90**           utilities for managing ice thickness distribution
|    **icepack_kinds.F90**         basic definitions of reals, integers, etc.
|    **icepack_mechred.F90**       mechanical redistribution (ridging)
|    **icepack_meltpond_cesm.F90** CESM melt pond parameterization
|    **icepack_meltpond_lvl.F90**  level-ice melt pond parameterization
|    **icepack_meltpond_topo.F90** topo melt pond parameterization
|    **icepack_mushy_physics.F90** physics routines for mushy thermodynamics
|    **icepack_ocean.F90**         mixed layer ocean model
|    **icepack_orbital.F90**       orbital parameters for Delta-Eddington shortwave parameterization
|    **icepack_parameters.F90**    basic model parameters including physical and numerical constants requried for column package
|    **icepack_shortwave.F90**     shortwave and albedo parameterizations
|    **icepack_therm_0layer.F90**  zero-layer thermodynamics of :cite:`Semtner76`
|    **icepack_therm_bl99.F90**    multilayer thermodynamics of :cite:`Bitz99`
|    **icepack_therm_itd.F90**     thermodynamic changes mostly related to ice thickness distribution
|    **icepack_therm_mushy.F90**   mushy-theory thermodynamics of :cite:`Turner13`
|    **icepack_therm_shared.F90**  code shared by all thermodynamics parameterizations
|    **icepack_therm_vertical.F90**  vertical growth rates and fluxes
|    **icepack_tracers.F90**       tracer information
|    **icepack_warnings.F90**      utilities for writing warning and error messages
|    **icepack_wavefracspec.F90**  wave impact on sea ice
|    **icepack_zbgc.F90**          driver for ice biogeochemistry and brine tracer motion
|    **icepack_zbgc_shared.F90**   parameters and shared code for biogeochemistry and brine height
|    **icepack_zsalinity.F90**     vertical salinity parameterization of :cite:`Jeffery11`


Coding Standard
------------------------------------

The column physics is a library that solves the sea ice column physics on a 
timestep by timestep and gridpoint by gridpoint basis.  It consists of Fortran routines with 
input and output arguments.  The model state is saved in the host model.  There is no 
communication between gridcells so the underlying implementation
supports no parallelization.  It however can be called in parallel by a driver
that is running on multiple pes with a decomposed grid.

The column physics does not have a time manager.  Calendaring is expected to be
dealt with by the host model.  The column physics does not read any forcing data,
that is passed into the column physics though interfaces.  In fact, 
there are no direct IO capabilities in the column physics.  That is to say, the
column physics does not open files to read or write.  The column physics is able to write 
data via several different routines that specifically have a fortran unit number as an input
argument.  In addition, there is a warning and abort package (see section :ref:`aborts`) that
provides the column package with the ability to store log output.  That output can
be queried by the host model or it can be written directly via a specific routine.
The warning package also provides access to an abort flag that the host model can
query after each call to check for successful completion of the column physics package.

All column physics public interfaces and public data are defined in the **icepack_intfc.F90**
file (see section :ref:`calling`).  Internal column physics settings should all be accessible through interfaces.
The internal constants, parameters, and tracer settings have init (set), query (get), and
write interfaces that provides access to internal column physics settings.  The host model
should not have to use "use" statements to access any of the column physics data outside
of what is provided through the icepack_intfc module.  
The public column physics interfaces use optional arguments where it makes sense and
there is an ongoing effort to extend the optional arguments supported.  It's strongly recommended
that calls to the icepack interfaces be done with keyword=value arguments.  All icepack arguments
support this method.

Overall, columnphysics changes in the Icepack model should include the following

  * All modules should have the following set at the top

    .. code-block:: fortran

       implicit none
       private

  * Any public module interfaces or data should be explicitly specified

  * All subroutines and functions should define the subname character parameter statement to match the interface name like

    .. code-block:: fortran

       character(len=*),parameter :: subname='(lateral_melt_bgc)'

  * All interfaces that are public outside the Icepack columnphysics should include autodocument_start and autodocument_end comment lines with appropriate syntax and location.  If any interfaces are added or updated, then the internal documentation should be updated via

    .. code-block:: bash

       ./icepack.setup --docintfc

    See also :ref:`docintfc` for more information about the docintfc option.

  * The icepack_warnings package should be used to cache log messages and set the abort flag.  To add a log message, use icepack_warnings_add like

    .. code-block:: fortran

       call icepack_warnings_add(subname//' algorithm did not converge')

    To formally set the abort flag, use

    .. code-block:: fortran

       call icepack_warnings_setabort(.true.,__FILE__,__LINE__)

    See also :ref:`aborts` for more information about how the external calling program will write those message and check whether Icepack aborted.

  * Every interface call within the columnphysics should be followed by

    .. code-block:: fortran

       if (icepack_warnings_aborted(subname)) return

    to support errors backing up the call tree to the external program

  * Variables defined in icepack_kinds, icepack_tracers, icepack_parameters, and icepack_orbital should be accessed within Icepack by Fortran use statements.  It's also possible to access some of those variables thru methods that query for the value, but this tends to be a little more cumbersome, so Fortran use statements are recommended within columnphysics.  From the icepack driver or other external programs, the columnphysics variables should ALWAYS be access thru the interface methods and icepack_intfc (see also :ref:`calling`).

  * Optional arguments are encouraged in the public Icepack interfaces but should generally be avoided in interfaces within the columnphysics.  There are several reasons for taking this approach.  There is a desire to support backwards compatible Icepack public interfaces as much as possible, so optional arguments will be used for some future extensions.  There is also a desire to allow users to pass only the data thru the Icepack interfaces that is needed.  To support optional tracers and features, optional arguments are needed.  Within the internal columnphysics calling tree, optional arguments are discouraged because they tend to add complexity to deep calling trees and often lead to implementations with many calls to the same interface that only vary by which arguments are passed.  In the long term, that approach is not sustainable.  As a result, a scheme has been developed to support optional arguments in the public interfaces while minimizing optional arguments within the columphysics.  Within the columnphysics, we suggest optional arguments available thru the public interfaces should generally be treated as follows

    * Check whether optional arguments are passed and create temporary data to store the values

    * The temporary data should be locally name l_${argument_name}

    * The temporary data should be allocated at runtime if it's not a scalar based on the size of the incoming argument

    * The optional argument values should be copied into the temporary data

    * The temporary data should be passed thru other columnphysics subroutines

    * The temporary data should be deallocated at the end of the method if it was allocated

    * The temporary data should be copied back to the argument if the argument intent is out or inout

    * If optional arguments are not passed, temporary data should be created of size 1 with values of c0, and they should be passed thru other columnphysics subroutines

    * A logical can be instantiated and passed down the columnphysics interface to manage any logic related to whether valid or fake data is being passed down the calling tree.  See **closing_flag** and **iso_flag** within the columnphysics as examples.  There may also be externally set logicals that can be used to control how the optional features are handles.  See **tr_iso** within the columnphysics as an example.

    * An example of how this might look is

      .. code-block:: fortran

         subroutine icepack_example_interface(arg1, arg2, ...)
         real (kind=dbl_kind), intent(inout) :: arg1
         real (kind=dbl_kind), optional, dimension(:), intent(inout) :: arg2
         !
         real (kind=dbl_kind), allocatable, dimension(:) :: l_arg2
         logical :: arg2_flag

         character(len=*), parameter :: subname = '(icepack_example_interface)'

         if (present(arg2)) then
            arg2_flag = .true.
            allocate(l_arg2(size(arg2)))
            l_arg2 = arg2
         else
            arg2_flag = .false.
            allocate(l_arg2(1))
            l_arg2 = c0
         endif

         ...

         call some_columnphysics_subroutine(arg1, l_arg2, arg2_flag, ...)

         ...

         if (present(arg2)) then
            arg2 = l_arg2
         endif
         deallocate(l_arg2)

         return
         end subroutine

