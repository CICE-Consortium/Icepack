:tocdepth: 3 

.. _adddiag:

Adding diagnostics
==================

Icepack produces separate ASCII (text) log output for four cells, each with a different initial condition (full ITD, slab ice, ice free, land) designated by the variable ``n`` here. Each of the diagnostic files contains the state information for that cell. The procedure for adding diagnostic variables to the output is outlined here.

#. For non-BGC variables, edit **icedrv\_diagnostics.F90**:

   -  If the variable is already defined within the code, then add it to a "use" statement in the subroutine
      ``runtime_diags``.

   -  Note that if the variable is not readily accessible through a use statement, then a global variable may need to
      be defined. This might be in **icedrv\_state.F90** or **icedrv\_flux.F90** for example.

   -  Additionally, if the variable is a derived quantity, then the variables needed to calculate the new quantity
      may need to be added to a use statement. For example, see how ``hiavg`` and ``hsavg`` are computed.

   -  If the variable is a scalar, then follow the example of ``aice`` or ``hiavg``, copying the write statement to
      an appropriate place in the output list, and editing as needed. The format "900" is appropriate for most scalars.
      The following example adds snow melt (``melts``).

    .. code-block:: fortran

       use icedrv_flux, only: melts

       write(nu_diag_out+n-1,900) 'snow melt (m)         = ',melts(n) ! snow melt

   -  If the variable is an array, then you can compute the mean value (e.g. ``hiavg``) or print the array values (e.g. ``fiso_evap``).
      This may requires adding the array sizes and a counter for the loop(s). E.g. to print
      the category ice area, ``aicen`` over ncat thickness categories:

    .. code-block:: fortran

       use icedrv_domain_size, only: ncat

       use icedrv_state, only: aicen

       ! local variables

       integer (kind=int_kind) :: &
          n, nc, k

       do nc = 1,ncat
          write(nu_diag_out+n-1,901) 'Category ice area =       ',aicen(n,nc),nc ! category ice area
       enddo

   -  If the variable is a tracer, then in addition to the variable trcr or trcrn, you will need the tracer
      index (e.g. ``nt_Tsfc``). 

   -  In some cases, a new format statement might be needed.

#. For BGC variables, edit **icedrv\_diagnostics\_bgc.F90**:

   -  If the variable is already defined within the code, then add it to a "use" statement in the subroutine
      ``hbrine_diags`` or ``bgc_diags`` and follow a similar procedure for state variables as above.

   -  Note that the BGC needs to be activated and the particular tracer turned on. 

In general, try to format the output statements to line up with the surrounding print messages. This may require a couple of tries to get it to compile and run.

