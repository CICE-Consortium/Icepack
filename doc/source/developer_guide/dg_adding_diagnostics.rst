:tocdepth: 3 

.. _adddiag:

Adding diagnostics
==================

Icepack only produces ASCII (text) log output for four points (full model with ITD, an initially ice free point, a land point, and a case with a slab ocean) designated by the variable ``n`` here. Each of these files contains the state information for that point. Sometimes additional variables are required in this output. The procedure for adding diagnostic variables is outlined here.

#. For non-BGC variables, one should edit **icedrv\_diagnostics.F90**:

   -  If the variable is already defined within the code, then add it to a "use" statement in the subroutine
      ``runtime_diags``.

   -  Note that if the variable is not readily accessible through a use statement, then a global variable needs to
      be defined. This might be in **icedrv\_state.F90** or **icedrv\_flux.F90** for example.

   -  Additionally, if the variable is a derived quantity, then one needs to include the variables from a use statement
      to calculate the new quantity. For example, see how ``hiavg`` and ``hsavg`` are computed.

   -  If the variable is just a scalar, then follow the example of "aice". Copy the write statement for Qa to
      a place in the output list where it is most appropriate. The format "900" is appropriate for most scalars.
      Edit the copied statement to be the variable you want. The following example adds snow-melt (melts).

    .. code-block:: fortran

       use icedrv_flux, only: melts

       write(nu_diag_out+n-1,900) 'snow melt             = ',melts(n)! snow melt

   -  If the variable is an array, say depending on ncat, then follow the example of fiso_evap. This just requires
      adding a loop for the print statement. Make sure ncat and a counter, nc are available. Say for example, 
      the category ice area, aicen. The variable ncat is the number of subgridscale categories.

    .. code-block:: fortran

       use icedrv_domain_size, only: ncat

       use icedrv_state, only: aicen

       ! local variables

       integer (kind=int_kind) :: &
          n, nc, k

       do nc = 1,ncat
          write(nu_diag_out+n-1,901) 'Category ice area =       ',aicen(n,nc),nc ! category ice area
       enddo

   -  If the variable is a tracer, then in addition to the variable trcr or trcrn, you will need to have the tracer
      index available. Here, you can look at the example of nt_Tsfc. 

   -  In some cases, a new format statement might be required if 900 or 901 are not correct.

#. For BGC variables, one should edit **icedrv\_diagnostics\_bgc.F90**:

   -  If the variable is already defined within the code, then add it to a "use" statement in the subroutine
      ``hbrine_diags`` or ``bgc_diags`` or ``zsal_diags``. The similar procedure for state variables is used here.

   -  Note that the BGC needs to be activated and and the particular tracer turned on. 

In general, try to format the output statements to line up with the surrounding print messages. This may require a couple of tries to get it to compile and run.

