:tocdepth: 3

.. _library:

Using Icepack in other models
=================================

This section documents how to use Icepack in other models.

.. _liboverview:

Overview
----------------

Icepack is a column physics package designed to be used in other broader sea ice models, such as
CICE, SIS, or even in ocean models.  
Icepack includes options for simulating sea ice thermodynamics, mechanical redistribution 
(ridging) and associated area and thickness changes. In addition, the model supports a number of 
tracers, including thickness, enthalpy, ice age, first-year ice area, deformed ice area and 
volume, melt ponds, and biogeochemistry.

Icepack is called on a grid point by grid point basis.  All data is passed in and out of the model
via subroutine interfaces.  Fortran "use" statements are not encouraged for accessing data inside
the Icepack model.

Icepack does not contain any parallelization or I/O.  The driver of Icepack is expected to support
those features.  Icepack can be called concurrently across multiple MPI tasks.  Icepack should also
be thread safe.

.. _initialization:

Icepack Initialization
----------------------

The subroutine icepack_configure should be called before any other icepack interfaces are called.
This subroutine initializes the abort flag and a few other important defaults.  We recommend that
call be implemented as::

      call icepack_configure()  ! initialize icepack
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call my_abort_method()

The 2nd and 3rd line above are described further in :ref:`aborts`.


.. _aborts:

Error Messages and Aborts
--------------------------------

Icepack does not understand the I/O (file units) or computing environment (MPI, etc).  It provides an
interface that allows the driver to write error messsages and check for an abort flag.  If Icepack
fails, it will make error messages available thru that interface and it will set an abort flag
that can be queried by the driver.
To best use those features, it's recommended that after every icepack interface call, the user
add the following::

      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call my_abort_method()

icepack_warnings_flush is a public interface in icepack that writes any warning or error messages
generated in icepack to the driver file unit number defined by nu_diag.  
The function icepack_warnings_aborted queries the internal icepack abort flag and
returns true if icepack generated an abort error.  
my_abort_method represents method that stops the driver model from 
running.  That interface or command is driver dependent.

.. _callingseq:

Calling Sequence
-----------------------

TBD

.. _interfaces:

Public Interfaces
---------------------

Below are a list of public icepack interfaces.

These interfaces are extracted directly from the icepack source code using the script
**doc/generate_interfaces.sh**.  There is documentation about how to use the script
in the script.

.. include:: interfaces.rst

