:tocdepth: 3

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

Icepack does not generally contain any parallelization or I/O.  The driver of Icepack is 
expected to support
those features.  Icepack can be called concurrently across multiple MPI tasks.  Icepack should also
be thread safe.

