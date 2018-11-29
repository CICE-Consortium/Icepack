:tocdepth: 3

.. _about:

About Icepack
=============

Modern sea ice models have evolved into highly complex collections of physical parameterizations and
infrastructural elements to support various configurations and computational approaches.  In particular,
numerical models may now be implemented for unstructured grids, requiring new approaches for referencing
information in neighboring grid cells and communication information across grid elements.  However, a
large portion of the physics in sea ice models can be described in a vertical column, without reference
to neighboring grid cells.

The column physics package of the sea ice model CICE, "Icepack", is maintained by the
CICE Consortium. This code includes several options for simulating sea ice
thermodynamics, mechanical redistribution (ridging) and associated area and thickness
changes.  In addition, the model supports a number of tracers, including
thickness, enthalpy, ice age, first-year ice area, deformed ice area and volume,
melt ponds, and biogeochemistry.

Icepack is implemented in CICE as a git submodule.
The purpose of Icepack is to provide the column physics model as a separate
library for use in other host models such as CICE.
Development and testing of CICE and Icepack may be done together,
but the repositories are independent.
This document describes the Icepack model. The Icepack code is 
available from https://github.com/CICE-Consortium/Icepack.

Icepack consists of three independent parts, the column physics code,
the icepack driver that supports stand-alone testing of the column physics code, and the
icepack scripts that build and test the Icepack model.  
The column physics is called from a host (driver) model
on a gridpoint by gridpoint basis.  Each gridpoint is independent
and the host model stores and passes the model state and forcing to
the column physics.

Major changes with each Icepack release (https://github.com/CICE-Consortium/Icepack/releases) 
will be detailed with the included release notes. Enhancements and bug fixes made to 
Icepack since the last numbered release can be found on the Icepack wiki
(https://github.com/CICE-Consortium/Icepack/wiki/Icepack-Recent-changes).
**Please cite any use of the Icepack code.** More information can be found at :ref:`citing`. 

This document uses the following text conventions: Variable names used in 
the code are ``typewritten``. Subroutine names are given in *italic*. File 
and directory names are in **boldface**. A comprehensive :ref:`index`, 
including glossary of symbols with many of their values, appears at the 
end of this guide.