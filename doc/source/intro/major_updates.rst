:tocdepth: 3

.. _updates:


Major Icepack updates
============================================

Modern sea ice models have evolved into highly complex collections of physical parameterizations and
infrastructural elements to support various configurations and computational approaches.  In particular,
numerical models may now be implemented for unstructured grids, requiring new approaches for referencing
information in neighboring grid cells and communication information across grid elements.  However, a
large portion of the physics in sea ice models can be described in a vertical column, without reference
to neighboring grid cells.  This part of the CICE model has been separated into its own modular software
package, Icepack.  The column physics code was separated from CICE version 5.1.2 by removing all references to
the horizontal grid and other infrastructural CICE elements (e.g. MPI tasks, calendar).  

To allow the column physics to be developed and maintained as a software package independent of CICE,
a simplified driver was created along with a full test suite and scripts for building and running
the code.  Icepack includes the simplified driver and scripts for configuring various tests of the 
column physics code in columnphysics/.

Enhancements and bug fixes made to Icepack since the last numbered release can be found on the
Icepack wiki https://github.com/CICE-Consortium/Icepack/wiki/Recent-changes. Major changes with
each Icepack release (found here: https://github.com/CICE-Consortium/Icepack/releases) will
be included as release notes.
