**********************
Introduction - Icepack
**********************

The column physics of the sea ice model CICE, "Icepack", is maintained by the
CICE Consortium. This code includes several options for simulating sea ice
thermodynamics, mechanical redistribution (ridging) and associated area and thickness
changes.  In addition, the model supports a number of tracers, including
thickness, enthalpy, ice age, first-year ice area, deformed ice area and volume,
melt ponds, and biogeochemistry.

Icepack is implemented in CICE as a git submodule.
Development and testing of CICE and Icepack may be done together,
but development within the repositories is independent.

This document uses the following text conventions:
Variable names used in the code are ``typewritten``.
Subroutine names are given in *italic*.
File and directory names are in **boldface**.
A comprehensive :ref:`index`, including glossary of symbols with many of their values, appears
at the end of this guide.

.. _quickstart:

Quick Start
===========

Download the model from the CICE-Consortium repository, 
    https://github.com/CICE-Consortium/Icepack

Instructions for working in github with Icepack (and CICE) can be
found in the `CICE Git and Workflow Guide <https://docs.google.com/document/d/1rR6WAvZQT9iAMUp-m_HZ06AUCCI19mguFialsMCYs9o>`_.

From your main Icepack directory, execute

``> ./icepack.create.case -c ~/mycase1 -m testmachine``

``> cd ~/mycase1``

``> ./icepack.build``

``> ./icepack.submit``

Note that it is necessary to have your computer set up for testmachine. Currently there are working ports for 
NCAR yellowstone and cheyenne, AFRL thunder, NavyDSRC gordon and conrad, and LANL's wolf machines. To add 
your own machine follow the instructions here (CHECK: LINK). 


Major Icepack updates since CICE v5.1.2
============================================

This model release is Icepack version 1.0.

The column physics code was separated from CICE version 5.1.2 by removing all references to
the horizontal grid and other infrastructural CICE elements (e.g. MPI tasks, calendar).  

- A simplified driver was developed for Icepack, for testing purposes. 
- Additional tests for the column physics are now available.
- This release includes the full vertical biogeochemistry code.

Acknowledgments and Copyright
=============================

Acknowledgements
----------------

This work has been completed through the CICE Consortium and its members with funding 
through the 
Department of Energy,
Department of Defense (Navy),
Department of Commerce (NOAA),
National Science Foundation
and Environment and Climate Change Canada.
Special thanks are due to the following people:

-  Elizabeth Hunke, Nicole Jeffery, Adrian Turner and Chris Newman at Los Alamos National Laboratory
 
-  David Bailey, Alice DuVivier and Marika Holland at the National Center for Atmospheric Research

-  Rick Allard and Matt Turner at the Naval Research Laboratory, Stennis Space Center,

-  Andrew Roberts of the Naval Postgraduate School,

-  Jean-Francois Lemieux and Frederic Dupont of Environment and Climate Change Canada,

-  Tony Craig and his supporters at the National Center for Atmospheric Research, the Naval Postgraduate School...CHECK,

-  Cecilia Bitz at the University of Washington, for her column forcing data,

-  and many others who contributed to previous versions of CICE.

Copyright
-----------
© Copyright 2017, Los Alamos National Security LLC. All rights reserved. 
This software was produced under U.S. Government contract 
DE-AC52-06NA25396 for Los Alamos National Laboratory (LANL), which is
operated by Los Alamos National Security, LLC for the U.S. Department
of Energy. The U.S. Government has rights to use, reproduce, and distribute
this software. NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC
MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE
OF THIS SOFTWARE. If software is modified to produce derivative works, such
modified software should be clearly marked, so as not to confuse it with the
version available from LANL. 

Additionally, redistribution and use in source and binary forms, with or
without modification, are permitted provided that the following conditions
are met:

- Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

- Neither the name of Los Alamos National Security, LLC, Los Alamos National Laboratory, LANL, the U.S. Government, nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT
NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL
SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


