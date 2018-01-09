:tocdepth: 3

.. _updates:


Major Icepack updates since CICE v5.1.2
============================================

This model release is Icepack version 1.0.

The column physics code was separated from CICE version 5.1.2 by removing all references to
the horizontal grid and other infrastructural CICE elements (e.g. MPI tasks, calendar).  

- A simplified driver was developed for Icepack, for testing purposes. 
- Additional tests for the column physics are now available.
- This release includes the full vertical biogeochemistry code.
- The ice velocity can optionally be included in the calculation of wind stress.
