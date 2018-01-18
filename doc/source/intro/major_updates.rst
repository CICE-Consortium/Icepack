:tocdepth: 3

.. _updates:


Major Icepack updates since CICE v5.1.2
============================================

This model release is Icepack version 1.0.

The column physics code was separated from CICE version 5.1.2 by removing all references to
the horizontal grid and other infrastructural CICE elements (e.g. MPI tasks, calendar).  

To allow the column physics to be developed and maintained as a software package independent of CICE,
a simplified driver was created along with a full test suite and scripts for building and running
the code.

Additional enhancements:
- This release includes the full vertical biogeochemistry code, including particulate iron, humic material, modal aerosols, proteins, stationary and mobile phases
- Biogeochemistry can now feed back on sea ice physics through the shortwave formulation
- The ice velocity can optionally be included in the calculation of wind stress
- A warning package captures diagnostic and error information from within the column physics, for printing by the driver

Bug fixes:
- Use net shortwave sum instead of cosine of the zenith angle to limit shortwave calculation for low-/no-light conditions
- Correct roundoff errors in the delta-Eddington shortwave calculation
- Define interface temperature in the brine height parameterization
- Provide flexibility for tracking frazil ice salt/water in the ice or ocean model component
- When ponds are very thin, ignore them in the radiation calculation
- Ensure fractions of snow, ponds and bare ice sum to 1
- Wrap pond tracers in conditional blocks in mushy thermodynamics and form drag parameterization
- Miscellaneous bug fixes for biogeochemistry and brine height tracer (hbrine)

