:tocdepth: 3

.. _adding:

Other things
=============


.. _debugger:

Running with a Debugger
-------------------------

Availability and usage of interactive debuggers varies across machines.  Contact your 
system administrator for additional information about what’s available on your system.  
To run with an interactive debugger, the following general steps should be taken.

- Setup a case
- Modify the env file and Macros file to add appropriate modules and compiler/ linker flags
- Build the model
- Get interactive hardware resources as needed
- Open a csh shell
- Source the env.${machine} file
- Source cice.settings
- Change directories to the run directory
- Manually launch the executable thru the debugger


