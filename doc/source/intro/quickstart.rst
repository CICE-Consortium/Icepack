:tocdepth: 3


.. _quickstart:

Quick Start
===========

Download the model from the CICE-Consortium repository, 
    https://github.com/CICE-Consortium/Icepack

Instructions for working in github with Icepack (and CICE) can be
found in the `CICE Git and Workflow Guide <https://github.com/CICE-Consortium/About-Us/wiki/Git-Workflow-Guidance>`_.

From your main Icepack directory, execute::

  ./icepack.setup -c ~/mycase1 -m testmachine
  cd ~/mycase1
  ./icepack.build
  ./icepack.submit

``testmachine`` is a generic machine name included with the icepack scripts.
The local machine name will have to be substituted for ``testmachine`` and
there are working ports for several different machines.  However, it may be necessary
to port the model to a new machine.  See :ref:`porting` for 
more information about how to port. See :ref:`scripts` for more information about 
how to use the icepack.setup and icepack.submit scripts.

Please cite any use of the Icepack code. More information can be found at :ref:`citing`.

