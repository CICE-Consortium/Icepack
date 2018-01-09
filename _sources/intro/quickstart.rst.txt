:tocdepth: 3


.. _quickstart:

Quick Start
===========

Download the model from the CICE-Consortium repository, 
    https://github.com/CICE-Consortium/Icepack

Instructions for working in github with Icepack (and CICE) can be
found in the `CICE Git and Workflow Guide <https://docs.google.com/document/d/1rR6WAvZQT9iAMUp-m_HZ06AUCCI19mguFialsMCYs9o>`_.

From your main Icepack directory, execute::

  ./icepack.create.case -c ~/mycase1 -m testmachine
  cd ~/mycase1
  ./icepack.build
  ./icepack.submit

``testmachine`` is a generic machine name included with the icepack scripts.
The local machine name will have to be substituted for ``testmachine`` and
there are working ports for several different machines.  However, it may be necessary
to port the model to a new machine.  See :ref:`porting` for 
more information about how to port and :ref:`scripts` for more information about 
how to use the icepack.create.case script.

