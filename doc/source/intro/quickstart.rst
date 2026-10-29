:tocdepth: 3


.. _quickstart:

Quick Start
===========

Download the model from the CICE-Consortium repository, 
    https://github.com/CICE-Consortium/Icepack

Instructions for working in github with Icepack (and CICE) can be
found in the `CICE Git and Workflow Guide <https://github.com/CICE-Consortium/About-Us/wiki/Git-Workflow-Guidance>`_.

On an HPC (high-performance computing) system that Icepack has already been ported to, from your main Icepack directory, execute::

  ./icepack.setup -c ~/mycase1 -m MACHINE -e ENVIRONMENT
  cd ~/mycase1
  ./icepack.build
  ./icepack.submit

where ``MACHINE`` is the machine name and ``ENVIRONMENT`` is the compiler environment,
for example ``-m derecho -e intel``.
The supported combinations are the **env.[machine]_[environment]** files in
**configuration/scripts/machines/**, and ``./icepack.setup --help`` lists the setup options.
See :ref:`scripts` for more information about how to use the icepack.setup and icepack.submit scripts.

To run Icepack on a laptop or personal computer, use the conda configuration
(``-m conda -e macos`` or ``-m conda -e linux``), see :ref:`laptops`.
To run on another HPC system, the model has to be ported first, see :ref:`porting`.
In both cases the input data must also be downloaded, see :ref:`force`.

Please cite any use of the Icepack code. More information can be found at :ref:`citing`.

