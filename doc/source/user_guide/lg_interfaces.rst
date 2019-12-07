:tocdepth: 3

.. _docintfc:

Public Interfaces
---------------------

Below are a list of public icepack interfaces.

These interfaces are extracted directly from the icepack source code using the script
``doc/generate_interfaces.sh``.  That script updates rst files in the
doc directory tree which are then incorporated into the sphinx documentation.
There is information about how ``generate_interfaces.sh`` parses
the source code in a comment section in that script.  In addition, 
executing ``icepack.setup --docintfc`` will also run the generate_interfaces 
script as noted in :ref:`case_options`.  
Once ``generate_interfaces`` is executed, the user
still has to add and commit the changes to the documentation manually.  A typical workflow
would be::

    ./icepack.setup --docintfc
    git add doc/source/user_guide/interfaces.rst
    git commit -m "update public interface documentation"

If the script is run, but no interfaces have changed, there should be no changes to the
documentation files.

.. include:: interfaces.include

