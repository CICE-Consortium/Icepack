:tocdepth: 3

.. _docintfc:

Public Interfaces
---------------------

Below are a list of public icepack interfaces.

The documentation for these interfaces is extracted directly from the icepack source code using the script
``doc/generate_interfaces.sh``.  That script updates the rst file ``interfaces.include`` in
the ``doc/source/user_guide directory``.  That file is part of the internal documentation.
There is information about how ``generate_interfaces.sh`` parses
the source code in a comment section in that script.  Executing ``icepack.setup --docintfc`` will 
run the generate_interfaces script as noted in :ref:`case_options`.  
Once ``generate_interfaces`` is executed, the user
still has to git add, commit, and push the changes to the documentation manually.  A typical workflow
would be::

    # verify all public interfaces in the columnphysics have appropriate autodocument comment line
    #  there should be a "!autodocument_start ${interface_name}" at the begining of the interface
    #  there should be a "!autodocument_end" at the end of the declaration of the interface arguments
    ./icepack.setup --docintfc
    git add doc/source/user_guide/interfaces.include
    git commit -m "update public interface documentation"

.. include:: interfaces.include

