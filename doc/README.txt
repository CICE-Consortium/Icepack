THE ICEPACK DOCUMENTATION SYSTEM
================================

This README file contains a brief summary of how to add to or change documentation for Icepack. Detailed information on GitHub at:

https://github.com/CICE-Consortium/CICE/wiki/Working-with-Sphinx-documentation

That documentation includes an expanation of Sphinx, which we use to document CICE, as well as detailed information of how to checkout a version of Icepack, make changes to the documentation, and place a pull request in GitHub to commit your changes to the Icepack repository of the CICE-Consortium.

This README file details only the information you need to compile the documentation in your local version of Icepack.  

Steps for Modifying Documentation
=================================

Installing Sphinx:
------------------

This must be done once on each platform. See Sphinx or Installing Sphinx for details. On a mac laptop, execute the following on the command line:

$ sudo pip install --ignore-installed sphinx
$ sudo pip install --ignore-installed sphinxcontrib-bibtex

Other platforms may require other steps. The CICE Consortium uses the following software to get successful Sphinx HTML builds including references:

-python 2.7.11

-Sphinx (1.6.3)

-sphinx-rtd-theme (0.1.9)

-sphinxcontrib-bibtex (0.3.5)

-sphinxcontrib-websupport (1.0.1)

You will need to use the CICE Consortium's conf.py file, which is found under /doc/source/conf.py in the repository.

To use linked references within the HTML you will need to have the sphinxcontrib-bibtex package as well as the zreferences.rst and master_list.bib files located in /doc/source/ in the master repository. The list of references in master_list.bib is currently ordered sequentially from oldest to newest and alphabetically within a given year. To add references for your documentation, edit the master_list.bib file using the Articles and Books as examples for your addition(s). Please follow the format for ordering the date/alphabetization as well as including a URL with the document's DOI.

Model sandbox and documentation:
--------------------------------

Follow the general CICE-Consortium Git Workflow and Developer's guide to clone the repository and create your personal fork for model modifications. Whenever you modify the model you should update documentation. You can update the documentation on the same branch of your fork on which you test code, or you can create a separate branch called documentation to test only the RST and HTML documentation.

Editing RST files:
------------------

Open the RST file using a text editor and make the changes necessary. Note that from the User's Guide documentation (see link above) there is a hyperlink called "Show Source" on the left hand column that will show you the RST source code for the HTML you are viewing. This is a good way to see the syntax for tables, equations, linking references, labeling tables or figures, and correctly identifying documentation sections or subsections.

Move into the /doc/ directory of your sandbox. Then do the following:

 $ make clean 
This gets rid of old HTML files.

 $ make html
This builds HTML into /build/html/ directory. It will also give you errors if there is a problem with the build that will help you figure out how you need to modify your RST files for a successful HTML build.

 $ open /build/html/FILE.html 
Open the HTML on your browser for testing.

For further information, see https://github.com/CICE-Consortium/CICE/wiki/Working-with-Sphinx-documentation#installing-sphinx, including information on editing RST files.



