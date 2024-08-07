:tocdepth: 3

.. _tutorial:

Icepack Tutorial
=================



Learning Goals
----------------

In this activity you will clone the Icepack model code from the Consortium GitHub repository, create a branch, add a new tracer to Icepack, and run standalone Icepack simulations. You will also make namelist changes and code modifications for experiments and make some basic plots.

Notes:

* Command line text is shown in highlighted boxes.
* When there is the <X> syntax, you need to fill in your personal information (e.g. a URL or username) for that command but without the angle brackets. Your GitHub and local computer usernames may not be the same, so check which you need to use.
* There is a lot of documentation.

  * Icepack User Guide, https://cice-consortium-icepack.readthedocs.io/en/latest/index.html
  * CICE-Consortium GitHub Usage Guide, https://github.com/CICE-Consortium/About-Us/wiki/Git-Workflow-Guide
  * CICE and Icepack Resources, https://github.com/CICE-Consortium/About-Us/wiki/Resource-Index


Create An Icepack Branch, Port, and Run Initial Case
------------------------------------------------------

You should fork the Icepack repository and create a new branch in your fork for development.  This is documented in https://github.com/CICE-Consortium/About-Us/wiki/Git-Workflow-Guide.  But basically,

* Create a github account for yourself
* Fork the Consortium Icepack repository using the "Fork" feature at https://github.com/CICE-Consortium/Icepack
* Clone your repository onto your machine
* Create a branch
* Checkout the branch
* Port the model and verify you can build and run the model

Many of these steps are done once per user or machine.  You may be able to leverage the conda port to build and run the model, you may be working on a supported machine, or you may have to port the model to your machine.  Porting, setting up cases, building, running, and running test suites are all documented in the Icepack User Guide, https://cice-consortium-icepack.readthedocs.io/en/latest/index.html.  There are many ways to setup and run the model.

**Note:** The workflow guide is oriented toward setting up CICE rather than Icepack, but the same workflow applies to Icepack standalone.

To summarize the steps in greater detail **assuming use of conda in a Mac or Linux environment**.

* Create a github account for yourself if you don't have one already (done once per user)

* Fork the Consortium Icepack repository, go to https://github.com/CICE-Consortium/Icepack and click on the fork button (done once per user)

* Clone Icepack, sync the main branch, and create a new branch from the main branch.  This will create a branch based on the lastest version of main::

    mkdir ~/icepack-dirs
    cd ~/icepack-dirs
    git clone https://github.com/<github-user>/Icepack
    cd Icepack
    git status        (Branch should be main)
    git remote --v    (Check the origin and NO upstream)
    git remote add upstream https://github.com/CICE-Consortium/Icepack
    git remote --v    (Check upstream has been added)
    git pull upstream main
    git push origin main
    git branch <branchname>
    git checkout <branchname>
    git status        (Branch should be <branchname>)

* Setup local env and download input datasets (done once per machine)::

    mkdir -p ~/icepack-dirs/runs ~/icepack-dirs/input ~/icepack-dirs/baseline
    cd ~/icepack-dirs/input
    curl -O https://zenodo.org/records/3728287/files/Icepack_data-20200326.tar.gz
    tar -xzf Icepack_data-20200326.tar.gz

* Setup the conda environment (done once per machine), see :ref:`laptops`::

    curl -L https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -o ~/Downloads/miniconda.sh
    bash ~/Downloads/miniconda.sh
    source $HOME/miniconda3/bin/activate
    conda init bash
    conda config --set auto_activate_base false
    source $HOME/miniconda3/bin/activate
    cd ~/icepack-dirs/Icepack
    conda env create -f configuration/scripts/machines/environment.yml
    conda activate icepack 

**NOTE:**  If you have a Windows machine, we recommend using the Ubuntu Linux application, https://ubuntu.com/desktop/wsl.

Again, there are many options for setting up the model on hardware, see the Icepack User Guide for more details.

* Set Up an Icepack Simulation.  See :ref:`quickstart` and :ref:`running_icepack`::

    cd ~/Icepack
    ./icepack.setup --case icepack_test0 --mach conda --env macos
    cd ~/icepack-dirs/cases/icepack_test0
    ./icepack.build
    ./icepack.submit

Several env variables are defined in **icepack.settings** and the Icepack namelist file is **icepack_in**.  Output files are copied from the run directory to a logs directory under the case.  If the run is successful, you will see the message “ICEPACK COMPLETED SUCCESSFULLY” in the icepack run log file. Note that this job runs quickly - you are running a column model with four grid cells!

Look at the output.  Go to the ICE_RUNDIR (defined in **icepack.settings**). A successful model integration will create ice_diag.* files and a file in the “restart” directory called “iced.yyyy-mm-dd-sssss” where yyyy-mm-dd-sssss is a model date/time stamp. The Icepack documentation has more information about :ref:`history`.

* Plot some output using the timeseries script provided, see :ref:`testplotting`. The conda icepack environment must be activated::

    cd $ICE_RUNDIR
    conda activate icepack
    ${ICE_SANDBOX}/configurations/scripts/tests/timeseries.csh ice_diag.full_ITD

**Note:** that you can run the plotting script on any of the four ice_diag.* files.  The .png files are created in the ICE_RUNDIR directory. View the png files.

* Questions to think about while looking at the output.

  * What time period does an out-of-the-box run cover? 
  * What are the differences between the full_ITD plots and the icefree plots (or any other combination of the ice_diag.* output files)? Which fields are the same? Which are different? Why would this be?
  * What happens to ice area and ice thickness around October 1, 2015? Why do you see this signal?
  * How does your output compare to the sample output provided for this release?


Modify the Configuration or Code
------------------------------------

* Set up a longer Run.  Modify ``npt`` in icepack_in.  ``npt`` defines the number of timesteps to run.  Details about namelist options are in the documentation (:ref:`case_settings`).

* Modify a physics option.  Change the thermodynamics option from ktherm=2 to ktherm=1 in **icepack_in**, and set sw_redist=.true.  The intent here is to change the namelist option for the current experiment in the case directory.  What is different compared to your first run?  What happens if sw_redist = .false. with ktherm = 1?  Why?

* Undo your latest **icepack_in** changes

* Change a Parameter in the Fortran Code.  Edit **icepack_mechred.F90** and set

    ``fsnowrdg = c1    , & ! snow fraction that survives in ridging``.  

  Rebuild the code before running.  What is different about this run?  What do you think the fsnowrdg parameter is doing here?

* Revert your latest code changes::

    cd ~/Icepack
    git status
    git checkout columnphysics/icepack_mechred.F90
    git status

.. _tutorialfluff:

Add a New Tracer
--------------------------------------

In this exercise, add a new tracer associated with fluffballs.
Call the tracer fluff and make it depend on ice area.  Follow the step-by-step instructions in the :ref:`addtrcr` documentation.
Once you have implemented the model changes, be sure to add fluffballs output to the standard output diagnostics and turn on the
fluff tracer.  Then update the timeseries plotting script to plot the fluffballs values over time.

* First, set the initial value, physics, sources, and sinks of fluff to zero and make sure fluff values remain zero throughout the run

* Add some constant atmospheric forcing and review results

* Change the dependency to ice volume, how do the results change?

* Modify the physics to create some physics processes, see isotopes or aerosols for some ideas

**NOTE:** The file, **doc/source/tutorial/fluff.diff** in the Icepack repository, demonstrates code differences for this fluffball activity as implemented in a version of Icepack from July, 2024.  These code differences may not be directly applicable to other code versions, but they provide an example of the typical code modifications required to add the tracer, fluff.
