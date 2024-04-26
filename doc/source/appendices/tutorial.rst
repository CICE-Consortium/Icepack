:tocdepth: 3

.. _tutorial:

Icepack Tutorial
=================



Learning Goals
----------------

In this activity you will clone the Icepack model code from the Consortium GitHub repository to run standalone Icepack simulations. You will also make namelist changes and code modifications for experiments and make some basic plots. If you run into issues, contact dbailey@ucar.edu.

Notes:

* Command line text is shown in highlighted boxes.
* When there is the <X> syntax, you need to fill in your personal information (e.g. a URL or username) for that command but without the angle brackets. Your GitHub and local computer usernames may not be the same, so check which you need to use.
* There is a lot of documentation.

  * Icepack User Guide, https://cice-consortium-icepack.readthedocs.io/en/latest/index.html
  * CICE-Consortium GitHub Usage Guide, https://github.com/CICE-Consortium/About-Us/wiki/Git-Workflow-Guide
  * CICE and Icepack Resources, https://github.com/CICE-Consortium/About-Us/wiki/Resource-Index


Github One-time Configuration
----------------------------------

You need to have your own GitHub account before you can start the following activities, and you should have already forked the Icepack repository.
For information about how to set up a GitHub account for the Icepack repository, see the Consortium documentation here, https://github.com/CICE-Consortium/About-Us/wiki/Git-Workflow-Guide.  The Consortium recommends that you keep your fork’s main branch in sync with the Consortium version and that you always work on branches.  This is all documented in the Git-Workflow-Guide linked above. 

Note: 

* The workflow guide is oriented toward setting up CICE rather than Icepack, but the same workflow applies to Icepack standalone.  Icepack can be set up and run as an independent model following the same workflow.


Clone Icepack
-------------------

Clone your Icepack repository fork (use the URL from your fork) to a local sandbox::

  mkdir ~/icepack-dirs
  cd ~/icepack-dirs
  git clone https://github.com/<github-user>/Icepack

If you have completed this correctly there should be an “Icepack” directory in the icepack-dirs directory. This is the “sandbox” we will be working in locally on your machine.

Move to the Icepack directory and check which branch you are using. This should be main::

  cd Icepack
  git status

Take a minute to orient yourself to the big picture structure of the directories and files in Icepack. The documentation has information about the Icepack :ref:`dirstructure`.

Make sure your main is up to date and create a branch. You can also update your fork directly in github by clicking the Sync fork button. If your code is already up to date, you can skip this step::

  git remote --v  (Check the origin and NO upstream)
  git remote add upstream https://github.com/CICE-Consortium/Icepack
  git remote --v  (Check upstream has been added)
  git pull upstream main
  git push origin main
  git branch <branchname>
  git checkout <branchname>


Conda and Laptop One-time Configuration
------------------------------------------

To build and run Icepack on your laptop, you need to install software via conda.  Instructions on how to do that can be found in the Icepack user guide, :ref:`laptops`.  If you have a Windows machine, we recommend using the Ubuntu Linux application, https://ubuntu.com/desktop/wsl.  Make sure to follow the instructions for installing miniconda. If your laptop has a conda environment already installed, you will still need to activate the icepack environment, and you may need to do so using the recommended miniconda distribution. Return here after completing section :ref:`laptops` in the documentation.  After installing miniconda, the main steps are::

  cd ~/icepack-dirs/Icepack
  conda env create -f configuration/scripts/machines/environment.yml
  conda activate icepack 

Before you can run Icepack, you have to set up a directory structure and download the input and forcing datasets::

  mkdir -p ~/icepack-dirs/runs ~/icepack-dirs/input ~/icepack-dirs/baseline
  cd ~/icepack-dirs/input
  curl -O https://zenodo.org/records/3728287/files/Icepack_data-20200326.tar.gz
  tar -xzf Icepack_data-20200326.tar.gz

You can also run Icepack on an external machine that is supported by the Consortium or to which you have ported the code. In this case, you do not need to port to your laptop.


Set Up an Icepack Simulation
-----------------------------

Use the online Icepack documentation and in particular the :ref:`quickstart` and :ref:`running_icepack` sections as guidance and for details on the command line settings::

  cd ~/icepack-dirs
  mkdir cases
  cd ~/Icepack
  ./icepack.setup --case ~/icepack-dirs/cases/icepack_test0 --mach <machine> --env <myenv> 

Notes:

* If you are doing this in the conda environment, the machine is “conda”.
* Similarly, the <myenv> variable is set to the compiler on your machine. For the conda environment, this is “macos” or “linux”.

The setup script creates a case consistent with the machine and other defined settings under ~/icepack-dirs/cases/ with the name you selected (icepack_test0). The case directory will contain build and run scripts, a namelist file, and other necessary files. Once the case is set up any of these files can be manually edited to refine the desired configuration.

Move to the new case directory and examine the settings::

  cd ~/icepack-dirs/cases/icepack_test0

Open the **icepack.settings** file and look at it briefly. Note the ICE_CASEDIR (it should match this directory) and the ICE_RUNDIR (where the model will be run and output created). Now look at the default namelist settings in **icepack_in**.

Build the code::

  ./icepack.build

The build script basically runs gmake under the covers, but there are a number of other tasks that are handled by the script to make the build more robust.  If the build is successful you will see the message “COMPILE SUCCESSFUL” at the bottom of the screen. You can also check the README.case file to check the status.

Submit the job. The submit script just submits the run scripts. Look at both **icepack.run** and **icepack.submit** files to see more details. The out-of-the-box run has default settings for the physics and other options. You can have a look at **icepack_in** and **icepack.settings** to review those settings. Then::

  ./icepack.submit

If the run is successful, you will see the message “ICEPACK COMPLETED SUCCESSFULLY” in the icepack run log file. Note that this job runs quickly - you are running a column model with four grid cells!

Look at the output!  Go to the ICE_RUNDIR where output was created. A successful model integration will create ice_diag.* files and a file in the “restart” directory called “iced.2016-01-01-00000”. The Icepack documentation has more information about :ref:`history`.

Follow the documentation to create some plots of the output using the tools provided with Icepack (:ref:`testplotting`). The conda icepack environment must be activated, if it isn’t already::
 
  cd ~/icepack-dirs/Icepack/configuration/scripts/tests/
  conda activate icepack
  ./timeseries.csh ~/icepack-dirs/runs/icepack_test0/ice_diag.full_ITD

Note that you can run the plotting script on any of the four ice_diag.* files.  The .png files are created in the ICE_RUNDIR directory. Open the files::

  cd ~/icepack-dirs/runs/icepack_test0/
  open <figurename>.png

Or use your file browser to navigate to the directory and double click on the images.

Questions to think about while looking at the output.

* What time period does an out-of-the-box run cover? 
* What are the differences between the full_ITD plots and the icefree plots (or any other combination of the ice_diag.* output files)? Which fields are the same? Which are different? Why would this be?
* What happens to ice area and ice thickness around October 1, 2015? Why do you see this signal?
* How does your output compare to the sample output provided for this release? (hint: see the wiki!)

Take a step back and think about all the directories and files you have created. The Icepack “sandbox” was cloned from GitHub and has the actual Icepack code.

* There is a particular case directory for building and launching the code, and some output (e.g. job log) are copied.
* There is a particular run directory for each case. This is where the model is run and big files are found.


Set Up a Longer Run
---------------------

Once you have had success with the previous step, you can run another, longer experiment to practice some basic changes for Icepack. Go back to your Icepack directory::

  cd ~/icepack-dirs/Icepack/

You need to set up a new out-of-the-box case (icepack_test1)::

  ./icepack.setup --case ~/icepack-dirs/cases/icepack_test1 --mach <machine> --env <myenv>

Go into the cases/icepack_test1 directory, and build the case.
Change the following namelist settings in **icepack_in**,

  npt = 8760

How long is this setting the model to run?  Change this to run for 10 years (hint: The timestep is one hour, and there are 24 steps per day, and 365 days per year).

Details about namelist options are in the documentation (:ref:`case_settings`).

Submit the job. Check the output and think about the following:

* Over what dates did the model run this time?
* What date would the model restart from?


Modify a physics option
---------------------------

Set up another case::

  ./icepack.setup --case ~/icepack-dirs/cases/icepack_test2 --mach <machine> --env <myenv>

Build the code.

Change the thermodynamics option from ktherm = 2 to ktherm = 1 in **icepack_in**, and set sw_redist = .true.  The intent here is to change the namelist option for the current experiment in the case directory.  Think about what would happen if you changed **icepack_in** in the source code before creating the case instead (hint: this experiment should work the same, but what about future experiments?).

Submit the job. Have a look at the output.

* What is different compared to your first run?
* What happens if sw_redist = .false. with ktherm = 1?  Why?


Change a Parameter in the Fortran Code
-----------------------------------------

Set up another case::

  ./icepack.setup --case ~/icepack-dirs/cases/icepack_test3 --mach <machine> --env <myenv>

Change to the source code directory::

  cd columnphysics

Edit icepack_mechred.F90 to change the line

  fsnowrdg = p5    , & ! snow fraction that survives in ridging

to

  fsnowrdg = c1    , & ! snow fraction that survives in ridging

Build the code and submit the job.

* What is different about this run?
* What do you think the fsnowrdg parameter is doing here?

Revert your code changes::

  cd ~/Icepack
  git status
  git checkout columnphysics/icepack_mechred.F90
  git status

