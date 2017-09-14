
Numerical implementation
========================

Icepack is written in FORTRAN90 and runs on platforms using UNIX, LINUX,
and other operating systems. The code is parallelized via grid
decomposition with MPI or OpenMP threads and includes some optimizations
for vector architectures.

.. _dirstructure:

Directory structure
-------------------

Need to update for Icepack.

*old stuff*
The present code distribution includes make files, several scripts and
some input files. The main directory is **cice/**, and a run directory
(**rundir/**) is created upon initial execution of the script
**comp\_ice**. One year of atmospheric forcing data is also available
from the code distribution web site (see the **README** file for
details).

*incomplete former list*

 history and restart modules

**ice\_history\_write.F90**
    subroutines with  output

**ice\_restart.F90**
    read/write   restart files


Grid, boundary conditions and masks
-----------------------------------

*Needs update for Icepack*

.. _bio-grid:

Bio-grid
~~~~~~~~

The bio-grid is a vertical grid used for solving the brine height
variable :math:`h_b`. In the future, it will also be used for
discretizing the vertical transport equations of biogeochemical tracers.
The bio-grid is a non-dimensional vertical grid which takes the value
zero at :math:`h_b` and one at the ice–ocean interface. The number of
grid levels is specified during compilation in **comp\_ice** by setting
the variable `NBGCLYR` equal to an integer (:math:`n_b`) .

Ice tracers and microstructural properties defined on the bio-grid are
referenced in two ways: as `bgrid` :math:`=n_b+2` points and as
igrid\ :math:`=n_b+1` points. For both bgrid and igrid, the first and
last points reference :math:`h_b` and the ice–ocean interface,
respectively, and so take the values :math:`0` and :math:`1`,
respectively. For bgrid, the interior points :math:`[2, n_b+1]` are
spaced at :math:`1/n_b` intervals beginning with `bgrid(2)` :math:` =
1/(2n_b)`. The `igrid` interior points :math:`[2, n_b]` are also
equidistant with the same spacing, but physically coincide with points
midway between those of `bgrid`.

Column configuration
~~~~~~~~~~~~~~~~~~~~

A column modeling capability is available. Because of the boundary
conditions and other spatial assumptions in the model, this is not a
single column, but a small array of columns (minimum grid size is 5x5).
However, the code is set up so that only the single, central column is
used (all other columns are designated as land). The column is located
near Barrow (71.35N, 156.5W). Options for choosing the column
configuration are given in **comp\_ice** (choose `RES col`) and in the
namelist file, **input\_templates/col/ice\_in**. Here, `istep0` and the
initial conditions are set such that the run begins September 1 with no
ice. The grid type is rectangular, dynamics are turned off (`kdyn` = 0) and
one processor is used.

History variables available for column output are ice and snow
temperature, `Tinz` and `Tsnz`. These variables also include thickness
category as a fourth dimension.

.. _init:

Initialization and coupling
---------------------------

The ice model’s parameters and variables are initialized in several
steps. Many constants and physical parameters are set in
**ice\_constants.F90**. Namelist variables (:ref:`tabnamelist`),
whose values can be altered at run time, are handled in *input\_data*
and other initialization routines. These variables are given default
values in the code, which may then be changed when the input file
**ice\_in** is read. Other physical constants, numerical parameters, and
variables are first set in initialization routines for each ice model
component or module. Then, if the ice model is being restarted from a
previous run, core variables are read and reinitialized in
*restartfile*, while tracer variables needed for specific configurations
are read in separate restart routines associated with each tracer or
specialized parameterization. Finally, albedo and other quantities
dependent on the initial ice state are set. Some of these parameters
will be described in more detail in :ref:`tabnamelist`.

The restart files supplied with the code release include the core
variables on the default configuration, that is, with seven vertical
layers and the ice thickness distribution defined by `kcatbound` = 0.
Restart information for some tracers is also included in the  restart
files.

Three namelist variables control model initialization, `ice\_ic`, `runtype`,
and `restart`, as described in :ref:`tab-ic`. It is possible to do an
initial run from a file **filename** in two ways: (1) set runtype =
‘initial’, restart = true and ice\_ic = **filename**, or (2) runtype =
‘continue’ and pointer\_file = **./restart/ice.restart\_file** where
**./restart/ice.restart\_file** contains the line
“./restart/[filename]". The first option is convenient when repeatedly
starting from a given file when subsequent restart files have been
written. With this arrangement, the tracer restart flags can be set to
true or false, depending on whether the tracer restart data exist. With
the second option, tracer restart flags are set to ‘continue’ for all
active tracers.

An additional namelist option, `restart\_ext` specifies whether halo cells
are included in the restart files. This option is useful for tripole and
regional grids, but can not be used with PIO.

MPI is initialized in *init\_communicate* for both coupled and
stand-alone MPI runs. The ice component communicates with a flux coupler
or other climate components via external routiines that handle the
variables listed in :ref:`tab-flux-cpl`. For stand-alone runs,
routines in **ice\_forcing.F90** read and interpolate data from files,
and are intended merely to provide guidance for the user to write his or
her own routines. Whether the code is to be run in stand-alone or
coupled mode is determined at compile time, as described below.

:ref:`tab-ic` : *Ice initial state resulting from combinations of*
`ice\_ic`, `runtype` and `restart`. :math:`^a`\ *If false, restart is reset to
true.* :math:`^b`\ *restart is reset to false.* :math:`^c`\ ice\_ic *is
reset to ‘none.’*

.. _tab-ic:

.. table:: Table 4

   +----------------+--------------------------+--------------------------------------+----------------------------------------+
   | ice\_ic        |                          |                                      |                                        |
   +================+==========================+======================================+========================================+
   |                | initial/false            | initial/true                         | continue/true (or false\ :math:`^a`)   |
   +----------------+--------------------------+--------------------------------------+----------------------------------------+
   | none           | no ice                   | no ice\ :math:`^b`                   | restart using **pointer\_file**        |
   +----------------+--------------------------+--------------------------------------+----------------------------------------+
   | default        | SST/latitude dependent   | SST/latitude dependent\ :math:`^b`   | restart using **pointer\_file**        |
   +----------------+--------------------------+--------------------------------------+----------------------------------------+
   | **filename**   | no ice\ :math:`^c`       | start from **filename**              | restart using **pointer\_file**        |
   +----------------+--------------------------+--------------------------------------+----------------------------------------+

.. _parameters:

Choosing an appropriate time step
---------------------------------

Transport in thickness space imposes a restraint on the time
step, given by the ice growth/melt rate and the smallest range of
thickness among the categories,
:math:`\Delta t<\min(\Delta H)/2\max(f)`, where :math:`\Delta H` is the
distance between category boundaries and :math:`f` is the thermodynamic
growth rate. For the 5-category ice thickness distribution used as the
default in this distribution, this is not a stringent limitation:
:math:`\Delta t < 19.4` hr, assuming :math:`\max(f) = 40` cm/day.


Model output
------------

.. _history:

Diagnostic files
~~~~~~~~~~~~~~~~

*Needs specifics for icepack*


Restart files
~~~~~~~~~~~~~

CICE now provides restart data in binary unformatted or  formats, via
the `IO\_TYPE` flag in **comp\_ice** and namelist variable
`restart\_format`. Restart and history files must use the same format. As
with the history output, there is also an option for writing parallel
restart files using PIO.

The restart files created by CICE contain all of the variables needed
for a full, exact restart. The filename begins with the character string
‘iced.’, and the restart dump frequency is given by the namelist
variables `dumpfreq` and `dumpfreq\_n`. The pointer to the filename from
which the restart data is to be read for a continuation run is set in
`pointer\_file`. The code assumes that auxiliary binary tracer restart
files will be identified using the same pointer and file name prefix,
but with an additional character string in the file name that is
associated with each tracer set. All variables are included in  restart
files.

Additional namelist flags provide further control of restart behavior.
`dump\_last` = true causes a set of restart files to be written at the end
of a run when it is otherwise not scheduled to occur. The flag
`use\_restart\_time` enables the user to choose to use the model date
provided in the restart files. If `use\_restart\_time` = false then the
initial model date stamp is determined from the namelist parameters.
lcdf64 = true sets 64-bit  output, allowing larger file sizes with
version 3.

Routines for gathering, scattering and (unformatted) reading and writing
of the “extended" global grid, including the physical domain and ghost
(halo) cells around the outer edges, allow exact restarts on regional
grids with open boundary conditions, and they will also simplify
restarts on the various tripole grids. They are accessed by setting
`restart\_ext` = true in namelist. Extended grid restarts are not
available when using PIO; in this case extra halo update calls fill
ghost cells for tripole grids (do not use PIO for regional grids).

Two restart files are included with the CICE v5 code distribution, for
the gx3 and gx1 grids. The were created using the default model
configuration (settings as in **comp\_ice** and **ice\_in**), but
initialized with no ice. The gx3 case was run for 1 year using the 1997
forcing data provided with the code. The gx1 case was run for 20 years,
so that the date of restart in the file is 1978-01-01. Note that the
restart dates provided in the restart files can be overridden using the
namelist variables `use\_restart\_time`, `year\_init` and `istep0`. The
forcing time can also be overridden using `fyear\_init`.

Execution procedures
--------------------

* point to appropriate info online?*


Adding things
-------------

.. _addtrcr:

Tracers
~~~~~~~

*point to workflow online?*

Each optional tracer has its own module, **ice\_[tracer].F90**, which
also contains as much of the additional tracer code as possible, and for
backward compatibility of binary restart files, each new tracer has its
own binary restart file. We recommend that the logical namelist variable
`tr\_[tracer]` be used for all calls involving the new tracer outside of
**ice\_[tracer].F90**, in case other users do not want to use that
tracer.

A number of optional tracers are available in the code, including ice
age, first-year ice area, melt pond area and volume, brine height,
aerosols, and level ice area and volume (from which ridged ice
quantities are derived). Salinity, enthalpies, age, aerosols, level-ice
volume, brine height and most melt pond quantities are volume-weighted
tracers, while first-year area, pond area, level-ice area and all of the
biogeochemistry tracers in this release are area-weighted tracers. In
the absence of sources and sinks, the total mass of a volume-weighted
tracer such as aerosol (kg) is conserved under transport in horizontal
and thickness space (the mass in a given grid cell will change), whereas
the aerosol concentration (kg/m) is unchanged following the motion, and
in particular, the concentration is unchanged when there is surface or
basal melting. The proper units for a volume-weighted mass tracer in the
tracer array are kg/m.

In several places in the code, tracer computations must be performed on
the conserved “tracer volume" rather than the tracer itself; for
example, the conserved quantity is :math:`h_{pnd}a_{pnd}a_{lvl}a_{i}`,
not :math:`h_{pnd}`. Conserved quantities are thus computed according to
the tracer dependencies, and code must be included to account for new
dependencies (e.g., :math:`a_{lvl}` and :math:`a_{pnd}` in
**ice\_itd.F90** and **ice\_mechred.F90**).

To add a tracer, follow these steps using one of the existing tracers as
a pattern.

#. **ice\_domain\_size.F90**: increase `max\_ntrcr` (can also add option
   to **comp\_ice** and **bld/Macros.\***)

#. **ice\_state.F90**: declare `nt\_[tracer]` and `tr\_[tracer]`

#. **ice\_[tracer].F90**: create initialization, physics, restart
   routines

#. **ice\_fileunits.F90**: add new dump and restart file units

#. **ice\_init.F90**: (some of this may be done in **ice\_[tracer].F90**
   instead)

   -  add new module and `tr\_[tracer]` to list of used modules and
      variables

   -  add logical namelist variable `tr\_[tracer]`

   -  initialize namelist variable

   -  broadcast namelist variable

   -  print namelist variable to diagnostic output file

   -  increment number of tracers in use based on namelist input (`ntrcr`)

   -  define tracer types (`trcr\_depend` = 0 for ice area tracers, 1 for
      ice volume, 2 for snow volume, 2+nt\_[tracer] for dependence on
      other tracers)

#. **ice\_itd.F90**, **ice\_mechred.F90**: Account for new dependencies
   if needed.

#. **CICE\_InitMod.F90**: initialize tracer (includes reading restart
   file)

#. **CICE\_RunMod.F90**, **ice\_step\_mod.F90**:

   -  call routine to write tracer restart data

   -  call physics routines in **ice\_[tracer].F90** (often called from
      **ice\_step\_mod.F90**)

#. **ice\_restart.F90**: define restart variables (for binary,  and PIO)

#. **ice\_history\_[tracer].F90**: add history variables

#. **ice\_in**: add namelist variables to *tracer\_nml* and
   *icefields\_nml*

#. If strict conservation is necessary, add diagnostics as noted for
   topo ponds in Section :ref:`ponds`.

