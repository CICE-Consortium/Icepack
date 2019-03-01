:tocdepth: 3

.. _implementation:

Implementation
========================

Icepack is written in FORTRAN90 and runs on platforms using UNIX, LINUX,
and other operating systems. The code is not parallelized. (CHANGE IF OPENMP IS IMPLEMENTED)

Icepack consists of the sea ice column physics code, contained in the 
**columnphysics/** directory, and a **configuration/** directory that includes
a driver for testing the column physics and a set of scripts for configuring the tests.
Icepack is designed such that the column physics code may be used by a host sea ice
model without direct reference to the driver or scripts, although these may be consulted for 
guidance when coupling the column physics code to the host sea ice model 
(`CICE <https://github.com/CICE-Consortium/CICE>`_ may also be useful for this.)  Information
about the interface between the column physics and the driver or host sea ice model is
located in the :ref:`init` section.

.. _dirstructure:

Directory structure
-------------------

The present code distribution includes source code for the column physics,
source code for the driver, and the scripts.  Forcing data is available from the ftp site.
The directory structure of Icepack is as follows.  All columnphysics filename have a prefix
of icepack\_ and all driver files are prefixed with icedrv\_*.

**LICENSE.pdf**
  license for using and sharing the code

**DistributionPolicy.pdf**
  policy for using and sharing the code

**README.md**
  basic information and pointers

**columnphysics/**
  columnphysics source code, see :ref:`dev_colphys`

**configuration/scripts/**
  support scripts, see :ref:`dev_scripts`

**configuration/driver/**
  icepack driver code, see :ref:`dev_driver`

**doc/**
  documentation

**icepack.setup**
  main icepack script for creating cases

A case (compile) directory is created upon initial execution of the script 
**icepack.setup** at the user-specified location provided after the -c flag. 
Executing the command ``./icepack.setup -h`` provides helpful information for 
this tool.

.. _grids:

Grid and boundary conditions 
----------------------------

The driver configures a collection of grid cells on which the column physics code 
will be run. This "horizontal" grid is a vector of length ``nx``, with a minimum length 
of 4.   
The grid vector is initialized with different sea ice conditions, such as open 
water, a uniform slab of ice, a multi-year ice thickness distribution with snow, 
and land. For simplicity, the same forcing values are applied to all grid cells. 

Icepack includes two vertical grids.  The basic vertical grid contains 
``nilyr`` equally spaced grid cells.  
History variables available for column output are ice and snow
temperature, ``Tinz`` and ``Tsnz``. These variables also include thickness
category as a fourth dimension.

In addition, there is a bio-grid that 
can be more finely resolved and includes additional nodes for boundary conditions.
It is used for solving the brine height variable :math:`h_b` and for
discretizing the vertical transport equations of biogeochemical tracers.
The bio-grid is a non-dimensional vertical grid which takes the value
zero at :math:`h_b` and one at the ice–ocean interface. The number of
grid levels is specified during compilation by setting
the variable ``NBGCLYR`` equal to an integer (:math:`n_b`) .

Ice tracers and microstructural properties defined on the bio-grid are
referenced in two ways: as ``bgrid`` :math:`=n_b+2` points and as
igrid\ :math:`=n_b+1` points. For both bgrid and igrid, the first and
last points reference :math:`h_b` and the ice–ocean interface,
respectively, and so take the values :math:`0` and :math:`1`,
respectively. For bgrid, the interior points :math:`[2, n_b+1]` are
spaced at :math:`1/n_b` intervals beginning with `bgrid(2)` = 
:math:`1/(2n_b)`. The ``igrid`` interior points :math:`[2, n_b]` are also
equidistant with the same spacing, but physically coincide with points
midway between those of ``bgrid``.

.. _init:

Initialization and Forcing
--------------------------

Icepack’s parameters and variables are initialized in several
steps. Many constants and physical parameters are set in
**icepack\_parameters.F90**. In the current driver implementation,
a namelist file is read to setup the model.
Namelist values are given default
values in the code, which may then be changed when the input file
**icepack\_in** is read. Other physical constants, numerical parameters, and
variables are first set in initialization routines for each ice model
component or module. Then, if the ice model is being restarted from a
previous run, core variables are read and reinitialized in
*restartfile*, while tracer variables needed for specific configurations
are read in separate restart routines associated with each tracer or
specialized parameterization. Finally, albedo and other quantities
dependent on the initial ice state are set. Some of these parameters
will be described in more detail in the :ref:`tabnamelist`.

Two namelist variables control model initialization, ``ice_ic``
and ``restart``.  Setting ``ice_ic`` = 'default' causes the model to run using
initial values set in the code.  To start
from a file **filename**, set 
``restart`` = .true. and ``ice_ic`` = **filename**.  When restarting using the Icepack
driver, for simplicity the tracers are assumed to be set the same way (on/off) as in the
run that created the restart file; i.e. that the restart file contains exactly the 
information needed for the new run.  CICE is more flexible in this regard.

For stand-alone runs,
routines in **icedrv\_forcing.F90** read and interpolate data from files,
and are intended merely for testing, although they can also provide guidance for 
the user to write his or her own routines. 

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


.. _history:

Model output
------------

History output from Icepack is not currently supported in the Icepack driver, except
in restart files.
The sea ice model `CICE <https://github.com/CICE-Consortium/CICE>`_ provides extensive 
options for model output, including many derived output variables.

Diagnostic files
~~~~~~~~~~~~~~~~

Icepack writes diagnostic information for each grid cell as a separate file, 
**ice\_diag.\***, identified by the initial ice state of the grid cell (no ice, slab, land, etc).


Restart files
~~~~~~~~~~~~~

Icepack provides restart data in binary unformatted format. The restart files 
created by the Icepack driver contain all of the variables needed
for a full, exact restart. The filename begins with the character string
‘iced.’ and is placed in the directory specified by the namelist variable
``restart_dir``. The restart dump frequency is given by the namelist
variable ``dumpfreq``. The namelist variable ``ice_ic`` contains the
pointer to the filename from which the restart data is to be read and 
the namelist option ``restart`` must be set to ``.true.`` to use the file.
``dump_last`` namelist can also be set to true to trigger restarts automatically
at then end of runs.

.. _bgc-hist:

Biogeochemistry History Fields
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

History output is not provided with Icepack.  This documentation
indicates what is available for output and is implemented in CICE.

Table :ref:`tab-bio-history` lists the
biogeochemical tracer history flags along with a short description and
the variable or variables saved. Not listed are flags appended with
\_ai, i.e. f\_fbio\_ai. These fields are identical to their counterpart.
i.e. f\_fbio, except they are averaged by ice area.

.. _tab-bio-history:

.. csv-table:: *Biogeochemical History variables*
   :header: "History Flag", "Definition", "Variable(s)", "Units"
   :widths: 10, 25, 20, 10

   "f\_faero\_atm", "atmospheric aerosol deposition flux", "faero\_atm", "kg m\ :math:`^{-2}` s\ :math:`^{-1}`"
   "f\_faero\_ocn", "aerosol flux from ice to ocean", "faero\_ocn", "kg m\ :math:`^{-2}` s\ :math:`^{-1}`"
   "f\_aero", "aerosol mass (snow and ice ssl and int)", "aerosnossl, aerosnoint,aeroicessl, aeroiceint", "kg/kg"
   "f\_fbio", "biological ice to ocean flux", "fN, fDOC, fNit, fAm,fDON,fFep\ :math:`^a`, fFed\ :math:`^a`, fSil,fhum, fPON, fDMSPd,fDMS, fDMSPp, fzaero", "mmol m\ :math:`^{-2}` s\ :math:`^{-1}`"
   "f\_zaero", "bulk z-aerosol mass fraction", "zaero", "kg/kg"
   "f\_bgc\_S", "bulk z-salinity", "bgc\_S", "ppt"
   "f\_bgc\_N", "bulk algal N concentration", "bgc\_N", "mmol m\ :math:`^{-3}`"
   "f\_bgc\_C", "bulk algal C concentration", "bgc\_C", "mmol m\ :math:`^{-3}`"
   "f\_bgc\_DOC", "bulk DOC concentration", "bgc\_DOC", "mmol m\ :math:`^{-3}`"
   "f\_bgc\_DON", "bulk DON concentration", "bgc\_DON", "mmol m\ :math:`^{-3}`"
   "f\_bgc\_DIC", "bulk DIC concentration", "bgc\_DIC", "mmol m\ :math:`^{-3}`"
   "f\_bgc\_chl", "bulk algal chlorophyll concentration", "bgc\_chl", "mg chl m\ :math:`^{-3}`"
   "f\_bgc\_Nit", "bulk nitrate concentration", "bgc\_Nit", "mmol m\ :math:`^{-3}`"
   "f\_bgc\_Am", "bulk ammonium concentration", "bgc\_Am", "mmol m\ :math:`^{-3}`"
   "f\_bgc\_Sil", "bulk silicate concentration", "bgc\_Sil", "mmol m\ :math:`^{-3}`"
   "f\_bgc\_DMSPp", "bulk particulate DMSP concentration", "bgc\_DMSPp", "mmol m\ :math:`^{-3}`"
   "f\_bgc\_DMSPd", "bulk dissolved DMSP concentration", "bgc\_DMSPd", "mmol m\ :math:`^{-3}`"
   "f\_bgc\_DMS", "bulk DMS concentration", "bgc\_DMS", "mmol m\ :math:`^{-3}`"
   "f\_bgc\_Fe", "bulk dissolved and particulate iron conc.", "bgc\_Fed, bgc\_Fep", ":math:`\mu\,`\ mol m\ :math:`^{-3}`"
   "f\_bgc\_hum", "bulk humic matter concentration", "bgc\_hum", "mmol m\ :math:`^{-3}`"
   "f\_bgc\_PON", "bulk passive mobile tracer conc.", "bgc\_PON", "mmol m\ :math:`^{-3}`"
   "f\_upNO", "Total algal :math:`{\mbox{NO$_3$}}` uptake rate", "upNO", "mmol m\ :math:`^{-2}` d\ :math:`^{-1}`"
   "f\_upNH", "Total algal :math:`{\mbox{NH$_4$}}` uptake rate", "upNH", "mmol m\ :math:`^{-2}` d\ :math:`^{-1}`"
   "f\_bgc\_ml", "upper ocean tracer concentrations", "ml\_N, ml\_DOC, ml\_Nit,ml\_Am, ml\_DON, ml\_Fep\ :math:`^b`,ml\_Fed\ :math:`^b`, ml\_Sil, ml\_hum, ml\_PON,ml\_DMS, ml\_DMSPd, ml\_DMSPp", "mmol m\ :math:`^{-3}`"
   "f\_bTin", "ice temperature on the bio grid", "bTizn", ":math:`^o`\ C"
   "f\_bphi", "ice porosity on the bio grid", "bphizn", "%"
   "f\_iDin", "brine eddy diffusivity on the interface bio grid", "iDin", "m\ :math:`^{2}` d\ :math:`^{-1}`"
   "f\_iki", "ice permeability on the interface bio grid", "ikin", "mm\ :math:`^{2}`"
   "f\_fbri", "ratio of brine tracer height to ice thickness", "fbrine", "1"
   "f\_hbri", "brine tracer height", "hbrine", "m"
   "f\_zfswin", "internal ice PAR on the interface bio grid", "zfswin", "W m\ :math:`^{-2}`"
   "f\_bionet", "brine height integrated tracer concentration", "algalN\_net, algalC\_net, chl\_net, pFe\ :math:`^c`\ \_net, dFe\ :math:`^c`\ \_net, Sil\_net, Nit\_net, Am\_net, hum\_net, PON\_net, DMS\_net, DMSPd\_net, DMSPp\_net, DOC\_net, zaero\_net, DON\_net", "mmol m\ :math:`^{-2}`"
   "f\_biosnow", snow integrated tracer concentration", "algalN\_snow, algalC\_snow,chl\_snow, pFe\ :math:`^c`\ \_snow, dFe\ :math:`^c`\ \_snow,Sil\_snow, Nit\_snow, Am\_snow, hum\_snow, PON\_snow, DMS\_snow, DMSPd\_snow, DMSPp\_snow, DOC\_snow, zaero\_snow, DON\_snow", "mmol m\ :math:`^{-2}`"
   "f\_grownet", "Net specific algal growth rate", "grow\_net", "m d\ :math:`^{-1}`"
   "f\_PPnet", "Net primary production", "PP\_net", "mgC m\ :math:`^{-2}` d\ :math:`^{-1}`"
   "f\_algalpeak", "interface bio grid level of peak chla", "peak\_loc", "1"
   "f\_zbgc\_frac", "mobile fraction of tracer", "algalN\_frac, chl\_frac, pFe\_frac,dFe\_frac, Sil\_frac, Nit\_frac,Am\_frac, hum\_frac, PON\_frac,DMS\_frac, DMSPd\_frac, DMSPp\_frac,DOC\_frac, zaero\_frac, DON\_frac", "1"


:math:`^a` units are :math:`\mu`\ mol m\ :math:`^{-2}` s\ :math:`^{-1}`

:math:`^b` units are :math:`\mu`\ mol m\ :math:`^{-3}`

:math:`^c` units are :math:`\mu`\ mol m\ :math:`^{-2}`
