:tocdepth: 3

.. _case_settings:

Case Settings, Model Namelist, and CPPs
====================================================

There are two important files that define the case, **icepack.settings** and 
**icepack_in**.  **icepack.settings** is a list of env variables that define many
values used to setup, build and run the case.  **icepack_in** is the input namelist file
for the icepack driver.  Variables in both files are described below.  In addition,
the first table documents available C Preprocessor Macros.

.. _tabcpps:

Table of C Preprocessor (CPP) Macros
---------------------------------------------------

The Icepack model supports a few C Preprocessor (CPP) Macros.  These
can be turned on during compilation to activate different pieces of source
code.  The main purpose is to introduce build-time code modifications to
include or exclude certain libraries or Fortran language features, in part to
support CICE or other applications.  More information
can be found in :ref:`cicecpps`.  The following CPPs are available.

.. csv-table:: **CPP Macros**
   :header: "CPP name", "description"
   :widths: 15, 60

   "",""
   "**General Macros**", ""
   "NO_I8", "Converts ``integer*8`` to ``integer*4``."
   "NO_R16", "Converts ``real*16`` to ``real*8``."
   "NO_SNICARHC", "Does not compile hardcoded (HC) 5 band snicar tables tables needed by ``shortwave=dEdd_snicar_ad``. May reduce compile time."
   "USE_NETCDF", "Turns on netcdf capabilities in Icepack.  By default and generally, Icepack does not need netcdf."
   "",""
   "**Application Macros**", ""
   "CESMCOUPLED", "Turns on code changes for the CESM coupled application                          "
   "CICE_IN_NEMO", "Turns on code changes for coupling in the NEMO ocean model"


.. _tabsettings:

Table of Icepack Settings
--------------------------

The **icepack.settings** file contains a number of environment variables that define
configuration, file system, run, and build settings.  Several variables are set
by the **icepack.setup** script.  This file is created on a case by case basis and
can be modified as needed.

.. csv-table:: **Icepack settings**
   :header: "variable", "options/format", "description", "default value"
   :widths: 15, 15, 25, 20

   "ICE_CASENAME", "string", "case name", "set by icepack.setup"
   "ICE_SANDBOX", "string", "sandbox directory", "set by icepack.setup"
   "ICE_MACHINE", "string", "machine name", "set by icepack.setup"
   "ICE_ENVNAME", "string", "compilation environment", "set by icepack.setup"
   "ICE_MACHCOMP", "string", "machine_environment name", "set by icepack.setup"
   "ICE_SCRIPTS", "string", "scripts directory", "set by icepack.setup"
   "ICE_CASEDIR", "string", "case directory", "set by icepack.setup"
   "ICE_RUNDIR", "string", "run directory", "set by icepack.setup"
   "ICE_OBJDIR", "string", "compile directory", "${ICE_RUNDIR}/compile"
   "ICE_RSTDIR", "string", "unused", "${ICE_RUNDIR}/restart"
   "ICE_HSTDIR", "string", "unused", "${ICE_RUNDIR}/history"
   "ICE_LOGDIR", "string", "log directory", "${ICE_CASEDIR}/logs"
   "ICE_RSTPFILE", "string", "unused", "undefined"
   "ICE_DRVOPT", "string", "unused", "icepack"
   "ICE_IOTYPE", "none,netcdf", "IO options", "none"
   "ICE_CLEANBUILD", "true,false", "automatically clean before building", "true"
   "ICE_CPPDEFS", "string", "user defined preprocessor macros for build", "null"
   "ICE_QUIETMODE", "true, false", "reduce build output to the screen", "false"
   "ICE_GRID", "col", "grid", "col"
   "ICE_NXGLOB", "4", "number of gridcells", "4"
   "ICE_NTASKS", "1", "number of tasks, must be set to 1", "1"
   "ICE_NTHRDS", "1", "number of threads per task, must be set to 1", "1"
   "ICE_TEST", "string", "test setting if using a test", "set by icepack.setup"
   "ICE_TESTNAME", "string", "test name if using a test", "set by icepack.setup"
   "ICE_BASELINE", "string", "baseline directory name, associated with icepack.setup -bd", "set by icepack.setup"
   "ICE_BASEGEN", "string", "baseline directory name for regression generation, associated with icepack.setup -bg ", "set by icepack.setup"
   "ICE_BASECOM", "string", "baseline directory name for regression comparison, associated with icepack.setup -bc ", "set by icepack.setup"
   "ICE_BFBCOMP", "string", "location of case for comparison, associated with icepack.setup -td", "set by icepack.setup"
   "ICE_SPVAL", "string", "unused", "UnDeFiNeD"
   "ICE_RUNLENGTH", "string", "batch run length default", "  00:10:00"
   "ICE_ACCOUNT", "string", "batch account number", "set by icepack.setup or by default"
   "ICE_QUEUE", "string", "batch queue name", "set by icepack.setup or by default"
   "ICE_THREADED", "true,false", "force threading in compile, will always compile threaded if NTHRDS is gt 1", "false"
   "NICELYR", "integer", "number of vertical layers in the ice", "7"
   "NSNWLYR", "integer", "number of vertical layers in the snow", "1"
   "NICECAT", "integer", "number of ice thickness categories", "5"
   "NFSDCAT", "integer", "number of floe size categories", "12"
   "TRAGE", "0,1", "ice age tracer", "1"
   "TRFY", "0,1", "first-year ice area tracer", "1"
   "TRLVL", "0,1", "deformed ice tracer", "1"
   "TRPND", "0,1", "melt pond tracer", "1"
   "NTRAERO", "integer", "number of aerosol tracers", "1"
   "NTRISO", "integer", "number of water isotope tracers", "1"
   "TRBRI", "0,1", "brine height tracer", "0"
   "TRZS", "", "DEPRECATED", ""
   "TRBGCS", "0,1", "skeletal layer tracer, needs TRBGCZ=0", "0"
   "TRBGCZ", "0,1", "zbgc tracers, needs TRBGCS=0 and TRBRI=1", "0"
   "NBGCLYR", "integer", "number of zbgc layers", "1"
   "TRZAERO", "0-6", "number of z aerosol tracers", "0"
   "TRALG", "0,1,2,3", "number of algal tracers", "0"
   "TRDOC", "0,1,2", "number of dissolved organic carbon", "0"
   "TRDIC", "0,1", "number of dissolved inorganic carbon", "0"
   "TRDON", "0,1", "number of dissolved organic nitrogen", "0"
   "TRFEP", "0,1", "number of particulate iron tracers", "0"
   "TRFED", "0,1", "number of dissolved iron tracers", "0"
   "ICE_SNICARHC", "true,false", "include hardcoded (HC) snicar tables", "false"
   "ICE_BLDDEBUG", "true,false", "turn on compile debug flags", "false"
   "ICE_COVERAGE", "true,false", "turn on code coverage flags", "false"


.. _tabnamelist:

Tables of Namelist Options
----------------------------

The Icepack driver reads a namelist input file, **icepack_in**, consisting of several namelist groups.  The tables below
summarize the different groups and the variables in each group.  The variables are organized alphabetically 
and the default values listed are the values defined in the source code.  Those values will be 
used unless overridden by the Icepack namelist file, **icepack_in**.  The source code default values as listed 
in the table are not necessarily the recommended production values.

setup_nml
~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: **setup_nml namelist options**
   :header: "variable", "options/format", "description", "default value"
   :widths: 15, 15, 30, 15 

   "", "", "", ""
   "``conserv_check``", "logical", "check conservation", "``.false.``"
   "``cpl_bgc``", "logical", "couple bgc thru driver", "``.false.``"
   "``days_per_year``", "integer", "number of days in a model year", "365"
   "``diagfreq``", "integer", "frequency of diagnostic output in timesteps", "24"
   "``diag_file``", "string", "diagnostic output filename", "'ice_diag'"
   "``dumpfreq``", "``d``", "write restart every ``dumpfreq_n`` days", "``y``"
   "", "``m``", "write restart every ``dumpfreq_n`` months", ""
   "", "``y``", "write restart every ``dumpfreq_n`` years", ""
   "``dump_last``", "true/false", "write restart at end of run", "false"
   "``dt``", "seconds", "thermodynamics time step length", "3600."
   "``hbar_init_itd``", "real", "initial modal ice thickness for itd-initialized grid cell in m", "3.0"
   "``hi_init_slab``", "real", "initial ice thickness for slab-initialized grid cell in m", "2.0"
   "``history_format``", "``cdf``", "history file output in netcdf format", "``none``"
   "","``none``","no history output",""
   "``hsno_init_itd``", "real", "initial snow depth for itd-initialized grid cell in m", "0.25"
   "``hsno_init_slab``", "real", "initial snow depth for slab-initialized grid cell in m", "0.0"
   "``ice_ic``", "``default``", "latitude and sst dependent initial condition", "``default``"
   "", "``none``", "no ice", ""
   "", "'path/file'", "restart file name", ""
   "``istep0``", "integer", "initial time step number", "0"
   "``ndtd``", "integer", "number of dynamics/advection/ridging/steps per thermo timestep", "1"
   "``npt``", "integer", "total number of time steps to take", "99999"
   "``restart``", "logical", "initialize using restart file", "``.false.``"
   "``restart_dir``", "string", "path to restart directory", "'./'"
   "``restart_file``", "string", "output file prefix for restart dump", "'iced'"
   "``restart_format``", "``bin``", "restart file output in binary format", "``bin``"
   "","``cdf``","restart file output in netcdf format",""
   "``sst_init``", "real", "initial ocean mixed layer temperature in C", "-1.8"
   "``use_leap_years``", "logical", "include leap days", "``.false.``"
   "``year_init``", "integer", "the initial year if not using restart", "0"
   "", "", "", ""

grid_nml
~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: **grid_nml namelist options**
   :header: "variable", "options/format", "description", "default value"
   :widths: 15, 15, 30, 15 

   "", "", "", ""
   "``kcatbound``", "``-1``", "single category formulation", "1"
   "", "``0``", "old formulation", ""
   "", "``1``", "new formulation with round numbers", ""
   "", "``2``", "WMO standard categories", ""
   "", "``3``", "asymptotic scheme", ""
   "", "", "", ""

tracer_nml
~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: **tracer_nml namelist options**
   :header: "variable", "options/format", "description", "default value"
   :widths: 15, 15, 30, 15 

   "", "", "", ""
   "``tr_aero``", "logical", "aerosols", "``.false.``"
   "``tr_fsd``", "logical", "floe size distribution", "``.false.``"
   "``tr_FY``", "logical", "first-year ice area", "``.false.``"
   "``tr_iage``", "logical", "ice age", "``.false.``"
   "``tr_iso``", "logical", "isotopes", "``.false.``"
   "``tr_lvl``", "logical", "level ice area and volume", "``.false.``"
   "``tr_pond_lvl``", "logical", "level-ice melt ponds", "``.false.``"
   "``tr_pond_topo``", "logical", "topo melt ponds", "``.false.``"
   "``tr_snow``", "logical", "advanced snow physics", "``.false.``"
   "", "", "", ""

thermo_nml
~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: **thermo_nml namelist options**
   :header: "variable", "options/format", "description", "default value"
   :widths: 15, 15, 30, 15 

   "", "", "", ""
   "``a_rapid_mode``", "real", "brine channel diameter in m", "0.5e-3"
   "``aspect_rapid_mode``", "real", "brine convection aspect ratio", "1.0"
   "``conduct``", "``bubbly``", "conductivity scheme :cite:`Pringle07`", "``bubbly``"
   "", "``MU71``", "conductivity :cite:`Maykut71`", ""
   "``dSdt_slow_mode``", "real", "slow drainage strength parameter m/s/K", "-1.5e-7"
   "``floediam``", "real", "effective floe diameter for lateral melt in m", "300.0"
   "``hfrazilmin``", "real", "min thickness of new frazil ice in m", "0.05"
   "``hi_min``", "real", "minimum ice thickness allowed for thermo in m", "0.01"
   "``kitd``", "``0``", "delta function ITD approximation", "1"
   "", "``1``", "linear remapping ITD approximation", ""
   "``ksno``", "real", "snow thermal conductivity", "0.3"
   "``ktherm``", "``-1``", "thermodynamic model disabled", "1"
   "", "``1``", "Bitz and Lipscomb thermodynamic model", ""
   "", "``2``", "mushy-layer thermodynamic model", ""
   "``phi_c_slow_mode``", ":math:`0<\phi_c < 1`", "critical liquid fraction", "0.05"
   "``phi_i_mushy``", ":math:`0<\phi_i < 1`", "solid fraction at lower boundary", "0.85"
   "``Rac_rapid_mode``", "real", "critical Rayleigh number", "10.0"
   "``Tliquidus_max``", "real", "maximum liquidus temperature of mush (C)", "0.0"
   "``tscale_pnd_drain``", "real", "mushy pond macroscopic drainage timescale in days", "10."
   "", "", "", ""


dynamics_nml
~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: **dynamics_nml namelist options**
   :header: "variable", "options/format", "description", "default value"
   :widths: 15, 15, 30, 15 

   "", "", "", ""
   "``Cf``", "real", "ratio of ridging work to PE change in ridging", "17.0"
   "``kstrength``", "``0``", "ice strength formulation :cite:`Hibler79`", "1"
   "", "``1``", "ice strength formulation :cite:`Rothrock75`", ""
   "``krdg_partic``", "``0``", "old ridging participation function", "1"
   "", "``1``", "new ridging participation function", ""
   "``krdg_redist``", "``0``", "old ridging redistribution function", "1"
   "", "``1``", "new ridging redistribution function", ""
   "``mu_rdg``", "real", "e-folding scale of ridged ice for ``krdg_partic`` = 1 in m^0.5", "3.0"
   "", "", "", ""

shortwave_nml
~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: **shortwave_nml namelist options**
   :header: "variable", "options/format", "description", "default value"
   :widths: 15, 15, 30, 15 

   "", "", "", ""
   "``ahmax``", "real", "albedo is constant above this thickness in meters", "0.3"
   "``albedo_type``", "``ccsm3``", "NCAR CCSM3 albedo implementation", "``ccsm3``"
   "", "``constant``", "four constant albedos", ""
   "``albicei``", ":math:`0<\alpha <1`", "near infrared ice albedo for thicker ice", "0.36"
   "``albicev``", ":math:`0<\alpha <1`", "visible ice albedo for thicker ice", "0.78"
   "``albsnowi``", ":math:`0<\alpha <1`", "near infrared, cold snow albedo", "0.70"
   "``albsnowv``", ":math:`0<\alpha <1`", "visible, cold snow albedo", "0.98"
   "``dT_mlt``", "real", ":math:`\Delta` temperature per :math:`\Delta` snow grain radius", "1.5"
   "``kalg``", "real", "absorption coefficient for algae", "0.6"
   "``rsnw_mlt``", "real", "maximum melting snow grain radius", "1500."
   "``R_ice``", "real", "tuning parameter for sea ice albedo from Delta-Eddington shortwave", "0.0"
   "``R_pnd``", "real", "tuning parameter for ponded sea ice albedo from Delta-Eddington shortwave", "0.0"
   "``R_snw``", "real", "tuning parameter for snow (broadband albedo) from Delta-Eddington shortwave", "1.5"
   "``shortwave``", "``ccsm3``", "NCAR CCSM3 shortwave distribution method", "``dEdd``"
   "", "``dEdd``", "Delta-Eddington method (3-band)", ""
   "", "``dEdd_snicar_ad``", "Delta-Eddington method with 5-band snow", ""
   "``snw_ssp_table``", "``snicar``", "lookup table for `dEdd_snicar_ad`", "test"
   "", "``test``", "reduced lookup table for `dEdd_snicar_ad` testing", ""
   "``sw_dtemp``", "real", "temperature from melt for sw_redist", "0.02"
   "``sw_frac``", "real", "fraction of shortwave redistribution", "0.9"
   "``sw_redist``", "logical", "shortwave redistribution", ".false."
   "", "", "", ""

ponds_nml
~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: **ponds_nml namelist options**
   :header: "variable", "options/format", "description", "default value"
   :widths: 15, 15, 30, 15 

   "", "", "", ""
   "``dpscale``", "real", "scaling factor for flushing in permeable ice (ktherm=1)", "1.e-3"
   "``frzpnd``", "``cesm``", "CESM pond refreezing forumulation", "``cesm``"
   "", "``hlid``", "Stefan refreezing with pond ice thickness", ""
   "``hp1``", "real", "critical ice lid thickness for topo ponds in m", "0.01"
   "``hs0``", "real", "snow depth of transition to bare sea ice in m", "0.03"
   "``hs1``", "real", "snow depth of transition to pond ice in m", "0.03"
   "``pndaspect``", "real", "aspect ratio of pond changes (depth:area)", "0.8"
   "``rfracmax``", ":math:`0 \le r_{max} \le 1`", "maximum melt water added to ponds", "0.85"
   "``rfracmin``", ":math:`0 \le r_{min} \le 1`", "minimum melt water added to ponds", "0.15"
   "", "", "", ""

snow_nml
~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: **snow_nml namelist options**
   :header: "variable", "options/format", "description", "default value"
   :widths: 15, 15, 30, 15 

   "", "", "", ""
   "``drhosdwind``", "real", "wind compaction factor for snow", "27.3"
   "``rhosmin``", "real", "minimum snow density", "100.0"
   "``rhosmax``", "real", "maximum snow density", "450.0"
   "``rhosnew``", "real", "new snow density", "100.0"
   "``rsnw_fall``", "real", "radius of new snow (um)", "54.526"
   "``rsnw_tmax``", "real", "maximum snow radius (um)", "1500.0"
   "``snw_aging_table``", "test", "snow aging lookup table", "test"
   "", "snicar", "(not available in Icepack)", ""
   "``snwgrain``",  "logical", "snow grain metamorphosis", ".false."
   "``snwlvlfac``", "real", "fraction increase in bulk snow redistribution", "0.3"
   "``snwredist``", "``snwITDrdg``", "snow redistribution using ITD/ridges", "none"
   "", "``bulk``", "bulk snow redistribution", ""
   "", "``none``", "no snow redistribution", ""
   "``use_smliq_pnd``", "logical", "use liquid in snow for ponds", ".false."
   "``windmin``",  "real", "minimum wind speed to compact snow", "10.0"
   "", "", "", ""

forcing_nml
~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: **forcing_nml namelist options**
   :header: "variable", "options/format", "description", "default value"
   :widths: 15, 15, 30, 15 

   "", "", "", ""
   "``atmbndy``", "string", "bulk transfer coefficients", "``similarity``"
   "", "``similarity``", "stability-based boundary layer", ""
   "", "``constant``", "constant-based boundary layer", ""
   "", "``mixed``", "stability-based, but constant for sensible+latent heatfluxes", ""
   "``atmiter_conv``", "real", "convergence criteria for ustar", "0.0"
   "``atm_data_file``", "string", "file containing atmospheric data", "' '"
   "``atm_data_format``", "``bin``", "read direct access binary forcing files", "``bin``"
   "", "``nc``", "read netCDF forcing files", ""
   "``atm_data_type``", "``clim``", "monthly climatology (see :ref:`force`)", "``default``"
   "", "``CFS``", "CFS model output  (see :ref:`force`)", ""
   "", "``default``", "constant values defined in the code", ""
   "", "``ISPOL``", "ISPOL experiment data  (see :ref:`force`)", ""
   "", "``MDF``", "Merged Data File (MDF) formatted forcing data  (see :ref:`init`)", ""
   "", "``NICE``", "N-ICE experiment data  (see :ref:`force`)", ""
   "``bgc_data_file``", "string", "file containing biogeochemistry data", "' '"
   "``bgc_data_format``", "``bin``", "read direct access binary forcing files", "``bin``"
   "``bgc_data_type``", "``clim``", "bgc climatological data", "``default``"
   "", "``default``", "constant values defined in the code", ""
   "", "``ncar``", "POP ocean forcing data", ""
   "``calc_strair``", "``.false.``", "read wind stress and speed from files", "``.true.``"
   "", "``.true.``", "calculate wind stress and speed", ""
   "``calc_Tsfc``", "logical", "calculate surface temperature", "``.true.``"
   "``congel_freeze``", "``one-step``", "immediately freeze congelation ice", "``two-step``"
   "", "``two-step``", "delayed freezing of congelation ice", ""
   "``cpl_frazil``", "``external``", "frazil water/salt fluxes are handled outside of Icepack", "``fresh_ice_correction``"
   "", "``fresh_ice_correction``", "correct fresh-ice frazil water/salt fluxes for mushy physics", ""
   "", "``internal``", "send full frazil water/salt fluxes for mushy physics", ""
   "``data_dir``", "string", "path to forcing data directory", "' '"
   "``default_season``", "``summer``", "forcing initial summer values", "``winter``"
   "", "``winter``", "forcing initial winter values", ""
   "``emissivity``", "real", "emissivity of snow and ice", "0.985"
   "``fbot_xfer_type``", "``Cdn_ocn``", "variabler ocean heat transfer coefficient scheme", "``constant``"
   "", "``constant``", "constant ocean heat transfer coefficient", ""
   "``formdrag``", "logical", "calculate form drag", "``.false.``"
   "``fyear_init``", "integer", "first year of atmospheric forcing data", "1998"
   "``highfreq``", "logical", "high-frequency atmo coupling", "``.false.``"
   "``hmix_fixed``", "real", "constant ocean mixed layer depth in m", "20.0"
   "``lateral_flux_type``", "``uniform_ice``", "flux ice with identical properties into the cell when closing (Icepack only)", ""
   "", "``none``", "advect open water into the cell when closing (Icepack only)", ""
   "``ice_data_file``", "string", "file containing ice opening, closing data", "' '"
   "``l_mpond_fresh``", "``.false.``", "release pond water immediately to ocean", "``.false.``"
   "", "``true``", "retain (topo) pond water until ponds drain", ""
   "``natmiter``", "integer", "number of atmo boundary layer iterations", "5"
   "``oceanmixed_ice``", "logical", "active ocean mixed layer calculation", "``.false.``"
   "``ocn_data_file``", "string", "file containing ocean data", "' ' "
   "``ocn_data_format``", "``bin``", "read direct access binary forcing files", "``bin``"
   "", "``nc``", "read netCDF forcing files", ""
   "``ocn_data_type``", "``default``", "constant values defined in the code", "``default``"
   "", "``ISPOL``", "ISPOL experiment data  (see :ref:`force`)", ""
   "", "``MDF``", "Merged Data File (MDF) formatted forcing data  (see :ref:`init`)", ""
   "", "``NICE``", "N-ICE experiment data  (see :ref:`force`)", ""
   "", "``SHEBA``", "Opening/closing dataset from SHEBA", ""
   "``precalc_forc``", "logical", "average/interpolate forcing data on initialization", "``.false.``"
   "``precip_units``", "``mks``", "liquid precipitation data units", "``mks``"
   "", "``mm_per_month``", "", ""
   "", "``mm_per_sec``", "(same as MKS units)", ""
   "", "``m_per_sec``", "", ""
   "``qdp_fixed``", "real", "constant oceanic heat flux convergence W/m^2", "0.0"
   "``restore_ocn``", "logical", "restore sst to data", "``.false.``"
   "``saltflux_option``", "``constant``","salt flux is referenced to a constant salinity","``constant``"
   "","``prognostic``","use actual sea ice bulk salinity in flux"
   "``sss_fixed``", "real", "constant ocean mixed layer salinity in ppt", "34.0"
   "``tfrz_option``","``constant``", "constant ocean freezing temperature (Tocnfrz)","``mushy``"
   "", "``linear_salt``", "linear function of salinity (ktherm=1)", ""
   "", "``minus1p8``", "constant ocean freezing temperature (:math:`-1.8^{\circ} C`)", ""
   "", "``mushy``", "matches mushy-layer thermo (ktherm=2)", ""
   "``trestore``", "integer", "sst restoring time scale (days)", "90"
   "``update_ocn_f``", "``.false.``", "do not include frazil water/salt fluxes in ocn fluxes", "``.false.``"
   "", "``true``", "include frazil water/salt fluxes in ocn fluxes", ""
   "``ustar_min``", "real", "minimum value of ocean friction velocity in m/s", "0.005"
   "``wave_spec_type``", "``constant``", "wave data file is provided, sea surface height generated using constant phase (1 iteration of wave fracture)", "``none``"
   "", "``none``", "no wave data provided, no wave-ice interactions (not recommended when using the FSD)", ""
   "", "``profile``", "no wave data file is provided, use fixed dummy wave spectrum, for testing, sea surface height generated using constant phase (1 iteration of wave fracture)", ""
   "", "``random``", "wave data file is provided, sea surface height generated using random number (multiple iterations of wave fracture)", ""
   "``ycycle``", "integer", "number of years in forcing data cycle", "1"
   "", "", "", ""
   
* = If Icepack is run stand-alone and wave_spec_type is not set to none, then a fixed wave spectrum is defined in the code to use for testing. As with other input data, this spectrum should not be used for production runs or publications.

zbgc_nml
~~~~~~~~~~~~~~~~~~~~~~~~~

.. _tab-bio-tracers-namelist2:

.. csv-table:: **zbgc_nml namelist options**
   :header: "variable", "options/format", "description", "default value"
   :widths: 15, 15, 30, 15

   "", "", "", ""
   "``algaltype_diatoms``", "real", "mobility type between stationary and mobile algal diatoms", "0.0"
   "``algaltype_phaeo``", "real", "mobility type between stationary and mobile algal phaeocystis", "0.0"
   "``algaltype_sp``", "real", "mobility type between stationary and mobile small plankton", "0.0"
   "``algal_vel``", "real", "maximum speed of algae (m/s) :cite:`Lavoie05`", "1.0e-7"
   "``alpha2max_low_diatoms``", "real", "light limitation diatoms 1/(W/m^2)", "0.3"
   "``alpha2max_low_phaeo``", "real", "light limitation phaeocystis 1/(W/m^2)", "0.17"
   "``alpha2max_low_sp``", "real", "light limitation small plankton 1/(W/m^2)", "0.2"
   "``ammoniumtype``", "real", "mobility type between stationary and mobile ammonium", "0.0"
   "``beta2max_diatoms``", "real", "light inhibition diatoms 1/(W/m^2)", "0.001"
   "``beta2max_phaeo``", "real", "light inhibition phaeocystis 1/(W/m^2)", "0.04"
   "``beta2max_sp``", "real", "light inhibition small plankton 1/(W/m^2)", "0.001"
   "``bgc_flux_type``", "``constant``", "constant ice–ocean flux velocity", "``Jin2006``"
   "", "``Jin2006``", "ice–ocean flux velocity of :cite:`Jin06`", ""
   "``chlabs_diatoms``", "real", "chl absorbtion diatoms 1/m/(mg/m^3)", "0.03"
   "``chlabs_phaeo``", "real", "chl absorbtion phaeocystis 1/m/(mg/m^3)", "0.05"
   "``chlabs_sp``", "real", "chl absorbtion small plankton 1/m/(mg/m^3)", "0.01"
   "``dEdd_algae``", "logical", "include ice bgc and/or aerosols in radiative transfer", "``.false.``"
   "``dmspdtype``", "real", "mobility type between stationary and mobile dmspd", "0.0"
   "``dmspptype``", "real", "mobility type between stationary and mobile dmspp", "0.5"
   "``doctype_l``", "real", "mobility type between stationary and mobile doc lipids", "0.0"
   "``doctype_s``", "real", "mobility type between stationary and mobile doc saccharids", "0.0"
   "``dontype_protein``", "real", "mobility type between stationary and mobile don proteins", "0.0"
   "``dustFe_sol``", "real", "solubility fraction", "0.005"
   "``fedtype_1``", "real", "mobility type between stationary and mobile fed lipids", "0.0"
   "``feptype_1``", "real", "mobility type between stationary and mobile fep lipids", "0.5"
   "``frazil_scav``", "real", "increase in initial bio bracer from ocean scavenging", "0.8"
   "``fr_dFe``", "real", "fraction of remineralized iron from algae", "1.0"
   "``fr_graze_diatoms``", "real", "fraction grazed diatoms", "0.19"
   "``fr_graze_e``", "real", "fraction of assimilation excreted", "0.5"
   "``fr_graze_phaeo``", "real", "fraction grazed phaeocystis", "0.19"
   "``fr_graze_s``", "real", "fraction of grazing spilled or slopped", "0.5"
   "``fr_graze_sp``", "real", "fraction grazed small plankton", "0.19"
   "``fr_mort2min``", "real", "fractionation of mortality to Am", "0.9"
   "``fr_resp``", "real", "fraction of algal growth lost due to respiration", "0.05"
   "``fr_resp_s``", "real", "DMSPd fraction of respiration loss as DMSPd", "0.9"
   "``fsal``", "real", "salinity limitation ppt", "1.0"
   "``F_abs_chl_diatoms``", "real", "scales absorbed radiation for dEdd chl diatoms", "2.0"
   "``F_abs_chl_phaeo``", "real", "scales absorbed radiation for dEdd chl phaeocystis", "5.0"
   "``F_abs_chl_sp``", "real", "scales absorbed radiation for dEdd small plankton", "4.0"
   "``f_doc_l``", "real", "fraction of mortality to DOC lipids", "0.5"
   "``f_doc_s``", "real", "fraction of mortality to DOC saccharides", "0.5"
   "``f_don_Am_protein``", "real", "fraction of remineralized DON to ammonium", "1.0"
   "``f_don_protein``", "real", "fraction of spilled grazing to proteins", "0.6"
   "``f_exude_l``", "real", "fraction of exudation to DOC lipids", "1.0"
   "``f_exude_s``", "real", "fraction of exudation to DOC saccharids", "1.0"
   "``grid_o``", "real", "z biology length scale for bottom flux (m)", "0.006"
   "``grid_oS``", "real", "DEPRECATED", ""
   "``grow_Tdep_diatoms``", "real", "temperature dependence growth diatoms per degC", "0.063"
   "``grow_Tdep_phaeo``", "real", "temperature dependence growth phaeocystis per degC", "0.063"
   "``grow_Tdep_sp``", "real", "temperature dependence growth small plankton per degC", "0.063"
   "``humtype``", "real", "mobility type between stationary and mobile hum", "0.0"
   "``initbio_frac``", "real", "fraction of ocean trcr concentration in bio tracers", "1.0"
   "``K_Am_diatoms``", "real", "ammonium half saturation diatoms mmol/m^3", "0.3"
   "``K_Am_phaeo``", "real", "ammonium half saturation phaeocystis mmol/m^3", "0.3"
   "``K_Am_sp``", "real", "ammonium half saturation small plankton mmol/m^3", "0.3"
   "``k_bac_l``", "real", "Bacterial degredation of DOC lipids per day", "0.03"
   "``k_bac_s``", "real", "Bacterial degredation of DOC saccharids per day", "0.03"
   "``k_exude_diatoms``", "real", "algal exudation diatoms per day", "0.0"
   "``k_exude_phaeo``", "real", "algal exudation phaeocystis per day", "0.0"
   "``k_exude_sp``", "real", "algal exudation small plankton per day", "0.0"
   "``K_Fe_diatoms``", "real", "iron half saturation diatoms nM", "1.0"
   "``K_Fe_phaeo``", "real", "iron half saturation phaeocystis nM", "0.1"
   "``K_Fe_sp``", "real", "iron half saturation small plankton nM", "0.2"
   "``k_nitrif``", "real", "nitrification rate per day", "0.046"
   "``K_Nit_diatoms``", "real", "nitrate half saturation diatoms mmol/m^3", "1.0"
   "``K_Nit_phaeo``", "real", "nitrate half saturation phaeocystis mmol/m^3", "1.0"
   "``K_Nit_sp``", "real", "nitrate half saturation small plankton mmol/m^3", "1.0"
   "``K_Sil_diatoms``", "real", "silicate half saturation diatoms mmol/m^3", "4.0"
   "``K_Sil_phaeo``", "real", "silicate half saturation phaeocystis mmol/m^3", "0.0"
   "``K_Sil_sp``", "real", "silicate half saturation small plankton mmol/m^3", "0.0"
   "``kn_bac_protein``", "real", "bacterial degradation of DON per day", "0.2"
   "``l_sk``", "real", "characteristic diffusive scale in m", "2.0"
   "``l_skS``", "real", "DEPRECATED", ""
   "``max_dfe_doc1``", "real", "max ratio of dFe to saccharides in the ice in nm Fe / muM C", "0.2"
   "``max_loss``", "real", "restrict uptake to percent of remaining value", "0.9"
   "``modal_aero``", "logical", "modal aerosols", "``.false.``"
   "``mort_pre_diatoms``", "real", "mortality diatoms per day", "0.007"
   "``mort_pre_phaeo``", "real", "mortality phaeocystis per day", "0.007"
   "``mort_pre_sp``", "real", "mortality small plankton per day", "0.007"
   "``mort_Tdep_diatoms``", "real", "temperature dependence of mortality diatoms per degC", "0.03"
   "``mort_Tdep_phaeo``", "real", "temperature dependence of mortality phaeocystis per degC", "0.03"
   "``mort_Tdep_sp``", "real", "temperature dependence of mortality small plankton per degC", "0.03"
   "``mu_max_diatoms``", "real", "maximum growth rate diatoms per day", "1.44"
   "``mu_max_phaeo``", "real", "maximum growth rate phaeocystis per day", "0.63"
   "``mu_max_sp``", "real", "maximum growth rate small plankton per day", "0.41"
   "``nitratetype``", "real", "mobility type between stationary and mobile nitrate", "-1.0"
   "``op_dep_min``", "real", "light attenuates for optical depths exceeding min", "0.1"
   "``phi_snow``", "real", "snow porosity for brine height tracer (compute from snow density if negative)", "-1.0"
   "``ratio_chl2N_diatoms``", "real", "algal chl to N in mg/mmol diatoms", "2.1"
   "``ratio_chl2N_phaeo``", "real", "algal chl to N in mg/mmol phaeocystis", "0.84"
   "``ratio_chl2N_sp``", "real", "algal chl to N in mg/mmol small plankton", "1.1"
   "``ratio_C2N_diatoms``", "real", "algal C to N in mol/mol diatoms", "7.0"
   "``ratio_C2N_phaeo``", "real", "algal C to N in mol/mol phaeocystis", "7.0"
   "``ratio_C2N_proteins``", "real", "algal C to N in mol/mol proteins", "5.0"
   "``ratio_C2N_sp``", "real", "algal C to N in mol/mol small plankton", "7.0"
   "``ratio_Fe2C_diatoms``", "real", "algal Fe to C in mmol/mol diatoms", "0.0033"
   "``ratio_Fe2C_phaeo``", "real", "algal Fe to C in mmol/mol phaeocystis", "0.1"
   "``ratio_Fe2C_sp``", "real", "algal Fe to C in mmol/mol small plankton", "0.0033"
   "``ratio_Fe2N_diatoms``", "real", "algal Fe to N in mmol/mol diatoms", "0.023"
   "``ratio_Fe2N_phaeo``", "real", "algal Fe to N in mmol/mol phaeocystis", "0.7"
   "``ratio_Fe2N_sp``", "real", "algal Fe to N in mmol/mol small plankton", "0.023"
   "``ratio_Fe2DOC_s``", "real", "Fe to C of DON saccharids nmol/umol", "0.1"
   "``ratio_Fe2DOC_l``", "real", "Fe to C of DOC lipids nmol/umol", "0.033"
   "``ratio_Fe2DON``", "real", "Fe to C of DON nmol/umol", "0.023"
   "``ratio_Si2N_diatoms``", "real", "algal Si to N in mol/mol diatoms", "1.8"
   "``ratio_Si2N_phaeo``", "real", "algal Si to N in mol/mol phaeocystis", "0.0"
   "``ratio_Si2N_sp``", "real", "algal Si to N in mol/mol small plankton", "0.0"
   "``ratio_S2N_diatoms``", "real", "algal S to N in mol/mol diatoms", "0.03"
   "``ratio_S2N_phaeo``", "real", "algal S to N in mol/mol phaeocystis", "0.03"
   "``ratio_S2N_sp``", "real", "algal S to N in mol/mol small plankton", "0.03"
   "``restore_bgc``", "logical", "DEPRECATED", "``.false.``"
   "``R_dFe2dust``", "real", "g/g :cite:`Tagliabue09`", "0.035"
   "``scale_bgc``", "logical", "initialize bgc by scaling with salinity", "``.false.``"
   "``silicatetype``", "real", "mobility type between stationary and mobile silicate", "-1.0"
   "``skl_bgc``", "logical", "DEPRECATED: skeletal layer biogeochemistry", "``.false.``"
   "``solve_zbgc``", "logical", "solve z-biogeochemistry reactions", "``.false.``"
   "``solve_zsal``", "logical", "DEPRECATED", "``.false.``"
   "``tau_max``", "real", "long time mobile to stationary exchanges (s)", "604800.0"
   "``tau_min``", "real", "rapid module to stationary exchanges (s)", "3600.0"
   "``tr_bgc_Am``", "logical", "ammonium tracer", "``.false.``"
   "``tr_bgc_C``", "logical", "algal carbon tracer", "``.false.``"
   "``tr_bgc_chl``", "logical", "algal chlorophyll tracer", "``.false.``"
   "``tr_bgc_DMS``", "logical", "DMS tracer", "``.false.``"
   "``tr_bgc_DON``", "logical", "DON tracer", "``.false.``"
   "``tr_bgc_Fe``", "logical", "iron tracer", "``.false.``"
   "``tr_bgc_hum``", "logical", "refractory DOC", "``.false.``"
   "``tr_bgc_Nit``", "logical", "nitrate tracer", "``.false.``"
   "``tr_bgc_PON``", "logical", "Non-reactive nitrate tracer", "``.false.``"
   "``tr_bgc_Sil``", "logical", "silicate tracer", "``.false.``"
   "``tr_brine``", "logical", "brine height tracer", "``.false.``"
   "``tr_zaero``", "logical", "vertical aerosol tracers", "``.false.``"
   "``t_iron_conv``", "real", "desorption loss pFe to dFe in days", "3065."
   "``t_sk_conv``", "real", "Stefels conversion time in days", "5.0"
   "``t_sk_ox``", "real", "DMS oxidation time in days", "12.0"
   "``T_max``", "real", "maximum brine temperature degC", "0.0"
   "``y_sk_DMS``", "real", "fraction conversion given high yield", "0.7"
   "``zaerotype_bc1``", "real", "mobility type between stationary and mobile zaero bc1", "-1.0"
   "``zaerotype_bc2``", "real", "mobility type between stationary and mobile zaero bc2", "-1.0"
   "``zaerotype_dust1``", "real", "mobility type between stationary and mobile zaero dust1", "-1.0"
   "``zaerotype_dust2``", "real", "mobility type between stationary and mobile zaero dust2", "-1.0"
   "``zaerotype_dust3``", "real", "mobility type between stationary and mobile zaero dust3", "-1.0"
   "``zaerotype_dust4``", "real", "mobility type between stationary and mobile zaero dust4", "-1.0"
   "``z_tracers``", "logical", "allows vertically resolved bgc and/or z-aerosol tracers", "``.false.``"
   "", "", "", ""


.. commented out below
..   "``dbug``", "true/false", "if true, write extra diagnostics", "``.false.``"
..   "``atm_data_format``", "``nc``", "read  atmo forcing files", ""
..   "", "``bin``", "read direct access, binary files", ""
..   "", "``NICE``", "N-ICE experiment data", ""
..   "", "``NICE``", "N-ICE experiment data", ""
..   "", "``NICE``", "N-ICE experiment data", ""
..   "``grid_o_t``", "real", "z biology for top flux", "5.0"
..   "``restart_bgc``", "logical", "restart tracer values from file", "``.false.``"
..   "``restart_hbrine``", "logical", "", "``.false.``"
..   "``solve_zsal``", "logical", "update salinity tracer profile", "``.false.``"
..   "TRZS", "0,1", "zsalinity tracer, needs TRBRI=1", "0"
  
.. _tuning:

BGC Parameter Arrays
------------------------


Biogeochemical parameter arrays are specified in the code from namelist options in
**icepack\_in**. Table :ref:`tab-bio-tracers2` provides a list of parameter arrays
used in the reaction equations, their representation in the code, a
short description of each and the namelist parameters (Table :ref:`tab-bio-tracers-namelist2`)
used in their definition.

.. _tab-bio-tracers2:

.. csv-table:: *Biogeochemical Reaction Parameters*
   :header: "Text Variable", "Variable in code", "Description", "Value", "units"
   :widths: 7, 20, 15, 15, 15

   ":math:`f_{graze}`", "fr\_graze(1:3)", "fraction of growth grazed", "(``fr_graze_diatoms``, ``fr_graze_sp``, ``fr_graze_phaeo``)", "1"
   ":math:`m_{pre}`", "mort\_pre(1:3)", "maximum mortality rate", "(``mort_pre_diatoms``, ``mort_pre_sp``, ``mort_pre_phaeo``)", "day\ :math:`^{-1}`"
   ":math:`m_{T}`", "mort\_Tdep(1:3)", "mortality temperature decay", "(``mort_Tdep_diatoms``, ``mort_Tdep_sp``, ``mort_Tdep_phaeo``)", ":math:`^o`\ C\ :math:`^{-1}`"
   ":math:`f_{cg}`", "f\_doc(1:2)", "fraction of mortality to :math:`{\mbox{DOC}}`", "(``f_doc_s``, ``f_doc_l``)", "1"
   ":math:`R_{c:n}`", "R\_C2N(1:3)", "algal carbon to nitrogen ratio", "(``ratio_C2N_diatoms``, ``ratio_C2N_sp``, ``ratio_C2N_phaeo``)", "mol/mol"
   ":math:`k_{cb}`", "k\_bac(1:2)", "bacterial degradation of DOC", "(``k_bac_s``, ``k_bac_l``)", "day\ :math:`^{-1}`"
   ":math:`R_{fe:n}`", "R\_Fe2N(1:3)", "algal Fe to N ratio", "(``ratio_Fe2N_diatoms``, ``ratio_Fe2N_sp``, ``ratio_Fe2N_phaeo``)", "mmol/mol"
   ":math:`R_{s:n}`", "R\_S2N(1:3)", "algal S to N ratio", "(``ratio_S2N_diatoms``, ``ratio_S2N_sp``, ``ratio_S2N_phaeo``)", "mol/mol"
   ":math:`K_{{\mbox{NO$_3$}}}`", "K\_Nit(1:3)", ":math:`{\mbox{NO$_3$}}` half saturation constant", "(``K_Nit_diatoms``, ``K_Nit_sp``, ``K_Nit_phaeo``)", "mmol/m\ :math:`^{3}`"
   ":math:`K_{{\mbox{NH$_4$}}}`", "K\_Am(1:3)", ":math:`{\mbox{NH$_4$}}` half saturation constant", "(``K_Am_diatoms``, ``K_Am_sp``, ``K_Am_phaeo``)", "mmol/m\ :math:`^{-3}`"
   ":math:`K_{{\mbox{SiO$_3$}}}`", "K\_Sil(1:3)", "silicate half saturation constant", "(``K_Sil_diatoms``, ``K_Sil_sp``, ``K_Sil_phaeo``)", "mmol/m\ :math:`^{-3}`"
   ":math:`K_{{\mbox{fed}}}`", "K\_Fe(1:3)", "iron half saturation constant", "(``K_Fe_diatoms``, ``K_Fe_sp``, ``K_Fe_phaeo``)", ":math:`\mu`\ mol/m\ :math:`^{-3}`"
   ":math:`chlabs`", "chlabs(1:3)", "light absorption length per chla conc.", "(``chlabs_diatoms``, ``chlabs_sp``, ``chlabs_phaeo``)", "1\ :math:`/`\ m\ :math:`/`\ (mg\ :math:`/`\ m\ :math:`^{3}`)"
   
   ":math:`\alpha`", "alpha2max\_low(1:3)", "light limitation factor", "(``alpha2max_low_diatoms``, ``alpha2max_low_sp``, ``alpha2max_low_phaeo``)", "m\ :math:`^2`/W"
   ":math:`\beta`", "beta2max(1:3)", "light inhibition factor", "(``beta2max_diatoms``, ``beta2max_sp``, ``beta2max_phaeo``)", "m\ :math:`^2`/W"
   ":math:`\mu_{max}`", "mu\_max(1:3)", "maximum algal growth rate", "(``mu_max_diatoms``, ``mu_max_sp``, ``mu_max_phaeo``)", "day\ :math:`^{-1}`"
   ":math:`\mu_T`", "grow\_Tdep(1:3)", "temperature growth factor", "(``grow_Tdep_diatoms``, ``grow_Tdep_sp``, ``grow_Tdep_phaeo``)", "day\ :math:`^{-1}`"
   ":math:`R_{si:n}`", "R\_Si2N(1:3)", "algal silicate to nitrogen", "(``ratio_Si2N_diatoms``, ``ratio_Si2N_sp``, ``ratio_Si2N_phaeo``)", "mol/mol"
