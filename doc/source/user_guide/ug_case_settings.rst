:tocdepth: 3

.. _case_settings:

Case Settings
=====================

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
   "ICE_IOTYPE", "string", "unused", "none"
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
   "TRZS", "0,1", "zsalinity tracer, needs TRBRI=1", "0"
   "TRBGCS", "0,1", "skeletal layer tracer, needs TRBGCZ=0", "0"
   "TRBGCZ", "0,1", "zbgc tracers, needs TRBGCS=0 and TRBRI=1", "0"
   "NBGCLYR", "integer", "number of zbgc layers", "7"
   "TRZAERO", "0-6", "number of z aerosol tracers", "0"
   "TRALG", "0,1,2,3", "number of algal tracers", "0"
   "TRDOC", "0,1,2,3", "number of dissolved organic carbon", "0"
   "TRDIC", "0,1", "number of dissolved inorganic carbon", "0"
   "TRDON", "0,1", "number of dissolved organic nitrogen", "0"
   "TRFEP", "0,1,2", "number of particulate iron tracers", "0"
   "TRFED", "0,1,2", "number of dissolved iron tracers", "0"
   "ICE_BLDDEBUG", "true,false", "turn on compile debug flags", "false"
   "ICE_COVERAGE", "true,false", "turn on code coverage flags", "false"


.. _tabnamelist:

Table of Namelist Inputs
--------------------------

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
   "``ice_ic``", "``default``", "latitude and sst dependent initial condition", "``default``"
   "", "``none``", "no ice", ""
   "", "'path/file'", "restart file name", ""
   "``istep0``", "integer", "initial time step number", "0"
   "``ndtd``", "integer", "number of dynamics/advection/ridging/steps per thermo timestep", "1"
   "``npt``", "integer", "total number of time steps to take", "99999"
   "``restart``", "logical", "initialize using restart file", "``.false.``"
   "``restart_dir``", "string", "path to restart directory", "'./'"
   "``restart_file``", "string", "output file prefix for restart dump", "'iced'"
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
   "``tr_pond_cesm``", "logical", "CESM melt ponds", "``.false.``"
   "``tr_pond_lvl``", "logical", "level-ice melt ponds", "``.false.``"
   "``tr_pond_topo``", "logical", "topo melt ponds", "``.false.``"
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
   "``kitd``", "``0``", "delta function ITD approximation", "1"
   "", "``1``", "linear remapping ITD approximation", ""
   "``ksno``", "real", "snow thermal conductivity", "0.3"
   "``ktherm``", "``-1``", "thermodynamic model disabled", "1"
   "", "``0``", "zero-layer thermodynamic model", ""
   "", "``1``", "Bitz and Lipscomb thermodynamic model", ""
   "", "``2``", "mushy-layer thermodynamic model", ""
   "``phi_c_slow_mode``", ":math:`0<\phi_c < 1`", "critical liquid fraction", "0.05"
   "``phi_i_mushy``", ":math:`0<\phi_i < 1`", "solid fraction at lower boundary", "0.85"
   "``Rac_rapid_mode``", "real", "critical Rayleigh number", "10.0"
   "``sw_redist``", "logical", "shortwave redistribution", ".false."
   "``sw_frac``", "real", "fraction of shortwave redistribution", "0.9"
   "``sw_dtemp``", "real", "temperature from melt for sw_redist", "0.02"
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
   "``albedo_type``", "`ccsm3``", "NCAR CCSM3 albedo implementation", "``ccsm3``"
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
   "", "``dEdd``", "Delta-Eddington method", ""
   "", "", "", ""

ponds_nml
~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: **ponds_nml namelist options**
   :header: "variable", "options/format", "description", "default value"
   :widths: 15, 15, 30, 15 

   "", "", "", ""
   "``dpscale``", "real", "time scale for flushing in permeable ice", "1.0"
   "``frzpnd``", "``cesm``", "CESM pond refreezing forumulation", "``cesm``"
   "", "``hlid``", "Stefan refreezing with pond ice thickness", ""
   "``hp1``", "real", "critical ice lid thickness for topo ponds in m", "0.01"
   "``hs0``", "real", "snow depth of transition to bare sea ice in m", "0.03"
   "``hs1``", "real", "snow depth of transition to pond ice in m", "0.03"
   "``pndaspect``", "real", "aspect ratio of pond changes (depth:area)", "0.8"
   "``rfracmax``", ":math:`0 \le r_{max} \le 1`", "maximum melt water added to ponds", "0.85"
   "``rfracmin``", ":math:`0 \le r_{min} \le 1`", "minimum melt water added to ponds", "0.15"
   "", "", "", ""

forcing_nml
~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: **forcing_nml namelist options**
   :header: "variable", "options/format", "description", "default value"
   :widths: 15, 15, 30, 15 

   "", "", "", ""
   "``atmbndy``", "``constant``", "bulk transfer coefficients", "``default``"
   "", "``default``", "stability-based boundary layer", ""
   "``atmiter_conv``", "real", "convergence criteria for ustar", "0.0"
   "``atm_data_file``", "string", "file containing atmospheric data", "' '"
   "``atm_data_format``", "``bin``", "read direct access binary forcing files", "``bin``"
   "``atm_data_type``", "``clim``", "monthly climatology (see :ref:`force`)", "``default``"
   "", "``CFS``", "CFS model output  (see :ref:`force`)", ""
   "", "``default``", "constant values defined in the code", ""
   "", "``ISPOL``", "ISPOL experiment data  (see :ref:`force`)", ""
   "", "``NICE``", "N-ICE experiment data  (see :ref:`force`)", ""
   "``bgc_data_file``", "string", "file containing biogeochemistry data", "' '"
   "``bgc_data_format``", "``bin``", "read direct access binary forcing files", "``bin``"
   "``bgc_data_type``", "``clim``", "bgc climatological data", "``default``"
   "", "``default``", "constant values defined in the code", ""
   "", "``ncar``", "POP ocean forcing data", ""
   "``calc_strair``", "``.false.``", "read wind stress and speed from files", "``.true.``"
   "", "``.true.``", "calculate wind stress and speed", ""
   "``calc_Tsfc``", "logical", "calculate surface temperature", "``.true.``"
   "``data_dir``", "string", "path to forcing data directory", "' '"
   "``default_season``", "``summer``", "forcing initial summer values", "``winter``"
   "", "``winter``", "forcing initial winter values", ""
   "``emissivity``", "real", "emissivity of snow and ice", "0.95"
   "``fbot_xfer_type``", "``Cdn_ocn``", "variabler ocean heat transfer coefficient scheme", "``constant``"
   "", "``constant``", "constant ocean heat transfer coefficient", ""
   "``formdrag``", "logical", "calculate form drag", "``.false.``"
   "``fyear_init``", "integer", "first year of atmospheric forcing data", "1998"
   "``highfreq``", "logical", "high-frequency atmo coupling", "``.false.``"
   "``ice_data_file``", "string", "file containing ice opening, closing data", "' '"
   "``l_mpond_fresh``", "``.false.``", "release pond water immediately to ocean", "``.false.``"
   "", "``true``", "retain (topo) pond water until ponds drain", ""
   "``natmiter``", "integer", "number of atmo boundary layer iterations", "5"
   "``oceanmixed_ice``", "logical", "active ocean mixed layer calculation", "``.false.``"
   "``ocn_data_file``", "string", "file containing ocean data", "' ' "
   "``ocn_data_format``", "``bin``", "read direct access binary forcing files", "``bin``"
   "``ocn_data_type``", "``default``", "constant values defined in the code", "``default``"
   "", "``ISPOL``", "ISPOL experiment data  (see :ref:`force`)", ""
   "", "``NICE``", "N-ICE experiment data  (see :ref:`force`)", ""
   "", "``SHEBA``", "Opening/closing dataset from SHEBA", ""
   "``precip_units``", "``mks``", "liquid precipitation data units", "``mks``"
   "", "``mm_per_month``", "", ""
   "", "``mm_per_sec``", "(same as MKS units)", ""
   "", "``m_per_sec``", "", ""
   "``restore_ocn``", "logical", "restore sst to data", "``.false.``"
   "``tfrz_option``", "``linear_salt``", "linear functino of salinity (ktherm=1)", "``mushy``"
   "", "``minus1p8``", "constant ocean freezing temperature (:math:`-1.8^{\circ} C`)", ""
   "", "``mushy``", "matches mushy-layer thermo (ktherm=2)", ""
   "``trestore``", "integer", "sst restoring time scale (days)", "90"
   "``update_ocn_f``", "``.false.``", "do not include frazil water/salt fluxes in ocn fluxes", "``.false.``"
   "", "``true``", "include frazil water/salt fluxes in ocn fluxes", ""
   "``ustar_min``", "real", "minimum value of ocean friction velocity in m/s", "0.005"
   "``wave_spec_type``", "``constant``", "wave data file is provided, constant wave spectrum, for testing", "``none``"
   "", "``none``", "no wave data provided, no wave-ice interactions", ""
   "", "``profile``", "no wave data file is provided, use fixed dummy wave spectrum, for testing", ""
   "", "``random``", "wave data file is provided, wave spectrum generated using random number", ""
   "``ycycle``", "integer", "number of years in forcing data cycle", "1"
   "", "", "", ""

zbgc_nml
~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table:: **zbgc_nml namelist options**
   :header: "variable", "options/format", "description", "default value"
   :widths: 15, 15, 30, 15 

   "", "", "", ""
   "``algaltype_diatoms``", "real", "mobility type between stationary and mobile algal diatoms", "0.0"
   "``algaltype_phaeo``", "real", "mobility type between stationary and mobile algal phaeocystis", "0.5"
   "``algaltype_sp``", "real", "mobility type between stationary and mobile small plankton", "0.5"
   "``algal_vel``", "real", ":cite:`Lavoie05`", "1.11e-8"
   "``alpha2max_low_diatoms``", "real", "light limitation diatoms 1/(W/m^2)", "0.8"
   "``alpha2max_low_phaeo``", "real", "light limitation phaeocystis 1/(W/m^2)", "0.67"
   "``alpha2max_low_sp``", "real", "light limitation small plankton 1/(W/m^2)", "0.67"
   "``ammoniumtype``", "real", "mobility type between stationary and mobile ammonium", "1.0"
   "``beta2max_diatoms``", "real", "light inhibition diatoms 1/(W/m^2)", "0.18"
   "``beta2max_phaeo``", "real", "light inhibition phaeocystis 1/(W/m^2)", "0.01"
   "``beta2max_sp``", "real", "light inhibition small plankton 1/(W/m^2)", "0.0025"
   "``bgc_data_type``", "``clim``", "bgc climatological data", "``default``"
   "", "``default``", "constant values defined in the code", ""
   "", "``ncar``", "POP ocean forcing data", ""
   "``bgc_flux_type``", "``constant``", "constant ice–ocean flux velocity", "``Jin2006``"
   "", "``Jin2006``", "ice–ocean flux velocity of :cite:`Jin06`", ""
   "``chlabs_diatoms``", "real", "chl absorbtion diatoms 1/m/(mg/m^3)", "0.03"
   "``chlabs_phaeo``", "real", "chl absorbtion phaeocystis 1/m/(mg/m^3)", "0.05"
   "``chlabs_sp``", "real", "chl absorbtion small plankton 1/m/(mg/m^3)", "0.01"
   "``dEdd_algae``", "logical", "", "``.false.``"
   "``dmspdtype``", "real", "mobility type between stationary and mobile dmspd", "-1.0"
   "``dmspptype``", "real", "mobility type between stationary and mobile dmspp", "0.5"
   "``doctype_l``", "real", "mobility type between stationary and mobile doc lipids", "0.5"
   "``doctype_s``", "real", "mobility type between stationary and mobile doc saccharids", "0.5"
   "``dontype_protein``", "real", "mobility type between stationary and mobile don proteins", "0.5"
   "``dustFe_sol``", "real", "solubility fraction", "0.005"
   "``fedtype_1``", "real", "mobility type between stationary and mobile fed lipids", "0.5"
   "``feptype_1``", "real", "mobility type between stationary and mobile fep lipids", "0.5"
   "``frazil_scav``", "real", "increase in initial bio bracer from ocean scavenging", "1.0"
   "``fr_dFe``", "real", "fraction of remineralized nitrogen in units of algal iron", "0.3"
   "``fr_graze_diatoms``", "real", "fraction grazed diatoms", "0.01"
   "``fr_graze_e``", "real", "fraction of assimilation excreted", "0.5"
   "``fr_graze_phaeo``", "real", "fraction grazed phaeocystis", "0.1"
   "``fr_graze_s``", "real", "fraction of grazing spilled or slopped", "0.5"
   "``fr_graze_sp``", "real", "fraction grazed small plankton", "0.1"
   "``fr_mort2min``", "real", "fractionation of mortality to Am", "0.5"
   "``fr_resp``", "real", "frac of algal growth lost due to respiration", "0.05"
   "``fr_resp_s``", "real", "DMSPd fraction of respiration loss as DMSPd", "0.75"
   "``fsal``", "real", "salinity limitation ppt", "1.0"
   "``F_abs_chl_diatoms``", "real", "scales absorbed radiation for dEdd chl diatoms", "2.0"
   "``F_abs_chl_phaeo``", "real", "scales absorbed radiation for dEdd chl phaeocystis", "5.0"
   "``F_abs_chl_sp``", "real", "scales absorbed radiation for dEdd small plankton", "4.0"
   "``f_doc_l``", "real", "fraction of mortality to DOC lipids", "0.4"
   "``f_doc_s``", "real", "fraction of mortality to DOC saccharides", "0.4"
   "``f_don_Am_protein``", "real", "fraction of remineralized DON to ammonium", "0.25"
   "``f_don_protein``", "real", "fraction of spilled grazing to proteins", "0.6"
   "``f_exude_l``", "real", "fraction of exudation to DOC lipids", "1.0"
   "``f_exude_s``", "real", "fraction of exudation to DOC saccharids", "1.0"
   "``grid_o``", "real", "z biology for bottom flux", "5.0"
   "``grid_oS``", "real", "z salinity for bottom flux", "5.0"
   "``grow_Tdep_diatoms``", "real", "temperature dependence growth diatoms per degC", "0.06"
   "``grow_Tdep_phaeo``", "real", "temperature dependence growth phaeocystis per degC", "0.06"
   "``grow_Tdep_sp``", "real", "temperature dependence growth small plankton per degC", "0.06"
   "``humtype``", "real", "mobility type between stationary and mobile hum", "1.0"
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
   "``k_nitrif``", "real", "nitrification rate per day", "0.0"
   "``K_Nit_diatoms``", "real", "nitrate half saturation diatoms mmol/m^3", "1.0"
   "``K_Nit_phaeo``", "real", "nitrate half saturation phaeocystis mmol/m^3", "1.0"
   "``K_Nit_sp``", "real", "nitrate half saturation small plankton mmol/m^3", "1.0"
   "``K_Sil_diatoms``", "real", "silicate half saturation diatoms mmol/m^3", "4.0"
   "``K_Sil_phaeo``", "real", "silicate half saturation phaeocystis mmol/m^3", "0.0"
   "``K_Sil_sp``", "real", "silicate half saturation small plankton mmol/m^3", "0.0"
   "``kn_bac_protein``", "real", "bacterial degradation of DON per day", "0.03"
   "``l_sk``", "real", "characteristic diffusive scale in m", "7.0"
   "``l_skS``", "real", "z salinity characteristic diffusive scale in m", "7.0"
   "``max_dfe_doc1``", "real", "max ratio of dFe to saccharides in the ice in nm Fe / muM C", "0.2"
   "``max_loss``", "real", "restrict uptake to percent of remaining value", "0.9"
   "``modal_aero``", "logical", "modal aersols", "``.false.``"
   "``mort_pre_diatoms``", "real", "mortality diatoms", "0.007"
   "``mort_pre_phaeo``", "real", "mortality phaeocystis", "0.007"
   "``mort_pre_sp``", "real", "mortality small plankton", "0.007"
   "``mort_Tdep_diatoms``", "real", "temperature dependence of mortality diatoms per degC", "0.03"
   "``mort_Tdep_phaeo``", "real", "temperature dependence of mortality phaeocystis per degC", "0.03"
   "``mort_Tdep_sp``", "real", "temperature dependence of mortality small plankton per degC", "0.03"
   "``mu_max_diatoms``", "real", "maximum growth rate diatoms per day", "1.2"
   "``mu_max_phaeo``", "real", "maximum growth rate phaeocystis per day", "0.851"
   "``mu_max_sp``", "real", "maximum growth rate small plankton per day", "0.851"
   "``nitratetype``", "real", "mobility type between stationary and mobile nitrate", "-1.0"
   "``op_dep_min``", "real", "light attenuates for optical depths exceeding min", "0.1"
   "``phi_snow``", "real", "snow porosity for brine height tracer", "0.5"
   "``ratio_chl2N_diatoms``", "real", "algal chl to N in mg/mmol diatoms", "2.1"
   "``ratio_chl2N_phaeo``", "real", "algal chl to N in mg/mmol phaeocystis", "0.84"
   "``ratio_chl2N_sp``", "real", "algal chl to N in mg/mmol small plankton", "1.1"
   "``ratio_C2N_diatoms``", "real", "algal C to N in mol/mol diatoms", "7.0"
   "``ratio_C2N_phaeo``", "real", "algal C to N in mol/mol phaeocystis", "7.0"
   "``ratio_C2N_proteins``", "real", "algal C to N in mol/mol proteins", "7.0"
   "``ratio_C2N_sp``", "real", "algal C to N in mol/mol small plankton", "7.0"
   "``ratio_Fe2C_diatoms``", "real", "algal Fe to C in umol/mol diatoms", "0.0033"
   "``ratio_Fe2C_phaeo``", "real", "algal Fe to C in umol/mol phaeocystis", "1.0"
   "``ratio_Fe2C_sp``", "real", "algal Fe to C in umol/mol small plankton", "0.0033"
   "``ratio_Fe2N_diatoms``", "real", "algal Fe to N in umol/mol diatoms", "0.23"
   "``ratio_Fe2N_phaeo``", "real", "algal Fe to N in umol/mol phaeocystis", "0.7"
   "``ratio_Fe2N_sp``", "real", "algal Fe to N in umol/mol small plankton", "0.23"
   "``ratio_Fe2DOC_s``", "real", "Fe to C of DON saccharids nmol/umol", "1.0"
   "``ratio_Fe2DOC_l``", "real", "Fe to C of DOC lipids nmol/umol", "0.033"
   "``ratio_Fe2DON``", "real", "Fe to C of DON nmol/umol", "0.023"
   "``ratio_Si2N_diatoms``", "real", "algal Si to N in mol/mol diatoms", "1.8"
   "``ratio_Si2N_phaeo``", "real", "algal Si to N in mol/mol phaeocystis", "0.0"
   "``ratio_Si2N_sp``", "real", "algal Si to N in mol/mol small plankton", "0.0"
   "``ratio_S2N_diatoms``", "real", "algal S to N in mol/mol diatoms", "0.03"
   "``ratio_S2N_phaeo``", "real", "algal S to N in mol/mol phaeocystis", "0.03"
   "``ratio_S2N_sp``", "real", "algal S to N in mol/mol small plankton", "0.03"
   "``restore_bgc``", "logical", "restore bgc to data", "``.false.``"
   "``R_dFe2dust``", "real", "g/g :cite:`Tagliabue09`", "0.035"
   "``scale_bgc``", "logical", "", "``.false.``"
   "``silicatetype``", "real", "mobility type between stationary and mobile silicate", "-1.0"
   "``skl_bgc``", "logical", "biogeochemistry", "``.false.``"
   "``solve_zbgc``", "logical", "", "``.false.``"
   "``solve_zsal``", "logical", "update salinity tracer profile", "``.false.``"
   "``tau_max``", "real", "long time mobile to stationary exchanges", "1.73e-5"
   "``tau_min``", "real", "rapid module to stationary exchanges", "5200."
   "``tr_bgc_Am``", "logical", "ammonium tracer", "``.false.``"
   "``tr_bgc_C``", "logical", "algal carbon tracer", "``.false.``"
   "``tr_bgc_chl``", "logical", "algal chlorophyll tracer", "``.false.``"
   "``tr_bgc_DMS``", "logical", "DMS tracer", "``.false.``"
   "``tr_bgc_DON``", "logical", "DON tracer", "``.false.``"
   "``tr_bgc_Fe``", "logical", "iron tracer", "``.false.``"
   "``tr_bgc_hum``", "logical", "", "``.false.``"
   "``tr_bgc_Nit``", "logical", "", "``.false.``"
   "``tr_bgc_PON``", "logical", "PON tracer", "``.false.``"
   "``tr_bgc_Sil``", "logical", "silicate tracer", "``.false.``"
   "``tr_brine``", "logical", "brine height tracer", "``.false.``"
   "``tr_zaero``", "logical", "vertical aerosol tracers", "``.false.``"
   "``t_iron_conv``", "real", "desorption loss pFe to dFe in days", "3065."
   "``t_sk_conv``", "real", "Stefels conversion time in days", "3.0"
   "``t_sk_ox``", "real", "DMS oxidation time in days", "10.0"
   "``T_max``", "real", "maximum temperature degC", "0.0"
   "``y_sk_DMS``", "real", "fraction conversion given high yield", "0.5"
   "``zaerotype_bc1``", "real", "mobility type between stationary and mobile zaero bc1", "1.0"
   "``zaerotype_bc2``", "real", "mobility type between stationary and mobile zaero bc2", "1.0"
   "``zaerotype_dust1``", "real", "mobility type between stationary and mobile zaero dust1", "1.0"
   "``zaerotype_dust2``", "real", "mobility type between stationary and mobile zaero dust2", "1.0"
   "``zaerotype_dust3``", "real", "mobility type between stationary and mobile zaero dust3", "1.0"
   "``zaerotype_dust4``", "real", "mobility type between stationary and mobile zaero dust4", "1.0"
   "``z_tracers``", "logical", "", "``.false.``"
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
..   "``restart_zsal``", "logical", "", "``.false.``"

* = If Icepack is run stand-alone and wave_spec_type is not set to none, then a fixed wave spectrum is defined in the code to use for testing. As with other input data, this spectrum should not be used for production runs or publications.
  
.. _tuning:

BGC Tuning Parameters
------------------------

Biogeochemical tuning parameters are specified as namelist options in
**icepack\_in**. Table :ref:`tab-bio-tracers2` provides a list of parameters
used in the reaction equations, their representation in the code, a
short description of each and the default values. Please keep in mind
that there has only been minimal tuning of the model.

.. _tab-bio-tracers2:

.. csv-table:: *Biogeochemical Reaction Parameters*
   :header: "Text Variable", "Variable in code", "Description", "Value", "units"
   :widths: 7, 20, 15, 15, 15

   ":math:`f_{graze}`", "fr\_graze(1:3)", "fraction of growth grazed", "0, 0.1, 0.1", "1"
   ":math:`f_{res}`", "fr\_resp", "fraction of growth respired", "0.05", "1"
   ":math:`l_{max}`", "max\_loss", "maximum tracer loss fraction", "0.9", "1"
   ":math:`m_{pre}`", "mort\_pre(1:3)", "maximum mortality rate", "0.007, 0.007, 0.007", "day\ :math:`^{-1}`"
   ":math:`m_{T}`", "mort\_Tdep(1:3)", "mortality temperature decay", "0.03, 0.03, 0.03", ":math:`^o`\ C\ :math:`^{-1}`"
   ":math:`T_{max}`", "T\_max", "maximum brine temperature", "0", ":math:`^o`\ C"
   ":math:`k_{nitr}`", "k\_nitrif", "nitrification rate", "0", "day\ :math:`^{-1}`"
   ":math:`f_{ng}`", "fr\_graze\_e", "fraction of grazing excreted", "0.5", "1"
   ":math:`f_{gs}`", "fr\_graze\_s", "fraction of grazing spilled", "0.5", "1"
   ":math:`f_{nm}`", "fr\_mort2min", "fraction of mortality to :math:`{\mbox{NH$_4$}}`", "0.5", "1"
   ":math:`f_{dg}`", "f\_don", "frac. spilled grazing to :math:`{\mbox{DON}}`", "0.6", "1"
   ":math:`k_{nb}`", "kn\_bac :math:`^a`", "bacterial degradation of :math:`{\mbox{DON}}`", "0.03", "day\ :math:`^{-1}`"
   ":math:`f_{cg}`", "f\_doc(1:3)", "fraction of mortality to :math:`{\mbox{DOC}}`", "0.4, 0.4, 0.2 ", "1"
   ":math:`R_{c:n}^c`", "R\_C2N(1:3)", "algal carbon to nitrogen ratio", "7.0, 7.0, 7.0", "mol/mol"
   ":math:`k_{cb}`", "k\_bac1:3\ :math:`^a`", "bacterial degradation of DOC", "0.03, 0.03, 0.03", "day\ :math:`^{-1}`"
   ":math:`\tau_{fe}`", "t\_iron\_conv", "conversion time pFe :math:`\leftrightarrow` dFe", "3065.0 ", "day"
   ":math:`r^{max}_{fed:doc}`", "max\_dfe\_doc1", "max ratio of dFe to saccharids", "0.1852", "nM Fe\ :math:`/\mu`\ M C"
   ":math:`f_{fa}`", "fr\_dFe  ", "fraction of remin. N to dFe", "0.3", "1"
   ":math:`R_{fe:n}`", "R\_Fe2N(1:3)", "algal Fe to N ratio", "0.023, 0.023, 0.7", "mmol/mol"
   ":math:`R_{s:n}`", "R\_S2N(1:3)", "algal S to N ratio", "0.03, 0.03, 0.03", "mol/mol"
   ":math:`f_{sr}`", "fr\_resp\_s", "resp. loss as DMSPd", "0.75", "1"
   ":math:`\tau_{dmsp}`", "t\_sk\_conv", "Stefels rate", "3.0", "day"
   ":math:`\tau_{dms}`", "t\_sk\_ox", "DMS oxidation rate", "10.0", "day"
   ":math:`y_{dms}`", "y\_sk\_DMS", "yield for DMS conversion", "0.5", "1"
   ":math:`K_{{\mbox{NO$_3$}}}`", "K\_Nit(1:3)", ":math:`{\mbox{NO$_3$}}` half saturation constant", "1,1,1", "mmol/m\ :math:`^{3}`"
   ":math:`K_{{\mbox{NH$_4$}}}`", "K\_Am(1:3)", ":math:`{\mbox{NH$_4$}}` half saturation constant", "0.3, 0.3, 0.3", "mmol/m\ :math:`^{-3}`"
   ":math:`K_{{\mbox{SiO$_3$}}}`", "K\_Sil(1:3)", "silicate half saturation constant", "4.0, 0, 0", "mmol/m\ :math:`^{-3}`"
   ":math:`K_{{\mbox{fed}}}`", "K\_Fe(1:3)", "iron half saturation constant", "1.0, 0.2, 0.1", ":math:`\mu`\ mol/m\ :math:`^{-3}`"
   ":math:`op_{min}`", "op\_dep\_min", "boundary for light attenuation", "0.1", "1"
   ":math:`chlabs`", "chlabs(1:3)", "light absorption length per chla conc.", "0.03, 0.01, 0.05", "1\ :math:`/`\ m\ :math:`/`\ (mg\ :math:`/`\ m\ :math:`^{3}`)"
   ":math:`\alpha`", "alpha2max\_low(1:3)", "light limitation factor", "0.25, 0.25, 0.25", "m\ :math:`^2`/W"
   ":math:`\beta`", "beta2max(1:3)", "light inhibition factor", "0.018, 0.0025, 0.01", "m\ :math:`^2`/W"
   ":math:`\mu_{max}`", "mu\_max(1:3)", "maximum algal growth rate", "1.44, 0.851, 0.851", "day\ :math:`^{-1}`"
   ":math:`\mu_T`", "grow\_Tdep(1:3)", "temperature growth factor", "0.06, 0.06, 0.06", "day\ :math:`^{-1}`"
   ":math:`f_{sal}`", "fsal", "salinity growth factor", "1", "1"
   ":math:`R_{si:n}`", "R\_Si2N(1:3)", "algal silicate to nitrogen", "1.8, 0, 0", "mol/mol"

:math:`^a` only (1:2) of DOC and DOC parameters have physical meaning
