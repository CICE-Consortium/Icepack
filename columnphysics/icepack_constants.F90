!  SVN:$Id: icepack_constants.F90 1226 2017-05-22 22:45:03Z tcraig $
!=======================================================================
!
! This module defines a variety of physical and numerical constants
! used in the column package
!
! author Elizabeth C. Hunke, LANL

      module icepack_constants

      use icepack_kinds
      use icepack_warnings, only: warnstr, icepack_warnings_add
      use icepack_warnings, only: icepack_warnings_setabort, icepack_warnings_aborted

      implicit none
      private

      public :: icepack_init_constants
      public :: icepack_query_constants
      public :: icepack_write_constants
      public :: icepack_recompute_constants

      !-----------------------------------------------------------------
      ! parameter constants
      !-----------------------------------------------------------------

      integer (kind=int_kind), parameter, public :: & 
         nspint = 3                ! number of solar spectral intervals
                    
      real (kind=dbl_kind), parameter, public :: &
         secday = 86400.0_dbl_kind ,&! seconds in calendar day
         c0   = 0.0_dbl_kind, &
         c1   = 1.0_dbl_kind, &
         c1p5 = 1.5_dbl_kind, &
         c2   = 2.0_dbl_kind, &
         c3   = 3.0_dbl_kind, &
         c4   = 4.0_dbl_kind, &
         c5   = 5.0_dbl_kind, &
         c6   = 6.0_dbl_kind, &
         c8   = 8.0_dbl_kind, &
         c10  = 10.0_dbl_kind, &
         c15  = 15.0_dbl_kind, &
         c16  = 16.0_dbl_kind, &
         c20  = 20.0_dbl_kind, &
         c25  = 25.0_dbl_kind, &
         c100 = 100.0_dbl_kind, &
         c1000= 1000.0_dbl_kind, &
         p001 = 0.001_dbl_kind, &
         p01  = 0.01_dbl_kind, &
         p1   = 0.1_dbl_kind, &
         p2   = 0.2_dbl_kind, &
         p4   = 0.4_dbl_kind, &
         p5   = 0.5_dbl_kind, &
         p6   = 0.6_dbl_kind, &
         p05  = 0.05_dbl_kind, &
         p15  = 0.15_dbl_kind, &
         p25  = 0.25_dbl_kind, &
         p75  = 0.75_dbl_kind, &
         p333 = c1/c3, &
         p666 = c2/c3, &
         spval_const= -1.0e36_dbl_kind

      !-----------------------------------------------------------------
      ! derived physical constants
      !    Lfresh = Lsub-Lvap     ,&! latent heat of melting of fresh ice (J/kg)
      !    cprho  = cp_ocn*rhow   ,&! for ocean mixed layer (J kg / K m^3)
      !    Cp     = 0.5_dbl_kind*gravit*(rhow-rhoi)*rhoi/rhow ,&! proport const for PE 
      !-----------------------------------------------------------------

      real (kind=dbl_kind), public :: &
         pih    = spval_const, &
         piq    = spval_const, &
         pi2    = spval_const, &
         Lfresh = spval_const, &! latent heat of melting of fresh ice (J/kg)
         cprho  = spval_const, &! for ocean mixed layer (J kg / K m^3)
         Cp     = spval_const   ! proport const for PE 

      !-----------------------------------------------------------------
      ! settable physical constants
      !-----------------------------------------------------------------

      real (kind=dbl_kind), public :: &
         rhos      = 330.0_dbl_kind   ,&! density of snow (kg/m^3)
         rhoi      = 917.0_dbl_kind   ,&! density of ice (kg/m^3)
         rhow      = 1026.0_dbl_kind  ,&! density of seawater (kg/m^3)
         cp_air    = 1005.0_dbl_kind  ,&! specific heat of air (J/kg/K)
         ! (Briegleb JGR 97 11475-11485  July 1992)
         emissivity= 0.95_dbl_kind    ,&! emissivity of snow and ice
         cp_ice    = 2106._dbl_kind   ,&! specific heat of fresh ice (J/kg/K)
         cp_ocn    = 4218._dbl_kind   ,&! specific heat of ocn    (J/kg/K)
                                        ! freshwater value needed for enthalpy
         depressT  = 0.054_dbl_kind   ,&! Tf:brine salinity ratio (C/ppt)
         dragio    = 0.00536_dbl_kind ,&! ice-ocn drag coefficient
         albocn    = 0.06_dbl_kind    ,&! ocean albedo
         gravit    = 9.80616_dbl_kind ,&! gravitational acceleration (m/s^2)
         viscosity_dyn = 1.79e-3_dbl_kind, & ! dynamic viscosity of brine (kg/m/s)
         Tocnfrz   = -1.8_dbl_kind    ,&! freezing temp of seawater (C),
                                        ! used as Tsfcn for open water
         rhofresh  = 1000.0_dbl_kind  ,&! density of fresh water (kg/m^3)
         zvir      = 0.606_dbl_kind   ,&! rh2o/rair - 1.0
         vonkar    = 0.4_dbl_kind     ,&! von Karman constant
         cp_wv     = 1.81e3_dbl_kind  ,&! specific heat of water vapor (J/kg/K)
         stefan_boltzmann = 567.0e-10_dbl_kind,&!  W/m^2/K^4
         Tffresh   = 273.15_dbl_kind  ,&! freezing temp of fresh ice (K)
         Lsub      = 2.835e6_dbl_kind ,&! latent heat, sublimation freshwater (J/kg)
         Lvap      = 2.501e6_dbl_kind ,&! latent heat, vaporization freshwater (J/kg)
         Timelt    = 0.0_dbl_kind     ,&! melting temperature, ice top surface  (C)
         Tsmelt    = 0.0_dbl_kind     ,&! melting temperature, snow top surface (C)
         ice_ref_salinity = 4._dbl_kind ,&! (ppt)
         ! ocn_ref_salinity = 34.7_dbl_kind,&! (ppt)
         iceruf   = 0.0005_dbl_kind   ,&! ice surface roughness (m)

         ! for ice strength
         Cf       = 17._dbl_kind      ,&! ratio of ridging work to PE change in ridging 
         Pstar    = 2.75e4_dbl_kind   ,&! constant in Hibler strength formula 
                                        ! (kstrength = 0) 
         Cstar    = 20._dbl_kind      ,&! constant in Hibler strength formula 
                                        ! (kstrength = 0) 

         ! (Ebert, Schramm and Curry JGR 100 15965-15975 Aug 1995)
         kappav = 1.4_dbl_kind ,&! vis extnctn coef in ice, wvlngth<700nm (1/m)
         !kappan = 17.6_dbl_kind,&! vis extnctn coef in ice, wvlngth<700nm (1/m)

         ! kice is not used for mushy thermo
         kice   = 2.03_dbl_kind  ,&! thermal conductivity of fresh ice(W/m/deg)
         ! kseaice is used only for zero-layer thermo
         kseaice= 2.00_dbl_kind  ,&! thermal conductivity of sea ice (W/m/deg)
                                   ! (used in zero layer thermodynamics option)
         ksno   = 0.30_dbl_kind  ,&! thermal conductivity of snow  (W/m/deg)
         zref   = 10._dbl_kind   ,&! reference height for stability (m)
         hs_min = 1.e-4_dbl_kind ,&! min snow thickness for computing zTsn (m)
         snowpatch = 0.02_dbl_kind, &  ! parameter for fractional snow area (m)
         rhosi     = 940.0_dbl_kind, & ! average sea ice density
                                       ! Cox and Weeks, 1982: 919-974 kg/m^2
         sk_l      = 0.03_dbl_kind, &  ! skeletal layer thickness (m)

         ! from parameters
         saltmax = 3.2_dbl_kind  , & ! max salinity at ice base for BL99 (ppt)
         ! phi_init and dSin0_frazil are used for mushy thermo, ktherm=2
         phi_init = 0.75_dbl_kind, & ! initial liquid fraction of frazil
         min_salin = p1          , & ! threshold for brine pocket treatment
         salt_loss =0.4_dbl_kind , & ! fraction of salt retained in zsalinity
         min_bgc  = 0.01_dbl_kind, & ! fraction of ocean bgc concentration in surface melt
         dSin0_frazil = c3,        & ! bulk salinity reduction of newly formed frazil
         hi_ssl = 0.050_dbl_kind,  & ! ice surface scattering layer thickness (m)
         hs_ssl = 0.040_dbl_kind     ! snow surface scattering layer thickness (m)

      ! weights for albedos 
      ! 4 Jan 2007 BPB  Following are appropriate for complete cloud
      ! in a summer polar atmosphere with 1.5m bare sea ice surface:
      ! .636/.364 vis/nir with only 0.5% direct for each band.
      real (kind=dbl_kind), public :: &           ! currently used only
         awtvdr = 0.00318_dbl_kind, &! visible, direct  ! for history and
         awtidr = 0.00182_dbl_kind, &! near IR, direct  ! diagnostics
         awtvdf = 0.63282_dbl_kind, &! visible, diffuse
         awtidf = 0.36218_dbl_kind   ! near IR, diffuse

      real (kind=dbl_kind), public :: &
         qqqice  = 11637800._dbl_kind   ,&! for qsat over ice
         TTTice  = 5897.8_dbl_kind      ,&! for qsat over ice
         qqqocn  = 627572.4_dbl_kind    ,&! for qsat over ocn
         TTTocn  = 5107.4_dbl_kind        ! for qsat over ocn
    
      real (kind=dbl_kind), public :: &
        puny   = 1.0e-11_dbl_kind, &
        bignum = 1.0e+30_dbl_kind, &
        pi     = 3.14159265358979323846_dbl_kind

!=======================================================================
      contains
!=======================================================================

      subroutine icepack_init_constants( &
         rhos_in, rhoi_in, rhow_in, cp_air_in, emissivity_in, &
         cp_ice_in, cp_ocn_in, &
         depressT_in, dragio_in, albocn_in, gravit_in, viscosity_dyn_in, &
         Tocnfrz_in, rhofresh_in, zvir_in, vonkar_in, cp_wv_in, &
         stefan_boltzmann_in, ice_ref_salinity_in, &
         Tffresh_in, Lsub_in, Lvap_in, Timelt_in, Tsmelt_in, &
         iceruf_in, Cf_in, Pstar_in, Cstar_in, kappav_in, &
         kice_in, kseaice_in, ksno_in, &
         zref_in, hs_min_in, snowpatch_in, rhosi_in, sk_l_in, &
         saltmax_in, phi_init_in, min_salin_in, salt_loss_in, &
         min_bgc_in, dSin0_frazil_in, hi_ssl_in, hs_ssl_in, &
         awtvdr_in, awtidr_in, awtvdf_in, awtidf_in, &
         qqqice_in, TTTice_in, qqqocn_in, TTTocn_in, &
         puny_in, bignum_in, pi_in )

      real (kind=dbl_kind), intent(in), optional :: &
         rhos_in,       & ! density of snow (kg/m^3)
         rhoi_in,       & ! density of ice (kg/m^3)
         rhow_in,       & ! density of seawater (kg/m^3)
         cp_air_in,     & ! specific heat of air (J/kg/K)
         emissivity_in, & ! emissivity of snow and ice
         cp_ice_in,     & ! specific heat of fresh ice (J/kg/K)
         cp_ocn_in,     & ! specific heat of ocn    (J/kg/K)
         depressT_in,   & ! Tf:brine salinity ratio (C/ppt)
         dragio_in,     & ! ice-ocn drag coefficient
         albocn_in,     & ! ocean albedo
         gravit_in,     & ! gravitational acceleration (m/s^2)
         viscosity_dyn_in, & ! dynamic viscosity of brine (kg/m/s)
         Tocnfrz_in,    & ! freezing temp of seawater (C)
         rhofresh_in,   & ! density of fresh water (kg/m^3)
         zvir_in,       & ! rh2o/rair - 1.0
         vonkar_in,     & ! von Karman constant
         cp_wv_in,      & ! specific heat of water vapor (J/kg/K)
         stefan_boltzmann_in, & !  W/m^2/K^4
         Tffresh_in,    & ! freezing temp of fresh ice (K)
         Lsub_in,       & ! latent heat, sublimation freshwater (J/kg)
         Lvap_in,       & ! latent heat, vaporization freshwater (J/kg)
         Timelt_in,     & ! melting temperature, ice top surface  (C)
         Tsmelt_in,     & ! melting temperature, snow top surface (C)
         ice_ref_salinity_in, & ! (ppt)
         iceruf_in,     & ! ice surface roughness (m)
         Cf_in,         & ! ratio of ridging work to PE change in ridging 
         Pstar_in,      & ! constant in Hibler strength formula 
         Cstar_in,      & ! constant in Hibler strength formula 
         kappav_in,     & ! vis extnctn coef in ice, wvlngth<700nm (1/m)
         kice_in,       & ! thermal conductivity of fresh ice(W/m/deg)
         kseaice_in,    & ! thermal conductivity of sea ice (W/m/deg)
         ksno_in,       & ! thermal conductivity of snow  (W/m/deg)
         zref_in,       & ! reference height for stability (m)
         hs_min_in,     & ! min snow thickness for computing zTsn (m)
         snowpatch_in,  & ! parameter for fractional snow area (m)
         rhosi_in,      & ! average sea ice density (kg/m2)
         sk_l_in,       & ! skeletal layer thickness (m)
         saltmax_in,    & ! max salinity at ice base for BL99 (ppt)
         phi_init_in,   & ! initial liquid fraction of frazil
         min_salin_in,  & ! threshold for brine pocket treatment
         salt_loss_in,  & ! fraction of salt retained in zsalinity
         min_bgc_in,    & ! fraction of ocean bgc concentration in surface melt
         dSin0_frazil_in, & ! bulk salinity reduction of newly formed frazil
         hi_ssl_in,     & ! ice surface scattering layer thickness (m)
         hs_ssl_in,     & ! visible, direct 
         awtvdr_in,     & ! visible, direct  ! for history and
         awtidr_in,     & ! near IR, direct  ! diagnostics
         awtvdf_in,     & ! visible, diffuse
         awtidf_in,     & ! near IR, diffuse
         qqqice_in,     & ! for qsat over ice
         TTTice_in,     & ! for qsat over ice
         qqqocn_in,     & ! for qsat over ocn
         TTTocn_in,     & ! for qsat over ocn
         puny_in,       & !
         bignum_in,     & !
         pi_in            !

      character(len=*),parameter :: subname='(icepack_init_constants)'

         if (present(rhos_in))       rhos = rhos_in
         if (present(rhoi_in))       rhoi = rhoi_in
         if (present(rhow_in))       rhow = rhow_in
         if (present(cp_air_in))     cp_air = cp_air_in
         if (present(emissivity_in)) emissivity = emissivity_in
         if (present(cp_ice_in))     cp_ice = cp_ice_in
         if (present(cp_ocn_in))     cp_ocn = cp_ocn_in
         if (present(depressT_in))   depressT = depressT_in
         if (present(dragio_in))     dragio = dragio_in
         if (present(albocn_in))     albocn = albocn_in
         if (present(gravit_in))     gravit = gravit_in
         if (present(viscosity_dyn_in)) viscosity_dyn = viscosity_dyn_in
         if (present(Tocnfrz_in))    Tocnfrz = Tocnfrz_in
         if (present(rhofresh_in))   rhofresh = rhofresh_in
         if (present(zvir_in))       zvir   = zvir_in
         if (present(vonkar_in))     vonkar = vonkar_in
         if (present(cp_wv_in))      cp_wv  = cp_wv_in
         if (present(stefan_boltzmann_in)) stefan_boltzmann = stefan_boltzmann_in
         if (present(Tffresh_in))    Tffresh = Tffresh_in
         if (present(Lsub_in))       Lsub = Lsub_in
         if (present(Lvap_in))       Lvap = Lvap_in
         if (present(Timelt_in))     Timelt = Timelt_in
         if (present(Tsmelt_in))     Tsmelt = Tsmelt_in
         if (present(ice_ref_salinity_in)) ice_ref_salinity = ice_ref_salinity_in
         if (present(iceruf_in))     iceruf = iceruf_in
         if (present(Cf_in))         Cf     = Cf_in
         if (present(Pstar_in))      Pstar  = Pstar_in
         if (present(Cstar_in))      Cstar  = Cstar_in
         if (present(kappav_in))     kappav = kappav_in
         if (present(kice_in))       kice   = kice_in
         if (present(kseaice_in))    kseaice = kseaice_in
         if (present(ksno_in))       ksno   = ksno_in
         if (present(zref_in))       zref   = zref_in
         if (present(hs_min_in))     hs_min = hs_min_in
         if (present(snowpatch_in))  snowpatch = snowpatch_in
         if (present(rhosi_in))      rhosi  = rhosi_in
         if (present(sk_l_in))       sk_l   = sk_l_in
         if (present(saltmax_in))    saltmax = saltmax_in
         if (present(phi_init_in))   phi_init = phi_init_in
         if (present(min_salin_in))  min_salin = min_salin_in
         if (present(salt_loss_in))  salt_loss = salt_loss_in
         if (present(min_bgc_in))    min_bgc = min_bgc_in
         if (present(dSin0_frazil_in)) dSin0_frazil = dSin0_frazil_in
         if (present(hi_ssl_in))     hi_ssl = hi_ssl_in
         if (present(hs_ssl_in))     hs_ssl = hs_ssl_in
         if (present(awtvdr_in))     awtvdr = awtvdr_in
         if (present(awtidr_in))     awtidr = awtidr_in
         if (present(awtvdf_in))     awtvdf = awtvdf_in
         if (present(awtidf_in))     awtidf = awtidf_in
         if (present(qqqice_in))     qqqice = qqqice_in
         if (present(TTTice_in))     TTTice = TTTice_in
         if (present(qqqocn_in))     qqqocn = qqqocn_in
         if (present(TTTocn_in))     TTTocn = TTTocn_in
         if (present(puny_in))       puny   = puny_in
         if (present(bignum_in))     bignum = bignum_in
         if (present(pi_in))         pi     = pi_in

         call icepack_recompute_constants()
         if (icepack_warnings_aborted(subname)) return

      end subroutine icepack_init_constants

!=======================================================================

      subroutine icepack_query_constants( &
         rhos_out, rhoi_out, rhow_out, cp_air_out, emissivity_out, &
         cp_ice_out, cp_ocn_out, &
         depressT_out, dragio_out, albocn_out, gravit_out, viscosity_dyn_out, &
         Tocnfrz_out, rhofresh_out, zvir_out, vonkar_out, cp_wv_out, &
         stefan_boltzmann_out, ice_ref_salinity_out, &
         Tffresh_out, Lsub_out, Lvap_out, Timelt_out, Tsmelt_out, &
         iceruf_out, Cf_out, Pstar_out, Cstar_out, kappav_out, &
         kice_out, kseaice_out, ksno_out, &
         zref_out, hs_min_out, snowpatch_out, rhosi_out, sk_l_out, &
         saltmax_out, phi_init_out, min_salin_out, salt_loss_out, &
         min_bgc_out, dSin0_frazil_out, hi_ssl_out, hs_ssl_out, &
         awtvdr_out, awtidr_out, awtvdf_out, awtidf_out, &
         qqqice_out, TTTice_out, qqqocn_out, TTTocn_out, &
         Lfresh_out, cprho_out, Cp_out, &
         puny_out, bignum_out, pi_out, &
         secday_out, c0_out, c1_out, c1p5_out, c2_out, c3_out, c4_out, &
         c5_out, c6_out, c8_out, c10_out, c15_out, c16_out, c20_out, &
         c25_out, c100_out, c1000_out, p001_out, p01_out, p1_out, &
         p2_out, p4_out, p5_out, p6_out, p05_out, p15_out, p25_out, p75_out, &
         p333_out, p666_out, spval_const_out, pih_out, piq_out, pi2_out)


      real (kind=dbl_kind), intent(out), optional :: &
         rhos_out,       & ! density of snow (kg/m^3)
         rhoi_out,       & ! density of ice (kg/m^3)
         rhow_out,       & ! density of seawater (kg/m^3)
         cp_air_out,     & ! specific heat of air (J/kg/K)
         emissivity_out, & ! emissivity of snow and ice
         cp_ice_out,     & ! specific heat of fresh ice (J/kg/K)
         cp_ocn_out,     & ! specific heat of ocn    (J/kg/K)
         depressT_out,   & ! Tf:brine salinity ratio (C/ppt)
         dragio_out,     & ! ice-ocn drag coefficient
         albocn_out,     & ! ocean albedo
         gravit_out,     & ! gravitational acceleration (m/s^2)
         viscosity_dyn_out, & ! dynamic viscosity of brine (kg/m/s)
         Tocnfrz_out,    & ! freezing temp of seawater (C)
         rhofresh_out,   & ! density of fresh water (kg/m^3)
         zvir_out,       & ! rh2o/rair - 1.0
         vonkar_out,     & ! von Karman constant
         cp_wv_out,      & ! specific heat of water vapor (J/kg/K)
         stefan_boltzmann_out, & !  W/m^2/K^4
         Tffresh_out,    & ! freezing temp of fresh ice (K)
         Lsub_out,       & ! latent heat, sublimation freshwater (J/kg)
         Lvap_out,       & ! latent heat, vaporization freshwater (J/kg)
         Timelt_out,     & ! melting temperature, ice top surface  (C)
         Tsmelt_out,     & ! melting temperature, snow top surface (C)
         ice_ref_salinity_out, & ! (ppt)
         iceruf_out,     & ! ice surface roughness (m)
         Cf_out,         & ! ratio of ridging work to PE change in ridging 
         Pstar_out,      & ! constant in Hibler strength formula 
         Cstar_out,      & ! constant in Hibler strength formula 
         kappav_out,     & ! vis extnctn coef in ice, wvlngth<700nm (1/m)
         kice_out,       & ! thermal conductivity of fresh ice(W/m/deg)
         kseaice_out,    & ! thermal conductivity of sea ice (W/m/deg)
         ksno_out,       & ! thermal conductivity of snow  (W/m/deg)
         zref_out,       & ! reference height for stability (m)
         hs_min_out,     & ! min snow thickness for computing zTsn (m)
         snowpatch_out,  & ! parameter for fractional snow area (m)
         rhosi_out,      & ! average sea ice density (kg/m2)
         sk_l_out,       & ! skeletal layer thickness (m)
         saltmax_out,    & ! max salinity at ice base for BL99 (ppt)
         phi_init_out,   & ! initial liquid fraction of frazil
         min_salin_out,  & ! threshold for brine pocket treatment
         salt_loss_out,  & ! fraction of salt retained in zsalinity
         min_bgc_out,   & ! fraction of ocean bgc concentration in surface melt
         dSin0_frazil_out, & ! bulk salinity reduction of newly formed frazil
         hi_ssl_out,     & ! ice surface scattering layer thickness (m)
         hs_ssl_out,     & ! visible, direct 
         awtvdr_out,     & ! visible, direct  ! for history and
         awtidr_out,     & ! near IR, direct  ! diagnostics
         awtvdf_out,     & ! visible, diffuse
         awtidf_out,     & ! near IR, diffuse
         qqqice_out,     & ! for qsat over ice
         TTTice_out,     & ! for qsat over ice
         qqqocn_out,     & ! for qsat over ocn
         TTTocn_out,     & ! for qsat over ocn
         Lfresh_out,     & ! latent heat of melting of fresh ice (J/kg)
         cprho_out,      & ! for ocean mixed layer (J kg / K m^3)
         Cp_out,         & ! proport const for PE 
         puny_out,       & !
         bignum_out,     & !
         pi_out,         & !
         secday_out, c0_out, c1_out, c1p5_out, c2_out, c3_out, c4_out, &
         c5_out, c6_out, c8_out, c10_out, c15_out, c16_out, c20_out, &
         c25_out, c100_out, c1000_out, p001_out, p01_out, p1_out, &
         p2_out, p4_out, p5_out, p6_out, p05_out, p15_out, p25_out, p75_out, &
         p333_out, p666_out, spval_const_out, pih_out, piq_out, pi2_out

      character(len=*),parameter :: subname='(icepack_query_constants)'

         if (present(rhos_out))       rhos_out = rhos
         if (present(rhoi_out))       rhoi_out = rhoi
         if (present(rhow_out))       rhow_out = rhow
         if (present(cp_air_out))     cp_air_out = cp_air
         if (present(emissivity_out)) emissivity_out = emissivity
         if (present(cp_ice_out))     cp_ice_out = cp_ice
         if (present(cp_ocn_out))     cp_ocn_out = cp_ocn
         if (present(depressT_out))   depressT_out = depressT
         if (present(dragio_out))     dragio_out = dragio
         if (present(albocn_out))     albocn_out = albocn
         if (present(gravit_out))     gravit_out = gravit
         if (present(viscosity_dyn_out)) viscosity_dyn_out = viscosity_dyn
         if (present(Tocnfrz_out))    Tocnfrz_out = Tocnfrz
         if (present(rhofresh_out))   rhofresh_out = rhofresh
         if (present(zvir_out))       zvir_out   = zvir
         if (present(vonkar_out))     vonkar_out = vonkar
         if (present(cp_wv_out))      cp_wv_out  = cp_wv
         if (present(stefan_boltzmann_out)) stefan_boltzmann_out = stefan_boltzmann
         if (present(Tffresh_out))    Tffresh_out = Tffresh
         if (present(Lsub_out))       Lsub_out = Lsub
         if (present(Lvap_out))       Lvap_out = Lvap
         if (present(Timelt_out))     Timelt_out = Timelt
         if (present(Tsmelt_out))     Tsmelt_out = Tsmelt
         if (present(ice_ref_salinity_out)) ice_ref_salinity_out = ice_ref_salinity
         if (present(iceruf_out))     iceruf_out = iceruf
         if (present(Cf_out))         Cf_out     = Cf
         if (present(Pstar_out))      Pstar_out  = Pstar
         if (present(Cstar_out))      Cstar_out  = Cstar
         if (present(kappav_out))     kappav_out = kappav
         if (present(kice_out))       kice_out   = kice
         if (present(kseaice_out))    kseaice_out = kseaice
         if (present(ksno_out))       ksno_out   = ksno
         if (present(zref_out))       zref_out   = zref
         if (present(hs_min_out))     hs_min_out = hs_min
         if (present(snowpatch_out))  snowpatch_out = snowpatch
         if (present(rhosi_out))      rhosi_out  = rhosi
         if (present(sk_l_out))       sk_l_out   = sk_l
         if (present(saltmax_out))    saltmax_out = saltmax
         if (present(phi_init_out))   phi_init_out = phi_init
         if (present(min_salin_out))  min_salin_out = min_salin
         if (present(salt_loss_out))  salt_loss_out = salt_loss
         if (present(min_bgc_out))    min_bgc_out = min_bgc
         if (present(dSin0_frazil_out)) dSin0_frazil_out = dSin0_frazil
         if (present(hi_ssl_out))     hi_ssl_out = hi_ssl
         if (present(hs_ssl_out))     hs_ssl_out = hs_ssl
         if (present(awtvdr_out))     awtvdr_out = awtvdr
         if (present(awtidr_out))     awtidr_out = awtidr
         if (present(awtvdf_out))     awtvdf_out = awtvdf
         if (present(awtidf_out))     awtidf_out = awtidf
         if (present(qqqice_out))     qqqice_out = qqqice
         if (present(TTTice_out))     TTTice_out = TTTice
         if (present(qqqocn_out))     qqqocn_out = qqqocn
         if (present(TTTocn_out))     TTTocn_out = TTTocn
         if (present(Lfresh_out))     Lfresh_out = Lfresh
         if (present(cprho_out))      cprho_out = cprho
         if (present(Cp_out))         Cp_out = Cp
         if (present(puny_out))       puny_out   = puny
         if (present(bignum_out))     bignum_out = bignum
         if (present(pi_out))         pi_out     = pi

         if (present(secday_out)) secday_out = secday
         if (present(c0_out))   c0_out = c0
         if (present(c1_out))   c1_out = c1
         if (present(c1p5_out)) c1p5_out = c1p5
         if (present(c2_out))   c2_out = c2
         if (present(c3_out))   c3_out = c3
         if (present(c4_out))   c4_out = c4
         if (present(c5_out))   c5_out = c5
         if (present(c6_out))   c6_out = c6
         if (present(c8_out))   c8_out = c8
         if (present(c10_out))  c10_out = c10
         if (present(c15_out))  c15_out = c15
         if (present(c16_out))  c16_out = c16
         if (present(c20_out))  c20_out = c20
         if (present(c25_out))  c25_out = c25
         if (present(c100_out)) c100_out = c100
         if (present(c1000_out)) c1000_out = c1000
         if (present(p001_out)) p001_out = p001
         if (present(p01_out))  p01_out = p01
         if (present(p1_out))   p1_out = p1
         if (present(p2_out))   p2_out = p2
         if (present(p4_out))   p4_out = p4
         if (present(p5_out))   p5_out = p5
         if (present(p6_out))   p6_out = p6
         if (present(p05_out))  p05_out = p05
         if (present(p15_out))  p15_out = p15
         if (present(p25_out))  p25_out = p25
         if (present(p75_out))  p75_out = p75
         if (present(p333_out)) p333_out = p333
         if (present(p666_out)) p666_out = p666
         if (present(spval_const_out)) spval_const_out = spval_const
         if (present(pih_out))  pih_out = pih
         if (present(piq_out))  piq_out = piq
         if (present(pi2_out))  pi2_out = pi2

         call icepack_recompute_constants()
         if (icepack_warnings_aborted(subname)) return

      end subroutine icepack_query_constants

!=======================================================================

      subroutine icepack_write_constants(iounit)

      integer (kind=int_kind), intent(in) :: &
         iounit           ! file unit number

      character(len=*),parameter :: subname='(icepack_write_constants)'

         write(iounit,*) subname
         write(iounit,*) "  rhos   = ",rhos
         write(iounit,*) "  rhoi   = ",rhoi
         write(iounit,*) "  rhow   = ",rhow
         write(iounit,*) "  cp_air = ",cp_air
         write(iounit,*) "  emissivity = ",emissivity
         write(iounit,*) "  cp_ice = ",cp_ice
         write(iounit,*) "  cp_ocn = ",cp_ocn
         write(iounit,*) "  depressT = ",depressT
         write(iounit,*) "  dragio = ",dragio
         write(iounit,*) "  albocn = ",albocn
         write(iounit,*) "  gravit = ",gravit
         write(iounit,*) "  viscosity_dyn = ",viscosity_dyn
         write(iounit,*) "  Tocnfrz = ",Tocnfrz
         write(iounit,*) "  rhofresh = ",rhofresh
         write(iounit,*) "  zvir   = ",zvir
         write(iounit,*) "  vonkar = ",vonkar
         write(iounit,*) "  cp_wv  = ",cp_wv
         write(iounit,*) "  stefan_boltzmann = ",stefan_boltzmann
         write(iounit,*) "  Tffresh = ",Tffresh
         write(iounit,*) "  Lsub   = ",Lsub
         write(iounit,*) "  Lvap   = ",Lvap
         write(iounit,*) "  Timelt = ",Timelt
         write(iounit,*) "  Tsmelt = ",Tsmelt
         write(iounit,*) "  ice_ref_salinity = ",ice_ref_salinity
         write(iounit,*) "  iceruf = ",iceruf
         write(iounit,*) "  Cf     = ",Cf
         write(iounit,*) "  Pstar  = ",Pstar
         write(iounit,*) "  Cstar  = ",Cstar
         write(iounit,*) "  kappav = ",kappav
         write(iounit,*) "  kice   = ",kice
         write(iounit,*) "  kseaice = ",kseaice
         write(iounit,*) "  ksno   = ",ksno
         write(iounit,*) "  zref   = ",zref
         write(iounit,*) "  hs_min = ",hs_min
         write(iounit,*) "  snowpatch = ",snowpatch
         write(iounit,*) "  rhosi  = ",rhosi
         write(iounit,*) "  sk_l   = ",sk_l
         write(iounit,*) "  saltmax   = ",saltmax
         write(iounit,*) "  phi_init  = ",phi_init
         write(iounit,*) "  min_salin = ",min_salin
         write(iounit,*) "  salt_loss = ",salt_loss
         write(iounit,*) "  min_bgc   = ",min_bgc
         write(iounit,*) "  dSin0_frazil = ",dSin0_frazil
         write(iounit,*) "  hi_ssl = ",hi_ssl
         write(iounit,*) "  hs_ssl = ",hs_ssl
         write(iounit,*) "  awtvdr = ",awtvdr
         write(iounit,*) "  awtidr = ",awtidr
         write(iounit,*) "  awtvdf = ",awtvdf
         write(iounit,*) "  awtidf = ",awtidf
         write(iounit,*) "  qqqice = ",qqqice
         write(iounit,*) "  TTTice = ",TTTice
         write(iounit,*) "  qqqocn = ",qqqocn
         write(iounit,*) "  TTTocn = ",TTTocn
         write(iounit,*) "  puny   = ",puny
         write(iounit,*) "  bignum = ",bignum
         write(iounit,*) "  pi     = ",pi
         write(iounit,*) "  pih    = ",pih
         write(iounit,*) "  piq    = ",piq
         write(iounit,*) "  pi2    = ",pi2
         write(iounit,*) "  Lfresh = ",Lfresh
         write(iounit,*) "  cprho  = ",cprho
         write(iounit,*) "  Cp     = ",Cp

      end subroutine icepack_write_constants

!=======================================================================

      subroutine icepack_recompute_constants()

      character(len=*),parameter :: subname='(icepack_recompute_constants)'

        cprho  = cp_ocn*rhow
        Lfresh = Lsub-Lvap
        Cp     = 0.5_dbl_kind*gravit*(rhow-rhoi)*rhoi/rhow
        pih    = p5*pi
        piq    = p5*p5*pi
        pi2    = c2*pi

      end subroutine icepack_recompute_constants

!=======================================================================

      end module icepack_constants

!=======================================================================
