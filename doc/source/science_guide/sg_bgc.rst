:tocdepth: 3

.. _ice-bgc:

Biogeochemistry
===============

Aerosols
--------

Basic Aerosols
~~~~~~~~~~~~~~

Aerosols may be deposited on the ice and gradually work their way
through it until the ice melts and they are passed into the ocean. They
are defined as ice and snow volume tracers (Eq. 15 and 16 in CICE.v5
documentation), with the snow and ice each having two tracers for each
aerosol species, one in the surface scattering layer (delta-Eddington
SSL) and one in the snow or ice interior below the SSL.

Rather than updating aerosols for each change to ice/snow thickness due
to evaporation, melting, snow-ice formation, etc., during the
thermodynamics calculation, these changes are deduced from the
diagnostic variables (melts, meltb, snoice, etc) in
**icepack\_aerosol.F90**. Three processes change the volume of ice or snow
but do not change the total amount of aerosol, thus causing the aerosol
concentration (the value of the tracer itself) to increase: evaporation,
snow deposition and basal ice growth. Basal and lateral melting remove
all aerosols in the melted portion. Surface ice and snow melt leave a
significant fraction of the aerosols behind, but they do "scavenge" a
fraction of them given by the parameter ``kscav`` = [0.03, 0.2, 0.02, 0.02,
0.01, 0.01] (only the first 3 are used in CESM, for their 3 aerosol
species). Scavenging also applies to snow-ice formation. When sea ice
ridges, a fraction of the snow on the ridging ice is thrown into the
ocean, and any aerosols in that fraction are also lost to the ocean.

As upper SSL or interior layers disappear from the snow or ice, aerosols
are transferred to the next lower layer, or into the ocean when no ice
remains. The atmospheric flux ``faero_atm`` contains the rates of aerosol
deposition for each species, while ``faero_ocn`` has the rate at which the
aerosols are transferred to the ocean.

The aerosol tracer flag ``tr_aero`` must be set to true in **icepack\_in**, and
the number of aerosol species is set in **icepack.settings**; CESM uses 3.

Z-Aerosols
~~~~~~~~~~

An alternate scheme for aerosols in sea ice is available using
the brine motion based transport scheme of the biogeochemical tracers.
All vertically resolved biogeochemical tracers (z-tracers), including
aerosols, have the potential to be atmospherically deposited onto the
snow or ice, scavenged during snow melt, and passed into the brine. The
mobile fraction (discussed in :ref:`mobile-and-stationary`) is
then transported via brine drainage processes
(Eq. :eq:`mobile-transport`) while a
stationary fraction (discussed in :ref:`mobile-and-stationary`)
adheres to the ice crystals. Snow deposition and the process of
scavenging aerosols during snow melt is consistent with the basic
aerosol scheme, though parameters have been generalized to accomodate
potential atmospheric deposition for all z-tracers. For an example, see
the scavenging parameter ``kscavz`` for z-tracers defined in
**icepack\_zbgc\_shared.F90**.

Within the snow, z-tracers are defined as concentrations in the snow
surface layer (:math:`h_{ssl}`) and the snow interior
(:math:`h_s-h_{ssl}`). The total snow content of z-tracers per ice area
per grid cell area, :math:`C_{snow}` is

.. math::
   C_{snow} = C_{ssl}h_{ssl} + C_{sint}(h_{s}-h_{ssl})

One major difference in how the two schemes model snow aerosol transport
is that the fraction scavenged from snow melt in the z-tracer scheme is
not immediately fluxed into the ocean, but rather, enters the ice as a
source of low salinity but potentially tracer-rich brine. The snow melt
source is included as a surface flux condition in **icepack\_algae.F90**.

All the z-aerosols are nonreactive with the exception of the dust
aerosols. We assume that a small fraction of the dust flux into the ice
has soluble iron (``dustFe_sol`` in **icepack\_in**) and so is
passed to the dissolved iron tracer. The remaining dust passes through
the ice without reactions.

To use z-aerosols, ``tr_zaero`` must be set to true in **icepack\_in**, and the
number of z-aerosol species is set in **icepack.settings**, ``TRZAERO``. Note, the
basic tracers ``tr_aero`` must be false and ``NTRAERO`` in **icepack.settings**
should be 0. In addition, z-tracers and the brine height tracer must
also be active. These are set in **icepack\_in** with ``tr_brine`` and
``z_tracer`` set to true. In addition, to turn on the radiative coupling
between the aerosols and the Delta-Eddington radiative scheme, ``shortwave``
must equal ’dEdd’ and ``dEdd_algae`` must be true in **icepack\_in**.

.. _brine-ht:

Brine height
------------

The brine height, :math:`h_b`, is the distance from the ice-ocean
interface to the brine surface. When ``tr_brine`` is set true in
**icepack\_in** and ``TRBRI`` is set equal to 1 in **icepack.settings**, the brine
surface can move relative to the ice surface. Physically, this occurs
when the ice is permeable and there is a nonzero pressure head: the
difference between the brine height and the equilibrium sea surface.
Brine height motion is computed in **icepack\_brine.F90** from thermodynamic
variables and the ice microstructural state
deduced from internal bulk salinities and temperature. This tracer is
required for the transport of vertically resolved biogeochemical tracers
and is closely coupled to the z-salinity prognostic salinity model.

Vertical transport processes are, generally, a result of the brine
motion. Therefore the vertical transport equations for biogeochemical
tracers will be defined only where brine is present. This region, from
the ice-ocean interface to the brine height, defines the domain of the
vertical bio-grid. The resolution of the bio-grid is specified in
**icepack.settings** by setting the variable ``NBGCLYR``. A detailed description of
the bio-grid is given in section :ref:`grids`. The ice
microstructural state, determined in **icepack\_brine.F90**, is computed
from sea ice salinities and temperatures linearly interpolated to the
bio-grid. When :math:`h_b > h_i`, the upper surface brine is assumed to
have the same temperature as the ice surface.

Brine height is transported horizontally as the fraction
:math:`f_{bri} = h_b/h_i`, a volume conserved tracer. Note that unlike the sea ice porosity, brine height
fraction may be greater than 1 when :math:`h_b > h_i`.

Changes to :math:`h_b` occur from ice and snow melt, ice bottom boundary
changes, and from pressure adjustments. The computation of :math:`h_b`
at :math:`t+\Delta
t` is a two step process. First, :math:`h_b` is updated from changes in
ice and snow thickness, ie.

.. math::
   h_b' = h_b(t) + \Delta h_b|_{h_i,h_s} .
   :label: hb_thickness_changes

Second, pressure driven adjustments arising from meltwater flushing and
snow loading are applied to :math:`h'_b`. Brine flow due to pressure
forces are governed by Darcy’s equation 

.. math::
   w = -\frac{\Pi^* \bar{\rho} g}{\mu}\frac{h_p}{h_i}.
   :label: Darcy1

The vertical component of the net permeability tensor :math:`\Pi^*` is
computed as

.. math::
   \Pi^* = \left(\frac{1}{h}\sum_{i=1}^N{\frac{\Delta
         z_i}{\Pi_i}}\right)^{-1}
   :label: netPi1

where the sea ice is composed of :math:`N` vertical layers with
:math:`i`\ th layer thickness :math:`\Delta z_i` and permeability
:math:`\Pi_i`. The average sea ice density is :math:`\bar{\rho}`
specified in **icepack\_zbgc\_shared.F90**. The hydraulic head is
:math:`h_p = h_b - h_{sl}` where :math:`h_{sl}` is the sea level given
by

.. math::
   h_{sl} = \frac{\bar{\rho}}{\rho_w}h_i + \frac{\rho_s}{\rho_w}h_s .
   :label: hsl

Assuming constant :math:`h_i` and :math:`h_s` during Darcy flow, the
rate of change of :math:`h_b` is

.. math::
   \frac{\partial h_b}{\partial t} = -w h_p
   :label: h_p

where :math:`w_o = \Pi^* \bar{\rho}
g/(h_i\mu\phi_{top})` and :math:`\phi_{top}` is the upper surface
porosity. When the Darcy flow is downward into the ice
(:math:`w_o < 0`), then :math:`\phi_{top}` equals the sea ice porosity
in the uppermost layer. However, when the flow is upwards into the snow,
then :math:`\phi_{top}` equals the snow porosity phi\_snow specified in
**icepack\_in**. If a negative number is specified for phi\_snow, then the
default value is used: phi\_snow :math:`=1 - \rho_s/\rho_w`.

Since :math:`h_{sl}` remains relatively unchanged during Darcy flow,
:eq:`h_p` has the approximate solution

.. math::
   \begin{aligned}
   h_b(t+\Delta t) \approx h_{sl}(t+\Delta t) +  [h'_b - h_{sl}(t+\Delta t)]\exp\left\{-w \Delta t\right\}.\end{aligned}
   :label: brine_height

The contribution :math:`\Delta h_b|_{h_i,h_s}` arises from snow and ice
melt and bottom ice changes. Since the ice and brine bottom boundaries
coincide, changes in the ice bottom from growth or melt,
:math:`(\Delta h_i)_{bot}`, equal the bottom brine boundary changes. The
surface contribution from ice and snow melt, however, is opposite in
sign. The ice contribution is as follows. If :math:`h_i > h_b` and the
ice surface is melting, ie. :math:`(\Delta h_i)_{top} <
0`), then meltwater increases the brine height:

.. math::
   \begin{aligned}
   \left(\Delta h_b\right)_{top} = \frac{\rho_i}{\rho_o} \cdot \left\{ \begin{array}{ll}
    -(\Delta h_i)_{top} &  \mbox{if }
    |(\Delta h_i)_{top}| < h_i-h_b  \\
   h_i-h_b & \mbox{otherwise.}   \end{array} \right.  \end{aligned}
   :label: delta-hb

For snow melt (:math:`\Delta h_s < 0`), it is assumed that all snow
meltwater contributes a source of surface brine. The total change from
snow melt and ice thickness changes is

.. math::
   \Delta h_b|_{h_i,h_s} = \left( \Delta
   h_b\right)_{top} -\left(\Delta h_i\right)_{bot} -\frac{\rho_s}{\rho_o}\Delta h_s.
   :label: dzdt_meltwater

The above brine height calculation is used only when :math:`h_i` and
:math:`h_b` exceed a minimum thickness, thinS, specified in
**icepack\_zbgc\_shared.F90**. Otherwise

.. math::
   h_b(t+\Delta t) = h_b(t) + \Delta h_i
   :label: thinbrine1

provided that :math:`|h_{sl}-h_b| \leq 0.001`. This formulation ensures
small Darcy velocities when :math:`h_b` first exceeds thinS.


Both the volume fraction :math:`f_{bri}` and the area-weighted brine
height :math:`h_b` are available for output.

.. math:: 
   {{\sum f_{bri} v_i} \over {\sum v_i}},
   :label: volume-frac

while ``hbri`` is comparable to hi (:math:`h_i`)

.. math:: 
   {{\sum f_{bri} h_i a_i} \over {\sum a_i}},
   :label: volume-frac2

where the sums are taken over thickness categories.

Sea ice ecosystem
-----------------

There are two options for modeling biogeochemistry in sea ice: 1) a
skeletal layer or bottom layer model that assumes biology
and biological molecules are restricted to a single layer at the base of
the sea ice; and 2) a vertically resolved model (zbgc) that allows for
biogeochemical processes throughout the ice column. The two models may
be run with the same suite of biogeochemical tracers and use the same
module **algal\_dyn** in **icepack\_algae.F90** to determine the biochemical
reaction terms for the tracers at each vertical grid level. In the case
of the skeletal-layer model this is a single layer, while for zbgc there are
``NBGCLYR``\ :math:`+1` vertical layers. The primary difference between the
two schemes is in the vertical transport assumptions for each
biogeochemical tracer. This includes the parameterizations of fluxes
between ocean and ice.

In order to run with the skeletal-layer model, the code must be built with the
following options in **icepack.settings**:

::

    setenv TRBGCS 1   # set to 1 for skeletal layer tracers
    setenv TRBGCZ 0   # set to 1 for zbgc tracers

For zbgc with 8 vertical layers:

::

    setenv TRBRI  1   # set to 1 for brine height tracer
    setenv TRBGCS 0   # set to 1 for skeletal layer tracers
    setenv TRBGCZ 1   # set to 1 for zbgc tracers
    setenv NBGCLYR 7  # number of zbgc layers 

There are also environmental variables in **icepack.settings** that, in part,
specify the complexity of the ecosystem and are used for both zbgc and
the skeletal-layer model. These are 1) ``TRALG``, the number of algal species; 2)
``TRDOC``, the number of dissolved organic carbon groups, 3) ``TRDIC``, the
number of dissolved inorganic carbon groups (this is currently not yet
implemented and should be set to 0); 4) ``TRDON``, the number of dissolved
organic nitrogen groups, 5) ``TRFEP``, the number of particulate iron
groups; and 6) ``TRFED``, the number of dissolved iron groups. The current
version of **algal\_dyn** biochemistry has parameters for up to 3 algal
species (diatoms, small phytoplankton and *Phaeocystis* sp,
respectively), 2 DOC tracers (polysaccharids and lipids, respectively),
0 DIC tracers, 1 DON tracer (proteins/amino acids), 1 particulate iron
tracer and 1 dissolved iron tracer. Note, for tracers with multiple
species/groups, the order is important. For example, specifying
``TRALG`` = 1 will compute reaction terms using parameters
specific to ice diatoms.  However, many of these parameters can be modified in **icepack\_in**. 

The complexity of the algal ecosystem must be specified in both
**icepack.settings** during the build and in the namelist, **icepack\_in**. The
procedure is equivalent for both the skeletal-layer model and zbgc. The namelist
specification is described in detail in section :ref:`zbgc`

Biogeochemical upper ocean concentrations are initialized in the
subroutine **icepack\_init\_ocean\_conc** in **icepack\_zbgc.F90** unless
coupled to the ocean biogeochemistry. Silicate and nitrate may be read
from a file. This option is specified in the namelist by setting the
variables ``bgc_data_type`` to ``ISPOL``  or ``NICE``. The location of 
forcing files is specified in ``data_dir`` and the filename is also in
namelist, ``bgc_data_file``.


Skeletal Layer BGC
~~~~~~~~~~~~~~~~~~

In the skeletal layer model, biogeochemical processing is modelled as a
single layer of reactive tracers attached to the sea ice bottom.
Optional settings are available via the *zbgc\_nml* namelist in
**icepack\_in**. In particular, ``skl_bgc`` must be true and ``z_tracers`` and
``solve_zbgc`` must both be false.

Skeletal tracers :math:`T_b` are ice area conserved and follow the
horizontal transport Equation :eq:`itd-transport`. For each
horizontal grid point, local biogeochemical tracer equations are solved
in **icepack\_algae.F90**. There are two types of ice-ocean tracer flux
formulations: 1) ‘Jin2006’ modeled after the growth rate dependent
piston velocity and 2) ‘constant’ modeled after a constant piston
velocity. The formulation is specified in **icepack\_in** by setting
``bgc_flux_type`` equal to ‘Jin2006’ or ‘constant’.

In addition to horizontal advection and transport among thickness
categories, biogeochemical tracers (:math:`T_b` where
:math:`b = 1,\ldots, N_b`) satisfy a set of local coupled equations of
the form

.. math::
   \frac{d T_b}{dt} = w_b \frac{\Delta T_b}{\Delta z} +  R_b({T_j : j = 1,\ldots,N_b})
   :label: bgc_Tracer

where :math:`R_b` represents the nonlinear biochemical reaction terms
(described in section :ref:`reactions`) and :math:`\Delta z` is a length
scale representing the molecular sublayer of the ice-ocean interface.
Its value is absorbed in the piston velocity parameters. The piston
velocity :math:`w_b` depends on the particular tracer and the flux
formulation.

For ‘Jin2006’, the piston velocity is a function of ice growth and melt
rates. All tracers (algae included) flux with the same piston velocity
during ice growth, :math:`dh/dt > 0`:

.. math::
   \begin{aligned}
   w_b  & =  & - p_g\left|m_1 + m_2 \frac{dh}{dt} - m_3
     \left(\frac{dh}{dt} \right)^2\right|\end{aligned}
   :label: pwJin_growth

with parameters :math:`m_1`, :math:`m_2`, :math:`m_3` and :math:`p_g`
defined in **skl\_biogeochemistry** in **icepack\_algae.F90**. For ice melt,
:math:`dh/dt < 0`, all tracers with the exception of ice algae flux with

.. math::
   \begin{aligned}
   w_b  & =  & p_m\left|m_2 \frac{dh}{dt} - m_3
       \left(\frac{dh}{dt}  \right)^2\right| \end{aligned}
   :label: pwJin_melt

with :math:`p_m` defined in **skl\_biogeochemistry**. The ‘Jin2006’
formulation also requires that for both expressions,
:math:`|w_b| \leq 0.9
h_{sk}/\Delta t`. The concentration difference at the ice-ocean boundary
for each tracer, :math:`\Delta
T_b`, depends on the sign of :math:`w_b`. For growing ice,
:math:`w_b <0`, :math:`\Delta T_b  = T_b/h_{sk} - T_{io}`, where
:math:`T_{io}` is the ocean concentration of tracer :math:`i`. For
melting ice, :math:`w_b > 0`, :math:`\Delta T_b = T_b/h_{sk}`.

In ‘Jin2006’, the algal tracer (:math:`N_a`) responds to ice melt in the
same manner as the other tracers :eq:`pwJin_melt`. However, this is
not the case for ice growth. Unlike dissolved nutrients, algae are able
to cling to the ice matrix and resist expulsion during desalination. For
this reason, algal tracers do not flux between ice and ocean during ice
growth unless the ice algal brine concentration is less than the ocean
algal concentration (:math:`N_o`). Then the ocean seeds the sea ice
concentration according to

.. math::
   \begin{aligned}
   w_b \frac{\Delta N_a}{\Delta z} = \frac{N_oh_{sk}/\phi_{sk} -
     N_a}{\Delta t}\end{aligned}
   :label: seed2

The ‘constant’ formulation uses a fixed piston velocity (PVc) for
positive ice growth rates for all tracers except :math:`N_a`. As in
‘Jin2006’, congelation ice growth seeds the sea ice algal population
according to :eq:`seed2` when :math:`N_a < N_o
h_{sk}/\phi_{sk}`. For bottom ice melt, all tracers follow the
prescription

.. math::
   \begin{aligned}
    w_b \frac{\Delta T_b}{\Delta z} & = &  \left\{ \begin{array}{ll}
      T_b   |dh_i/dt|/h_{sk} \ \ \ \ \ &   \mbox{if }
    |dh_i/dt|\Delta t/h_{sk} < 1  \\
   T_b/\Delta t & \mbox{otherwise.}   \end{array} \right. \end{aligned} 
   :label: constant_melt

A detailed description of the biogeochemistry reaction terms is given in
section :ref:`reactions`.


.. _zbgc:

Vertical BGC (''zbgc'')
~~~~~~~~~~~~~~~~~~~~~~~

In order to solve for the vertically resolved biogeochemistry, several
flags in **icepack\_in** must be true: a) ``tr_brine``, b) ``z_tracers``, and c)
``solve_zbgc``.

-  ``tr_brine`` = true turns on the dynamic brine height tracer,
   :math:`h_b`, which defines the vertical domain of the biogeochemical
   tracers. z-Tracer horizontal transport is conserved on ice
   volume\ :math:`\times`\ brine height fraction.

-  ``z_tracers`` = true indicates use of vertically resolved
   biogeochemical and z-aerosol tracers. This flag alone turns on the
   vertical transport scheme but not the biochemistry.

-  ``solve_zbgc`` = true turns on the biochemistry for the vertically
   resolved tracers and automatically turns on the algal nitrogen tracer
   flag tr\_bgc\_N. If false, ``tr_bgc_N`` is set false and any other
   biogeochemical tracers in use are transported as passive tracers.
   This is appropriate for the black carbon and dust aerosols specified
   by ``tr_zaero`` true.

In addition, a halodynamics scheme must also be used. The default
thermo-halodynamics is mushy layer ``ktherm`` set to 2. An alternative uses
the Bitz and Lipscomb thermodynamics ``ktherm`` set to 1 and ``solve_zsal``
true (referred to as "zsalinity").

With the above flags, the default biochemistry is a simple
algal-nitrate system: ``tr_bgc_N`` and ``tr_bgc_Nit`` are true. Options
exist in **icepack\_in** to use a more complicated ecosystem which includes up
to three algal classes, two DOC groups, one DON pool, limitation by
nitrate, silicate and dissolved iron, sulfur chemistry plus refractory
humic material.

The **icepack\_in** namelist options are described in the :ref:`tabnamelist`.


Vertically resolved z-tracers are brine- volume conserved and thus depend
on both the ice volume and the brine height fraction tracer
(:math:`v_{in}f_b`). These tracers follow the conservation equations for
multiply dependent tracers (see, for example Equation :eq:`transport-apnd-lvl` where :math:`a_{pnd}` is a tracer on :math:`a_{lvl}a_{i}`)  

The following sections describe the vertical
transport equation for mobile tracers, the partitioning of tracers into
mobile and stationary fractions and the biochemical reaction equations.
The vertical bio-grid is described in the :ref:`grids` section.

.. _mobile-and-stationary:

*Mobile and stationary phases*
``````````````````````````````
Purely mobile tracers are tracers which move with the brine and thus, in
the absence of biochemical reactions, evolve like salinity. For vertical
tracer transport of purely mobile tracers, the flux conserved quantity
is the bulk tracer concentration multiplied by the ice thickness, i.e.
:math:`C = h\phi
[c]`, where :math:`h` is the ice thickness, :math:`\phi` is the
porosity, and :math:`[c]` is the tracer concentration in the brine.
:math:`\phi`, :math:`[c]` and :math:`C` are defined on the interface bio
grid (igrid):

.. math::
   \mbox{igrid}(k) = \Delta (k-1) \ \ \ \mbox{for }k = 1:n_b+1 \ \ \mbox{and }\Delta = 1/n_b.

The biogeochemical module solves the following equation:

.. math::
   \begin{aligned}
   \frac{\partial C}{\partial t} & =& \frac{\partial }{\partial x}\left\{
   \left( \frac{v}{h} + \frac{w_f}{h\phi} -
     \frac{\tilde{D}}{h^2\phi^2}\frac{\partial \phi}{\partial x} \right) C
   + \frac{\tilde{D}}{h^2\phi}\frac{\partial C}{\partial x} 
   \right\} + h\phi R([c])\end{aligned}
   :label: mobile-transport

where :math:`D_{in} = \tilde{D}/h^2 =  (D + \phi D_m)/h^2` and
:math:`R([c])` is the nonlinear biogeochemical interaction term (see
:cite:`Jeffery11`).

The solution to :eq:`mobile-transport` is flux-corrected and
positive definite. This is accomplished using a finite element Galerkin
discretization. Details are in :ref:`tracer-numerics`.

In addition to purely mobile tracers, some tracers may also adsorb or
otherwise adhere to the ice crystals. These tracers exist in both the
mobile and stationary phases. In this case, their total brine
concentration is a sum :math:`c_m + c_s` where :math:`c_m` is the mobile
fraction transported by equation :eq:`mobile-transport` and :math:`c_s`
is fixed vertically in the ice matrix. The algae are an exception,
however. We assume that algae in the stationary phase resist brine
motion, but rather than being fixed vertically, these tracers maintain
their relative position in the ice. Algae that adhere to the ice
interior (bottom, surface), remain in the ice interior (bottom, surface)
until release to the mobile phase.

In order to model the transfer between these fractions, we assume that
tracers adhere (are retained) to the crystals with a time-constant of
:math:`\tau_{ret}`, and release with a time constant :math:`\tau_{rel}`,
i.e.

.. math::
   \begin{aligned}
   \frac{\partial c_m}{\partial t} & = & -\frac{c_m}{\tau_{ret}} + \frac{c_s}{\tau_{rel}} \\
   \frac{\partial c_s}{\partial t} & = &\frac{c_m}{\tau_{ret}} - \frac{c_s}{\tau_{rel}}\end{aligned}

We use the exponential form of these equations:

.. math::
   \begin{aligned}
   c_m^{t+dt} & = & c_m^t\exp\left(-\frac{dt}{\tau_{ret}}\right) +
   c^t_s\left(1-\exp\left[-\frac{dt}{\tau_{rel}}\right]\right) \end{aligned}

.. math::
   \begin{aligned}
   c_s^{t+dt} & = & c_s^t\exp\left(-\frac{dt}{\tau_{rel}}\right) +
   c_m^t\left(1-\exp\left[-\frac{dt}{\tau_{ret}}\right]\right) \end{aligned}

The time constants are functions of the ice growth and melt rates
(:math:`dh/dt`). All tracers except algal nitrogen diatoms follow the
simple case: when :math:`dh/dt \geq 0`, then
:math:`\tau_{rel} \rightarrow \infty` and :math:`\tau_{ret}` is finite.
For :math:`dh/dt < 0`, then :math:`\tau_{ret} \rightarrow \infty` and
:math:`\tau_{rel}` is finite. In other words, ice growth promotes
transitions to the stationary phase and ice melt enables transitions to
the mobile phase.

The exception is the diatom pool. We assume that diatoms, the first
algal nitrogen group, can actively maintain their relative position
within the ice, i.e. bottom (interior, upper) algae remain in the bottom
(interior, upper) ice, unless melt rates exceed a threshold. The
namelist parameter ``algal_vel`` sets this threshold.

The variable ``bgc_tracer_type`` determines the mobile to stationary
transition timescales for each z-tracer. It is multi-dimensional with a
value for each z-tracer. For ``bgc_tracer_type``(k) equal to -1, the kth
tracer remains solely in the mobile phase. For ``bgc_tracer_type``
equal to 1, the tracer has maximal rates in the retention phase and
minimal in the release. For ``bgc_tracer_type`` equal to 0, the tracer
has maximal rates in the release phase and minimal in the retention.
Finally for ``bgc_tracer_type`` equal to 0.5, minimum timescales are
used for both transitions. Table :ref:`tab-phases` summarizes the
transition types. The tracer types are: ``algaltype_diatoms``,
``algaltype_sp`` (small plankton), ``algaltype_phaeo`` (*phaeocystis*),
``nitratetype``, ``ammoniumtype``, ``silicatetype``, ``dmspptype``, 
``dmspdtype``, ``humtype``,
``doctype_s`` (saccharids), ``doctype_l`` (lipids), ``dontype_protein``,
``fedtype_1``, ``feptype_1``, ``zaerotype_bc1`` (black carbon class 1),
``zaerotype_bc2`` (black carbon class 2), and four dust classes,
``zaerotype_dustj``, where j takes values 1 to 4. These may be modified to
increase or decrease retention. Another option is to alter the minimum
``tau_min`` and maximum ``tau_max`` timescales which would impact all the
z-tracers.

.. _tab-phases:

.. table:: *Types of Mobile and Stationary Transitions*

   +-----------------+--------------------+--------------------+------------------------------+
   | bgc_tracer_type | :math:`\tau_{ret}` | :math:`\tau_{rel}` |        Description           |
   +=================+====================+====================+==============================+
   |     -1.0        | :math:`\infty`     |         0          | entirely in the mobile phase |
   +-----------------+--------------------+--------------------+------------------------------+
   |      0.0        |       min          |        max         |     retention dominated      |
   +-----------------+--------------------+--------------------+------------------------------+
   |      1.0        |       max          |        min         |      release dominated       |
   +-----------------+--------------------+--------------------+------------------------------+
   |      0.5        |       min          |        min         |  equal but rapid exchange    |
   +-----------------+--------------------+--------------------+------------------------------+
   |      2.0        |       max          |        max         |  equal but slow exchange     |
   +-----------------+--------------------+--------------------+------------------------------+

The fraction of a given tracer in the mobile phase is independent of ice
depth and stored in the tracer variable zbgc\_frac. The horizontal
transport of this tracer is conserved on brine volume and so is
dependent on two tracers: brine height fraction (:math:`f_b`) and ice
volume (:math:`v_{in}`). The conservation equations are given by

.. math::
   {\partial\over\partial t} (f_{b}v_{in}) + \nabla \cdot (f_{b}v_{in} {\bf u}) = 0.

The tracer, ``zbgc_frac``, is initialized to 1 during new ice formation,
because all z-tracers are initially in the purely mobile phase.
Similarly, as the ice melts, z-tracers return to the mobile phase. Very
large release timescales will prevent this transition and could result
in an unphysically large accumulation during the melt season.

.. _tracer-numerics:

*Flux-corrected, positive definite transport scheme*
````````````````````````````````````````````````````

Numerical solution of the vertical tracer transport equation is
accomplished using the finite element Galerkin discretization. Multiply
:eq:`mobile-transport` by "w" and integrate by parts

.. math::
   \begin{aligned}
   & & \int_{h}\left[ w\frac{\partial C}{\partial t} -   \frac{\partial
       w}{\partial x} \left(-\left[\frac{v}{h} + \frac{w_f}{h\phi}\right]C + \frac{D_{in}}{\phi^2}\frac{\partial
         \phi}{\partial x}C  -  \frac{D_{in}}{\phi}\frac{\partial C}{\partial
         x} \right) \right]dx \nonumber \\
   & + &  w\left.\left(
       -\left[\frac{1}{h}\frac{dh_b}{dt}+  \frac{w_f}{h\phi}\right]C + \frac{D_{in}}{\phi^2}\frac{\partial \phi}{\partial
       x}C -\frac{D_{in}}{\phi}\frac{\partial C}{\partial
       x}\right)\right|_{bottom} + w\left[\frac{1}{h}\frac{dh_t}{dt} +
   \frac{w_f}{h\phi}\right]C|_{top}  =  0\end{aligned}

The bottom boundary condition indicated by :math:`|_{bottom}` satisfies

.. math::
   \begin{aligned}
   -w\left.\left(
       -\left[\frac{1}{h}\frac{dh_b}{dt}+  \frac{w_f}{h\phi}\right]C + \frac{D_{in}}{\phi^2}\frac{\partial \phi}{\partial
       x}C -\frac{D_{in}}{\phi}\frac{\partial C}{\partial
       x}\right)\right|_{bottom} & = & \nonumber \\
    w\left[\frac{1}{h}\frac{dh_b}{dt} +
   \frac{w_f}{h \phi_{N+1}}\right](C_{N+2} \mbox{ or }C_{N+1}) -
   w\frac{D_{in}}{\phi_{N+1}(\Delta h+g_o)}\left(C_{N+1} - C_{N+2}\right) & & \end{aligned}

where :math:`C_{N+2} = h\phi_{N+1}[c]_{ocean}` and :math:`w = 1` at the
bottom boundary and the top. The component :math:`C_{N+2}` or
:math:`C_{N+1}` depends on the sign of the advection boundary term. If
:math:`dh_b + w_f/\phi > 0` then use :math:`C_{N+2}` otherwise
:math:`C_{N+1}`.

Define basis functions as linear piecewise, with two nodes (boundary
nodes) in each element. Then for :math:`i > 1` and :math:`i < N+1`

.. math::
   \begin{aligned}
   w_i(x) & = & \left\{ \begin{array}{llll}
   0 & x < x_{i-1} \\
   (x - x_{i-1})/\Delta & x_{i-1}< x \leq x_{i} \\
   1 - (x-x_i)/\Delta & x_i \leq x < x_{i+1} \\
   0, & x \geq x_{i+1} 
   \end{array}
   \right.\end{aligned}

For :math:`i=1`

.. math::
   \begin{aligned}
   w_1(x) & = & \left\{ \begin{array}{ll}
   1 - x/\Delta & x < x_{2} \\
   0, & x \geq x_{2} 
   \end{array}
   \right.\end{aligned}

and :math:`i = N+1`

.. math::
   \begin{aligned}
   w_{N+1}(x) & = & \left\{ \begin{array}{ll}
   0, & x < x_{N} \\
   (x-x_{N})/\Delta & x \geq x_{N}
   \end{array}
   \right.\end{aligned}

Now assume a form

.. math::
   C_h =  \sum_j^{N+1} c_j w_j

Then

.. math::
   \begin{aligned}
   \int_h C_h dx & = & c_1\int_0^{x_{2}}\left(1-\frac{x}{\Delta}\right)dx
     +  c_{N+1}\int_{x_N}^{x_{N+1}}\frac{x-x_{N}}{\Delta}dx  \nonumber \\
   & + &
     \sum_{j=2}^{N}c_j\left\{\int_{j-1}^{j}\frac{x-x_{j-1}}{\Delta}dx +
       \int_{j}^{j+1}\left[1 - \frac{(x-x_j)}{\Delta}\right]dx\right\}
     \nonumber \\
   & = & \Delta \left[\frac{c_1}{2} + \frac{c_{N+1}}{2} + \sum_{j = 2}^{N}c_{j}\right]\end{aligned}

Now this approximate solution form is substituted into the variational
equation with :math:`w = w_h \in \{w_j\}`

.. math::
   \begin{aligned}
   0 &= & \int_{h}\left[ w_h\frac{\partial C_h}{\partial t} -   \frac{\partial
       w_h}{\partial x} \left(\left[-\frac{v}{h} - \frac{w_f}{h\phi}+ \frac{D_{in}}{\phi^2}\frac{\partial
         \phi}{\partial x}\right]C_h  -  \frac{D_{in}}{\phi}\frac{\partial C_h}{\partial
         x} \right) \right]dx \nonumber \\
   & + &  w_h\left.\left(
       -\left[\frac{1}{h}\frac{dh_b}{dt}+  \frac{w_f}{h\phi}\right]C_h + \frac{D_{in}}{\phi^2}\frac{\partial \phi}{\partial
       x}C -\frac{D_{in}}{\phi}\frac{\partial C_h}{\partial
       x}\right)\right|_{bottom} + w_h\left[\frac{1}{h}\frac{dh_t}{dt} +
   \frac{w_f}{h\phi}\right]C_h|_{top} \end{aligned}

The result is a linear matrix equation

.. math::
   M_{jk}\frac{\partial C_k(t)}{\partial t} = [K_{jk}+S_{jk}] C_k(t) + q_{in}

where

.. math::
   \begin{aligned}
   M_{jk} & = & \int_h w_j(x)w_k(x)dx \nonumber \\
   K_{jk} & = & \left[-\frac{v}{h} - \frac{w_f}{h\phi}+ \frac{D_{in}}{\phi^2}\frac{\partial
         \phi}{\partial x}\right]\int_h \frac{\partial w_j}{\partial x}
     w_k dx \nonumber \\
   &-&
   \left. w_j\left(-\left[\frac{v}{h} + \frac{w_f}{h\phi_k}\right]w_k +
       \frac{D_{in}}{\phi^2}\frac{\partial \phi_k}{\partial x} w_k - \frac{D_{in}}{\phi_k}\frac{\partial w_k}{\partial
         x}\right)\right|_{bot} \nonumber \\
   & = & -V_k\int_h \frac{\partial w_j}{\partial x} w_k dx -
   \left. w_j\left(-V_kw_k - \frac{D_{in}}{\phi_k}\frac{\partial w_k}{\partial
         x}\right)\right|_{bot} \nonumber \\
   S_{jk} & = & -\frac{D_{in}}{\phi_k}\int_h \frac{\partial w_j}{\partial x} \cdot
   \frac{\partial w_k}{\partial x}dx \nonumber \\
   q_{in} & = & -V C_{t} w_j(x)|_{t}\end{aligned}

and :math:`C_{N+2} = h\phi_{N+1}[c]_{ocean}`.

For the top condition :math:`q_{in}` is applied to the upper value
:math:`C_2` when :math:`VC_t < 0`, i.e. :math:`q_{in}` is a source.

Compute the :math:`M_{jk}` integrals:

.. math::
   \begin{aligned}
   M_{jj} & = & \int_{x_{j-1}}^{x_j}\frac{(x- x_{j-1})^2}{\Delta^2}dx +
   \int_{x_{j}}^{x_{j+1}}\left[ 1-\frac{(x- x_{j})}{\Delta}\right]^2dx =
   \frac{2\Delta}{3} \ \ \ \mbox{for }1 < j < N+1 \nonumber \\
   M_{11} & = & \int_{x_{1}}^{x_{2}}\left[ 1-\frac{(x- x_{2})}{\Delta}\right]^2dx =
   \frac{\Delta}{3}   \nonumber \\
   M_{N+1,N+1} & = &\int_{x_{N}}^{x_{N+1}}\frac{(x- x_{N})^2}{\Delta^2}dx
   = \frac{\Delta}{3}\nonumber \end{aligned}

Off diagonal components:

.. math::
   \begin{aligned}
   M_{j,j+1} & = &  \int_{x_{j}}^{x_{j+1}}\left[1 - \frac{(x-
       x_{j})}{\Delta}\right]\left[ \frac{x-x_{j}}{\Delta}\right]dx =
   \frac{\Delta}{6} \ \ \ \mbox{for }j < N+1 \nonumber \\
   M_{j,j-1} & = &  \int_{x_{j-1}}^{x_{j}}\left[ \frac{x-x_{j-1}}{\Delta}\right]\left[1 - \frac{(x-
       x_{j-1})}{\Delta}\right]dx =
   \frac{\Delta}{6} \ \ \ \mbox{for }j > 1 \nonumber \\\end{aligned}

Compute the :math:`K_{jk}` integrals:

.. math::
   \begin{aligned}
   K_{jj} &=& k'_{jj}\left[ \int_{x_{j-1}}^{x_j} \frac{\partial w_j}{\partial x}
   w_j dx + \int_{x_{j}}^{x_{j+1}} \frac{\partial w_j}{\partial x}
   w_j dx \right] \nonumber \\
   &  = &  \frac{1}{2} + -\frac{1}{2} = 0  \ \ \ \mbox{for } 1 < j < N+1 \nonumber \\ 
   K_{11} & = &  -\frac{k'_{11}}{2} = \frac{1}{2}\left[\frac{v}{h} +
     \frac{w_f}{h\phi}\right] \nonumber \\
   K_{N+1,N+1}  & = & \frac{k'_{N+1,N+1}}{2} +\min\left[0, \left(\frac{1}{h}\frac{dh_b}{dt} +
   \frac{w_f}{h \phi_{N+1}}\right)\right] -
   \frac{D_{in}}{\phi_{N+1}(g_o/h)} \nonumber \\
   & = & \left[-\frac{v}{h} - \frac{w_f}{h\phi}+ \frac{D_{in}}{\phi^2}\frac{\partial
         \phi}{\partial x}\right]\frac{1}{2} +\min\left[0, \left(\frac{1}{h}\frac{dh_b}{dt} +
   \frac{w_f}{h \phi_{N+1}}\right)\right] -
   \frac{D_{in}}{\phi_{N+1}(g_o/h)} \end{aligned}

Off diagonal components:

.. math::
   \begin{aligned}
   K_{j(j+1)} &=& k'_{j(j+1)}\int_{x_{j}}^{x_{j+1}} \frac{\partial w_j}{\partial x}
   w_{j+1} dx  =
   -k'_{j(j+1)}\int_{x_{j}}^{x_{j+1}}\frac{(x-x_j)}{\Delta^2}dx \nonumber
   \\
   & = & -\frac{k'_{j(j+1)}}{\Delta^2}\frac{\Delta^2}{2} =
   -\frac{k'_{j(j+1)}}{2}  = p5*\left[\frac{v}{h} + \frac{w_f}{h\phi}- \frac{D_{in}}{\phi^2}\frac{\partial
         \phi}{\partial x}\right]_{(j+1)} \ \ \ \mbox{for }  j < N+1 \nonumber \\
   K_{j(j-1)} &=& k'_{j(j-1)}\int_{x_{j-1}}^{x_{j}} \frac{\partial w_j}{\partial x}
   w_{j-1} dx  =
   k'_{j(j-1)}\int_{x_{j-1}}^{x_{j}}\left[1 -
     \frac{(x-x_{j-1})}{\Delta^2}\right] dx \nonumber
   \\
   & = & \frac{k'_{j(j-1)}}{\Delta^2}\frac{\Delta^2}{2} =
   \frac{k'_{j(j-1)}}{2}  = -p5*\left[\frac{v}{h} + \frac{w_f}{h\phi}- \frac{D_{in}}{\phi^2}\frac{\partial
         \phi}{\partial x}\right]_{(j-1)} \ \ \ \ \ \mbox{for }  j > 1 \end{aligned}

For :math:`K_{N+1,N}`, there is a boundary contribution:

.. math::
   K_{N+1,N} = \frac{k'_{N+1(N)}}{2} - \frac{D_N}{\Delta \phi_{N}}

The bottom condition works if :math:`C_{bot} = h
\phi_{N+2}[c]_{ocean}`, :math:`\phi^2` is :math:`\phi_{N+1}\phi_{N+2}`
and

.. math::
   \begin{aligned}
   \left. \frac{\partial \phi}{\partial x}\right|_{bot} & = &
   \frac{\phi_{N+2} - \phi_{N}}{2\Delta};\end{aligned}

then the :math:`D_{N+1}/\phi_{N+1}/\Delta` cancels properly with the
porosity gradient. In general

.. math::
   \begin{aligned}
   \left. \frac{\partial \phi}{\partial x}\right|_{k} & = &
   \frac{\phi_{k+2} - \phi_{k}}{2\Delta}.\end{aligned}

When evaluating the integrals for the diffusion term, we will assume
that :math:`D/\phi` is constant in an element. For :math:`D_{in}/i\phi`
defined on interface points, :math:`D_1 = 0` and for :math:`j = 2,...,N`
:math:`D_j/b\phi_j = (D_{in}(j) + D_{in}(j+1))/(i\phi_j + i\phi_{j+1})`.
Then the above integrals will be modified as follows:

Compute the :math:`S_{jk}` integrals:

.. math::
   \begin{aligned}
   S_{jj} & = & - \left[\frac{D_{j-1}}{b\phi_{j-1}} \int_{x_{j-1}}^{x_j}\left( \frac{\partial
         w_j}{\partial x}\right)^2 dx + \frac{D_{j}}{b\phi_{j}} \int_{x_{j}}^{x_{j+1}}
     \left(\frac{\partial w_j}{\partial x}\right)^2 dx \right] \nonumber
   \\
   & = & -\frac{1}{\Delta}\left[\frac{D_{j-1}}{b\phi_{j-1}} + \frac{D_{j}}{b\phi{j}}\right]\ \ \mbox{for }1 < j < N+1 \nonumber \\
   S_{11} & = & \frac{s'_{11}}{\Delta}  = 0 \nonumber \\
   S_{N+1,N+1} & = & \frac{s'_{N+1,N+1}}{\Delta} = -\frac{(D_{N})}{b\phi_{N}\Delta}\end{aligned}

Compute the off-diagonal components of :math:`S_{jk}`:

.. math::
   \begin{aligned}
   S_{j(j+1)} & = & s'_{j(j+1)}\int_{x_j}^{x_{j+1}}\frac{\partial
     w_j}{\partial x} \frac{\partial w_{j+1}}{\partial x} dx =
   -\frac{s'_{j(j+1)}}{\Delta} =
   \frac{D_{j}}{b\phi_{j}\Delta} \ \ \ \mbox{for } j < N+1
   \nonumber \\
   S_{j(j-1)} & = & s'_{j(j-1)}\int_{x_{j-1}}^{x_{j}}\frac{\partial
     w_j}{\partial x} \frac{\partial w_{j-1}}{\partial x} dx =
   -\frac{s'_{j(j-1)}}{\Delta} =
   \frac{D_{j-1}}{b\phi_{j-1}} \ \ \ \mbox{for } j > 1\end{aligned}

We assume that :math:`D/\phi^2 \partial \phi/\partial x` is constant in
the element :math:`i`. If :math:`D/\phi_j` is constant, and
:math:`\partial \phi/\partial x` is constant then both the Darcy and
:math:`D` terms go as :math:`\phi^{-1}`. Then :math:`\phi = (\phi_j -
\phi_{j-1})(x-x_j)/\Delta + \phi_j` and :math:`m = (\phi_j -
\phi_{j-1})/\Delta` and :math:`b = \phi_j - mx_j`.

The first integral contribution to the Darcy term is:

.. math::
   \begin{aligned}
   K^1_{jj} & = &
   \frac{-1}{\Delta ^2}\left(\frac{w_f}{h}-\frac{D}{\phi}\frac{\partial
       \phi}{\partial x}\right)\int_{j-1}^{j}(x-x_{j-1})\frac{1}{mx
     + b}dx \nonumber \\
   & = &-\left(\frac{w_f}{h}-\frac{D}{\phi}\frac{\partial
       \phi}{\partial x}\right) \frac{1}{\Delta^2}\left[ \int_{j-1}^{j}\frac{x}{mx + b}dx - x_{j-1}\int_{j-1}^{j}\frac{1}{mx
     + b}dx  \right] \nonumber \\
   & = &- \left(\frac{w_f}{h}-\frac{D}{\phi}\frac{\partial
       \phi}{\partial x}\right) \frac{1}{\Delta^2}\left[ \frac{mx - b\log(b + mx)}{m^2} -
     x_{j-1}\frac{\log(b+mx)}{m}\right]^{x_j}_{x_{j-1}} \nonumber \\
   & = & -\left(\frac{w_f}{h}-\frac{D}{\phi}\frac{\partial
       \phi}{\partial x}\right)\frac{1}{\Delta_{\phi}}\left[ 1 + \log \left( \frac{\phi_j}{\phi_{j-1}} \right) -
     \frac{\phi_j}{\Delta_{\phi_j}}\log \left(
       \frac{\phi_j}{\phi_{j-1}} \right)\right] \nonumber \\
   & = & -\left(\frac{w_f}{h}-\frac{D}{\phi}\frac{\partial
       \phi}{\partial x}\right)\frac{1}{\Delta_{\phi}}\left[ 1 +  \frac{\phi_{j-1}}{\Delta_{\phi}}\log \left( \frac{\phi_j}{\phi_{j-1}} \right) \right]\end{aligned}

.. math::
   \begin{aligned}
   K^2_{jj} & = & \left(\frac{w_f}{h}-\frac{D}{\phi}\frac{\partial
       \phi}{\partial x}\right)\frac{1}{\Delta}\int_{x_{j}}^{x_{j+1}}\left[1 -
     \frac{(x-x_{j})}{\Delta}\right]\frac{1}{mx+b} dx \nonumber \\
   & = & \left(\frac{w_f}{h}-\frac{D}{\phi}\frac{\partial
       \phi}{\partial x}\right)\frac{1}{\Delta}\left[\frac{
       (b+m(x_j+\Delta))\log(b+mx)-mx}{\Delta
       m^2}\right]_{x_{j}}^{x_{j+1}} \nonumber \\
   & = & \left(\frac{w_f}{h}-\frac{D}{\phi}\frac{\partial
       \phi}{\partial x}\right)\frac{1}{\Delta_{\phi}}\left[ 1 -  \frac{\phi_{j+1}}{\Delta_{\phi}}\log \left( \frac{\phi_{j+1}}{\phi_{j}} \right) \right]\end{aligned}

Now :math:`m = (\phi_{j+1}-\phi_{j})/\Delta` and
:math:`b = \phi_{j+1} -
mx_{j+1}`.

Source terms :math:`q_{bot} = q_{N+1}` and :math:`q_{top} = q_{1}` (both
positive)

.. math::
   \begin{aligned}
   q_{bot} & = &\max\left[0, \left(\frac{1}{h}\frac{dh_b}{dt} +
   \frac{w_f}{h \phi_{N+1}}\right)\right]C|_{bot} +
   \frac{D_{in}}{\phi_{N+1}(g_o/h)}C|_{bot} \nonumber \\
     C|_{bot}&= & \phi_{N+1}[c]_{ocean}\end{aligned}

where :math:`g_o` is not zero.

.. math::
   \begin{aligned}
   q_{in}&  = & -\min\left[ 0, \left(\frac{1}{h}\frac{dh_t}{dt} +
   \frac{w_f}{h\phi}\right)C|_{top} \right]  \nonumber \\
   C|_{top}& = & h [c]_o\phi_{min}\end{aligned}

Calculating the low order solution: 

1) Find the lumped mass matrix
:math:`M_l = diag\{m_i\}`

.. math::
   \begin{aligned}
   m_j & = & \sum_i m_{ji} = m_{j(j+1)} + m_{j(j-1)} + m_{jj} \\
    & = & \frac{\Delta}{6} + \frac{\Delta}{6} + \frac{2\Delta}{3} =
    \Delta \ \ \ \mbox{for }1 < j < N+1 \\
   m_1 & = & m_{11} + m_{12} = \frac{\Delta}{3} + \frac{\Delta}{6} =
   \frac{\Delta}{2} \\
   m_{N+1} & = & m_{N+1,N} + m_{N+1,N+1} = \frac{\Delta}{6} + \frac{\Delta}{3} =
   \frac{\Delta}{2}\end{aligned}

2) Define artificial diffusion :math:`D_a`

.. math::
   \begin{aligned}
   d_{j,(j+1)} & = & \max\{ -k_{j(j+1)},0,-k_{(j+1)j}\} = d_{(j+1)j} \\
   d_{jj} & = & -\sum_{i\neq j} d_{ji}\end{aligned}

3) Add artificial diffusion to :math:`K`: :math:`L = K + D_a`. 

4) Solve for the low order predictor solution:

.. math::
   (M_l -\Delta t [L+S])C^{n+1} = M_l C^{n} + \Delta t q

Conservations terms for the low order solution are:

.. math::
   \begin{aligned}
   \int \left[C^{n+1} -C^{n}\right] w(x)dx & = & 
    \Delta \left[\frac{c^{n+1}_1-c^{n}_1}{2} +
      \frac{c^{n+1}_{N+1}-c^{n}_{N+1}}{2} + \sum_{j =
        2}^{N}(c^{n+1}_{j}-c^{n}_{j})\right] \nonumber \\
   &  = &   \Delta t \left[ q_{bot} +
   q_{in} + (K_{N+1,N+1}+ K_{N,N+1} )C^{n+1}_{N+1} + (K_{1,1} +
   K_{2,1})C^{n+1}_{1}\right]  \nonumber \end{aligned}

Now add the antidiffusive flux:  compute the F matrix using the low
order solution :math:`c^{n+1}`. Diagonal components are zero. For
:math:`i \neq j`

.. math::
   \begin{aligned}
   f_{ij} & = & m_{ij}\left[\frac{\Delta c_i}{\Delta t} - \frac{\Delta
       c_j}{\Delta t} + d_{ij}(c^{n+1}_i - c^{n+1}_j\right].\end{aligned}

.. _reactions:

Reaction Equations
~~~~~~~~~~~~~~~~~~

The biogeochemical reaction terms for each biogeochemical tracer (see
Table :ref:`tab-bio-tracer` for tracer definitions) are defined in
**icepack\_algae.F90** in the subroutine *algal\_dyn*. The same set of
equations is used for the bottom layer model (when ``skl_bgc`` is true) and
the multi-layer biogeochemical model (when ``z_tracers`` and ``solve_zbgc``
are true).

.. _tab-bio-tracer:

.. csv-table:: *Biogeochemical Tracers*
   :header: "Text Variable", "Variable in code", "flag", "Description", "units"
   :widths: 7, 15, 15, 15, 15

   "N (1)", "Nin(1)", "`tr_bgc_N`", "diatom", ":math:`mmol` :math:`N/m^3`"
   "N (2)", "Nin(2)", "`tr_bgc_N`", "small phytoplankton", ":math:`mmol` :math:`N/m^3`"
   "N (3)", "Nin(3)", "`tr_bgc_N`", "*Phaeocystis sp*", ":math:`mmol` :math:`N/m^3`"
   "DOC (1)", "DOCin(1)", "`tr_bgc_DOC`", "polysaccharids", ":math:`mmol` :math:`C/m^3`"
   "DOC (2)", "DOCin(2)", "`tr_bgc_DOC`", "lipids", ":math:`mmol` :math:`C/m^3`"
   "DON", "DONin(1)", "`tr_bgc_DON`", "proteins", ":math:`mmol` :math:`C/m^3`"
   "fed", "Fedin(1)", "`tr_bgc_Fe`", "dissolved iron", ":math:`\mu` :math:`Fe/m^3`"
   "fep", "Fepin(1)", "`tr_bgc_Fe`", "particulate iron", ":math:`\mu` :math:`Fe/m^3`"
   ":math:`{\mbox{NO$_3$}}`", "Nitin", "`tr_bgc_Nit`", ":math:`{\mbox{NO$_3$}}`", ":math:`mmol` :math:`N/m^3`"
   ":math:`{\mbox{NH$_4$}}`", "Amin", "`tr_bgc_Am`", ":math:`{\mbox{NH$_4$}}`", ":math:`mmol` :math:`N/m^3`"
   ":math:`{\mbox{SiO$_3$}}`", "Silin", "`tr_bgc_Sil`", ":math:`{\mbox{SiO$_2$}}`", ":math:`mmol` :math:`Si/m^3`"
   "DMSPp", "DMSPpin", "`tr_bgc_DMS`", "particulate DMSP", ":math:`mmol` :math:`S/m^3`"
   "DMSPd", "DMSPdin", "`tr_bgc_DMS`", "dissolved DMSP", ":math:`mmol` :math:`S/m^3`"
   "DMS", "DMSin", "`tr_bgc_DMS`", "DMS", ":math:`mmol` :math:`S/m^3`"
   "PON", "PON :math:`^a`", "`tr_bgc_PON`", "passive mobile tracer", ":math:`mmol` :math:`N/m^3`"
   "hum", "hum :math:`^{ab}`", "`tr_bgc_hum`", "passive sticky tracer", ":math:`mmol` :math:`/m^3`"
   "BC (1)", "zaero(1) :math:`^a`", "`tr_zaero`", "black carbon species 1", ":math:`kg` :math:`/m^3`"
   "BC (2)", "zaero(2) :math:`^a`", "`tr_zaero`", "black carbon species 2", ":math:`kg` :math:`/m^3`"
   "dust (1)", "zaero(3) :math:`^a`", "`tr_zaero`", "dust species 1", ":math:`kg` :math:`/m^3`"
   "dust (2)", "zaero(4) :math:`^a`", "`tr_zaero`", "dust species 2", ":math:`kg` :math:`/m^3`"
   "dust (3)", "zaero(5) :math:`^a`", "`tr_zaero`", "dust species 3", ":math:`kg` :math:`/m^3`"
   "dust (4)", "zaero(6) :math:`^a`", "`tr_zaero`", "dust species 4", ":math:`kg` :math:`/m^3`"

:math:`^a` not modified in *algal_dyn*

:math:`^b` may be in C or N units depending on the ocean concentration

The biochemical reaction term for each algal species has the form:

.. math::
   \Delta {\mbox{N}}/dt = R_{{\mbox{N}}} = \mu (1- f_{graze} - f_{res}) - M_{ort}

where :math:`\mu` is the algal growth rate, :math:`M_{ort}` is a
mortality loss, :math:`f_{graze}` is the fraction of algal growth that
is lost to predatory grazing, and :math:`f_{res}` is the fraction of
algal growth lost to respiration. Algal mortality is temperature
dependent and limited by a maximum loss rate fraction (:math:`l_{max}`):

.. math::
   M_{ort} = \min( l_{max}[{\mbox{N}}], m_{pre} \exp\{m_{T}(T-T_{max})\}[{\mbox{N}}])

Note, :math:`[\cdot]` denotes brine concentration.

Nitrate and ammonium reaction terms are given by

.. math::

   \begin{aligned}
   \Delta{\mbox{NO$_3$}}/dt & = & R_{{\mbox{NO$_3$}}} =  [{\mbox{NH$_4$}}] k_{nitr}- U^{tot}_{{\mbox{NO$_3$}}} \nonumber \\
   \Delta{\mbox{NH$_4$}}/dt & = & R_{{\mbox{NH$_4$}}} = -[{\mbox{NH$_4$}}] k_{nitr} -U^{tot}_{{\mbox{NH$_4$}}} +
   (f_{ng}f_{graze}(1-f_{gs})+f_{res})\mu^{tot} \nonumber \\
    &  +  & f_{nm} M_{ort}
   \nonumber \\
                        & = &  -[{\mbox{NH$_4$}}]k_{nitr} -U^{tot}_{{\mbox{NH$_4$}}} + N_{remin}\end{aligned}

where the uptake :math:`U^{tot}` and algal growth :math:`\mu^{tot}` are
accumulated totals for all algal species. :math:`k_{nitr}` is the
nitrification rate and :math:`f_{ng}` and :math:`f_{nm}` are the
fractions of grazing and algal mortality that are remineralized to
ammonium and :math:`f_{gs}` is the fraction of grazing spilled or lost.
Algal growth and nutrient uptake terms are described in more detail in
:ref:`growth-uptake`.

Dissolved organic nitrogen satisfies the equation

.. math::

   \begin{aligned}
   \Delta {\mbox{DON}}/dt & = & R_{{\mbox{DON}}} = f_{dg}f_{gs}f_{graze}\mu^{tot} - [{\mbox{DON}}]k_{nb}\end{aligned}

With a loss from bacterial degration (rate :math:`k_{nb}`) and a gain
from spilled grazing that does not enter the :math:`{\mbox{NH$_4$}}`
pool.

A term Z\ :math:`_{oo}` closes the nitrogen cycle by summing all the
excess nitrogen removed as zooplankton/bacteria in a timestep. This term
is not a true tracer, i.e. not advected horizontally with the ice
motion, but provides a diagnostic comparison of the amount of :math:`N`
removed biogeochemically from the ice
:math:`{\mbox{N}}`-:math:`{\mbox{NO$_3$}}`-:math:`{\mbox{NH$_4$}}`-:math:`{\mbox{DON}}`
cycle at each timestep.

.. math::

   \begin{aligned}
   \mbox{Z}_{oo} & = & [(1-f_{ng}(1-f_{gs}) - f_{dg}f_{gs}]f_{graze}\mu^{tot}dt + (1-f_{nm})M_{ort}dt  +
   [{\mbox{DON}}]k_{nb}dt \nonumber\end{aligned}

Dissolved organic carbon may be divided into polysaccharids and lipids.
Parameters are two dimensional (indicated by superscript :math:`i`) with
index 1 corresponding to polysaccharids and index 2 appropriate for
lipids. The :math:`{\mbox{DOC}}^i` equation is:

.. math::

   \begin{aligned}
   \Delta {\mbox{DOC}}^i/dt & = & R_{{\mbox{DOC}}} = f^i_{cg}f_{ng}\mu^{tot} + R^i_{c:n}M_{ort}-[{\mbox{DOC}}]k^i_{cb}\end{aligned}

Silicate has no biochemical source terms within the ice and is lost only
through algal uptake:

.. math::

   \begin{aligned}
   \Delta {\mbox{SiO$_3$}}/dt & = & R_{{\mbox{SiO$_3$}}} = -U_{{\mbox{SiO$_3$}}}^{tot}\end{aligned}

Dissolved iron has algal uptake and remineralization pathways. In
addition, :math:`{\mbox{fed}}` may be converted to or released from the
particulate iron pool depending on the dissolve iron
(:math:`{\mbox{fed}}`) to polysaccharid (:math:`{\mbox{DOC}}(1)`)
concentration ratio. If this ratio exceeds a maximum value
:math:`r^{max}_{fed:doc}` then the change in concentration for dissolved
and particulate iron is

.. math::

   \begin{aligned}
   \Delta_{fe}{\mbox{fed}}/dt & = & -[{\mbox{fed}}]/\tau_{fe} \nonumber \\
   \Delta_{fe}{\mbox{fep}}/dt & = & [{\mbox{fed}}]/\tau_{fe}\end{aligned}

For values less than :math:`r^{max}_{fed:doc}`

.. math::

   \begin{aligned}
   \Delta_{fe}{\mbox{fed}}/dt & = & [{\mbox{fep}}]/\tau_{fe} \nonumber \\
   \Delta_{fe}{\mbox{fep}}/dt & = & -[{\mbox{fep}}]/\tau_{fe}\end{aligned}

Very long timescales :math:`\tau_{fe}` will remove this source/sink
term. The default value is currently set at 3065 days to turn off this
dependency (any large number will do to turn it off). 
61-65 days is a more realistic option (Parekh et al., 2004).

The full equation for :math:`{\mbox{fed}}` including uptake and
remineralization is

.. math::

   \begin{aligned}
   \Delta {\mbox{fed}}/dt & = & R_{{\mbox{fed}}} = -U^{tot}_{{\mbox{fed}}} + f_{fa}R_{fe:n}N_{remin}
   + \Delta_{fe}{\mbox{fed}}/dt\end{aligned}

Particulate iron also includes a source term from algal mortality and
grazing that is not immediately bioavailable. The full equation for
:math:`{\mbox{fep}}` is

.. math::

   \begin{aligned}
   \Delta {\mbox{fep}}/dt & = & R_{{\mbox{fep}}} =  R_{fe:n}[\mbox{Z}_{oo}/dt + (1-f_{fa})]N_{remin}
   + \Delta_{fe}{\mbox{fep}}/dt\end{aligned}

The sulfur cycle includes :math:`{\mbox{DMS}}` and dissolved DMSP
(:math:`{\mbox{DMSPd}}`). Particulate DMSP is assumed to be proportional
to the algal concentration, i.e.
:math:`{\mbox{DMSPp}}= R^i_{s:n}{\mbox{N}}^i` for algal species
:math:`i`. For :math:`{\mbox{DMSP}}` and :math:`{\mbox{DMS}}`,

.. math::

   \begin{aligned}
   \Delta {\mbox{DMSPd}}/dt & = & R_{{\mbox{DMSPd}}} = R_{s:n}[ f_{sr}f_{res}\mu^{tot}
   +f_{nm}M_{ort} ] - [{\mbox{DMSPd}}]/\tau_{dmsp} \nonumber \\
   \Delta {\mbox{DMS}}/dt & = & R_{{\mbox{DMS}}} =  y_{dms}[{\mbox{DMSPd}}]/\tau_{dmsp} - [{\mbox{DMS}}]/\tau_{dms}\end{aligned}

See :ref:`tuning` for a more complete list and description of biogeochemical parameters.

.. _growth-uptake:

Algal growth and nutrient uptake
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Nutrient limitation terms are defined in the simplest ecosystem for
:math:`{\mbox{NO$_3$}}`. If the appropriate tracer flags are true, then
limitation terms may also be found for :math:`{\mbox{NH$_4$}}`,
:math:`{\mbox{SiO$_3$}}`, and :math:`{\mbox{fed}}`

.. math::
   \begin{aligned}
   {\mbox{NO$_3$}}_{lim} & = & \frac{[{\mbox{NO$_3$}}]}{[{\mbox{NO$_3$}}] + K_{{\mbox{NO$_3$}}}} \nonumber \\
   {\mbox{NH$_4$}}_{lim} & = & \frac{[{\mbox{NH$_4$}}]}{[{\mbox{NH$_4$}}] + K_{{\mbox{NH$_4$}}}}\nonumber \\
   N_{lim} & = &\mbox{min}(1,{\mbox{NO$_3$}}_{lim} + {\mbox{NH$_4$}}_{lim}) \nonumber \\
   {\mbox{SiO$_3$}}_{lim} & = & \frac{[{\mbox{SiO$_3$}}]}{[{\mbox{SiO$_3$}}] + K_{{\mbox{SiO$_3$}}}} \nonumber \\
   {\mbox{fed}}_{lim} & = & \frac{[{\mbox{fed}}]}{[{\mbox{fed}}] + K_{{\mbox{fed}}}} \end{aligned}

Light limitation :math:`L_{lim}` is defined in the following way:
:math:`I
_{sw}(z)` (in :math:`W/m^2`) is the shortwave radiation at the ice level
and the optical depth is proportional to the chlorophyll concentration,
:math:`op_{dep} =` ``chlabs`` [Chl*a*]. If ( :math:`op_{dep} > op_{min}`) then

.. math::
   I_{avg} = I_{sw}(1- \exp(-op_{dep}))/op_{dep}

otherwise :math:`I_{avg} = I_{sw}`.

.. math::
   L_{lim} = (1 - \exp(-\alpha I_{avg}))\exp(-\beta I_{avg})

The maximal algal growth rate before limitation is

.. math::
   \begin{aligned}
   \mu_o & = & \mu_{max}\exp(\mu_T\Delta T)f_{sal}[{\mbox{N}}] \\ 
   \mu' & = & min(L_{lim},N_{lim},{\mbox{SiO$_3$}}_{lim},{\mbox{fed}}_{lim}) \mu_o\end{aligned}

where :math:`\mu'` is the initial estimate of algal growth rate for a
given algal species and :math:`\Delta T` is the difference between the
local tempurature and the maximum (in this case
T\ :math:`_{max} = 0^o`\ C).

The initial estimate of the uptake rate for silicate and iron is

.. math::
   \begin{aligned}
   \tilde{U}_{{\mbox{SiO$_3$}}} & = & R_{si:n}\mu' \\
   \tilde{U}_{{\mbox{fed}}} & = & R_{fe:n}\mu' \end{aligned}

For nitrogen uptake, we assume that ammonium is preferentially acquired
by algae. To determine the nitrogen uptake needed for each algal growth
rate of :math:`\mu`, first determine the "potential" uptake rate of
ammonium:

.. math::
   U'_{{\mbox{NH$_4$}}} = {\mbox{NH$_4$}}_{lim}\mu_o

Then

.. math::
   \begin{aligned}
   \tilde{U}_{{\mbox{NH$_4$}}} & = & \min(\mu', U'_{{\mbox{NH$_4$}}}) \\
   \tilde{U}_{{\mbox{NO$_3$}}} & = & \mu' - \tilde{U}_{{\mbox{NH$_4$}}}\end{aligned}

We require that each rate not exceed a maximum loss rate
:math:`l_{max}/dt`. This is particularly important when multiple species
are present. In this case, the accumulated uptake rate for each nutrient
is found and the fraction (:math:`fU^i`) of uptake due to algal species
:math:`i` is saved. Then the total uptake rate is compared with the
maximum loss condition. For example, the net uptake of nitrate when
there are three algal species is

.. math::
   \tilde{U}^{tot}_{{\mbox{NO$_3$}}} = \sum^{3}_{i=1}\tilde{U}^i_{{\mbox{NO$_3$}}}\ \ \ .

Then the uptake fraction for species :math:`i` and the adjusted total
uptake is

.. math::
   \begin{aligned}
   fU^i_{{\mbox{NO$_3$}}} & = & \frac{\tilde{U}^i_{{\mbox{NO$_3$}}}}{\tilde{U}^{tot}_{{\mbox{NO$_3$}}}} \nonumber \\
   U^{tot}_{{\mbox{NO$_3$}}} & = & \min(\tilde{U}^{tot}_{{\mbox{NO$_3$}}}, l_{max}[{\mbox{NO$_3$}}]/dt)\end{aligned}

Now, for each algal species the nitrate uptake is

.. math::
   U^i_{{\mbox{NO$_3$}}} = fU^i_{{\mbox{NO$_3$}}} U^{tot}_{{\mbox{NO$_3$}}}

Similar expressions are found for all potentially limiting nutrients.
Then the true growth rate for each algal species :math:`i` is

.. math::
   \begin{aligned}
   \mu^i & = & \min(U^i_{{\mbox{SiO$_3$}}}/R_{si:n}, U^i_{{\mbox{NO$_3$}}} + U^i_{{\mbox{NH$_4$}}}, U^i_{{\mbox{fed}}}/R_{fe:n})\end{aligned}

Preferential ammonium uptake is assumed once again and the remaining
nitrogen is taken from the nitrate pool.





