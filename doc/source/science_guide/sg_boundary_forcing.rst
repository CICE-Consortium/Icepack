:tocdepth: 3

.. _boundary_forcing:

Atmosphere and ocean boundary forcing
=====================================

.. _tab-flux-cpl:

.. csv-table:: *External forcing data that are relevant to Icepack*
   :header: "Variable", "Description", "External Interactions"
   :widths: 10, 25, 25
     
   ":math:`z_o`", "Atmosphere level height", "From *atmosphere model* to *sea ice model*"
   ":math:`\vec{U}_a`", "Wind velocity", "From *atmosphere model* to *sea ice model*"
   ":math:`Q_a`", "Specific humidity", "From *atmosphere model* to *sea ice model*"
   ":math:`\rho_a`", "Air density", "From *atmosphere model* to *sea ice model*"
   ":math:`\Theta_a`", "Air potential temperature", "From *atmosphere model* to *sea ice model*"
   ":math:`T_a`", "Air temperature", "From *atmosphere model* to *sea ice model*"
   ":math:`F_{sw\downarrow}`", "Incoming shortwave radiation (4 bands)", "From *atmosphere model* to *sea ice model*"
   ":math:`F_{L\downarrow}`", "Incoming longwave radiation", "From *atmosphere model* to *sea ice model*"
   ":math:`F_{rain}`", "Rainfall rate", "From *atmosphere model* to *sea ice model*"
   ":math:`F_{snow}`", "Snowfall rate", "From *atmosphere model* to *sea ice model*"
   ":math:`F_{frzmlt}`", "Freezing/melting potential", "From *ocean model* to *sea ice model*"
   ":math:`T_w`", "Sea surface temperature", "From *ocean model* to *sea ice model*"
   ":math:`S`", "Sea surface salinity", "From *ocean model* to *sea ice model*"
   ":math:`\nabla H_o`", "Sea surface slope", "From *ocean model* via flux coupler to *sea ice model*"
   ":math:`\vec{U}_w`", "Surface ocean currents", "From *ocean model* to *sea ice model* (available in Icepack driver, not used directly in column physics)"
   ":math:`\vec{\tau}_a`", "Wind stress", "From *sea ice model* to *atmosphere model*"
   ":math:`F_s`", "Sensible heat flux", "From *sea ice model* to *atmosphere model*"
   ":math:`F_l`", "Latent heat flux", "From *sea ice model* to *atmosphere model*"
   ":math:`F_{L\uparrow}`", "Outgoing longwave radiation", "From *sea ice model* to *atmosphere model*"
   ":math:`F_{evap}`", "Evaporated water", "From *sea ice model* to *atmosphere model*"
   ":math:`\alpha`", "Surface albedo (4 bands)", "From *sea ice model* to *atmosphere model*"
   ":math:`T_{sfc}`", "Surface temperature", "From *sea ice model* to *atmosphere model*"
   ":math:`F_{sw\Downarrow}`", "Penetrating shortwave radiation", "From *sea ice model* to *ocean model*"
   ":math:`F_{water}`", "Fresh water flux", "From *sea ice model* to *ocean model*"
   ":math:`F_{hocn}`", "Net heat flux to ocean", "From *sea ice model* to *ocean model*"
   ":math:`F_{salt}`", "Salt flux", "From *sea ice model* to *ocean model*"
   ":math:`\vec{\tau}_w`", "Ice-ocean stress", "From *sea ice model* to *ocean model*"
   ":math:`F_{bio}`", "Biogeochemical fluxes", "From *sea ice model* to *ocean model*"
   ":math:`a_{i}`", "Ice fraction", "From *sea ice model* to both *ocean and atmosphere models*"
   ":math:`T^{ref}_{a}`", "2m reference temperature (diagnostic)", "From *sea ice model* to both *ocean and atmosphere models*"
   ":math:`Q^{ref}_{a}`", "2m reference humidity (diagnostic)", "From *sea ice model* to both *ocean and atmosphere models*"
   ":math:`F_{swabs}`", "Absorbed shortwave (diagnostic)", "From *sea ice model* to both *ocean and atmosphere models*"
   ":math:`E(f)`", "Wave spectrum as a function of frequency", "From *ocean surface wave model* to *sea ice model*"

The ice fraction :math:`a_i` (aice) is the total fractional ice
coverage of a grid cell. That is, in each cell,

.. math::
   \begin{array}{cl}
                  a_{i}=0 & \mbox{if there is no ice} \\ 
                  a_{i}=1 & \mbox{if there is no open water} \\ 
                  0<a_{i}<1 & \mbox{if there is both ice and open water,}
   \end{array}

where :math:`a_{i}` is the sum of fractional ice areas for each category
of ice. The ice fraction is used by the flux coupler to merge fluxes
from the sea ice model with fluxes from the other earth system components. For example,
the penetrating shortwave radiation flux, weighted by :math:`a_i`, is
combined with the net shortwave radiation flux through ice-free leads,
weighted by (:math:`1-a_i`), to obtain the net shortwave flux into the
ocean over the entire grid cell. The CESM flux coupler requires the fluxes to
be divided by the total ice area so that the ice and land models are
treated identically (land also may occupy less than 100% of an
atmospheric grid cell). These fluxes are "per unit ice area" rather than
"per unit grid cell area."

In some coupled climate models (for example, recent versions of the U.K.	
Hadley Centre model) the surface air temperature and fluxes are computed	
within the atmosphere model and are passed to CICE for use in the column physics. In this case the	
logical parameter ``calc_Tsfc`` in *ice_therm_vertical* is set to false.	
The fields ``fsurfn`` (the net surface heat flux from the atmosphere), ``flatn``	
(the surface latent heat flux), and ``fcondtopn`` (the conductive flux at	
the top surface) for each ice thickness category are copied or derived	
from the input coupler fluxes and are passed to the thermodynamic driver	
subroutine, *thermo_vertical*. At the end of the time step, the surface	
temperature and effective conductivity (i.e., thermal conductivity	
divided by thickness) of the top ice/snow layer in each category are	
returned to the atmosphere model via the coupler. Since the ice surface	
temperature is treated explicitly, the effective conductivity may need	
to be limited to ensure stability. As a result, accuracy may be	
significantly reduced, especially for thin ice or snow layers. A more	
stable and accurate procedure would be to compute the temperature	
profiles for both the atmosphere and ice, together with the surface	
fluxes, in a single implicit calculation. This was judged impractical,	
however, given that the atmosphere and sea ice models generally exist on	
different grids and/or processor sets.

.. _atmo:

Atmosphere
----------

The wind velocity, specific humidity, air density and potential
temperature at the given level height :math:`z_\circ` are used to
compute transfer coefficients used in formulas for the surface wind
stress and turbulent heat fluxes :math:`\vec\tau_a`, :math:`F_s`, and
:math:`F_l`, as described below. The sensible and latent heat fluxes,
:math:`F_s` and :math:`F_l`, along with shortwave and longwave
radiation, :math:`F_{sw\downarrow}`, :math:`F_{L\downarrow}`
and :math:`F_{L\uparrow}`, are included in the flux balance that
determines the ice or snow surface temperature when ``calc_Tsfc`` is true.
As described in the :ref:`thermo` section, these fluxes depend nonlinearly
on the ice surface temperature :math:`T_{sfc}`. The balance
equation is iterated until convergence, and the resulting fluxes and
:math:`T_{sfc}` are then passed to the flux coupler.

The snowfall precipitation rate (provided as liquid water equivalent and
converted by the ice model to snow depth) also contributes to the heat
and water mass budgets of the ice layer. Melt ponds generally form on
the ice surface in the Arctic and refreeze later in the fall, reducing
the total amount of fresh water that reaches the ocean and altering the
heat budget of the ice; this version includes two new melt pond
parameterizations. Rain and all melted snow end up in the ocean.

Wind stress and transfer coefficients for the
turbulent heat fluxes are computed in subroutine
*atmo\_boundary\_layer* following :cite:`Kauffman02`, with additions and changes as detailed in Appendix A of :cite:`Roberts15` for high frequency coupling (namelist variable ``highfreq``).
The resulting equations are provided here.

The wind stress and turbulent heat flux calculation accounts for both
stable and unstable atmosphere–ice boundary layers. Define the
"stability"

.. math::
   \Upsilon = {\kappa g z_\circ\over u^{*2}}
   \left({\Theta^*\over\Theta_a\left(1+0.606Q_a\right)}  +
   {Q^*\over 1/0.606 + Q_a}\right),
   :label: upsilon

where :math:`\kappa` is the von Karman constant, :math:`g` is
gravitational acceleration, and :math:`u^*`, :math:`\Theta^*` and
:math:`Q^*` are turbulent scales for velocity difference, temperature, and humidity,
respectively, given the ice velocity :math:`\vec{U}_i`:

.. math::
   \begin{aligned}
   u^*&=&c_u\;\textrm{max}\left(U_{\Delta\textrm{min}}, \left|\vec{U}_a - \vec{U}_i \right|\right), \\
   \Theta^*&=& c_\theta\left(\Theta_a-T_{sfc}\right), \\
   Q^*&=&c_q\left(Q_a-Q_{sfc}\right).
   \end{aligned}
   :label: stars

Within the :math:`u^*` expression, :math:`U_{\Delta\textrm{min}}` is the minimum allowable value of :math:`|\vec{U}_{a} - \vec{U}_{i}|` , which is set to of 0.5 m/s for high frequency coupling (``highfreq`` =.true.). 
When high frequency coupling is turned off (``highfreq`` =.false.), it is assumed in equation :eq:`stars` that:

.. math::
 \vec{U}_{a} - \vec{U}_{i} \approx  \vec{U}_{a} 
 :label: lowfreq 

and a higher threshold is taken for :math:`U_{\Delta\textrm{min}}` of 1m/s. Equation :eq:`lowfreq` is a poor assumption when resolving inertial oscillations in ice-ocean configurations where the ice velocity vector may make a complete rotation over a period of :math:`\ge` 11.96 hours, as discussed in :cite:`Roberts15`.
However, :eq:`lowfreq`  is acceptable for low frequency ice-ocean coupling on the order of a day or more, when transient ice-ocean Ekman transport is effectively filtered from the model solution.
For the :math:`\Theta^*` and :math:`Q^*` terms in :eq:`stars`, :math:`T_{sfc}` and :math:`Q_{sfc}` are the surface temperature and specific
humidity, respectively.  The latter is calculated by assuming a saturated
surface, as described in the :ref:`sfc-forcing` section.

Neglecting form drag, the exchange coefficients :math:`c_u`,
:math:`c_\theta` and :math:`c_q` are initialized as

.. math:: 
   \kappa\over \ln(z_{ref}/z_{ice})
   :label: kappa

and updated during a short iteration, as they depend upon the turbulent
scales. The number of iterations is set by the namelist variable
``natmiter``, nominally set to five but sometimes increased by users employing the ``highfreq`` option.
Here, :math:`z_{ref}` is a reference height of 10m and
:math:`z_{ice}` is the roughness length scale for the given
sea ice category. :math:`\Upsilon` is constrained to have magnitude less
than 10. Further, defining
:math:`\chi = \left(1-16\Upsilon\right)^{0.25}` and :math:`\chi \geq 1`,
the "integrated flux profiles" for momentum and stability in the
unstable (:math:`\Upsilon <0`) case are given by

.. math::
   \begin{aligned}
   \psi_m = &\mbox{}&2\ln\left[0.5(1+\chi)\right] +
            \ln\left[0.5(1+\chi^2)\right] -2\tan^{-1}\chi +
            {\pi\over 2}, \\
   \psi_s = &\mbox{}&2\ln\left[0.5(1+\chi^2)\right].\end{aligned}
   :label: psi1

In a departure from the parameterization used in
:cite:`Kauffman02`, we use profiles for the stable case
following :cite:`Jordan99`,

.. math::
   \psi_m = \psi_s = -\left[0.7\Upsilon + 0.75\left(\Upsilon-14.3\right)
            \exp\left(-0.35\Upsilon\right) + 10.7\right].
   :label: psi2

The coefficients are then updated as

.. math::
   \begin{aligned}
   c_u^\prime&=&{c_u\over 1+c_u\left(\lambda-\psi_m\right)/\kappa} \\
   c_\theta^\prime&=& {c_\theta\over 1+c_\theta\left(\lambda-\psi_s\right)/\kappa}\\
   c_q^\prime&=&c_\theta^\prime\end{aligned}
   :label: coeff1

where :math:`\lambda = \ln\left(z_\circ/z_{ref}\right)`. The
first iteration ends with new turbulent scales from
equations :eq:`stars`. After ``natmiter`` iterations the latent and sensible
heat flux coefficients are computed, along with the wind stress:

.. math::
   \begin{aligned}
   C_l&=&\rho_a \left(L_{vap}+L_{ice}\right) u^* c_q \\
   C_s&=&\rho_a c_p u^* c_\theta^* + 1 \\
   \vec{\tau}_a&=&{\rho_a (u^{*})^2 \left( \vec{U}_{a} - \vec{U}_{i} \right) \over  \left| \vec{U}_{a} - \vec{U}_{i} \right|}
   \end{aligned}
   :label: coeff2

where :math:`L_{vap}` and :math:`L_{ice}` are
latent heats of vaporization and fusion, :math:`\rho_a` is the density
of air and :math:`c_p` is its specific heat. Again following
:cite:`Jordan99`, we have added a constant to the sensible
heat flux coefficient in order to allow some heat to pass between the
atmosphere and the ice surface in stable, calm conditions. 
For the atmospheric stress term in :eq:`coeff2`, we make the assumption in :eq:`lowfreq` when ``highfreq`` =.false..

The atmospheric reference temperature :math:`T_a^{ref}` is computed from
:math:`T_a` and :math:`T_{sfc}` using the coefficients
:math:`c_u`, :math:`c_\theta` and :math:`c_q`. Although the sea ice
model does not use this quantity, it is convenient for the ice model to
perform this calculation. The atmospheric reference temperature is
returned to the flux coupler as a climate diagnostic. The same is true
for the reference humidity, :math:`Q_a^{ref}`.

Additional details about the latent and sensible heat fluxes and other
quantities referred to here can be found in
the :ref:`sfc-forcing` section.

.. _ocean:

Ocean
-----

New sea ice forms when the ocean temperature drops below its freezing
temperature. In the Bitz and Lipscomb thermodynamics,
:cite:`Bitz99` :math:`T_f=-\mu S`, where :math:`S` is the
seawater salinity and :math:`\mu=0.054^\circ`/ppt is the ratio of the
freezing temperature of brine to its salinity (linear liquidus
approximation). For the mushy thermodynamics, :math:`T_f` is given by a
piecewise linear liquidus relation. The ocean model calculates the new
ice formation; if the freezing/melting potential
:math:`F_{frzmlt}` is positive, its value represents a certain
amount of frazil ice that has formed in one or more layers of the ocean
and floated to the surface. (The ocean model assumes that the amount of
new ice implied by the freezing potential actually forms.)

If :math:`F_{frzmlt}` is negative, it is used to heat already
existing ice from below. In particular, the sea surface temperature and
salinity are used to compute an oceanic heat flux :math:`F_w`
(:math:`\left|F_w\right| \leq \left|F_{frzmlt}\right|`) which
is applied at the bottom of the ice. The portion of the melting
potential actually used to melt ice is returned to the coupler in
:math:`F_{hocn}`. The ocean model adjusts its own heat budget
with this quantity, assuming that the rest of the flux remained in the
ocean.

In addition to runoff from rain and melted snow, the fresh water flux
:math:`F_{water}` includes ice melt water from the top surface
and water frozen (a negative flux) or melted at the bottom surface of
the ice. This flux is computed as the net change of fresh water in the
ice and snow volume over the coupling time step, excluding frazil ice
formation and newly accumulated snow. Setting the namelist option
``update_ocn_f`` to true causes frazil ice to be included in the fresh
water and salt fluxes.

There is a flux of salt into the ocean under melting conditions, and a
(negative) flux when sea water is freezing. However, melting sea ice
ultimately freshens the top ocean layer, since the ocean is much more
saline than the ice. The ice model passes the net flux of salt
:math:`F_{salt}` to the flux coupler, based on the net change
in salt for ice in all categories. In the present configuration,
``ice_ref_salinity`` is used for computing the salt flux, although the ice
salinity used in the thermodynamic calculation has differing values in
the ice layers.

A fraction of the incoming shortwave :math:`F_{sw\Downarrow}`
penetrates the snow and ice layers and passes into the ocean, as
described in the :ref:`sfc-forcing` section.

A thermodynamic slab ocean mixed-layer parameterization is available 
in **icepack\_ocean.F90** and can be run in the full CICE configuration.
The turbulent fluxes are computed above the water surface using the same
parameterizations as for sea ice, but with parameters appropriate for
the ocean. The surface flux balance takes into account the turbulent
fluxes, oceanic heat fluxes from below the mixed layer, and shortwave
and longwave radiation, including that passing through the sea ice into
the ocean. If the resulting sea surface temperature falls below the
salinity-dependent freezing point, then new ice (frazil) forms.
Otherwise, heat is made available for melting the ice.

.. _formdrag:

Variable exchange coefficients
------------------------------

In the default configuration, atmospheric and oceanic neutral drag
coefficients (:math:`c_u` and :math:`c_w`) are assumed constant in time
and space. These constants are chosen to reflect friction associated
with an effective sea ice surface roughness at the ice–atmosphere and
ice–ocean interfaces. Sea ice (in both Arctic and Antarctic) contains
pressure ridges as well as floe and melt pond edges that act as discrete
obstructions to the flow of air or water past the ice, and are a source
of form drag. Following :cite:`Tsamados14` and based on
recent theoretical developments :cite:`Lupkes12,Lu11`, the
neutral drag coefficients can now be estimated from properties of the
ice cover such as ice concentration, vertical extent and area of the
ridges, freeboard and floe draft, and size of floes and melt ponds. The
new parameterization allows the drag coefficients to be coupled to the
sea ice state and therefore to evolve spatially and temporally. This
parameterization is contained in the subroutine *neutral\_drag\_coeffs*
and is accessed by setting ``formdrag`` = true in the namelist.
(Note:  see also :ref:`bugs`.)

Following :cite:`Tsamados14`, consider the general case of
fluid flow obstructed by N randomly oriented obstacles of height
:math:`H` and transverse length :math:`L_y`, distributed on a domain
surface area :math:`S_T`. Under the assumption of a logarithmic fluid
velocity profile, the general formulation of the form drag coefficient
can be expressed as

.. math:: 
   C_d=\frac{N c S_c^2 \gamma L_y  H}{2 S_T}\left[\frac{\ln(H/z_0)}{\ln(z_{ref}/z_0)}\right]^2,
   :label: formdrag

where :math:`z_0` is a roughness length parameter at the top or bottom
surface of the ice, :math:`\gamma` is a geometric factor, :math:`c` is
the resistance coefficient of a single obstacle, and :math:`S_c` is a
sheltering function that takes into account the shielding effect of the
obstacle,

.. math:: 
   S_{c}=\left(1-\exp(-s_l D/H)\right)^{1/2},
   :label: shelter

with :math:`D` the distance between two obstacles and :math:`s_l` an
attenuation parameter.

As in the original drag formulation in CICE (:ref:`atmo` and
:ref:`ocean` sections), :math:`c_u` and :math:`c_w` along with the transfer
coefficients for sensible heat, :math:`c_{\theta}`, and latent heat,
:math:`c_{q}`, are initialized to a situation corresponding to neutral
atmosphere–ice and ocean–ice boundary layers. The corresponding neutral
exchange coefficients are then replaced by coefficients that explicitly
account for form drag, expressed in terms of various contributions as

.. math::
   \tt{Cdn\_atm}  = \tt{Cdn\_atm\_rdg} + \tt{Cdn\_atm\_floe} + \tt{Cdn\_atm\_skin} + \tt{Cdn\_atm\_pond} ,
   :label: Cda

.. math::
   \tt{Cdn\_ocn}  =  \tt{Cdn\_ocn\_rdg} + \tt{Cdn\_ocn\_floe} + \tt{Cdn\_ocn\_skin}. 
   :label: Cdw

The contributions to form drag from ridges (and keels underneath the
ice), floe edges and melt pond edges can be expressed using the general
formulation of equation :eq:`formdrag` (see :cite:`Tsamados14` for
details). Individual terms in equation :eq:`Cdw` are fully described in
:cite:`Tsamados14`. Following :cite:`Arya75`
the skin drag coefficient is parametrized as

.. math:: 
   { \tt{Cdn\_(atm/ocn)\_skin}}=a_{i} \left(1-m_{(s/k)} \frac{H_{(s/k)}}{D_{(s/k)}}\right)c_{s(s/k)}, \mbox{       if  $\displaystyle\frac{H_{(s/k)}}{D_{(s/k)}}\ge\frac{1}{m_{(s/k)}}$,}
   :label: skindrag

where :math:`m_s` (:math:`m_k`) is a sheltering parameter that depends
on the average sail (keel) height, :math:`H_s` (:math:`H_k`), but is
often assumed constant, :math:`D_s` (:math:`D_k`) is the average
distance between sails (keels), and :math:`c_{ss}` (:math:`c_{sk}`) is
the unobstructed atmospheric (oceanic) skin drag that would be attained
in the absence of sails (keels) and with complete ice coverage,
:math:`a_{ice}=1`.

Calculation of equations :eq:`formdrag` – :eq:`skindrag` requires that small-scale geometrical
properties of the ice cover be related to average grid cell quantities
already computed in the sea ice model. These intermediate quantities are
briefly presented here and described in more detail in
:cite:`Tsamados14`. The sail height is given by

.. math:: 
   H_{s} = \displaystyle 2\frac{v_{rdg}}{a_{rdg}}\left(\frac{\alpha\tan \alpha_{k} R_d+\beta \tan \alpha_{s} R_h}{\phi_r\tan \alpha_{k} R_d+\phi_k \tan \alpha_{s} R_h^2}\right),
   :label: Hs

and the distance between sails\ 

.. math:: 
   D_{s} = \displaystyle 2 H_s\frac{a_{i}}{a_{rdg}} \left(\frac{\alpha}{\tan \alpha_s}+\frac{\beta}{\tan \alpha_k}\frac{R_h}{R_d}\right),
   :label: Ds

where :math:`0<\alpha<1` and :math:`0<\beta<1` are weight functions,
:math:`\alpha_{s}` and :math:`\alpha_{k}` are the sail and keel slope,
:math:`\phi_s` and :math:`\phi_k` are constant porosities for the sails
and keels, and we assume constant ratios for the average keel depth and
sail height (:math:`H_k/H_s=R_h`) and for the average distances between
keels and between sails (:math:`D_k/D_s=R_d`). With the assumption of
hydrostatic equilibrium, the effective ice plus snow freeboard is
:math:`H_{f}=\bar{h_i}(1-\rho_i/\rho_w)+\bar{h_s}(1-\rho_s/\rho_w)`,
where :math:`\rho_i`, :math:`\rho_w` and :math:`\rho_s` are
respectively the densities of sea ice, water and snow, :math:`\bar{h_i}`
is the mean ice thickness and :math:`\bar{h_s}` is the mean snow
thickness (means taken over the ice covered regions). For the melt pond
edge elevation we assume that the melt pond surface is at the same level
as the ocean surface surrounding the floes
:cite:`Flocco07,Flocco10,Flocco12` and use the simplification
:math:`H_p = H_f`. Finally to estimate the typical floe size
:math:`L_A`, distance between floes, :math:`D_F`, and melt pond size,
:math:`L_P` we use the parameterizations of :cite:`Lupkes12`
to relate these quantities to the ice and pond concentrations. All of
these intermediate quantities are available for output, along
with ``Cdn_atm``, ``Cdn_ocn`` and the ratio ``Cdn_atm_ratio_n`` between the
total atmospheric drag and the atmospheric neutral drag coefficient.

We assume that the total neutral drag coefficients are thickness
category independent, but through their dependance on the diagnostic
variables described above, they vary both spatially and temporally. The
total drag coefficients and heat transfer coefficients will also depend
on the type of stratification of the atmosphere and the ocean, and we
use the parameterization described in the :ref:`atmo` section that accounts
for both stable and unstable atmosphere–ice boundary layers. In contrast
to the neutral drag coefficients the stability effect of the atmospheric
boundary layer is calculated separately for each ice thickness category.

The transfer coefficient for oceanic heat flux to the bottom of the ice
may be varied based on form drag considerations by setting the namelist
variable ``fbot_xfer_type`` to ``Cdn_ocn``; this is recommended when using
the form drag parameterization. The default value of the transfer
coefficient is 0.006 (``fbot_xfer_type = ’constant’``).
