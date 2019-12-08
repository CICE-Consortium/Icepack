:tocdepth: 3

.. _model_comp:

Ice thickness distribution
=============================

The Arctic and Antarctic sea ice packs are mixtures of open water, thin
first-year ice, thicker multiyear ice, and thick pressure ridges. The
thermodynamic and dynamic properties of the ice pack depend on how much
ice lies in each thickness range. Thus the basic problem in sea ice
modeling is to describe the evolution of the ice thickness distribution
(ITD) in time and space.

The fundamental equation solved by CICE is :cite:`Thorndike75`:

.. math::
   \frac{\partial g}{\partial t} = -\nabla \cdot (g {\bf u}) 
    - \frac{\partial}{\partial h} (f g) + \psi - L,
   :label: transport-g

where :math:`{\bf u}` is the horizontal ice velocity,
:math:`\nabla = (\frac{\partial}{\partial x}, \frac{\partial}{\partial y})`,
:math:`f` is the rate of thermodynamic ice growth, :math:`\psi` is a
ridging redistribution function, 
:math:`L` is the lateral melt rate
and :math:`g` is the ice thickness
distribution function. We define :math:`g({\bf x},h,t)\,dh` as the
fractional area covered by ice in the thickness range :math:`(h,h+dh)`
at a given time and location.  Icepack represents all of the terms in this
equation except for the divergence (the first term on the right).

Equation :eq:`transport-g` is solved by partitioning the ice pack in
each grid cell into discrete thickness categories. The number of
categories can be set by the user, with a default value :math:`N_C = 5`.
(Five categories, plus open water, are generally sufficient to simulate
the annual cycles of ice thickness, ice strength, and surface fluxes
:cite:`Bitz01,Lipscomb01`.) Each category :math:`n` has
lower thickness bound :math:`H_{n-1}` and upper bound :math:`H_n`. The
lower bound of the thinnest ice category, :math:`H_0`, is set to zero.
The other boundaries are chosen with greater resolution for small
:math:`h`, since the properties of the ice pack are especially sensitive
to the amount of thin ice :cite:`Maykut82`. The continuous
function :math:`g(h)` is replaced by the discrete variable
:math:`a_{in}`, defined as the fractional area covered by ice in the
open water by :math:`a_{i0}`, giving :math:`\sum_{n=0}^{N_C} a_{in} = 1`
by definition.

Category boundaries are computed in *init\_itd* using one of several
formulas, summarized in Table :ref:`tab-itd`. 
Setting the namelist variable ``kcatbound`` equal to 0 or 1 gives lower 
thickness boundaries for any number of thickness categories :math:`N_C`.
Table :ref:`tab-itd` shows the boundary values for :math:`N_C` = 5 and linear remapping 
of the ice thickness distribution. A third option specifies the boundaries 
based on the World Meteorological Organization classification; the full WMO
thickness distribution is used if :math:`N_C` = 7; if :math:`N_C` = 5 or
6, some of the thinner categories are combined. The original formula
(``kcatbound`` = 0) is the default. Category boundaries differ from those
shown in Table :ref:`tab-itd` for the delta-function ITD. Users may
substitute their own preferred boundaries in *init\_itd*.

Table :ref:`tab-itd` shows lower boundary values for thickness categories, in meters, for 
the three distribution options (*``kcatbound``*) and linear remapping (*``kitd``* = 1). 
In the WMO case, the distribution used depends on the number of categories used.

.. _tab-itd:

.. table:: *Lower boundary values* 

   +----------------+------------+---------+--------+--------+--------+
   | distribution   | original   | round   |           WMO            |
   +================+============+=========+========+========+========+
   | ``kcatbound``  | 0          | 1       |            2             |
   +----------------+------------+---------+--------+--------+--------+
   | :math:`N_C`    | 5          | 5       | 5      | 6      | 7      |
   +----------------+------------+---------+--------+--------+--------+
   | categories     |             lower bound (m)                     |
   +----------------+------------+---------+--------+--------+--------+
   | 1              | 0.00       | 0.00    | 0.00   | 0.00   | 0.00   |
   +----------------+------------+---------+--------+--------+--------+
   | 2              | 0.64       | 0.60    | 0.30   | 0.15   | 0.10   |
   +----------------+------------+---------+--------+--------+--------+
   | 3              | 1.39       | 1.40    | 0.70   | 0.30   | 0.15   |
   +----------------+------------+---------+--------+--------+--------+
   | 4              | 2.47       | 2.40    | 1.20   | 0.70   | 0.30   |
   +----------------+------------+---------+--------+--------+--------+
   | 5              | 4.57       | 3.60    | 2.00   | 1.20   | 0.70   |
   +----------------+------------+---------+--------+--------+--------+
   | 6              |            |         |        | 2.00   | 1.20   |
   +----------------+------------+---------+--------+--------+--------+
   | 7              |            |         |        |        | 2.00   |
   +----------------+------------+---------+--------+--------+--------+
   
Joint floe size and thickness distribution
=============================================

Sizes of individual sea ice floes vary over an extremely broad range, from centimeters
to hundreds of kilometers. The floe size distribution (FSD) is a probability function that
characterizes this variability :cite:`Rothrock84`. An option to include a 
prognostic sea ice floe size distribution is available and used if ``tr_fsd`` is set to true. 
The scheme is based on the theoretical framework described in :cite:`Horvat15` for a
*joint* floe size and thickness distribution (FSTD), and was implemented by :cite:`Roach18`.

In this theory, individual floes are identified with a size :math:`r` and area :math:`x(r)`, where
:math:`x(r)=4\alpha r^2` for :math:`\alpha=0.66 < \pi/4` (:cite:`Rothrock84`). The probability 
distribution :math:`f(r,h) dr dh` is the fraction of grid surface area 
covered by ice with thickness between :math:`h` and :math:`h + dh` and lateral floe
size between :math:`r` and :math:`r + dr`. The FSTD integrates over all floe sizes and
ice thicknesses to unity; over all floe sizes to the ITD; and over all thicknesses to the FSD.

For implementation in CICE,  the continuous function :math:`f(r,h)` is replaced
with a product of two discrete variables: :math:`a_{in}` as defined above and :math:`F_{in,k}`. 
:math:`F_{in,k}` is the fraction of ice belonging to thickness category :math:`n` with lateral 
floe size belonging to floe size class :math:`k`, giving
:math:`\sum_{n=0}^{N_C}\sum_{k=0}^{N_f} a_{in} F_{in,k} = 1` and :math:`\sum_{k=0}^{N_f}  F_{in,k} = 1`.
:math:`F_{in,k}` is carried as an area-weighted tracer.

The FSD may be ignored when considering processes that only modify ice thickness
(eg. vertical thermodynamics), and the ITD can be ignored when considering processes that only modify floe sizes (eg. wave fracture). For processes that affect both the ITD and the FSD, (eg. lateral melt), 
both :math:`a_{in}` and :math:`F_{in,k}` are evolved.

The FSTD evolves subject to lateral growth, lateral melt, new ice growth, floe welding and 
wave fracture, as described in :cite:`Roach18` and with some modifications described in 
:cite:`Roach19`. The equation for time evolution of the FSTD is (:cite:`Horvat15`),

:math:`\frac{\partial f(r,h)}{\partial t} = - \nabla \cdot (f(r,h)\mathbf{v}) + \mathcal{L}_T + \mathcal{L}_M + \mathcal{L}_W`,

where the terms on the right hand side represent the effects of advection, thermodynamics, mechanical 
redistribution and wave fracture respectively. Floe sizes do not explicitly appear in the equations of sea ice motion and therefore the FSTD is advected as an area tracer. We also assume that mechanical redistribution of sea ice through ridging does not impact floe sizes. Thus it remains only to compute the thermodynamic and wave fracture tendencies.

Thermodynamic changes to the FSTD are given by 

:math:`\mathcal{L}_T(r,h)=-\nabla_{(r,h)} \cdot (f(r,h) \mathbf{G}) +\frac{2}{r}f(r,h)G_r  + 
\delta(r-r_{\text{min}})\delta(h-h_{\text{min}})\dot{A}_p + \beta_{\text{weld}}.`

The first two terms on the right-hand side represent growth and melt of existing floes 
in thickness and lateral size, at a rate :math:`\mathbf{G} = (G_r,G_h)`. The third 
term represents growth of new ice: new floes are created at a rate :math:`\dot{A}_p` 
in the smallest thickness category and a given lateral size category. If wave forcing 
is provided, the size of newly formed floes is determined via a tensile stress limitation 
arising from the wave field (:cite:`Shen01`,:cite:`Roach19`); otherwise, all floes 
are presumed to grow as pancakes in the smallest floe size category resolved. 
To allow for the joining of individual floes to one another, we represent
the welding together of floes in freezing conditions via the fourth term, 
:math:`\beta_{\text{weld}}`, using a coagulation equation.

To compute the impact of wave fracture of the FSD, given a local ocean surface wave 
spectrum is provided, we generate a realization of the sea surface height field, which 
is uniquely determined by the spectrum up to a phase. In :cite:`Horvat15` this phase is 
randomly chosen, and multiple realizations of the resulting surface height field are used to 
obtain convergent statistics. However this stochastic component would lead to a model that is 
not bit-for-bit reproducible. Users can choose in the namelist (via ``wave_spec_type``)
to run the model with the phase set to be constant to obtain bit-for-bit reproducibility, or
to include the random phase, or to exclude wave effects completely.

We calculate the number and length of fractures that would occur if waves enter a fully ice-covered 
region defined in one dimension in the direction of propagation, and then apply
the outcome proportionally to the ice-covered fraction in each grid cell. 
Assuming that sea ice flexes with the sea surface height field, strains are computed
on this sub-grid-scale 1D domain. If the strain between successive extrema exceeds
a critical value new floes are formed with diameters equal to the distance between 
the extrema.

Floe size categories are set in *init\_fsd\_bounds* using an exponential spacing, beginning at 0.5 m with the
largest size resolved set by choice of :math:`N_f` (``nfsd``), the number of floe size categories.  Icepack
currently supports ``nfsd = 1, 12, 16, 24``.  Although ``nfsd = 1`` tracks the same ice floe diameter as
is assumed when ``tr_fsd=false``, the processes acting on the floes differ.
It is assumed that the floe size lies at the midpoint of each floe size category.

If simulations begin without ice (``ice_init='none'``), the FSD can emerge without initialization. This
is the recommended initialization for studies on the FSD itself. If simulations begin with ice cover, 
some initial FSD must be prescribed in ``init_fsd``. The default (used for ``ice_init='default'``) 
is a simple relationship determined from point observations by :cite:`Perovich14`, but its basin-wide 
applicability has not been tested. In Icepack, ``ice_init='default'`` is selected for the slab
and the full ITD cells.




