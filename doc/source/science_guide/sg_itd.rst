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
============================

Sizes of individual sea ice floes vary over an extremely broad range, from centimeters
to hundreds of kilometers. The floe size distribution (FSD) is a probability function that
characterizes this variability (Rothrock & Thorndike, 1984). An option to include a 
prognostic sea ice floe size distribution is available and used if ``tr_fsd`` is set to true. 
The scheme is based on the theoretical framework described in Hovart & Tziperman (2015) for a
*joint* floe size and thickness distribution (FSTD), and was implemented by Roach et al. (2018).

In this theory, individual floes are identified with a size :math:`r` and area :math:`x(r)`, where
:math:`x(r)=4\alpha r^2` for :math:`\alpha=0.66 < \pi/4` (Rothrock & Thorndike, 1984). The probability 
distribution :math:`f(r,h) dr dh` is the fraction of grid surface area 
covered by ice with thickness between :math:`h` and :math:`h + dh` and lateral floe
size between :math:`r` and :math:`r + dr`. The FSTD integrates over all floe sizes and ice thicknesses to unity;
over all floe sizes to the ITD; and over all thicknesses to the FSD.

For implementation in CICE,  the continuous function :math:`f(r,h)` is replaced
with a product of two discrete variables: :math:`a_{in}` as defined above and :math:`F_{in,k}`. 
:math:`F_{in,k}` is the fraction of ice belonging to thickness category :math:`n` with lateral 
floe size belonging to floe size class :math:`k`, giving
:math:`\sum_{n=0}^{N_C}\sum_{k=0}^{N_f} a_{in} F_{in,k} = 1` and :math:`\sum_{k=0}^{N_f}  F_{in,k} = 1`.

This implementation means that the FSD can be ignored for processes that only modify the ITD (eg. vertical thermodynamics),
and the ITD can be ignored for processes that only modify the FSD (eg. wave fracture). For processes that affect both 
the ITD and the FSD, (eg. lateral melt), both :math:`a_{in}` and :math:`F_{in,k}` are evolved.

Floe size categories are set in *init\_fsd\_bounds* using an exponential spacing, beginning at 0.5 m with the
largest size resolved set by choice of :math:`N_f`, the number of floe size categories. It is assumed that 
the floe size lies at the midpoint of each floe size category.

When ``tr_fsd`` is set to true, the floe size distribution evolves subject to lateral growth, 
lateral melt, new ice growth, floe welding and wave fracture. These processes are described in XX. 
:math:`F_{in,k}` is carried as an area-weighted tracer. Floe sizes do not appear directly in any terms
in the momentum equation or constitutive law, and mechanical redistribution reduces the area fractions
of all floes equally.

If simulations begin without ice (``ice_ice='none'``), the FSD can emerge without initialization. This
is the recommended initialization for studies on the FSD itself. If simulations begin with ice cover, 
some initial FSD must be prescribed in ``init_fsd``. The default is a simple relationship determined 
from point observations by Perovich & Jones (2014), but its basin-wide applicability has not been tested.
