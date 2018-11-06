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

