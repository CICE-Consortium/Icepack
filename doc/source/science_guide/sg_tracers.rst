:tocdepth: 3

.. _tracers:

Tracers
=======

Numerous tracers are available with the column physics.  Several of these are 
required (surface temperature and thickness, salinity and enthalpy of ice and snow layers),
and many others are options.  For instance, there are tracers to track the age of the ice;
the area of first-year ice, fractions of ice area and volume that are level, from which
the amount of deformed ice can be calculated; pond area, pond volume and volume of ice covering ponds;
a prognostic floe size distribution; aerosols and numerous other biogeochemical tracers.
Most of these tracers are presented in later sections.  Here we describe the ice age 
tracers and how tracers may depend on other tracers, using the pond tracers as an 
example.

.. _ice-age:

Ice age
-------

The age of the ice, :math:`\tau_{age}`, is treated as an
ice-volume tracer (`trcr\_depend` = 1). It is initialized at 0 when ice
forms as frazil, and the ice ages the length of the timestep during each
timestep. Freezing directly onto the bottom of the ice does not affect
the age, nor does melting. Mechanical redistribution processes and
advection alter the age of ice in any given grid cell in a conservative
manner following changes in ice area. The sea ice age tracer is
validated in :cite:`Hunke09`.

Another age-related tracer, the area covered by first-year ice
:math:`a_{FY}`, is an area tracer (`trcr\_depend` = 0) that corresponds
more closely to satellite-derived ice age data for first-year ice than
does :math:`\tau_{age}`. It is re-initialized each year on 15
September (``yday`` = 259) in the northern hemisphere and 15 March (``yday`` =
75) in the southern hemisphere, in non-leap years. This tracer is
increased when new ice forms in open water, in subroutine
*add\_new\_ice* in **icepack\_therm\_itd.F90**. The first-year area tracer
is discussed in :cite:`Armour11`.

.. _pondtr:

Tracers that depend on other tracers 
------------------------------------

Tracers may be defined that depend on other tracers. Melt pond tracers
provide an example (these equations pertain to cesm and topo tracers;
level-ice tracers are similar with an extra factor of :math:`a_{lvl}`,
see Equations :eq:`transport-lvl`–:eq:`transport-ipnd-lvl`). Conservation
equations for pond area fraction :math:`a_{pnd}a_i` and pond volume
:math:`h_{pnd}a_{pnd}a_i`, given the ice velocity :math:`\bf u`, are

.. math::
   {\partial\over\partial t} (a_{pnd}a_{i}) + \nabla \cdot (a_{pnd}a_{i} {\bf u}) = 0,
   :label: transport-apnd

.. math::
   {\partial\over\partial t} (h_{pnd}a_{pnd}a_{i}) + \nabla \cdot (h_{pnd}a_{pnd}a_{i} {\bf u}) = 0.
   :label: transport-hpnd

These equations represent quantities within one thickness category;
all melt pond calculations are performed for each category, separately.
Equation :eq:`transport-hpnd` expresses conservation of melt pond
volume, but in this form highlights that the quantity tracked in the
code is the pond depth tracer :math:`h_{pnd}`, which depends on the pond
area tracer :math:`a_{pnd}`. Likewise, :math:`a_{pnd}` is a tracer on
ice area (Equation :eq:`transport-apnd`), which is a state variable, not a
tracer.

For a generic quantity :math:`q` that represents a mean value over the
ice fraction, :math:`q a_i` is the average value over the grid cell.
Thus for cesm or topo melt ponds, :math:`h_{pnd}` can be considered the
actual pond depth, :math:`h_{pnd}a_{pnd}` is the mean pond depth over
the sea ice, and :math:`h_{pnd}a_{pnd}a_i` is the mean pond depth over
the grid cell. These quantities are illustrated in Figure :ref:`fig-tracers`.
The graphic on the right illustrates the *grid cell* fraction of ponds or 
level ice as defined by the tracers. The chart on the left provides 
corresponding ice thickness and pond depth averages over the grid cell, 
sea ice and pond area fractions. 

.. _fig-tracers:

.. figure:: ./figures/tracergraphic.png
   :align: center
   :scale: 50%  

   *Melt pond tracer definitions*

Tracers may need to be modified for physical reasons outside of the
"core" module or subroutine describing their evolution. For example,
when new ice forms in open water, the new ice does not yet have ponds on
it. Likewise when sea ice deforms, we assume that pond water (and ice)
on the portion of ice that ridges is lost to the ocean.

When new ice is added to a grid cell, the *grid cell* total area of melt
ponds is preserved within each category gaining ice,
:math:`a_{pnd}^{t+\Delta t}a_{i}^{t+\Delta t} = a_{pnd}^{t}a_{i}^{t}`, 
or

.. math::
   a_{pnd}^{t+\Delta t}= {a_{pnd}^{t}a_{i}^{t} \over a_{i}^{t+\Delta t} }.
   :label: apnd

Similar calculations are performed for all tracer types, for example
tracer-on-tracer dependencies such as :math:`h_{pnd}`, when needed:

.. math:: 
   h_{pnd}^{t+\Delta t}= {h_{pnd}^{t}a_{pnd}^{t}a_{i}^{t} \over a_{pnd}^{t+\Delta t}a_{i}^{t+\Delta t} }.
   :label: hpnd

In this case (adding new ice), :math:`h_{pnd}` does not change because
:math:`a_{pnd}^{t+\Delta t}a_{i}^{t+\Delta t} = a_{pnd}^{t}a_{i}^{t}`.

When ice is transferred between two thickness categories, we conserve
the total pond area summed over categories :math:`n`,

.. math:: 
   \sum_n a_{pnd}^{t+\Delta t}(n)a_{i}^{t+\Delta t}(n) = \sum_n a_{pnd}^{t}(n)a_{i}^{t}(n).
   :label: apnd2

Thus,

.. math::
   a_{pnd}^{t+\Delta t}(m) &= {\sum_n a_{pnd}^{t}(n)a_{i}^{t}(n) - \sum_{n\ne m} a_{pnd}^{t+\Delta t}(n)a_{i}^{t+\Delta t}(n) \over a_i^{t+\Delta t}(m)  } \\
   &= {a_{pnd}^t(m)a_i^t(m) + \sum_{n\ne m} \Delta \left(a_{pnd}a_i\right)^{t+\Delta t} \over a_i^{t+\Delta t}(m)  }
   :label: apnd3

This is more complicated because of the :math:`\Delta` term on the
right-hand side, which is handled in subroutine *icepack\_compute\_tracers*. Such
tracer calculations are scattered throughout the code, wherever there
are changes to the ice thickness distribution.

Note that if a quantity such as :math:`a_{pnd}` becomes zero in a grid
cell’s thickness category, then all tracers that depend on it also
become zero. If a tracer should be conserved (e.g., aerosols and the
liquid water in topo ponds), additional code must be added to track
changes in the conserved quantity.

Tracer dependencies and conserved quantities associated with tracers are tracked
using the arrays ``trcr_depend``, which defines the type of dependency (area, volume, snow, etc),
``n_trcr_strata``, the number of underlying layers, ``nt_strata``, the indices of the underlying
layers, and ``trcr_base``, a mask that is one for the tracer dependency and zero otherwise.
These arrays are used to convert between the tracer values themselves and the conserved
forms.

More information about the melt pond schemes is in the
:ref:`ponds` section.
