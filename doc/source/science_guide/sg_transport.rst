:tocdepth: 3

.. _itd-trans:

Transport in thickness space
============================

Next we solve the equation for ice transport in thickness space due to
thermodynamic growth and melt,

.. math::
   \frac{\partial g}{\partial t} + \frac{\partial}{\partial h} (f g) = 0,
   :label: itd-transport

which is obtained from Equation :eq:`transport-g` by neglecting the first and
third terms on the right-hand side. We use the remapping method of
:cite:`Lipscomb01`, in which thickness categories are
represented as Lagrangian grid cells whose boundaries are projected
forward in time. The thickness distribution function :math:`g` is
approximated as a linear function of :math:`h` in each displaced
category and is then remapped onto the original thickness categories.
This method is numerically smooth and is not too diffusive. It can be
viewed as a 1D simplification of the 2D incremental remapping scheme
described above.

We first compute the displacement of category boundaries in thickness
space. Assume that at time :math:`m` the ice areas :math:`a_n^m` and
mean ice thicknesses :math:`h_n^m` are known for each thickness
category. (For now we omit the subscript :math:`i` that distinguishes
ice from snow.) We use a thermodynamic model (:ref:`thermo`)
to compute the new mean thicknesses :math:`h_n^{m+1}` at time
:math:`m+1`. The time step must be small enough that trajectories do not
cross; i.e., :math:`h_n^{m+1} < h_{n+1}^{m+1}` for each pair of adjacent
categories. The growth rate at :math:`h = h_n` is given by
:math:`f_n = (h_n^{m+1} - h_n^m) / \Delta t`. By linear interpolation we
estimate the growth rate :math:`F_n` at the upper category boundary
:math:`H_n`:

.. math:: 
   F_n = f_n + \frac{f_{n+1}-f_n}{h_{n+1}-h_n} \, (H_n - h_n).
   :label: growth-rate

If :math:`a_n` or :math:`a_{n+1} = 0`, :math:`F_n` is set to the growth
rate in the nonzero category, and if :math:`a_n = a_{n+1} = 0`, we set
:math:`F_n = 0`. The temporary displaced boundaries are given by

.. math:: 
   H_n^* = H_n + F_n \, \Delta t, \ n = 1 \ {\rm to} \ N-1
   :label: displ

The boundaries must not be displaced by more than one category to the
left or right; that is, we require :math:`H_{n-1} < H_n^* < H_{n+1}`.
Without this requirement we would need to do a general remapping rather
than an incremental remapping, at the cost of added complexity.

Next we construct :math:`g(h)` in the displaced thickness categories.
The ice areas in the displaced categories are :math:`a_n^{m+1} = a_n^m`,
since area is conserved following the motion in thickness space (i.e.,
during vertical ice growth or melting). The new ice volumes are
:math:`v_n^{m+1} = (a_n h_n)^{m+1} = a_n^m h_n^{m+1}`. For conciseness,
define :math:`H_L = H_{n-1}^*` and :math:`H_R = H_{n}^*` and drop the
time index :math:`m+1`. We wish to construct a continuous function
:math:`g(h)` within each category such that the total area and volume at
time :math:`m+1` are :math:`a_n` and :math:`v_n`, respectively:

.. math::
   \int_{H_L}^{H_R} g \, dh = a_n,
   :label: area-cons

.. math::
   \int_{H_L}^{H_R} h \, g \, dh = v_n.
   :label: volume-cons

The simplest polynomial that can satisfy both equations is a line. It
is convenient to change coordinates, writing
:math:`g(\eta) = g_1 \eta + g_0`, where :math:`\eta = h - H_L` and the
coefficients :math:`g_0` and :math:`g_1` are to be determined. Then
Equations :eq:`area-cons` and :eq:`volume-cons` can be written as

.. math:: 
   g_1 \frac{\eta_R^2}{2} + g_0 \eta_R = a_n,
  :label: g1

.. math:: 
   g_1 \frac{\eta_R^3}{3} + g_0 \frac{\eta_R^2}{2} = a_n \eta_n,
   :label: g1a

where :math:`\eta_R = H_R - H_L` and :math:`\eta_n = h_n - H_L`. These
equations have the solution

.. math::
   g_0 = \frac{6 a_n}{\eta_R^2} \left(\frac{2 \eta_R}{3} - \eta_n\right),
   :label: g0

.. math::
   g_1 = \frac{12 a_n}{\eta_R^3} \left(\eta_n - \frac{\eta_R}{2}\right).
   :label: g1b

Since :math:`g` is linear, its maximum and minimum values lie at the
boundaries, :math:`\eta = 0` and :math:`\eta_R`:

.. math::
   g(0)=\frac{6 a_n}{\eta_R^2} \, \left(\frac{2 \eta_R}{3} - \eta_n\right) = g_0,
   :label: gmin
 
.. math::
   g(\eta_R) = \frac{6 a_n}{\eta_R^2} \, \left(\eta_n - \frac{\eta_R}{3}\right).
   :label: gmax

Equation :eq:`gmin` implies that :math:`g(0) < 0` when
:math:`\eta_n > 2 \eta_R/3`, i.e., when :math:`h_n` lies in the right
third of the thickness range :math:`(H_L, H_R)`. Similarly, Equation :eq:`gmax`
implies that :math:`g(\eta_R) < 0` when :math:`\eta_n < \eta_R/3`, i.e.,
when :math:`h_n` is in the left third of the range. Since negative
values of :math:`g` are unphysical, a different solution is needed when
:math:`h_n` lies outside the central third of the thickness range. If
:math:`h_n` is in the left third of the range, we define a cutoff
thickness, :math:`H_C = 3 h_n - 2 H_L`, and set :math:`g = 0` between
:math:`H_C` and :math:`H_R`. Equations :eq:`g0` and :eq:`g1` are then
valid with :math:`\eta_R` redefined as :math:`H_C - H_L`. And if
:math:`h_n` is in the right third of the range, we define
:math:`H_C = 3 h_n - 2 H_R` and set :math:`g = 0` between :math:`H_L`
and :math:`H_C`. In this case, :eq:`g0` and :eq:`g1` apply with
:math:`\eta_R = H_R - H_C` and :math:`\eta_n = h_n - H_C`.

Figure :ref:`fig-gplot` illustrates the linear reconstruction of :math:`g`
for the simple cases :math:`H_L = 0`, :math:`H_R = 1`, :math:`a_n = 1`,
and :math:`h_n =` 0.2, 0.4, 0.6, and 0.8. Note that :math:`g` slopes
downward (:math:`g_1 < 0`) when :math:`h_n` is less than the midpoint
thickness, :math:`(H_L + H_R)/2 = 1/2`, and upward when :math:`h_n`
exceeds the midpoint thickness. For :math:`h_n = 0.2` and 0.8,
:math:`g = 0` over part of the range.

.. _fig-gplot:

.. figure:: ./figures/gplot.png
   :align: center
   :scale: 20%

   *Linear approximation of thickness distribution function*

Finally, we remap the thickness distribution to the original boundaries
by transferring area and volume between categories. We compute the ice
area :math:`\Delta a_n` and volume :math:`\Delta v_n` between each
original boundary :math:`H_n` and displaced boundary :math:`H_n^*`. If
:math:`H_n^* > H_n`, ice moves from category :math:`n` to :math:`n+1`.
The area and volume transferred are

.. math::
   \Delta a_n = \int_{H_n}^{H_n^*} g \, dh,
   :label: move-area

.. math::
   \Delta v_n = \int_{H_n}^{H_n^*} h \, g \, dh.
   :label: move-volume

If :math:`H_n^* < H_N`, ice area and volume are transferred from
category :math:`n+1` to :math:`n` using Equations :eq:`move-area` and
:eq:`move-volume` with the limits of integration reversed. To evaluate
the integrals we change coordinates from :math:`h` to
:math:`\eta = h - H_L`, where :math:`H_L` is the left limit of the range
over which :math:`g > 0`, and write :math:`g(\eta)` using Equations :eq:`g0` and
:eq:`g1`. In this way we obtain the new areas :math:`a_n` and volumes
:math:`v_n` between the original boundaries :math:`H_{n-1}` and
:math:`H_n` in each category. The new thicknesses,
:math:`h_n = v_n/a_n`, are guaranteed to lie in the range
:math:`(H_{n-1}, H_n)`. If :math:`g = 0` in the part of a category that
is remapped to a neighboring category, no ice is transferred.

Other conserved quantities are transferred in proportion to the ice
volume :math:`\Delta v_{in}`. For example, the transferred ice energy in
layer :math:`k` is
:math:`\Delta e_{ink} = e_{ink} (\Delta v_{in} / v_{in})`.

The left and right boundaries of the domain require special treatment.
If ice is growing in open water at a rate :math:`F_0`, the left boundary
:math:`H_0` is shifted to the right by :math:`F_0 \Delta t` before
:math:`g` is constructed in category 1, then reset to zero after the
remapping is complete. New ice is then added to the grid cell,
conserving area, volume, and energy. If ice cannot grow in open water
(because the ocean is too warm or the net surface energy flux is
downward), :math:`H_0` is fixed at zero, and the growth rate at the left
boundary is estimated as :math:`F_0 = f_1`. If :math:`F_0 < 0`, all ice
thinner than :math:`\Delta h_0 = -F_0 \Delta t` is assumed to have
melted, and the ice area in category 1 is reduced accordingly. The area
of new open water is

.. math:: 
   \Delta a_0 = \int_{0}^{\Delta h_0} g \, dh.
   :label: a0

The right boundary :math:`H_N` is not fixed but varies with
:math:`h_N`, the mean ice thickness in the thickest category. Given
:math:`h_N`, we set :math:`H_N = 3 h_N - 2 H_{N-1}`, which ensures that
:math:`g(h) > 0` for :math:`H_{N-1} < h < H_N` and :math:`g(h) = 0` for
:math:`h \geq H_N`. No ice crosses the right boundary. If the ice growth
or melt rates in a given grid cell are too large, the thickness
remapping scheme will not work. Instead, the thickness categories in
that grid cell are treated as delta functions following
:cite:`Bitz01`, and categories outside their prescribed
boundaries are merged with neighboring categories as needed. For time
steps of less than a day and category thickness ranges of 10 cm or more,
this simplification is needed rarely, if ever.

The linear remapping algorithm for thickness is not monotonic for
tracers, although significant errors rarely occur. Usually they appear
as snow temperatures (enthalpy) outside the physical range of values in
very small snow volumes. In this case we transfer the snow and its heat
and tracer contents to the ocean.
