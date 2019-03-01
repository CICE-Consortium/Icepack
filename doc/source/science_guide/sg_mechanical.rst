:tocdepth: 3

.. _mech-red:

Mechanical redistribution
=========================

The last term on the right-hand side of Equation :eq:`transport-g`
is :math:`\psi`, which describes the redistribution
of ice in thickness space due to ridging and other mechanical processes.
The mechanical redistribution scheme in Icepack is based on
:cite:`Thorndike75`, :cite:`Rothrock75`,
:cite:`Hibler80`, :cite:`Flato95`, and
:cite:`Lipscomb07`. This scheme converts thinner ice to thicker
ice and is applied after horizontal transport. When the ice is
converging, enough ice ridges to ensure that the ice area does not
exceed the grid cell area.

First we specify the participation function: the thickness distribution
:math:`a_P(h) = b(h) \, g(h)` of the ice participating in ridging. (We
use "ridging" as shorthand for all forms of mechanical redistribution,
including rafting.) The weighting function :math:`b(h)` favors ridging
of thin ice and closing of open water in preference to ridging of
thicker ice. There are two options for the form of :math:`b(h)`. If
``krdg_partic`` = 0 in the namelist, we follow :cite:`Thorndike75`
and set

.. math::
   b(h) = \left\{\begin{array}{ll}  
          \frac{2}{G^*}(1-\frac{G(h)}{G^*}) & \mbox{if $G(h)<G^*$} \\
                    0                       & \mbox{otherwise}   
                 \end{array}  \right.
   :label: partic-old-contin

where :math:`G(h)` is the fractional area covered by ice thinner than
:math:`h`, and :math:`G^*` is an empirical constant. Integrating
:math:`a_P(h)` between category boundaries :math:`H_{n-1}` and
:math:`H_n`, we obtain the mean value of :math:`a_P` in category
:math:`n`:

.. math::
   a_{Pn} = \frac{2}{G^*} (G_n - G_{n-1})
            \left( 1 - \frac{G_{n-1}+G_n}{2 G^*} \right),
   :label: partic-old-discrete

where :math:`a_{Pn}` is the ratio of the ice area ridging (or open
water area closing) in category :math:`n` to the total area ridging and
closing, and :math:`G_n` is the total fractional ice area in categories
0 to :math:`n`. Equation :eq:`partic-old-discrete` applies to
categories with :math:`G_n < G^*`. If :math:`G_{n-1} < G^* < G_n`, then
Equation :eq:`partic-old-discrete` is valid with :math:`G^*` replacing
:math:`G_n`, and if :math:`G_{n-1} > G^*`, then :math:`a_{Pn} = 0`. If
the open water fraction :math:`a_0 > G^*`, no ice can ridge, because
"ridging" simply reduces the area of open water. As in
:cite:`Thorndike75` we set :math:`G^* = 0.15`.

If the spatial resolution is too fine for a given time step
:math:`\Delta t`, the weighting function Equation :eq:`partic-old-contin` can
promote numerical instability. For :math:`\Delta t = \mbox{1 hour}`,
resolutions finer than :math:`\Delta x \sim \mbox{10 km}` are typically
unstable. The instability results from feedback between the ridging
scheme and the dynamics via the ice strength. If the strength changes
significantly on time scales less than :math:`\Delta t`, the
viscous-plastic solution of the momentum equation is inaccurate and
sometimes oscillatory. As a result, the fields of ice area, thickness,
velocity, strength, divergence, and shear can become noisy and
unphysical.

A more stable weighting function was suggested by
:cite:`Lipscomb07`:

.. math::
   b(h) = \frac{\exp[-G(h)/a^*]}
               {a^*[1-\exp(-1/a^*)]}
   :label: partic-new-contin

When integrated between category boundaries, Equation :eq:`partic-new-contin`
implies

.. math::
   a_{Pn} = \frac {\exp(-G_{n-1}/a^*) - \exp(-G_{n}/a^*)}
                  {1 - \exp(-1/a^*)}
   :label: partic-new-discrete

This weighting function is used if ``krdg_partic`` = 1 in the namelist.
From Equation :eq:`partic-new-contin`, the mean value of :math:`G` for ice
participating in ridging is :math:`a^*`, as compared to :math:`G^*/3`
for Equation :eq:`partic-old-contin`. For typical ice thickness distributions,
setting :math:`a^* = 0.05` with ``krdg_partic`` = 1 gives participation
fractions similar to those given by :math:`G^* = 0.15` with ``krdg_partic``
= 0. See :cite:`Lipscomb07` for a detailed comparison of these
two participation functions.

Thin ice is converted to thick, ridged ice in a way that reduces the
total ice area while conserving ice volume and internal energy. There
are two namelist options for redistributing ice among thickness
categories. If ``krdg_redist`` = 0, ridging ice of thickness :math:`h_n`
forms ridges whose area is distributed uniformly between
:math:`H_{\min} = 2 h_n` and :math:`H_{\max} = 2 \sqrt{H^* h_n}`, as in
:cite:`Hibler80`. The default value of :math:`H^*` is 25 m, as
in earlier versions of CICE. Observations suggest that
:math:`H^* = 50` m gives a better fit to first-year ridges
:cite:`Amundrud04`, although the lower value may be appropriate
for multiyear ridges :cite:`Flato95`. The ratio of the mean
ridge thickness to the thickness of ridging ice is
:math:`k_n = (H_{\min} + H_{\max}) / (2 h_n)`. If the area of category
:math:`n` is reduced by ridging at the rate :math:`r_n`, the area of
thicker categories grows simultaneously at the rate :math:`r_n/k_n`.
Thus the *net* rate of area loss due to ridging of ice in category
:math:`n` is :math:`r_n(1-1/k_n)`.

The ridged ice area and volume are apportioned among categories in the
thickness range :math:`(H_{\min}, H_{\max})`. The fraction of the new
ridge area in category :math:`m` is

.. math::
   f_m^{\mathrm{area}} = \frac{H_R - H_L} 
                              {H_{\max} - H_{\min}},
   :label: ridge-area-old

where :math:`H_L = \max(H_{m-1},H_{\min})` and
:math:`H_R= \min(H_m,H_{\max})`. The fraction of the ridge volume going
to category :math:`m` is

.. math::
   f_m^{\mathrm{vol}} = \frac{(H_R)^2 - (H_L)^2}
                             {(H_{\max})^2 - (H_{\min})^2}.
   :label: ridge-volume-old

This uniform redistribution function tends to produce too little ice in
the 3–5 m range and too much ice thicker than 10 m
:cite:`Amundrud04`. Observations show that the ITD of ridges is
better approximated by a negative exponential. Setting ``krdg_redist`` = 1
gives ridges with an exponential ITD :cite:`Lipscomb07`:

.. math::
   g_R(h) \propto \exp[-(h - H_{\min})/\lambda]
   :label: redist-new

for :math:`h \ge H_{\min}`, with :math:`g_R(h) = 0` for
:math:`h < H_{\min}`. Here, :math:`\lambda` is an empirical *e*-folding
scale and :math:`H_{\min}=2h_n` (where :math:`h_n` is the thickness of
ridging ice). We assume that :math:`\lambda = \mu h_n^{1/2}`, where
:math:`\mu` (mu\_rdg) is a tunable parameter with units . Thus the mean
ridge thickness increases in proportion to :math:`h_n^{1/2}`, as in
:cite:`Hibler80`. The value :math:`\mu = 4.0`  gives
:math:`\lambda` in the range 1–4 m for most ridged ice. Ice strengths
with :math:`\mu = 4.0`  and ``krdg_redist`` = 1 are roughly comparable to
the strengths with :math:`H^* = 50` m and ``krdg_redist`` = 0.

From Equation :eq:`redist-new` it can be shown that the fractional area going
to category :math:`m` as a result of ridging is

.. math::
   f_m^{\mathrm{area}} = \exp[-(H_{m-1} - H_{\min}) / \lambda] 
                        - \exp[-(H_m - H_{\min}) / \lambda].
   :label: ridge-area-new

The fractional volume going to category :math:`m` is

.. math::
   f_m^{\mathrm{vol}} = \frac{(H_{m-1}+\lambda) \exp[-(H_{m-1}-H_{\min})/\lambda]
                              - (H_m + \lambda) \exp[-(H_m - H_{\min}) / \lambda]}
                                {H_{min} + \lambda}.
   :label: ridge-volume-new

Equations :eq:`ridge-area-new` and :eq:`ridge-volume-new` replace
Equations :eq:`ridge-area-old` and :eq:`ridge-volume-old` when ``krdg_redist``
= 1.

Internal ice energy is transferred between categories in proportion to
ice volume. Snow volume and internal energy are transferred in the same
way, except that a fraction of the snow may be deposited in the ocean
instead of added to the new ridge.

The net area removed by ridging and closing is a function of the strain
rates. Let :math:`R_{\mathrm{net}}` be the net rate of area loss for the
ice pack (i.e., the rate of open water area closing, plus the net rate
of ice area loss due to ridging). Following :cite:`Flato95`,
:math:`R_{\mathrm{net}}` is given by

.. math::
   R_{\mathrm{net}} = \frac{C_s}{2}
                    (\Delta - |D_D|) - \min(D_D,0),
   :label: Rnet

where :math:`C_s` is the fraction of shear dissipation energy that
contributes to ridge-building, :math:`D_D` is the divergence, and
:math:`\Delta` is a function of the divergence and shear. These strain
rates are computed by the dynamics scheme. The default value of
:math:`C_s` is 0.25.

Next, define :math:`R_{\mathrm{tot}} = \sum_{n=0}^N r_n`. This rate is
related to :math:`R_{\mathrm{net}}` by

.. math::
   R_{\mathrm{net}} =
      \left[ a_{P0} + \sum_{n=1}^N a_{Pn}\left(1-{1\over k_n}\right)\right]
       R_{\mathrm{tot}}.
   :label: Rtot-Rnet

Given :math:`R_{\mathrm{net}}` from Equation :eq:`Rnet`, we
use Equation :eq:`Rtot-Rnet` to compute :math:`R_{\mathrm{tot}}`. Then the area
ridged in category :math:`n` is given by :math:`a_{rn} = r_n \Delta t`,
where :math:`r_n = a_{Pn} R_{\mathrm{tot}}`. The area of new ridges is
:math:`a_{rn} / k_n`, and the volume of new ridges is :math:`a_{rn} h_n`
(since volume is conserved during ridging). We remove the ridging ice
from category :math:`n` and use Equations :eq:`ridge-area-old`
and :eq:`ridge-volume-old` (or :eq:`ridge-area-new` and
:eq:`ridge-volume-new`) to redistribute the ice among thicker
categories.

Occasionally the ridging rate in thickness category :math:`n` may be
large enough to ridge the entire area :math:`a_n` during a time interval
less than :math:`\Delta t`. In this case :math:`R_{\mathrm{tot}}` is
reduced to the value that exactly ridges an area :math:`a_n` during
:math:`\Delta t`. After each ridging iteration, the total fractional ice
area :math:`a_i` is computed. If :math:`a_i > 1`, the ridging is
repeated with a value of :math:`R_{\mathrm{net}}` sufficient to yield
:math:`a_i = 1`.

Two tracers for tracking the ridged ice area and volume are available.
The actual tracers are for level (undeformed) ice area (`alvl`) and volume
(`vlvl`), which are easier to implement for a couple of reasons: (1) ice
ridged in a given thickness category is spread out among the rest of the
categories, making it more difficult (and expensive) to track than the
level ice remaining behind in the original category; (2) previously
ridged ice may ridge again, so that simply adding a volume of freshly
ridged ice to the volume of previously ridged ice in a grid cell may be
inappropriate. Although the code currently only tracks level ice
internally, both level ice and ridged ice are available for output.
They are simply related:

.. math::
   \begin{aligned}
   a_{lvl} + a_{rdg} &=& a_i, \\
   v_{lvl} + v_{rdg} &=& v_i.\end{aligned}
   :label: alvl

Level ice area fraction and volume increase with new ice formation and
decrease steadily via ridging processes. Without the formation of new
ice, level ice asymptotes to zero because we assume that both level ice
and ridged ice ridge, in proportion to their fractional areas in a grid
cell (in the spirit of the ridging calculation itself which does not
prefer level ice over previously ridged ice).

The ice strength :math:`P` may be computed in either of two ways. If the
namelist parameter ``kstrength`` = 0, we use the strength formula from
:cite:`Hibler79`:

.. math::
   P = P^* h \exp[-C(1-a_i)],
   :label: hib-strength

where :math:`P^* = 27,500 \, \mathrm {N/m}` and :math:`C = 20` are
empirical constants, and :math:`h` is the mean ice thickness.
Alternatively, setting ``kstrength`` = 1 gives an ice strength closely
related to the ridging scheme. Following
:cite:`Rothrock75`, the strength is assumed proportional
to the change in ice potential energy :math:`\Delta E_P` per unit area
of compressive deformation. Given uniform ridge ITDs (``krdg_redist`` = 0),
we have

.. math::
   P = C_f \, C_p \, \beta \sum_{n=1}^{N_C}
     \left[ -a_{Pn} \, h_n^2  + \frac{a_{Pn}}{k_n}
        \left( \frac{(H_n^{\max})^3 - (H_n^{\min})^3}
                    {3(H_n^{\max}-H_n^{\min})} \right) \right],
   :label: roth-strength0

where :math:`C_P = (g/2)(\rho_i/\rho_w)(\rho_w-\rho_i)`,
:math:`\beta =R_{\mathrm{tot}}/R_{\mathrm{net}} > 1`
from Equation :eq:`Rtot-Rnet`, and :math:`C_f` is an empirical parameter that
accounts for frictional energy dissipation. Following
:cite:`Flato95`, we set :math:`C_f = 17`. The first term in
the summation is the potential energy of ridging ice, and the second,
larger term is the potential energy of the resulting ridges. The factor
of :math:`\beta` is included because :math:`a_{Pn}` is normalized with
respect to the total area of ice ridging, not the net area removed.
Recall that more than one unit area of ice must be ridged to reduce the
net ice area by one unit. For exponential ridge ITDs (``krdg_redist`` = 1),
the ridge potential energy is modified:

.. math::
   P = C_f \, C_p \, \beta \sum_{n=1}^{N_C}
     \left[ -a_{Pn} \, h_n^2  + \frac{a_{Pn}}{k_n}
        \left( H_{\min}^2 + 2H_{\min}\lambda + 2 \lambda^2 \right) \right]
   :label: roth-strength1

The energy-based ice strength given by Equations :eq:`roth-strength0` or
:eq:`roth-strength1` is more physically realistic than the strength
given by Equation :eq:`hib-strength`. However, use of Equation :eq:`hib-strength` is
less likely to allow numerical instability at a given resolution and
time step. See :cite:`Lipscomb07` for more details.
