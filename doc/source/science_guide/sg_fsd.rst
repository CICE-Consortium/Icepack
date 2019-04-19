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
