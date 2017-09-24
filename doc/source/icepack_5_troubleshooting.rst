.. role:: math(raw)
   :format: html latex
..

Troubleshooting 
================

Check the FAQ: https://github.com/CICE-Consortium/Icepack/wiki

.. _setup:

Initial setup
-------------

.. _restarttrouble:

Restarts
--------

Known bugs
----------

#. With the old CCSM radiative scheme (``shortwave`` = ‘default’ or
   ‘ccsm3’), a sizable fraction (more than 10%) of the total shortwave
   radiation is absorbed at the surface but should be penetrating into
   the ice interior instead. This is due to use of the aggregated,
   effective albedo rather than the bare ice albedo when
   ``snowpatch`` :math:`< 1`.

Interpretation of albedos
-------------------------

The snow-and-ice albedo, ``albsni``, and diagnostic albedos ``albice``, ``albsno``,
and ``albpnd`` are merged over categories but not scaled (divided) by the
total ice area. (This is a change from CICE v4.1 for ``albsni``.) The latter
three history variables represent completely bare or completely snow- or
melt-pond-covered ice; that is, they do not take into account the snow
or melt pond fraction (``albsni`` does, as does the code itself during
thermodyamic computations). This is to facilitate comparison with
typical values in measurements or other albedo parameterizations. The
melt pond albedo ``albpnd`` is only computed for the Delta-Eddington
shortwave case.

With the Delta-Eddington parameterization, the albedo depends on the
cosine of the zenith angle (:math:`\cos\varphi`, ``coszen``) and is zero if
the sun is below the horizon (:math:`\cos\varphi < 0`). Therefore
time-averaged albedo fields would be low if a diurnal solar cycle is
used, because zero values would be included in the average for half of
each 24-hour period. To rectify this, a separate counter is used for the
averaging that is incremented only when :math:`\cos\varphi > 0`. The
albedos will still be zero in the dark, polar winter hemisphere.

Proliferating subprocess parameterizations
------------------------------------------

With the addition of several alternative parameterizations for sea ice
processes, a number of subprocesses now appear in multiple parts of the
code with differing descriptions. For instance, sea ice porosity and
permeability, along with associated flushing and flooding, are
calculated separately for mushy thermodynamics, topo and level-ice melt
ponds, and for the brine height tracer, each employing its own
equations. Likewise, the BL99 and mushy thermodynamics compute freeboard
and snow–ice formation differently, and the topo and level-ice melt pond
schemes both allow fresh ice to grow atop melt ponds, using slightly
different formulations for Stefan freezing. These various process
parameterizations will be compared and their subprocess descriptions
possibly unified in the future.
