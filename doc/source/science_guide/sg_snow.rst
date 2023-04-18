:tocdepth: 3

.. _snow:

Advanced snow physics
=====================

Once deposited, the character and distribution of snow on sea ice depend on re-transport (wind), melting/wetting, and metamorphism (chiefly producing low-conductivity depth hoar or snow-ice). Each of these processes affects the other, and they are crucial for the evolution of the sea ice pack :cite:`Sturm02`. In particular, Wind slab and depth hoar resist densification, and Wind slab may prevent snow from drifting after deposition. Snow drifts around ridges cover only 6% of the ice surface area and are about 30% deeper than other snow-covered areas, but they prevent seawater filled cracks around the ridges from freezing, with important biological consequences.

The standard model configuration includes a basic snow formulation describing the essential effects of snow on sea ice, such as its albedo, vertical conduction, and growth/melt processes. It also incorporates more detailed processes such as snow-ice formation due to flooding and snow infiltration by melt water, which may form melt ponds. Several potentially important processes are not included in the standard configuration, such as compaction and redistribution of snow by wind and their effects on the thermal balance and on effective roughness. Snow metamorphism due to temperature gradients and liquid water content also are not included.

Setting ``tr_snow`` = ``.true.`` activates advanced snow physics parameterizations that represent the following processes, each of which has its own namelist flag for flexible configuration:

1. Radiative effects of snow redistribution by wind with respect to ice topography, including snow loss to leads and snow compaction by wind

2. Coupling effects associated with snow saturation and pond formation.

3. Radiative effects of snow grain metamorphism (variable grain size)

Snow can be scoured from level ice, blowing into leads or piling up on ridges.
The presence of liquid water in snow, such as rain or melt water, changes the surface albedo dramatically. It also alters the conductivity of the snow pack. These effects are associated mainly with the formation of depth hoar (change in grain size).

The standard model configuration assumes that the snow depth is uniform across each ice thickness category within a grid cell for the vertical thermodynamic calculation. However, there are separate radiation calculations for bare ice, snow-covered ice, and pond-covered ice; snow and ponds interact through snow saturation levels. Redistributing the snow alters these radiative calculations.

.. _snow_redistribution:

Snow redistribution
-------------------

Because the thermodynamic schemes in CICE assume a uniform snow depth over each category, ignoring the fractions of level and deformed ice, effects of snow redistribution are included only via the delta-Eddington radiation scheme. The redistributed snow depth is used to determine the effective area of bare ice (for very small snow depths) and the effective area and depth of melt ponds over level ice. Once those areas are determined, the redistributed snow volume over them is known, from which the snow depth for the remaining snow-covered area can be computed and used for its radiation balance calculation.

Two basic approaches are available for snow redistribution by wind, ``snwredist`` = ``bulk``, for which a user-defined parameter :math:`p` (``snwlvlfac``) determines the ratio of snow on ridges to that on level ice, and ``snwITDrdg``, in which snow can be compacted by the wind or eroded and redeposited on other thickness categories. For both, nonlocal redistribution of snow (i.e., between grid cells) is neglected, assuming that the difference between snow mass blowing into a grid cell and that blowing out is negligible, but snow can be blown into nearby leads and open water.



.. _snow_bulk:

Bulk snow redistribution
~~~~~~~~~~~~~~~~~~~~~~~~

:cite:`Sturm02` noted that on average during the SHEBA experiment, snow near ridged ice was 30% deeper than snow on undeformed ice. Using this rule of thumb, we can reduce the amount of snow on level ice in the model by reducing the snowfall rate over the sea ice and assuming the removed snow volume passes into the ocean through leads, instantaneously. This approach takes into account the area of open water available, as in the original code, by employing a precipitation flux in units of kg m :math:`^{-2}` s :math:`^{-1}`, which accumulates snow only on the ice-covered area of the grid cell.
      
This approach affects the simulation in two ways: (1) the snow removed from the level ice area is deposited into leads, and (2) using the snow remaining on the level ice area to adjust the effective melt pond and bare ice areas. Case (1) affects both the radiative and thermodynamic calculations by reducing the total amount of snow on the ice. Case (2) affects the radiative calculation directly, by possibly exposing more bare ice or melt ponds, but it affects the thermodynamic (conduction) calculation only through the altered radiative absorption, since the snow is always assumed to be equally deep over both level and deformed ice for the thermodynamic calculation.

When ``snwredist`` = ``bulk``, snow loss to leads is accomplished simply by reducing the volume of snowfall reaching the ice:

.. math::
   f_{s}^\prime = {f_s \left[ a_{lvl} \left({p\over{1+p}}\right)\right]},

where :math:`f_s` is the snowfall rate, :math:`a_{lvl}` is the average level-ice tracer value, and primed quantities represent their modified values.

Snow is redistributed between level and ridged ice within a single thickness category by solving a pair of equations for the modified level- and ridged-ice snow depths in terms of the original snow depth:

.. math::
   h_{lvl}^\prime = {1\over {1+p(1-a_{lvl})}} h_{lvl}

.. math::
   h_{rdg}^\prime = {{1 + p}\over {1+p(1-a_{lvl})}} h_{lvl}.

In the shortwave module for level-ice ponds, we create a new variable :math:`h_{lvl}^\prime` (``hsnlvl``) for snow depth over the level ice, and replace ``hsn`` with ``hsnlvl`` for the snow infiltration calculation and for the calculation of snow depth over refrozen melt ponds.

.. _snow_windredist:

Snow redistribution and compaction by wind
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Following :cite:`Lecomte15`, when ``snwredist`` = ``snwITDrdg`` we parameterize the amount of snow lost into the ocean through leads or redistributed to other thickness categories by defining the redistribution function :math:`\Phi` for snow mass as the sum of an erosion rate :math:`\Phi_E` and a redeposition rate :math:`\Phi_R` for each category of thickness :math:`h_i`:

.. math::
   \Phi_E = \left({\partial m \over \partial t}\right)_{erosion} = -{\gamma \over \sigma_{ITD}} \left(V-V^*\right){\rho_{max} - \rho_s \over \rho_{max}}

where :math:`\rho_s` and :math:`\rho_{max}` are the effective snow density and the maximum snow density in the model, respectively. For now, we take :math:`\rho_s` to be the wind-compacted snow density computed at the end of the snow model time step.

:math:`\Phi_E \Delta t` represents the maximum snow mass per unit area that may be suspended from each category, subject to the total mass (per unit area) available on each category.

Erosion begins when the instantaneous wind speed :math:`V` exceeds the seasonal wind speed required to compact the snow to a density :math:`\rho_s`, :math:`V^* = (\rho_s - \beta)/\alpha`. :math:`\sigma_{ITD}` is the standard deviation of the ice thicknesses from the thickness distribution :math:`g` within the grid cell. :math:`\gamma` is a tuning coefficient for
the eroded mass, which :cite:`Lecomte15` set to :math:`10^{-5}` kg m :math:`^{-2}`. From :cite:`Lecomte13`, :math:`\rho_s = 44.6V^* + 174` kg m :math:`^{-3}` for seasonal mean wind speed :math:`V`, i.e. :math:`\alpha=174` kg m :math:`^{-3}` and :math:`\beta=44.6` kg s m :math:`^{-4}`.

In :cite:`Lecomte15`, the fraction of this suspended snow lost in leads is

.. math::
   f = \left(1-a_i\right) \exp\left({-\sigma_{ITD}\over\sigma_{ref}}\right),

where the scale factor :math:`\sigma_{ref}=1` m and :math:`a_i` is the total ice area fraction within the grid cell.
Thus, the snow mass that is redistribution on the ice (i.e., not lost in leads) is

.. math::
   \Phi_R \Delta t = a_i \left(1-f\right) \Phi_E \Delta t.

We extend this approach by using the level and ridged ice thicknesses to compute the standard deviation of ice thickness across all categories.  That is,

.. math::
   \sigma_{ITD}^2 = \sum_{n=1}^N a_{in} a_{lvln} \left(h_{ilvln}-\sum_{k=1}^N a_{ik}h_{ik}\right)^2 + a_{in} a_{rdgn} \left(h_{irdgn} - \sum_{k=1}^N a_{ik} h_{ik} \right)^2.

When considering snow over ridged and level ice for the redistribution, we reapportion the fraction of snow on level ice as :math:`a_{slvl} = 1-(1+p)a_{rdg}` and note that with the average expression

.. math::
   a_{slvl} = {\sum_{n=1}^N a_{in}\left(a_{lvln} - p a_{rdgn}\right)  \over \sum_{n=1}^N a_{in}}

a conservative redistribution of snow across thickness categories is (for each category :math:`n`)

.. math::
   \Phi_R(n) \Delta t = a_i \left(1-f\right) \left[a_{rdgn}\left(1+p\right) + a_{slvl} \right] \Phi_E \Delta t,

where :math:`p \le a_{lvln}/a_{rdgn}`.

The snow volume and energy state variables are updated in two steps, first for erosion of snow into suspension, then snow redeposition. When redepositing the snow, the snow energy is distributed among the snow layers affected by erosion, proportionally to the fraction of snow eroded. Finally, snow layer thicknesses are re-equalized, conserving snow energy. The fraction of suspended snow mass and energy lost in leads is added to the fresh water and heat fluxes for strict conservation.

High wind speeds compact the upper portion of a snow pack into "wind slab," a dense and more conductive medium that resists further drifting. An effective snow density is computed based on wind speed, which is then used to limit snow erosion of denser snow.

:cite:`Lecomte15` note that once snow is deposited, its density changes very little. During deposition, the density primarily falls into one of two types, wind slab for wind velocities greater than about 10 m/s, and loose snow for lighter winds. Their table 3 indicates densities for a variety of snow types. "Hard slab," deposited at :math:`V` = 13 m/s, has a density of :math:`\rho_s` = 403 kg m :math:`^{-3}` and "soft slab" is :math:`\rho_s` = 321 kg m :math:`^{-3}`, deposited at :math:`V` = 10 m/s. Linearly interpolating between these values, we have :math:`\rho_s = 27.3V + 47.7`.  The slope is an adjustable namelist parameter, ``drhosdwind``.
For simplicity, we assign a minimum snow density of :math:`\rho_s^{min}` = 100 kg m :math:`^{-3}` s (``rhosmin``)
and add to it the gradient associated with wind speed from :cite:`Lecomte15` for wind speeds greater than 10 m/s:  :math:`\rho_{s}^{new} = \rho_{s}^{min} + 27.3 \max \left(V-10, 0\right)`.  The minimum wind speed to compact snow ``windmin`` is adjustable, and the maximum snow density is also a namelist parameter, ``rhosmax``.
This density is merged with preexisting layer densities only if new snow falls. The thickness of the wind slab is the larger of the depth of newly fallen snow or the thickness of snow redeposited by the wind. Following :cite:`Sturm02`, density does not evolve further, other than by transport, unless additional snow falls at high enough wind speeds to compact the snow.
   
.. _snow_liquid:

Ice and liquid water mass in snow
---------------------------------

The advanced snow physics option calculates ice and liquid water mass and effective snow grain radius, enabling them to interact with the radiation calculation.  The mass of ice and liquid water in snow are implemented as tracers on snow volume layers and used for the snow grain metamorphism.
Together with snow volume, they also can be used to determine effective snow density as :math:`\rho_s^{eff} = (m_{ice}+m_{liq}) / h_s`. Note that :math:`m_{ice}+m_{liq}`  is the snow water equivalent (kg/m :math:`^2`).

Sources of :math:`m_{ice}` are snowfall, condensation, and freezing of liquid water within the snowpack; sinks are sublimation and melting. All of the sources and sinks of :math:`m_{ice}` are already computed in the code except for freezing of liquid water within the snow pack.

Sources of :math:`m_{liq}` are rain and snow melt; freezing of liquid water within the snowpack and runoff are sinks. Runoff and meltwater entering a snow layer (i.e., runoff from the layer above) are associated with vertical flow through the snow column. As in :cite:`Oleson10`, when the liquid water within a snow layer exceeds the layer's holding capacity, the excess water is added to the underlying layer, limited by the effective porosity of the layer. When ``use_smliq_pnd`` is true, the excess water is supplied to the melt pond parameterization, which puts a fraction of it into the pond volume and allows the rest to run off into the ocean.

The snow mass fractions of precipitation and old ice are saved for metamorphosing the snow grain radius.

Except for the topo melt pond scheme, melt water and heat in ponds (which may be hidden within a partially saturated snow pack) are "virtual" in the sense that they are provided to the ocean model component immediately upon melting, even though the effects of the liquid water continue to be tracked as if it were retained on the ice. Retaining that water and heat in the sea ice component alters the timing, location and magnitude of fresh water runoff events into the ocean. All melt pond schemes include the meltwater effects, regardless of whether the liquid water is virtual.  The advanced snow physics option allows the liquid water calculated by the snow metamorphism scheme to be used for melt pond calculations, replacing the snow melt and rainfall terms.

.. _snow_metamorphosis:

Metamorphosis of snow grains
----------------------------

When ``snwgrain`` = ``.true.``, dynamic, effective snow radius, a snow volume tracer, evolves analytically as a function of snow temperature, temperature gradient, and density for radiative calculations using the delta-Eddington radiation scheme.
Wet metamorphism changes both density (through volume change) and effective grain size; here we only consider changes in grain radius.
In the formation of depth hoar, dry snow kinetic metamorphism (TG metamorphism) also increases the snow grain radius.

The tracers :math:`m_{liq}` and :math:`m_{ice}` characterize the snow in each snow layer, for each ice category and horizontal grid cell. The model's meltpond volume covers a fraction of the grid cell and represents liquid in excess of :math:`m_{liq}`. The radiative effects of snow grain radius in the fraction of ice covered by pond volume are only calculated when the pond volume has not yet saturated the snow pack; otherwise, delta-Eddington transfer uses meltpond properties. Therefore, modelled changes in snow grain radii from metamorphism are designed specifically for the fraction without exposed (i.e. effective) melt ponds.

Following :cite:`Oleson10`, the new snow grain radius is computed as a weighted function of existing and new (freshly fallen, ``rsnw_fall``) snow grain radii, using parameters from a look-up table that depends on snow temperature, temperature gradient and (effective) density.  The maximum snow radius is a namelist option, ``rsnw_tmax``.



