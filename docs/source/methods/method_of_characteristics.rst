.. _method_of_characteristics:

=========================
Method of Characteristics
=========================

The method of characteristics (MOC) is a widely used technique for solving partial differential equations, including the Boltzmann form of the neutron transport equation [Askew]_. MOC is used to solve the transport equation in 2D by discretizing both polar and azimuthal angles and integrating the multi-group characteristic form of the equation for a particular azimuthal and polar angle quadrature. The following sections detail the derivation of the characteristic form of the 2D neutron transport equation solved in the OpenMOC method of characteristics scheme [Boyd]_.

:ref:`Section 2.1 <section-boltzmann-eqn>` introduces the Boltzmann form of the neutron transport equation parametrized in 6-dimensional phase space over position, angle and energy. The following several sections introduce the various approximations made to this equation:

* characteristic transformation (:ref:`Section 2.2 <characteristic-transformation>`)
* energy discretization (:ref:`Section 2.3 <multi-group-approximation>`)
* discrete ordinates approximation (:ref:`Section 2.4 <discrete-ordinates-approximation>`)
* isotropic scattering approximation (:ref:`Section 2.5 <isotropic-scattering-approximation>`)
* flat source approximation (:ref:`Section 2.6 <flat-source-region-approximation>`)
* constant cross-section approximation (:ref:`Section 2.7 <constant-xs-approximation>`)
* integrating factor solution (:ref:`Section 2.8 <integrating-factor-solution>`)
* track area approximation (:ref:`Section 2.9 <track-area-approximation>`)
* azimuthal planar projection (:ref:`Section 2.10 <projection-azimuthal-plane>`)

The final equations applied in OpenMOC to solve for the FSR source and scalar flux derived in the following sections are summarized below:

**The source in each flat source region**

.. math::
   \boxed{Q_{i,g} = \frac{1}{4\pi}\left(\displaystyle\sum\limits_{g'=1}^G \Sigma^S_{i,g'\rightarrow g}\Phi_{i,g'} + \frac{\chi_{i,g}}{k_{eff}}\displaystyle\sum\limits_{g'=1}^G\nu\Sigma^F_{i,g'}\Phi_{i,g'}\right)}

**Change in angular flux along a track segment**

.. math::
   \boxed{\Delta\Psi_{k,i,g,p} = \Psi_{k,g,p}(s') - \Psi_{k,g,p}(s'') = \left(\Psi_{k,g,p}(s') - \frac{Q_{i,g}}{\Sigma^T_{i,g}}\right)(1 - e^{-\tau_{k,i,g,p}})}

**Scalar flux in each flat source region**

.. math::
   \boxed{\Phi_{i,g} = \frac{4\pi}{\Sigma_{i,g}}\left[Q_{i,g} + \frac{1}{A_i}\displaystyle\sum\limits_{k\in A_{i}}\displaystyle\sum\limits_{p=1}^{P}\omega_{m(k)}\omega_{p}\omega_{k}\sin\theta_{p}\Delta\Psi_{k,i,g,p}\right]}



.. _section-boltzmann-eqn:

Introduction to the Boltzmann Equation
======================================

The Boltzmann form of the steady-state neutron transport equation is given by the following: 

.. math::
   :label: boltzmann-eqn

   \mathbf{\Omega} \cdot \nabla \Psi(\mathbf{r},\mathbf{\Omega},E) + \Sigma^T(\mathbf{r},E)\Psi(\mathbf{r},\mathbf{\Omega},E) = \int_{0}^{\infty} \mathrm{d}E' \int_{4\pi} \mathrm{d}\mathbf{\Omega'}\Sigma^S(\mathbf{r},{\mathbf{\Omega'}\rightarrow\mathbf{\Omega}},{E'\rightarrow E}) \Psi(\mathbf{r},\mathbf{\Omega'},E') + \frac{\chi(\mathbf{r},E)}{4\pi k_{eff}} \int_{0}^{\infty} \mathrm{d}E' \nu\Sigma^F(\mathbf{r},E') \int_{4\pi} \mathrm{d}\mathbf{\Omega'}\Psi(\mathbf{r},\mathbf{\Omega'},E')

Each of the variables in use is defined in :ref:`Table 1 <table-variables>`. This is a balance equation between neutrons lost to transport, lost to absorption, produced or lost from scattering and those produced from fission. It should be noted that this equation assumes isotropic emission from fission.

.. _table-variables:

=======================   ===========
Variable                  Description
=======================   ===========
:math:`\mathbf{r}`        Spatial position vector
:math:`\mathbf{\Omega}`   Angular direction vector
:math:`E`                 Neutron energy
:math:`\Psi`              Angular neutron flux
:math:`k_{eff}`           Effective neutron multiplication factor
:math:`\Sigma^T`          Neutron total cross-section
:math:`\Sigma^S`          Neutron scattering cross-section
:math:`\Sigma^F`          Neutron fission cross-section
:math:`\chi`              Energy spectrum for fission neutrons
:math:`\nu`               Average number of neutrons emitted per fission
=======================   ===========  

**Table 1**: Variables in the Boltzmann equation.

The first step is to simplify this equation by defining those quantities on the right hand side as the total neutron source :math:`Q(\mathbf{r},\mathbf{\Omega},E)`:

.. math::
   :label: integral-source

   Q(\mathbf{r},\mathbf{\Omega},E) = \int_{0}^{\infty} \mathrm{d}E' \int_{4\pi} \mathrm{d}\mathbf{\Omega'}\Sigma^S(\mathbf{r},{\mathbf{\Omega'}\rightarrow\mathbf{\Omega}},{E'\rightarrow E}) \Psi(\mathbf{r},\mathbf{\Omega'},E') + \frac{\chi(\mathbf{r},E)}{4\pi k_{eff}} \int_{0}^{\infty} \mathrm{d}E' \int_{4\pi} \mathrm{d}\mathbf{\Omega'} \nu\Sigma^F(\mathbf{r},E')\Psi(\mathbf{r},\mathbf{\Omega'},E')

The transport equation can now be more concisely written as follows:

.. math::
   :label: transport-with-source

   \mathbf{\Omega} \cdot \nabla \Psi(\mathbf{r},\mathbf{\Omega},E) + \Sigma^T(\mathbf{r},E)\Psi(\mathbf{r},\mathbf{\Omega},E) = Q(\mathbf{r},\mathbf{\Omega},E)


.. _characteristic-transformation:

The Characteristic Transformation
=================================

The characteristic form of the Boltzmann equation is found by a change of variables by parametrizing :math:`\mathbf{r}` with respect to some reference location :math:`\mathbf{r_0}`:

.. math::
   :label: characteristics-parametrization

   \mathbf{r} = (x(s), y(s)) = (x_0+s\mathbf{\Omega_x}, y_0+s\mathbf{\Omega_y}) = \mathbf{r_0}+s\mathbf{\Omega}

For any location :math:`\mathbf{r}` of interest, each angular direction vector :math:`\mathbf{\Omega'}` is matched to a corresponding reference location :math:`\mathbf{r_{0}'}` defined such that :math:`\mathbf{r} = \mathbf{r_{0}'} + s\mathbf{\Omega'}`. This parametrization for position may be substituted into the source and transport equations to obtain the following form for each:

.. math::
   :label: source-parametrization

   Q(\mathbf{r},\mathbf{\Omega},E) = \int_{0}^{\infty} \mathrm{d}E' \int_{4\pi} \mathrm{d}\mathbf{\Omega'}\Sigma^S(\mathbf{r_0'}+s\mathbf{\Omega'},{\mathbf{\Omega'}\rightarrow\mathbf{\Omega}},{E'\rightarrow E}) \Psi(\mathbf{r_0'}+s\mathbf{\Omega'},\mathbf{\Omega'},E') + \frac{\chi(\mathbf{r_0}+s\mathbf{\Omega},E)}{4\pi k_{eff}} \int_{0}^{\infty} \mathrm{d}E' \nu\Sigma^F(\mathbf{r_0}+s\mathbf{\Omega},E') \int_{4\pi} \mathrm{d}\mathbf{\Omega'} \Psi(\mathbf{r_0'}+s\mathbf{\Omega'},\mathbf{\Omega'},E')

.. math::
   :label: boltzmann-parametrization

   \mathbf{\Omega} \cdot \nabla \Psi(\mathbf{r_0}+s\mathbf{\Omega},\mathbf{\Omega},E) + \Sigma^T(\mathbf{r_0}+s\mathbf{\Omega},E)\Psi(\mathbf{r_0}+s\mathbf{\Omega},\mathbf{\Omega},E) = Q(\mathbf{r_0}+s\mathbf{\Omega},\mathbf{\Omega},E)

Applying the differential operator to the angular flux in :eq:`boltzmann-parametrization` leads to the characteristic form of the Boltzmann equation:

.. math::
   :label: boltzmann-differential

   \frac{d}{ds}\Psi(\mathbf{r_0}+s\mathbf{\Omega},\mathbf{\Omega},E) + \Sigma^T(\mathbf{r_0}+s\mathbf{\Omega},E)\Psi(\mathbf{r_0}+s\mathbf{\Omega},\mathbf{\Omega},E) = Q(\mathbf{r_0}+s\mathbf{\Omega},\mathbf{\Omega},E)

For brevity, the remainder of this section will assume the dependence of :math:`s` on the reference position :math:`\mathbf{r_0}` and :math:`\mathbf{\Omega}` and will simplify this as :math:`\mathbf{r_0} + s\mathbf{\Omega} \rightarrow s` such that the characteristic equation can be written as the following:

.. math::
   :label: simple-boltzmann-differential

   \frac{d}{ds}\Psi(s,\mathbf{\Omega},E) + \Sigma^T(s,E)\Psi(s,\mathbf{\Omega},E) = Q(s,\mathbf{\Omega},E)

This equation can be solved through the use of an integrating factor:

.. math::
   :label: integrating-factor

   e^{-\int_{0}^s\mathrm{d}s'\Sigma^T(s',E)}

The final analytical solution to the characteristic equation is therefore:

.. math::
   :label: moc-eqn

   \Psi(s,\mathbf{\Omega},E) = \Psi(\mathbf{r_{0}},\mathbf{\Omega},E)e^{-\int_{0}^s\mathrm{d}s'\Sigma^T(s',E)} + \int_0^s\mathrm{d}s''Q(s'',\mathbf{\Omega},E)e^{-\int_{s''}^s\mathrm{d}s'\Sigma^T(s',E)}


.. _multi-group-approximation:

The Multi-Group Energy Approximation
====================================

Equation :eq:`moc-eqn` is defined with :math:`\Psi`, :math:`Q` and :math:`\Sigma^T` as continuous functions of energy. The first approximation to numerically solve this equation is to discretize the energy domain into distinct *energy groups* :math:`g \in G = \{1, 2, ..., G\}` where group :math:`g` spans the continuous range of energies from :math:`E_{g}` to :math:`E_{g-1}`. This is otherwise known as the *multi-group approximation*. The multi-group form of the Boltzmann equation is presented below:

.. math::
   :label: boltzmann-multigroup

   \mathbf{\Omega} \cdot \nabla \Psi_g(s,\mathbf{\mathbf{\Omega}}) + \Sigma^T_{g}(s)\Psi_g(s,\mathbf{\Omega}) = Q_g(s,\mathbf{\Omega})

The characteristic form of the equation given in :eq:`simple-boltzmann-differential` can also be written in multi-group form:

.. math::
   :label: characteristic-multigroup

   \frac{d}{ds}\Psi_{g}(s,\mathbf{\Omega}) + \Sigma^T_{g}(s)\Psi_{g}(s,\mathbf{\Omega}) = Q_g(s,\mathbf{\Omega})

Likewise, the multi-group form of the neutron source :eq:`source-parametrization` is given by:

.. math::
   :label: source-multigroup

   Q_g(s,\mathbf{\Omega}) = \displaystyle\sum\limits_{g'=1}^G \int_{4\pi} \mathrm{d}\mathbf{\Omega'}\Sigma_{g'\rightarrow g}^S(s,{\mathbf{\Omega'}\rightarrow\mathbf{\Omega}}) \Psi_{g'}(s,\mathbf{\Omega'}) + \frac{\chi_{g}(s)}{4\pi k_{eff}} \displaystyle\sum\limits_{g'=1}^G \nu\Sigma_{g'}^F(s) \int_{4\pi} \mathrm{d}\mathbf{\Omega'} \Psi_{g'}(s,\mathbf{\Omega'})

It directly follows from :eq:`moc-eqn` and :eq:`simple-boltzmann-differential` that the solution to the multi-group characteristic neutron transport equation is the following:

.. math::
   :label: moc-multigroup

   \Psi_g(s,\mathbf{\Omega}) = \Psi_g(\mathbf{r_{0}},\mathbf{\Omega})e^{-\int_{0}^s\mathrm{d}s'\Sigma_g^T(s')} + \int_0^s\mathrm{d}s''Q_g(s'',\mathbf{\Omega})e^{-\int_{s''}^s\mathrm{d}s'\Sigma_g^T(s')}

Where both :eq:`moc-multigroup` and :eq:`source-multigroup` make use of the energy condensed cross-sections :math:`\Sigma^T`, :math:`\Sigma^F`, :math:`\Sigma^S`, and :math:`\chi`:

.. math::
   :label: condensed-total-xs

   \Sigma_{g}^T(s) = \frac{\int_{E_{g}}^{E_{g-1}}\mathrm{d}E'\Sigma^T(s,E')\Psi(s,\mathbf{\Omega},E')}{\int_{E_{g}}^{E_{g-1}}\mathrm{d}E'\Psi(s,\mathbf{\Omega},E')}

.. math::
   :label: condensed-fission-xs

   \Sigma_{g}^F(s) = \frac{\int_{E_{g}}^{E_{g-1}}\mathrm{d}E'\Sigma^F(s,E')\Psi(s,\mathbf{\Omega},E')}{\int_{E_{g}}^{E_{g-}1}\mathrm{d}E'\Psi(s,\mathbf{\Omega},E')}


.. math::
   :label: condensed-scatter-xs

   \Sigma_{g'\rightarrow g}^S(s,\mathbf{\Omega'}\rightarrow \mathbf{\Omega}) = \frac{\int_{E_{g'}}^{E_{g'-1}}\mathrm{d}E'\int_{E_{g}}^{E_{g-1}}\mathrm{d}E''\Sigma^S(s,\mathbf{\Omega'}\rightarrow \mathbf{\Omega},E'\rightarrow E'')\Psi(s,\mathbf{\Omega'},E')}{\int_{E_{g'}}^{E_{g'-1}}\mathrm{d}E'\Psi(s,\mathbf{\Omega'},E')}

.. math::
   :label: condensed-chi

   \chi_{g'\rightarrow g}(s) = \frac{\int_{E_{g'}}^{E_{g'-1}}\mathrm{d}E'\int_{E_{g}}^{E_{g-1}}\mathrm{d}E''\chi(s,E'\rightarrow E'')\nu\Sigma^F(s,\mathbf{\Omega},E')\Psi(s,\mathbf{\Omega'},E')}{\int_{E_{g'}}^{E_{g'-1}}\mathrm{d}E'\nu\Sigma^F(s,\mathbf{\Omega},E')\Psi(s,\mathbf{\Omega'},E')}

Although :eq:`condensed-chi` assumes a dependence of :math:`\chi` on both the energy of the neutron causing fission :math:`g'` and the fission emission energy group :math:`g`, the former is typically summed over to simplify the multi-group :math:`\chi` to the following approximation:

.. math::
   :label: condensed-chi-sum

   \chi_{g}(s) = \displaystyle\sum\limits_{g=1}^{G}\chi_{g'\rightarrow g}(s)


.. _discrete-ordinates-approximation:

The Discrete Ordinates Approximation
====================================

The discrete ordinates approximation is introduced to approximate the integral over the angular domain in the source :eq:`source-multigroup`. This is equivalent to applying quadrature rules to evaluate the integral over the angular flux using a weighted sum of fluxes at specific angles where weights :math:`w_{m}` are introduced for each of the quadrature points :math:`\mathbf{\Omega_{m}} \; \forall \; m \in \{1, ..., M\}`.

.. math::
   :label: moc-quadrature

   \Phi_{g}(s) = \int_{4\pi}\mathrm{d}\mathbf{\Omega'}\Psi_{g}(s,\mathbf{\Omega'}) \approx \displaystyle\sum\limits_{m=1}^{M}w_{m}\Psi_{g}(s,\mathbf{\Omega_{m}})

The integrated angular flux :math:`\Phi_{g}(s)` is termed the *scalar flux*. Substituting this approximation to the angular flux integral into :eq:`source-multigroup` leads to the following approximation to the source :math:`Q_{m,g}(s) \approx Q_{g}(s,\mathbf{\Omega_{m}})` at each quadrature point :math:`\mathbf{\Omega_{m}}`:

.. math::
   :label: source-angular-quadrature

   Q_{m,g}(s) = \displaystyle\sum\limits_{g'=1}^G \displaystyle\sum\limits_{m'=1}^{M}w_{m'}\Sigma_{g'\rightarrow g}^S(s,{\mathbf{\Omega_{m'}}\rightarrow\mathbf{\Omega_{m}}}) \Psi_{g'}(s,\mathbf{\Omega_{m'}}) + \frac{\chi_{g}(s)}{4\pi k_{eff}} \displaystyle\sum\limits_{g'=1}^G \displaystyle\sum\limits_{m'=1}^{M}w_{m'}\nu\Sigma_{g'}^F(s)\Psi_{g'}(s,\mathbf{\Omega_{m'}})

Substituting this approximation to the source into :eq:`moc-multigroup` one obtains the characteristic solution for the angular flux :math:`\Psi_{m,g}(s) \approx \Psi_g(s,\mathbf{\Omega_{m}})` at each quadrature point :math:`\mathbf{\Omega_{m}}`:

.. math::
   :label: angular-flux-angular-quadrature

   \Psi_{m,g}(s) = \Psi_{m,g}(\mathbf{r_{0}})e^{-\int_{0}^{s}\mathrm{d}s'\Sigma_g^T(s')} + \int_0^{s_{m}}\mathrm{d}s''Q_{m,g}(s'')e^{-\int_{s''}^{s}\mathrm{d}s'\Sigma_g^T(s')}

Equations :eq:`source-angular-quadrature` and :eq:`angular-flux-angular-quadrature` may be further decomposed into azimuthal and polar angle quadratures :math:`m \in \{1, 2, ..., M\}` and :math:`p \in \{1, 2, ..., P\}` with weights :math:`w_{m}` and :math:`w_{p}` for the azimuthal plane and axial dimension, respectively:

.. math::
   :label: source-azimuthal-polar

   Q_{m,p,g}(s) = \displaystyle\sum\limits_{g'=1}^G \displaystyle\sum\limits_{m'=1}^{M} \displaystyle\sum\limits_{p'=1}^{P} w_{m'}w_{p'}\Sigma_{g'\rightarrow g}^S(s,{\mathbf{\Omega_{m',p'}}\rightarrow\mathbf{\Omega_{m,p}}}) \Psi_{g'}(s,\mathbf{\Omega_{m',p'}}) + \frac{\chi_{g}(s)}{4\pi k_{eff}} \displaystyle\sum\limits_{g'=1}^G \displaystyle\sum\limits_{m'=1}^{M} \displaystyle\sum\limits_{p'=1}^{P} w_{m'}w_{p'} \nu\Sigma_{g'}^F(s)\Psi_{g'}(s,\mathbf{\Omega_{m',p'}})

.. math:: 
   :label: angular-flux-azimuthal-polar 

   \Psi_{m,p,g}(s) = \Psi_{m,p,g}(\mathbf{r_{0}})e^{-\int_{0}^s\mathrm{d}s'\Sigma_g^T(s')} + \int_0^s\mathrm{d}s''Q_{m,p,g}(s'')e^{-\int_{s''}^s\mathrm{d}s'\Sigma_g^T(s')}


.. _isotropic-scattering-approximation:

The Isotropic Scattering Approximation
======================================

An additional approximation that is made to simplify the evaluation of the source in :eq:`source-azimuthal-polar` is to assume that the scattering source is isotropic. This approximation allows the total source to be expressed solely in terms of the scalar flux:

.. math::
   :label: source-isotropic

   Q_{g}(s) = \frac{1}{4\pi}\left(\displaystyle\sum\limits_{g'=1}^G \Sigma^S_{g'\rightarrow g}(s)\Phi_{g'}(s) + \frac{\chi_{g}(s)}{k_{eff}}\displaystyle\sum\limits_{g'=1}^G\nu\Sigma^F_{g'}(s)\Phi_{g'}(s)\right)

The subscripts :math:`m` and :math:`p` for the azimuthal and polar angles, respectively, have been dropped from :math:`Q_{g}(s)` since they have been embedded in the integral over angular phase space to obtain the scalar flux :math:`\Phi_{g}(s)`.


.. _flat-source-region-approximation:

The Flat Source Region Approximation
====================================

Another common approximation for MOC is to assume that the source :math:`Q_g` is constant across discrete spatial cells termed *flat source regions* (FSRs). This implies that the source does not vary along a characteristic :math:`k` entering FSR :math:`i` at :math:`s'` and exiting at :math:`s''`: 

.. math::
   :label: flat-source

   Q_{i,g} = Q_{g}(s') = Q_{g}(s'') = Q_{g}(s) \;\;\; , \;\;\; s \in [s', s'']

.. _linear-source-region-approximation:

The Linear Source Region Approximation
======================================

A more accurate description of the source in spatial regions is to assume a linear variation. This is typically sufficient for the moderator in a PWR, when each channel is also cut in azimuthal source regions. The source then varies along each characteristic lines. The reader should refer themselves to [Ferrer]_ and [Gunow]_ for more details on the track-based linear source approximation and its implementation in OpenMOC.

.. math::
   :label: linear-source

   Q_{g}(s) = q_{t,g,0} + q_{t,g,1} (s - l_{t} / 2)

:math:`l_{t}` is the length of the segment considered, while :math:`q_{g,0}` and :math:`q_{g,1}` are track dependent coefficients that describe the source.

.. _constant-xs-approximation:

The Constant Cross-Section Approximation
========================================

In addition to the flat source approximation, it is assumed that the material properties are constant across each FSR. The area-averaged cross-sections for FSR :math:`i \in \{1, 2, ..., I\}` with area :math:`A_{i}` are defined as:

.. math::
   :label: area-averaged-total-xs

   \Sigma_{i,g}^{T} = \frac{\int_{\mathbf{r}\in A_{i}}\mathrm{d}\mathbf{r}\Sigma_{g}^T(\mathbf{r})\Phi_{g}(\mathbf{r})}{\int_{\mathbf{r}\in A_{i}}\mathrm{d}\mathbf{r}\Phi_{g}(\mathbf{r})}

.. math::
   :label: area-averaged-fission-xs

   \Sigma_{i,g}^{F} = \frac{\int_{\mathbf{r}\in A_{i}}\mathrm{d}\mathbf{r}\Sigma_{g}^F(\mathbf{r})\Phi_{g}(\mathbf{r})}{\int_{\mathbf{r}\in A_{i}}\mathrm{d}\mathbf{r}\Phi_{g}(\mathbf{r})}

.. math::
   :label: area-averaged-scatter-xs

   \Sigma_{i,g'\rightarrow g}^{S} = \frac{\int_{\mathbf{r}\in A_{i}}\mathrm{d}\mathbf{r}\Sigma_{g'\rightarrow g}^S(\mathbf{r})\Phi_{g'}(\mathbf{r})}{\int_{\mathbf{r}\in A_{i}}\mathrm{d}\mathbf{r}\Phi_{g'}(\mathbf{r})}

.. math::
   :label: area-averaged-chi

   \chi_{i,g} = \frac{\int_{\mathbf{r}\in A_{i}}\mathrm{d}\mathbf{r}\chi_{g}(\mathbf{r})}{\int_{\mathbf{r}\in A_{i}}\mathrm{d}\mathbf{r}}

The flat source term :math:`Q_{i,g}` for FSR :math:`i` with area :math:`A_i` is defined in terms of both fission and scattering from the area-averaged scalar flux :math:`\Phi_{g,i}` within the FSR:

.. math::
   :label: final-source

   Q_{i,g} = \frac{1}{4\pi}\left(\displaystyle\sum\limits_{g'=1}^G \Sigma^S_{i,g'\rightarrow g}\Phi_{i,g'} + \frac{\chi_{i,g}}{k_{eff}}\displaystyle\sum\limits_{g'=1}^G\nu\Sigma^F_{i,g'}\Phi_{i,g'}\right)

.. math::
   :label: area-averaged-scalar-flux

   \Phi_{i,g} = \frac{\int_{\mathbf{r}\in A_{i}}\mathrm{d}\mathbf{r}\Phi_{g}(\mathbf{r})}{\int_{\mathbf{r}\in A_{i}}\mathrm{d}\mathbf{r}}

The multi-group nuclear cross-sections for each FSR are an input to OpenMOC. As a result, the area-averaging integrals must be performed by some pre-processing method such as Monte Carlo.


.. _integrating-factor-solution:

The Integrating Factor Solution
===============================

Each chracteristic may be discretized into *segments* across individual FSRs. This approximation allows :eq:`angular-flux-azimuthal-polar` to be localized to a segment of characteristic :math:`k` across FSR :math:`i` from its entry point at :math:`s'` to exit point at :math:`s''`. By defining the integrating factor in terms of the optical length :math:`\tau_{k,i,g} = \Sigma^T_{i,g}(s''-s')` one may analytically evaluate the integrals in :eq:`angular-flux-azimuthal-polar` and express the outgoing flux along the characteristic as follows:

.. math::
   :label: angular-flux-fsr

   \Psi_{k,g}(s'') = \Psi_{k,g}(s')e^{-\tau_{k,i,g}} + \frac{Q_{i,g}}{\Sigma^T_{i,g}}(1 - e^{-\tau_{k,i,g}})

With minor algebraic rearrangement, the change in angular flux along the characteristic is given by the following:

.. math::
   :label: delta-angular-flux-fsr

   \Delta\Psi_{k,g} = \Psi_{k,g}(s') - \Psi_{k,g}(s'') = \left(\Psi_{k,g}(s') - \frac{Q_{i,g}}{\Sigma^T_{i,g}}\right)(1 - e^{-\tau_{k,i,g}})


.. _track-area-approximation:

The Track Area Approximation
============================

The key quantity remaining to be determined is the integral over area for the FSR area-averaged scalar flux :math:`\Phi_{g,i}` in :eq:`area-averaged-scalar-flux`. The track area approximation is used to compute this value numerically. 

First, define :math:`l_{k,i}=s''-s'` such that the average angular flux in FSR :math:`i` along characteristic :math:`k` is the following integral:

.. math::
   :label: avg-angular-flux-integral

   \overline{\Psi}_{k,i,g} = \frac{1}{l_{k,i}}\int_{s'}^{s''} \Psi_{k,i,g}(s) \mathrm{d}s

Upon evaluating the integral, the average angular flux along the characteristic can be reduced to the following algebraic expression:

.. math::
   :label: avg-angular-flux

   \overline{\Psi}_{k,i,g} = \frac{1}{l_{k,i}}\left[\frac{\Psi_{k,g}(s')}{\Sigma_{i,g}^T}(1 - e^{-\tau_{k,i,g}}) + \frac{l_{k,i}Q_{i,g}}{\Sigma_{i,g}^T}\left(1 - \frac{(1 - e^{-\tau_{k,i,g}})}{\tau_{k,i,g}}\right)\right]

Assuming a constant source and cross-sections in FSR :math:`i`, the value given for the average angular flux in :eq:`avg-angular-flux` is exact. In order to exactly compute the area-averaged scalar flux, the average angular flux from every characteristic crossing FSR :math:`i` must be taken into account. This is numerically intractable; hence, an appropriate subset :math:`K` of characteristics, henceforth known as *tracks*, is chosen and the integral over the area of the FSR is performed using quadrature rules with a weight :math:`w_{k}` for each track :math:`k \in K` crossing through the FSR :math:`k \in A_{i}`. The contribution :math:`\overline{\Psi}_{k,i,g}` of track :math:`k` with azimuthal and polar quadrature weights denoted by :math:`w_{m(k)}` and :math:`w_{p(k)}`, respectively, is then integrated to find the area-averaged scalar flux in FSR :math:`i` as follows:

.. math::
   :label: area-averaged-scalar-flux-quadrature

   \Phi_{i,g} = \frac{\int_{\mathbf{r}\in A_{i}}\mathrm{d}\mathbf{r}\int_{4\pi}\mathrm{d}\mathbf{\Omega}\Psi_{g}(\mathbf{r},\mathbf{\Omega})}{\int_{\mathbf{r}\in A_{i}}\mathrm{d}\mathbf{r}} \approx \frac{4\pi\displaystyle\sum\limits_{k\in A_{i}}w_{m(k)}w_{p(k)}w_{k}l_{k,i}\sin\theta_{p(k)}\overline{\Psi}_{k,i,g}}{\displaystyle\sum\limits_{k\in A_{i}}w_kl_{k,i}\sin\theta_{p(k)}}

In :eq:`area-averaged-scalar-flux-quadrature`, the angle :math:`\theta_{p(k)}` formed by characteristic :math:`k` with respect to the polar axis is introduced to project the length of the characteristic segment :math:`l_{k,i}` onto the azimuthal plane. In this application of quadrature to approximate an area integral, the weights can be thought of as the *effective width* of each track :math:`k`.

The denominator in :eq:`area-averaged-scalar-flux-quadrature` then simplifies to the area :math:`A_i`:

.. math::
   :label: avg-scalar-flux-quadrature

   \Phi_{i,g} \approx \frac{4\pi}{A_{i}}\displaystyle\sum\limits_{k\in A_{i}}w_{m(k)}w_{p(k)}w_{k}l_{k,i}\sin\theta_{p(k)}\overline{\Psi}_{k,i,g}

The scalar flux can be found in terms of average angular fluxes from each track by substituting the expression for the average angular flux from :eq:`avg-angular-flux` into :eq:`avg-scalar-flux-quadrature` and rearranging:

.. math::
   :label: avg-scalar-flux-v2

   \Phi_{i,g} = \frac{4\pi}{\Sigma_{i,g}}\left[Q_{i,g} + \frac{1}{A_i}\displaystyle\sum\limits_{k\in A_{i}}\omega_{m(k)}\omega_{p(k)}\omega_{k}\sin\theta_{p(k)}\left(\Psi_{k,i,g}(s') - \frac{Q_{i,g}}{\Sigma_{i,g}^T}\right)(1 - e^{-\tau_{k,i,g}})\right]

The final form for the scalar flux can be simplified in terms of the change in angular flux :math:`\Delta\Psi_{k,i,g}` along each track segment as defined in :eq:`delta-angular-flux-fsr`:

.. math::
   :label: avg-scalar-flux-delta-angular-flux

   \Phi_{i,g} = \frac{4\pi}{\Sigma_{i,g}}\left[Q_{i,g} + \frac{1}{A_i}\displaystyle\sum\limits_{k\in A_{i}}\omega_{m(k)}\omega_{p(k)}\omega_{k}\sin\theta_{p(k)}\Delta\Psi_{k,i,g}\right]


.. _projection-azimuthal-plane:

Projection from the Azimuthal Plane
===================================

The preceding sections used track segment lengths :math:`l_{k,i}` in 3D. In practice, the memory footprint for storing track segment data is greatly reduced if the polar angle quadrature is replicated for each azimuthal quadrature point. Such a quadrature allows for track segments to be stored in the 2D azimuthal plane and projected into 3D for each polar angle when necessary. The projection results in some minor changes to the equations presented in the previous sections. 

In what follows, each track segment length :math:`l_{k,i}` will be assumed to reside within the azimuthal plane. Likewise, the optical length :math:`\tau_{k,i,g} = \Sigma^T_{k,i,g}l_{k,i}` also resides in the azimuthal plane. For notational simplicity, the 3D projection of the track segment length for polar angle :math:`p` will be denoted by :math:`l_{k,i,p} = \frac{l_{k,i}}{\sin\theta_{p}}` and the optical length by :math:`\tau_{k,i,g,p} = \Sigma^T_{k,i,g}l_{k,i,p}`.

First, the polar angle must be accounted for in the expression for the track segment average angular flux to project the segment length into the polar dimension:

.. math::
   :label: avg-angular-flux-polar

   \overline{\Psi}_{k,i,g,p} = \frac{1}{l_{k,i,p}}\left[\frac{\Psi_{k,g,p}(s')}{\Sigma_{i,g}^T}(1 - \exp(-\tau_{k,i,g,p})) + \frac{l_{k,i,p}Q_{i,g}}{\Sigma_{i,g}^T}\left(1 - \frac{(1 - \exp(-\tau_{k,i,g,p}))}{\tau_{k,i,g,p}}\right)\right]

Next, :math:`\sin\theta_{p(k)}` is dropped and a summation over polar angles is incorporated into the area-averaged scalar flux in :eq:`avg-scalar-flux-quadrature`:

.. math::
   :label: avg-scalar-flux-quadrature-polar

   \Phi_{i,g} = \frac{4\pi}{A_i}\displaystyle\sum\limits_{k \in A_i}\displaystyle\sum\limits_{p=1}^{P}\omega_{m(k)}\omega_{p}\omega_{k}l_{k,i}\overline{\Psi}_{k,i,g,p}

The scalar flux can be found in terms of average angular fluxes from each track by substituting the expression for the average angular flux from :eq:`avg-angular-flux-polar` into :eq:`avg-scalar-flux-quadrature-polar` and rearranging:

.. math::
   :label: avg-scalar-flux-polar

   \Phi_{i,g} = \frac{4\pi}{\Sigma_{i,g}}\left[Q_{i,g} + \frac{1}{A_i}\displaystyle\sum\limits_{k\in A_{i}}\displaystyle\sum\limits_{p=1}^{P}\omega_{m(k)}\omega_{p}\omega_{k}\sin\theta_{p}\left(\Psi_{k,i,g,p}(s') - \frac{Q_{i,g}}{\Sigma_{i,g}^T}\right)(1 - e^{-\tau_{k,i,g,p}})\right]

The final form for the scalar flux can be simplified in terms of the change in angular flux :math:`\Delta\Psi_{k,i,g,p}` along each track segment as defined in :eq:`avg-scalar-flux-delta-angular-flux`:

.. math::
   :label: avg-scalar-flux-polar-final

   \Phi_{i,g} = \frac{4\pi}{\Sigma_{i,g}}\left[Q_{i,g} + \frac{1}{A_i}\displaystyle\sum\limits_{k\in A_{i}}\displaystyle\sum\limits_{p=1}^{P}\omega_{m(k)}\omega_{p}\omega_{k}\sin\theta_{p}\Delta\Psi_{k,i,g,p}\right]

This is the form of the transport equation solved by the MOC formulation used in OpenMOC.


References
==========

.. [Askew] J. Askew, "A Characteristics Formulation of the Neutron Transport Equation in Complicated Geometries." Technical Report AAEW-M 1108, UK Atomic Energy Establishment (1972).

.. [Boyd] W. Boyd, "Massively Parallel Algorithms for Method of Characteristics Neutral Particle Transport on Shared Memory Computer Architectures." M.S. Thesis, Massachusetts Institute of Technology (2014). 

.. [Ferrer] R. Ferrer and J. Rhodes, “A Linear Source Approximation Scheme for the Method of Characteristics,” volume 77, p. 119–136, 1981.

.. [Gunow] G. Gunow "Full Core 3D Neutron Transport Simulation Using the Method of Characteristics with Linear Sources", PhD Thesis, Massachusetts Institute of Technology (2018).

