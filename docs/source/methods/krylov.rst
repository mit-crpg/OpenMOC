.. _krylov:

=======================
Krylov Subspace Methods
=======================

The MOC equations are traditionally solved using a non-linear power iteration method.  This is not the only way to solve MOC, but many of the alternatives require the equations to be rewritten in a linear fashion.  To convert into a linear operator form, we need to essentially repeat the derivation shown in :ref:`Section 2 <method_of_characteristics>`, with a slight modification.  Starting with the Boltzmann transport equation:

.. math::
   \mathbf{\Omega} \cdot \nabla \Psi(\mathbf{r},\mathbf{\Omega},E) + \Sigma^T(\mathbf{r},E)\Psi(\mathbf{r},\mathbf{\Omega},E) = \int_{0}^{\infty} \mathrm{d}E' \int_{4\pi} \mathrm{d}\mathbf{\Omega'}\Sigma^S(\mathbf{r},{\mathbf{\Omega'}\rightarrow\mathbf{\Omega}},{E'\rightarrow E}) \Psi(\mathbf{r},\mathbf{\Omega'},E') + \frac{\chi(\mathbf{r},E)}{4\pi k_{eff}} \int_{0}^{\infty} \mathrm{d}E' \nu\Sigma^F(\mathbf{r},E') \int_{4\pi} \mathrm{d}\mathbf{\Omega'}\Psi(\mathbf{r},\mathbf{\Omega'},E')

Instead of replacing the right hand side with :math:`Q(\mathbf{r}, \mathbf{\Omega}, E)`, we replace it with a fission and a scatter source.  The idea is to factor out the eigenvalue.  The resulting source terms are:

**Scatter Source**

.. math::
   Q^S(\mathbf{r}, \mathbf{\Omega}, E) &= \int_{0}^{\infty} \mathrm{d}E' \int_{4\pi} \mathrm{d}\mathbf{\Omega'}\Sigma^S(\mathbf{r},{\mathbf{\Omega'}\rightarrow\mathbf{\Omega}},{E'\rightarrow E}) \Psi(\mathbf{r},\mathbf{\Omega'},E')

**Fission Source**

.. math::
   \frac{1}{k_\text{eff}} Q^F(\mathbf{r}, \mathbf{\Omega}, E) &= \frac{\chi(\mathbf{r},E)}{4\pi k_\text{eff}} \int_{0}^{\infty} \mathrm{d}E' \nu\Sigma^F(\mathbf{r},E') \int_{4\pi} \mathrm{d}\mathbf{\Omega'}\Psi(\mathbf{r},\mathbf{\Omega'},E')

Similar to before, the transport equation can be more concisely written as:

.. math::
   \mathbf{\Omega} \cdot \nabla \Psi(\mathbf{r},\mathbf{\Omega},E) + \Sigma^T(\mathbf{r},E)\Psi(\mathbf{r},\mathbf{\Omega},E) =  Q^S(\mathbf{r}, \mathbf{\Omega}, E) +\frac{1}{k_\text{eff}} Q^F(\mathbf{r}, \mathbf{\Omega}, E)

From here, we follow the exact same process all the way to the end.  The MOC sweep operator :math:`\mathcal{T}` (explicitly for a vacuum boundary condition, as non-vacuum conditions add another non-linearity) is given by the following equation:

.. math::
   \Phi_{i,g} = \mathcal{T}(Q_{i,g}) \Phi_{i,g} &= \frac{4\pi}{\Sigma_{i,g}}\left[Q_{i,g} + \frac{1}{A_i}\displaystyle\sum\limits_{k\in A_{i}}\omega_{m(k)}\omega_{p(k)}\omega_{k}\sin\theta_{p(k)}\Delta\Psi_{k,i,g}(Q_{i,g}) \right]

Since this is the linear Boltzmann transport equation, we can separate the sources (as was done earlier), solve for them separately, and add the results back together.  Each source will generate its own angular fluxes as well.

.. math::
   \Phi_{i,g} = \mathcal{T}(Q^S_{i,g}) \Phi_{i,g} + \frac{1}{k_\text{eff}}\mathcal{T}(Q^F_{i,g}) \Phi_{i,g}

One can now view the above equation in the light of linear algebra, as each of these new :math:`\mathcal{T}` operators no longer have any dependence on the eigenvalue, and are themselves linear with respect to the flux.  They can now be treated as linear operators and manipulated.

.. math::
   \Phi &= \mathcal{T}(Q^S) \Phi + \frac{1}{k_\text{eff}} \mathcal{T} (Q^F) \Phi
   
.. math::
   \Phi &= \mathbf{T}_S \Phi + \frac{1}{k_\text{eff}} \mathbf{T}_F \Phi

.. math::
   k_\text{eff}\left(\mathbf{I} - \mathbf{T}_S\right) \Phi &= \mathbf{T}_F \Phi

.. math::
   k_\text{eff} \Phi &= \left(\mathbf{I} - \mathbf{T}_S\right)^{-1}\mathbf{T}_F \Phi

As each term here is not defined as a matrix, but as a linear operator, any algorithm that operates on it must necessarily only use operator-vector products.  Krylov subspace algorithms, such as the generalized minimal residual method (GMRES) for the inversion, and implicitly restarted Arnoldi iteration (IRAM) for the eigenvalue problem itself, are ideal.  IRAM also gives the added benefit of being able to calculate multiple eigenmodes of the problem, which can be used in stability analysis and dominance ratio calculations.
