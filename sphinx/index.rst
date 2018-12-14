Neutron Star Cooling
====================

Background
----------

This discussion is based on [Yakovlev01ne]_, [Page04mc]_, and Dany's
original documentation
(http://www.astroscu.unam.mx/neutrones/NSCool/NSCool_Guide_1.1_CodeStructure.pdf).
We will use the metric

.. math::
   ds^2 = c^2 e^{2 \phi} dt^2 - \left( 1-\frac{2 G m}{r
   c^2}\right)^{-1} dr^2 - r^2 \left( d\theta^2 +
   \sin^2 \theta d \phi^2 \right)

where :math:`r` is the radius, :math:`t` is time,
:math:`m` is the enclosed gravtiational mass, and :math:`\phi` is the
gravitational potential. The two main equations are those of
energy balance, 

.. math::
   \frac{d\left(L e^{2 \phi}\right)}{dr} = 
   - \frac{4 \pi r^2 e^{\phi}}{\sqrt{1-2 G m/c^2 r}} 
   \left[ C_v \frac{dT}{dt} + e^{\phi} \left( Q_{\nu} - 
   Q_h\right) \right] \, ,

and energy transport, 

.. math::
   \frac{d\left(T e^{\phi}\right)}{dr} = - \frac{1}{\lambda}
   \cdot \frac{L e^{\phi}}{4 \pi r^2 \sqrt{1-2 G m/c^2 r}} \, .

where :math:`L` is the luminosity (heat flux transported by photons
through a sphere of radius :math:`r`) , :math:`T` is the temperature,
:math:`\lambda` is the thermal conductivity (sometimes denoted
:math:`\kappa`), :math:`n_B` is the baryon number density, :math:`C_V`
is the specific heat capacity at constant volume, :math:`Q_{\nu}` is
the neutrino luminosity and :math:`Q_h` is the total of all heating
processes. Denoting the observed temperature at infinity, :math:`{\cal
T} = e^{\phi} T` and luminosity at infinity :math:`{\cal L} = e^{2
\phi} L`, and a defining a new baryon number coordinate

.. math::
   d a = \frac{4 \pi r^2 n_B~dr}{\sqrt{1 - 2 G m/c^2 r}} \, ,

one gets new equations

.. math::
   \frac{d {\cal L}}{d a} = - \frac{C_v}{n_B} \frac{d {\cal T}}{d t} 
   - e^{2 \phi} \frac{Q_\nu - Q_h}{n_B}

or 

.. math::
   \frac{d {\cal T}}{d t} = - e^{2 \phi} \frac{Q_{\nu} - Q_h}{C_v}
   - \frac{n_B}{C_v} \frac{d {\cal L}}{d a}

and 

.. math::
   \frac{d {\cal T}}{d a} = -\frac{1}{\lambda} \frac{{\cal L}}
   {\left(4 \pi r^2\right)^2 n_B e^{\phi}}

or 

.. math::
   {\cal L} = - \lambda \left( 4 \pi r^2\right)^2 n_B e^{\phi} 
   \frac{d {\cal T}}{d a}

Miscellaneous documentation to be updated later
-----------------------------------------------

In the diffusion approximation, the thermal conductivity due to the
transport of photons, :math:`\lambda_{\gamma}` (in units of
:math:`\mathrm{cm}^2/\mathrm{s}`) is
      
.. math::
   \lambda_{\gamma} = \frac{4 a c T^3}{3 \kappa \rho^2 c_P}
   
where :math:`a=4 \sigma_{\mathrm{SB}}/c = 7.566 \times 10^{-15}
\mathrm{erg}~\mathrm{cm}^{-3}~\mathrm{K}^{-4}` is the radiation constant,
:math:`\kappa` is the opacity, :math:`\rho` is the mass density and
:math:`c_P` is the specific heat capacity at constant pressure in
:math:`\mathrm{erg}~\mathrm{g}^{-1} \mathrm{K}^{-1}`. The
total thermal conductivity, is thus the sum of the individual
contributions

.. math::
   \lambda = \lambda_{\gamma} + \sum {\lambda}_{\mathrm{other}}

where the other contributions are computed in the `conduct`
subroutine.

The inner boundary condition is :math:`L(r=0)=0`. The outer
conditions are

.. math::
   L(r_b) = 4 \pi R^2 \sigma_{SB} \left[T_e(T_b)\right]^4

and :math:`T_b = T(r_b)` where :math:`r_b` is the 
boundary radius (typically defined as the radius at
which the mass(?) density is equal to :math:`10^{10}~\mathrm{g}
/\mathrm{cm}^3` ).

Here, I list some relevant variables in the `NSCool` subroutine and
their definition in terms of the quantities above. The luminosity is
measured in units of the solar luminosity, :math:`L_{\odot} =
3.826\times10^{33}~\mathrm{erg}/\mathrm{s}`, and the following
quantities are

.. math::
   \mathtt{acd}=\frac{a c}{3 \kappa \rho}

.. math::
   \mathtt{fp} = \lambda_{\mathrm{other}} +
   \frac{4 a c}{3 \kappa \rho} \frac{n_B T^3}{L_{\odot}}

.. math::
   \mathrm{a2ephi} = (4 \pi r^2)^2 e^{\phi}

.. math::
   \mathtt{lum} = {\cal L} = - \mathtt{fp} \times \mathtt{a2ephin}
   \times \frac{d {\cal T}}{d a}

.. math::
   \mathtt{fq} = n_B/C_V L_{\odot}

.. math::
   \mathtt{fr} = e^{2 \phi} (Q_{\nu} - \mathrm{Heat})/C_V -
   e^{\phi}P/C_V + \left(\mathrm{log}~\rho_{\mathrm{new}} -
   \mathrm{log}~\rho_{\mathrm{old}} \right)/dt \times
   \mathtt{contraction}

Emissivity list
---------------

- qeebrem: electron-electron pair bremsstrahlung in the crust (1)
- qnpb: neutrino pair bremsstrahlung in the crust (2)
- qplasma: plasma neutrino (3)
- qsynch: synchrotron neutrino (4)
- qbubble: bubble neutrino (5)
- qpair: neutrino pair production (6)
- qphoto: neutrino photoproduction in the crust (7)
- qbrem_nn: nucleon bremsstrahlung (8)
- qmurca_nucl: nucleon modified Urca (9)
- qbrem_nucl: nucleon bremsstrahlung (10)
- qmurca_hyp: hyperon modified urca (11)
- qbrem_hyp: hyperon bremsstrahlung (12)
- qdurca_np: :math:`n-p` direct Urca (13)
- qdurca_lap: :math:`\Lambda-p` direct Urca (14)
- qdurca_smn: :math:`\Sigma^{-}-n` direct Urca (15)
- qdurca_smla: :math:`\Sigma^{-}-\Lambda` direct Urca (16)
- qdurca_sms0: :math:`\Sigma^{-}-\Sigma^{0}` direct Urca (17)
- qfast: kaon and pion neutrino emission (18)
- qdurca_q: quark direct Urca (19)
- qmurca_q: quark modified Urca (20)
- qpbf_n1s0: singlet neutron PBF (21)
- qpbf_p1s0: singlet proton PBF (22)
- qpbf_n3p2: triplet neutron PBF (23)
- qpbf_q: quark PBF (24)
   
.. toctree::
   :maxdepth: 2

   bib
   nscool_wrap
   tc
   emissivities
   wrap_funct

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

   
  
