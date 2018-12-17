/*
  -------------------------------------------------------------------
  
  Copyright (C) 2015-2018, Andrew W. Steiner
  
  This file is part of nscool_wrap.
  
  nscool_wrap is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.
  
  nscool_wrap is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with O2scl. If not, see <http://www.gnu.org/licenses/>.

  -------------------------------------------------------------------
*/
#ifndef NSCOOL_WRAP_H
#define NSCOOL_WRAP_H

#include <vector>

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/table.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>
#include <o2scl/eos_had_base.h>
#include <o2scl/nstar_cold.h>
#include <o2scl/fermion_eff.h>
#include <o2scl/permutation.h>
#include <o2scl/lib_settings.h>

#include "emissivities.h"
#include "tc.h"

// Forward definition 
class nscool_wrap;

// Global pointers to nscool class
std::vector<nscool_wrap *> nscool_wrap_ptrs;

// Declaration for FORTRAN-defined main cooling subroutine
extern "C" void nscool_(int *irank, int *iret, double *neebrem_logt,
			double *neebrem_nalpha, double *neebrem_n2,
			double *sf_lgtau1, double *sf_logtau2,
			double *sf_lgr, double *sf_lgr2);

// Declaration for Fortran spline function
extern "C" void spline_(double *X, double *Y, int *IN,
			double *YP1, double *YP2, double *Y2);

// Declaration for Fortran spline2 function
extern "C" void spline2_(double *x1a, double *x2a, double *ya,
			 int *im, int *in, double *y2a);

// Ublas vector typedef
typedef boost::numeric::ublas::vector<double> ubvector;

/** \brief The HHJ parameterized EOS
 */
class eos_had_hhj : public o2scl::eos_had_eden_base {

 public:

  double s;
  double gamma;
  double eps0;
  double S0;
  double n0;
    
  /** \brief The coefficient of the pressure modification
   */
  double delta_P_K;

  /** \brief The exponent of the pressure modification
   */
  double delta_P_Gamma;

  eos_had_hhj() {
    s=0.1;
    gamma=0.7;
    eps0=15.8/o2scl_const::hc_mev_fm;
    S0=32.0/o2scl_const::hc_mev_fm;
    n0=0.16;
    delta_P_K=0.0;
    delta_P_Gamma=0.0;
  }

  /** \brief Equation of state as a function of density
   */
  virtual int calc_e(o2scl::fermion &ln, o2scl::fermion &lp,
		     o2scl::thermo &lth) {

    double barn=ln.n+lp.n;
      
    double xp;
    if (barn<=0.0) {
      xp=0.0;
      lth.ed=0.0;
      ln.mu=0.0;
      lp.mu=0.0;
      lth.pr=0.0;
      return 0;
    } else {
      xp=lp.n/barn;
    }

    double u=barn/n0;
    double sym=S0*pow(u,gamma);
    double symp=S0*gamma*pow(u,gamma-1.0);
      
    lth.ed=ln.m*ln.n+lp.m*lp.n+barn*(eps0*u*(u-2.0-s)/(1.0+s*u)+
				     sym*pow((1.0-2.0*xp),2.0));
    ln.mu=ln.m+(lth.ed-ln.m*ln.n-lp.m*lp.n)/barn+u*
      (eps0*(u-2.0-s)/(1.0+s*u)+eps0*u/(1.0+s*u)-eps0*u*(u-2.0-s)*s/
       pow(1.0+s*u,2.0)+pow(1-2.0*xp,2.0)*symp)+xp*4.0*(1-2.0*xp)*sym;
    lp.mu=lp.m+ln.mu-ln.m-4.0*(1.0-2.0*xp)*sym;
    lth.pr=-lth.ed+ln.mu*ln.n+lp.mu*lp.n;

    ln.kf=cbrt(3.0*o2scl_const::pi2*ln.n);
    lp.kf=cbrt(3.0*o2scl_const::pi2*lp.n);

    if (delta_P_Gamma>0.0 && lth.ed>1.5) {
      lth.pr+=(delta_P_K*pow(lth.ed,delta_P_Gamma)-
	       delta_P_K*pow(1.5,delta_P_Gamma));
    }
    
    return 0;
  }

};

/** \brief Base wrapper for neutron star cooling

    <b>SVN benchmark commits:
    2160, 2552, 4036

    <b>General notes</b>
    
    The purpose for this (somewhat obtuse) wrapper of Dany's Fortran
    NSCool code is to allow the use of O2scl EOS routines and to
    implement some MPI calls while still allowing C++ extensions of
    the original functionality with inheritance and virtual functions.

    Basic usage is to instantiate the class and set the global
    pointer list <tt>nscool_wrap_ptrs</tt>, change the input EOS
    or structure tables, and then call \ref run() .

    <hr>
    <b>Grid definition</b>

    The temperature is defined at odd indices and the luminosities at
    even indices. The values \c rhocore, \c rhodrip, and \c rhomax are
    the densities at the indexes defined by \c icore, \c idrip, and \c
    imax. The value \c rhosurf defines the surface density. The
    parameter \c rhoenv defines the envelope. If \c rhoenv is smaller
    than \c rhosurf, then the envelope is ignored. The parameter \c
    icore defines how many zones will be in the core, and \c idec
    gives the number of points per decade in density in the crust.

    <hr>
    <b>Pairing specification</b>

    Parameter:
    - sfn1s0: neutron 1S0 gap model to be used 
    - sfn3p2: neutron 3P2 gap model to be used
    - sfp1s0: proton 1S0 gap model to be used
    - sfl1s0: lambda hyperon (\f$ \Lambda \f$) 1S0 gap model to be used
    - fn1s0: scaling factor for neutron 1S0 gap 
    - fn3p2: scaling factor for neutron 3P2 gap
    - fp1s0: scaling factor for proton 1S0 gap
    - fl1s0: scaling factor for lambda hyperon (\f$ \Lambda \f$) 1S0 gap

    sfn1s0: 
    - 1: SFB
    - 2: CCDK
    - 3: WAP
    - 4: GC
    - 5: GIPSF
    - 201: Ioffe 1NS
    - 202: Ioffe 2NS
    - 203: Ioffe 3NS
    - 150: Three-parameter Gaussian

    sfn3p2:
    - 1: HGRR
    - 2: AO
    - 3: AO M1
    - 4: T72
    - 5: T72 M1
    - 6: BCLL92
    - 7: EEHJO96 NR
    - 8: EEHJO96 R
    - 101: Gap "a"
    - 102: Gap "b"
    - 103: Gap "c"
    - 150: three-parameter Gaussian
    - 201: Ioffe 1NT
    - 202: Ioffe 2NT
    - 203: Ioffe 3NT
    - >1000: Uniform

    sfp1s0:
    - 1: CCY MS
    - 2: CCY PS
    - 3: T73
    - 4: NS
    - 5: AO
    - 6: BCLL92
    - 7: CCDK
    - 21: T72
    - 22: AWP 2
    - 23: AWP 3
    - 201: Ioffe 1P
    - 202: Ioffe 2P
    - 203: Ioffe 3P
    - >1000: Uniform
    - 150: three-parameter Gaussian

    <hr>
    <b>Core conductivity</b>

    \verbatim embed:rst
    Based on [Baiko01tc]_, [Baym69si]_, [Gnedin95tc]_, 
    [Shternin07eh]_, and [Flowers81tp]_.
    \endverbatim

    <hr>
    <b>Crust conductivity</b>

    \verbatim embed:rst
    Based on [Shternin06et]_, [Baiko95ta]_, [Potekhin99tp]_, ...
    \endverbatim

    <hr>
    <b>Opacities</b>

    The parameter \c iopacity gives no photon opacity and a value of 1
    gives normal photon opacity. The parameter \c Q_imp gives the
    impurity parameter for electron-impurity scattering.

    The parameter \c icon_core takes a value of 1 for the
    simple Flowers and Itoh formula
    \f[
    \lambda = 10^{23} 
    \left( \frac{k_{F,n}}{1.6~\mathrm{fm}^{-1}} \right) 
    \left( \frac{T}{10^8~\mathrm{K}} \right)
    \f]
    and a value of 2 uses the full calculation of Yakovlev et al. 

    ICON_CRUST :

    This will essentially distinguish between the Itoh et al. and
    Yakovlev et al. calculations, both in the liquid and the crystal
    phases:

    Gamma > Gammac:
    1: e-phonon from Itoh et al. + e-impurity from Yakovlev & Urpin.
    2: e-phonon from Baiko & Yakovlev + e-impurity from Yakovlev & Urpin.
    3: e-phonon from Gnedin et al. (2001: appendix) + e-impurity from 
    Yakovlev & Urpin.

    Gamma < Gammal:
    1: e-ion from Itoh et al.
    2: e-ion from Itoh et al.
    3: e-ion from Gnedin et al. (2001: appendix). Gammal <Gamma<Gammac:

    interpolate between the two previous cases (to avoid a
    discontinuity in \f$ \lambda \f$ in cases 1 & 2). [If you set
    Gammal = Gammac then, of course, there will be no interpolation !]

    If rho < 107 g cm-3 (“envelope”): none of the above, just use
    Potekhin et al. (1999). After all this the e-e scattering
    contribution (\ref Shternin06) is added.

    Gammac = gammacryst > gammaliq = Gammal are defined in the
    included file gamma_limits.inc.f
    
    <hr>
    <b>Effective masses</b>

    Dany's original version supported several different models for the
    core nucleon effective masses from the literature. In this
    version, the effective mass in the core must be provided in the
    core composition table. The neutron effective mass in the crust is
    determined by a simple function of the neutron Fermi momentum.

    <hr>
*/
class nscool_wrap {
  
 public:

  /** \brief Object for computing critical temperatures
   */
  tc atc;
  
  /** \brief Object for computing emissivities
   */
  emissivities emis;
  
  /** \brief Parameter for the envelope composition
   */
  double eta;
  
  /** \name Superfluid parameters
   */
  //@{
  /// Default 1 (SFB)
  int sfn1s0;
  /// Default 101 (Minimal gap "a")
  int sfn3p2;
  /// Default 3 (T73)
  int sfp1s0;
  /// Maximum critical temperature of neutron triplet superfluid
  double n3_tc;
  /// Fermi momentum at peak for neutron triplet superfluid
  double n3_kf;
  /// Fermi momentum width parameter for neutron triplet superfluid
  double n3_dk;
  /// Maximum critical temperature of proton singlet superfluid
  double p1_tc;
  /// Fermi momentum at peak for proton singlet superfluid
  double p1_kf;
  /// Fermi momentum width parameter for proton singlet superfluid
  double p1_dk;
  /// Maximum critical temperature of neutron singlet superfluid
  double n1_tc;
  /// Fermi momentum at peak for neutron singlet superfluid
  double n1_kf;
  /// Fermi momentum width parameter for neutron singlet superfluid
  double n1_dk;
  //@}

  /** \brief The NSCool debug parameter
   */
  int nscool_debug;

  /** \brief Fix the direct Urca process at a specified density
   */
  double fix_durca;
  /** \brief Direct Urca modulation parameter
   */
  double alpha_durca;
  /** \brief Direct Urca modulation parameter
   */
  double beta_durca;
  
  /// \name Main cooling curve output
  //@{
  /** \brief Time
   */
  std::vector<double> v_time;

  /** \brief Temperature
   */
  std::vector<double> v_tptr;

  /** \brief Photon luminosity
   */
  std::vector<double> v_lphot;

  /** \brief Neutrino luminosity
   */
  std::vector<double> v_lneut;

  /** \brief Heating
   */
  std::vector<double> v_lheat;
  //@}

  /// \name Data for the pair brehmsstrahlung rate
  //@{
  double pb_n2[56], pb_logt[56], pb_nalpha[56];
  //@}

  /// \name Data for superfluid suppression factor
  //@{
  double sf_lgtau1[35], sf_lgtau2[35], sf_lgr[1225], sf_lgr2[1225];
  //@}

  /// \name Hydrostatic input
  //@{
  /** \brief Crust EOS table

      This table should include (at least) the following columns:
      - <tt>"rho"</tt>: the energy density in units of 
      \f$ \mathrm{g}/\mathrm{cm}^3 \f$
      - <tt>"P"</tt>: the pressure in units of 
      \f$ \mathrm{dyne}/\mathrm{cm}^2 \f$
      - <tt>"n"</tt>: baryon number density in units of
      \f$ \mathrm{fm}^{-3} \f$
      - <tt>"A_cell"</tt>: Total number of nucleons in the W-S cell
      - <tt>"A_ion"</tt>: Total number of nucleons inside the nucleus
      - <tt>"Z"</tt>: Total number of protons inside the nucleus

      The table must be ordered so that the first row has the highest
      density and pressure and the last row has the smallest density
      and pressure. The table is limited to a maximum number of 500
      lines.

      This table is copied to the fortran arrays by \ref
      crust_comp() and \ref crust_eos() .
      
      The default crust is in <tt>crust_HZD_NV.o2</tt> .
  */
  o2scl::table<> nscool_crust;

  /** \brief Core EOS table

      This table should include (at least) the following columns:
      - <tt>"Rho"</tt>: the energy density in units of 
      \f$ \mathrm{g}/\mathrm{cm}^3 \f$
      - <tt>"nbar"</tt>: baryon number density in units of
      \f$ \mathrm{fm}^{-3} \f$
      - <tt>"Ye"</tt>: Number of electrons per baryon
      - <tt>"Ymu"</tt>: Number of muons per baryon
      - <tt>"Yn"</tt>: Number of neutrons per baryon
      - <tt>"Yp"</tt>: Number of protons per baryon
      - <tt>"Yla"</tt>: Number of Lambda hyperons per baryon
      - <tt>"Ysm"</tt>: Number of Sigma minus hyperons per baryon
      - <tt>"Ys0"</tt>: Number of Sigma zero hyperons per baryon
      - <tt>"Ysp"</tt>: Number of Sigma plus hyperons per baryon
      - <tt>"mstp"</tt>: Proton reduced effective mass
      - <tt>"mstn"</tt>: Neutron reduced effective mass
      - <tt>"mstla"</tt>: Lambda hyperon reduced effective mass
      - <tt>"mstsm"</tt>: Sigma minus hyperon reduced effective mass
      - <tt>"msts0"</tt>: Sigma zero hyperon reduced effective mass
      - <tt>"mstsp"</tt>: Sigma plus hyperon reduced effective mass

      The table must go from high energy densities to low energy
      densities. This table is copied to the fortran arrays by \ref
      core_comp() .

      In Dany's original EOS files, the core EOS has a second colum
      for the pressure in units of \f$ \mathrm{dyne}/\mathrm{cm}^2 \f$
      which is not used by the code.
      
      Hyperons support is in progress at the moemnt.

      The default core EOS is in <tt>core_APR.o2</tt> .
  */
  o2scl::table<> nscool_core;

  /** \brief Stellar profile table

      This table should include (at least) the following columns:
      - <tt>"radius"</tt>: the radial coordinate (in \f$ \mathrm{m} \f$ )
      - <tt>"n_baryon"</tt>: baryon number density in units of
      \f$ \mathrm{fm}^{-3} \f$
      - <tt>"density"</tt>: the energy density in units of 
      \f$ \mathrm{g}/\mathrm{cm}^3 \f$
      - <tt>"pressure"</tt>: the energy density in units of 
      \f$ \mathrm{dyne}/\mathrm{cm}^2 \f$
      - <tt>"emass"</tt>: the enclosed gravitational mass in units of
      \f$ \mathrm{M}_{\odot} \f$
      - <tt>"phi"</tt>: the gravitational potential

      The table must be ordered from lower radii (the core) to larger
      radii (the surface). This table is used in \ref star_struct() .
      
      Dany's original tables had the first column labeled "step"
      and a final column labeled "bmass" that were ignored by
      the code.

      The default stellar profile is in <tt>tov_APR_14.o2</tt> .
  */
  o2scl::table<> nscool_tov;
  //@}

  /** \brief Output temperature and luminosity profiles
   */
  o2scl::table3d tl_prof;
  
  /** \brief Times to output
      
      This defaults to a 25-point grid which gives 2 time points every
      decade at late times.
  */
  std::vector<double> time_print;

  /** \brief Initial temperature at drip density relative to 
      initial core temperature (default 0.8)
  */
  double T_fact_drip;

  /** \brief Initial temperature at surface relative to 
      initial core temperature (default 0.5)
  */
  double T_fact_surf;

  /** \brief Flag for profile output (default 0.0)
   */
  double ptemp;

  /** Brief Electron
   */
  o2scl::fermion electron;

  /** Brief Muon
   */
  o2scl::fermion muon;

  /** Brief Fermion thermodynamics
   */
  o2scl::fermion_eff fe; 
  
  nscool_wrap(std::string dir) {
    o2scl_hdf::hdf_file hf;
    std::string name;

    int mpi_rank=0, mpi_size=1;

#ifdef O2SCL_MPI
    // Get MPI rank, etc.
    MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);

    // Ensure that multiple threads aren't writing to the
    // filesystem at the same time
    int tag=0, buffer=0;
    if (mpi_size>1 && mpi_rank>=1) {
      MPI_Recv(&buffer,1,MPI_INT,mpi_rank-1,
	       tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
#endif

    if (false) {
      std::cout << "In nscool_wrap, rank " << mpi_rank
		<< " reading data files." << std::endl;
    }
    
    // Read default crust EOS
    hf.open(dir+"/crust_HZD_NV.o2");
    hdf_input(hf,nscool_crust,name);
    hf.close();
    
    // Read default core EOS
    hf.open(dir+"/core_APR.o2");
    hdf_input(hf,nscool_core,name);
    hf.close();
    
    // Read default stellar profile
    hf.open(dir+"/tov_APR_14.o2");
    hdf_input(hf,nscool_tov,name);
    hf.close();

#ifdef O2SCL_MPI
    if (mpi_size>1 && mpi_rank<mpi_size-1) {
      MPI_Send(&buffer,1,MPI_INT,mpi_rank+1,
	       tag,MPI_COMM_WORLD);
    }
#endif
    
    ptemp=0.0;
    
    // Setup time_print
    time_print.resize(25);
    time_print[0]=1.0e-10;
    time_print[1]=1.0e-4;
    time_print[2]=3.0e-4;
    time_print[3]=1.0e-3;
    time_print[4]=3.0e-3;
    time_print[5]=1.0e-2;
    time_print[6]=3.0e-2;
    time_print[7]=1.0e-1;
    time_print[8]=3.0e-1;
    time_print[9]=1.0e+0;
    time_print[10]=3.0e+0;
    time_print[11]=1.0e+1;
    time_print[12]=3.0e+1;
    time_print[13]=1.0e+2;
    time_print[14]=3.0e+2;
    time_print[15]=1.0e+3;
    time_print[16]=3.0e+3;
    time_print[17]=1.0e+4;
    time_print[18]=3.0e+4;
    time_print[19]=1.0e+5;
    time_print[20]=3.0e+5;
    time_print[21]=1.0e+6;
    time_print[22]=3.0e+6;
    time_print[23]=1.0e+7;
    time_print[24]=3.0e+7;
    T_fact_drip=0.8;
    T_fact_surf=0.5;
    eta=0.0;

    /*
      This is the SFB neutrino singlet gap, the T73 proton singlet
      gap, and neutrino triplet gap "a" from the minimal cooling
      paper.
     */
    sfn1s0=1;
    sfn3p2=101;
    sfp1s0=3;

    nscool_debug=0;

    fix_durca=0.0;
    alpha_durca=1.0e-8;
    beta_durca=1.0;

    main_out_it=20;
    
    // Read data for pair bremsstrahlung
    pair_brem_data(pb_logt,pb_nalpha);
    int pb_n=56;
    double pb_bound1=1.0e30;
    double pb_bound2=1.0e30;
    spline_(pb_logt,pb_nalpha,&pb_n,&pb_bound1,&pb_bound2,pb_n2);
    
    // Data for superfluid suppression
    sf_suppress_data(sf_lgtau1,sf_lgtau2,sf_lgr);
    int sf_n1=35;
    int sf_n2=35;
    spline2_(sf_lgtau1,sf_lgtau2,sf_lgr,&sf_n1,&sf_n2,sf_lgr2);
    
    // Lepton inits
    electron.init(o2scl::o2scl_settings.get_convert_units().convert
		  ("kg","1/fm",o2scl_mks::mass_electron),2.0);
    muon.init(o2scl::o2scl_settings.get_convert_units().convert
	      ("kg","1/fm",o2scl_mks::mass_muon),2.0);
  }

  void P_electron(double ne, double T, double *pre) {
    electron.n=ne;
    fe.calc_density(electron,T);
    *pre=electron.pr;
    return;
  }
  
  /** \brief Load the default star, APR with M=1.4
   */
  void default_star(std::string dir=".") {
    
    o2scl_hdf::hdf_file hf;
    std::string name;
    
    // Read default crust EOS
    hf.open(dir+"/crust_HZD_NV.o2");
    hdf_input(hf,nscool_crust,name);
    hf.close();
    
    // Read default core EOS
    hf.open(dir+"/core_APR.o2");
    hdf_input(hf,nscool_core,name);
    hf.close();
    
    // Read default stellar profile
    hf.open(dir+"/tov_APR_14.o2");
    hdf_input(hf,nscool_tov,name);
    hf.close();
    
    return;
  }
    
  /** \brief Compute the HHJ EOS
   */
  void hhj_eos(double mass) {

    eos_had_hhj hhj;
    o2scl::nstar_cold nc;
    o2scl::convert_units &cu=o2scl::o2scl_settings.get_convert_units();

    // ------------------------------------------------------------
    // Construct the EOS from an nstar_cold object
      
    // Automatically construct beta-equilibrium EOS
    nc.verbose=0;
    nc.def_tov.verbose=0;
    nc.set_eos(hhj);
    nc.err_nonconv=false;
    nc.calc_eos();
    std::shared_ptr<o2scl::table_units<> > nc_eos=nc.get_eos_results();
    nc.def_eos_tov.read_table(*nc_eos,"ed","pr","nb");

    // Prepare to compute effective masses
    o2scl::fermion n(o2scl::o2scl_settings.get_convert_units().convert
	      ("kg","1/fm",o2scl_mks::mass_neutron),2.0);
    o2scl::fermion p(o2scl::o2scl_settings.get_convert_units().convert
	      ("kg","1/fm",o2scl_mks::mass_proton),2.0);
    n.non_interacting=false;
    p.non_interacting=false;
    o2scl::thermo th;
      
    // ------------------------------------------------------------
    // Create a core table of the correct size with the correct units

    // Clear table for new columns
    nscool_core.clear();
    nscool_core.line_of_names("Rho Press nbar Ye Ymu Yn Yp mstp mstn");
    nscool_core.line_of_names("Yla Ysm Ys0 Ysp mstla mstsm msts0 mstsp");
      
    // Baryon density grid
    double nb_min=0.09;
    double nb_max=nc_eos->max("nb");
    double dnb=(nb_max-nb_min)/178.0;
      
    // Fill core table
    for(double nb=nb_max;nb>nb_min-dnb/10.0;nb-=dnb) {
      n.n=nc_eos->interp("nb",nb,"nn");
      p.n=nc_eos->interp("nb",nb,"np");
      hhj.calc_e(n,p,th);
      double line[17]={cu.convert("1/fm^4","g/cm^3",
				 nc_eos->interp("nb",nb,"ed")),
		      cu.convert("1/fm^4","dyne/cm^2",
				 nc_eos->interp("nb",nb,"pr")),nb,
		      nc_eos->interp("nb",nb,"ne")/nb,
		      0.0,
		      nc_eos->interp("nb",nb,"nn")/nb,
		      nc_eos->interp("nb",nb,"np")/nb,
		      n.ms/n.m,p.ms/p.m,0.0,0.0,0.0,0.0,
		      1.0,1.0,1.0,1.0};
      nscool_core.line_of_data(17,line);
    }
      
    // ------------------------------------------------------------
    // Use the nstar_cold object to construct the profile

    // Construct neutron star profile
    nc.def_tov.calc_gpot=true;
    nc.calc_nstar();
    std::shared_ptr<o2scl::table_units<> > nc_prof=nc.get_tov_results();
    std::cout << "M_max: " << nc_prof->max("gm") << std::endl;
    std::cout << "R(M_max): "
	 << nc_prof->get("r",nc_prof->lookup("gm",nc_prof->max("gm")))
	 << std::endl;

    double mass2=mass;
    if (mass2<0.0) mass2=mass+nc_prof->max("gm");
    std::cout << "Set mass to " << mass2 << std::endl;
    nc.fixed(mass2);
      
    // ------------------------------------------------------------
    // Create a table with the right size and the right units

    // Clear table for new columns
    nscool_tov.clear();
    nscool_tov.line_of_names(((std::string)"step radius n_baryon density ")+
			     "pressure emass phi bmass");
      
    // Radial grid
    double r_max=nc_prof->max("r");
      
    // Fill TOV table
    size_t ix=0;
    for(size_t ix=0;ix<153;ix++) {
      // A grid which focuses most of the points on the outer part
      // of the star
      double r=r_max*2.0*(1.0-pow(2.0,-((double)ix)/152.0));
      double line[8]={((double)ix),r*1.0e3,nc_prof->interp("r",r,"nb"),
		      cu.convert("1/fm^4","g/cm^3",
				 nc_prof->interp("r",r,"ed")),
		      cu.convert("1/fm^4","dyne/cm^2",
				 nc_prof->interp("r",r,"pr")),
		      nc_prof->interp("r",r,"gm"),
		      nc_prof->interp("r",r,"gp"),
		      nc_prof->interp("r",r,"bm")};
      nscool_tov.line_of_data(8,line);
    }

    return;
  }
  
  /// \name Functions called by the Fortran cooling code
  //@{
  /** \brief Specify initial temperature profile

      Called in <tt>NSCool.f</tt>.

      Uses \ref T_fact_surf and \ref T_fact_drip .
  */
  virtual void tptr_init(int ifteff, double tempini, double ephi_surf,
			 double ephi_drip, double ephi_core, double &tsurface,
			 double &tdrip, double &tcore, double &tb_acc0) {

    if (ifteff!=15) {
      if (tempini>0.0) {
	tsurface=T_fact_surf*ephi_surf*tempini;
	tdrip=T_fact_drip*ephi_drip*tempini;
	tcore=1.0*ephi_core*tempini;
      } else {
	tsurface=1.0e9;
	tdrip=2.0e10;
	tcore=1.0e11;
      }
    } else {
      tb_acc0*=ephi_surf;
      tsurface=tb_acc0;
      tdrip=tb_acc0;
      tcore=tb_acc0;
    }
    return;
  }
  
  /** \brief Compute effective temperature in envelope from
      boundary temperature

      This is a replacement for code originally in <tt>boundary.f</tt>
      and called in <tt>NSCool.f</tt>. The parameter \c Tb and the
      return value are both local temperatures (not redshifted).

      \verbatim embed:rst
      This function is based on Appendix A, section 3 in [Potekhin97it]_.
      \endverbatim
      
      The parameters <tt>bfield, Z, A, compactness, ifteff, istep, time,
      Ts1, Ts2, Rho, debug</tt> are currently unused (but might
      be used in future versions?). 
  */
  virtual double Teff(double Tb, int ifteff, double eta_arg, double bfield,
		      int istep, double time, double Ts1, double Ts2,
		      double Z, double A, double Rho, int debug,
		      double gs14, double compactness) {

    // The internal temperature in units of 10^9 K
    double Tb9=Tb/1.0e9;

    // T_{*} in units of 10^{6} K
    double Ts=sqrt(7.0e0*Tb9*sqrt(gs14));

    // zeta
    double z=Tb9-Ts/1.0e3;

    // The effective temperature to the fourth power for an Iron
    // envelope in units of 10^{6} K
    double t4_iron=gs14*(pow(7.0e0*z,2.25)+pow(z/3.0e0,1.25));

    // The effective temperature to the fourth power for a
    // fully accreted envelope in units of 10^{6} K
    double t4_wacc=gs14*pow(18.1e0*Tb9,2.42);

    // For a partially accreted envelope
    double t4_acc;
    if (eta_arg>1.0e-30) {
      double a=(1.2e0+pow(5.3e-6/eta_arg,0.38))*pow(Tb9,5.0/3.0);
      t4_acc=(a*t4_iron+t4_wacc)/(a+1.0e0);
    } else {
      t4_acc=t4_iron;
    }

    // Return the final effective temperature in Kelvin
    return pow(t4_acc,0.25)*1.0e6;
  }
	
  /** \brief Function for printing out iteration progress

      This function is only called if <tt>ptemp&gt;=1.0</tt> as
      specified in \ref cool_param() (the default). (This is slightly
      different than the original code which only outputs if
      <tt>ptemp=1.0</tt>.) This function's wrapper is called in
      <tt>NSCool.f</tt>. If <tt>ptemp&gt;=1.0</tt> then the
      temperature information is output to \ref tl_prof. If also
      <tt>ptemp&gt;=2.0</tt> then the temperature information is
      output to <tt>std::cout</tt>.
      
      The parameter \c time is the time in years, \c t_effective is
      the effective temperature at \f$ \infty \f$ in K (includes \f$
      \exp(\phi) \f$), \c imax is the maximum array index (always
      odd). The parameters \c w1 and \c w2 are the weighting factors
      for interpolation. The arrays \c otemp and \c temp are the
      previous and next temperature arrays, and \c olum and \c lum are
      the previous and new temperatures (both local values without any
      factors of \f$ \exp \phi \f$). The arrays \c rad and \c rrho are
      the radius (in cm) and mass(?) density (in \f$
      \mathrm{g}/\mathrm{cm}^3 \f$). The array \c ephi is \f$ \exp
      \phi \f$ and \c e2phi is \f$ \exp ( 2 \phi ) \f$ . The array \c
      dvol is the physical volume between \c i and \c i-1. (All of
      these arrays are zero-indexed in the original FORTRAN code).
      
  */
  virtual void print_temp(int istep, int itprint,
			  double time, double t_effective, int imax,
			  double w1, double w2, double *otemp, double *temp,
			  double *olum, double *lum, double *rad,
			  double *rrho, double *ephi, double *dvol,
			  double *e2phi, double *qnu, double *qeebrem,
			  double *qnpb, double *qplasma,
			  double *qsynch, double *qbubble, double *qpair,
			  double *qphoto, double *qbrem_nn,
			  double *qmurca_nucl, double *qbrem_nucl,
			  double *qmurca_hyp, double *qbrem_hyp,
			  double *qdurca_np, double *qdurca_lap,
			  double *qdurca_smn, double *qdurca_smla,
			  double *qdurca_sms0, double *qfast,
			  double *qdurca_q, double *qmurca_q, 
			  double *qpbf_n1s0, double *qpbf_p1s0,
			  double *qpbf_n3p2, double *qpbf_q) {

    // T_eff, as passed by the Fortran code in t_effective,
    // is not currently stored in tl_prof.o2. The variable
    // istep is also not currently stored anywhere.

    if (fabs(time_print[itprint]-time)/fabs(time)>1.0e-6) {
      O2SCL_ERR("Temperature print sanity check.",o2scl::exc_einval);
    }
    
    if (itprint==0) {
      tl_prof.clear();
      std::vector<double> r_grid;
      for(int i=imax;i>=1;i-=2) {
	r_grid.push_back(rad[imax+1-i]/1.0e5);
      }
      tl_prof.set_xy("r",r_grid.size(),r_grid,
		     "t",time_print.size(),time_print);
      tl_prof.line_of_names(((std::string)"rho ephi vol Tinf Linf qnu ")+
			    "qeebrem qnpb qplasma qsynch qbubble qpair "+
			    "qphoto qbrem_nn qmurca_nucl qbrem_nucl "+
			    "qmurca_hyp qbrem_hyp qdurca_np qdurca_lap "+
			    "qdurca_smn qdurca_smla qdurca_sms0 qfast "+
			    "qdurca_q qmurca_q qpbf_n1s0 qpbf_p1s0 "+
			    "qpbf_n3p2 qpbf_q qmax "+
			    "cv cv_n cv_p cv_e cv_m cv_la "+
			    "cv_sm cv_s0 cv_sp cv_q");
      tl_prof.add_constant("it_last",itprint);
    }
    tl_prof.set_constant("it_last",itprint);

    if (ptemp>=2.0) {
      std::cout << "Time: " << time << " years, T_eff: "
		<< t_effective << " K" << std::endl;
    }
    
    for(int i=imax;i>=1;i-=2) {
      int io2=i/2;
      double logtemp=w1*log(otemp[i])+w2*log(temp[i]);
      double temperature=exp(logtemp);
      double lumino=0.0;
      if (i!=1) {
	double loglum=w1*log(fabs(olum[i-1]))+w2*log(fabs(lum[i-1]));
	lumino=exp(loglum);
	if (lum[i-1]<0.0) lumino*=-1.0;
      }
      tl_prof.set(io2,itprint,"rho",rrho[i]);
      tl_prof.set(io2,itprint,"ephi",ephi[i]);
      tl_prof.set(io2,itprint,"vol",dvol[i]+dvol[i+1]);
      tl_prof.set(io2,itprint,"Tinf",temperature/ephi[i]);
      tl_prof.set(io2,itprint,"Linf",lumino/e2phi[i-1]);
      tl_prof.set(io2,itprint,"qnu",qnu[i]);
      tl_prof.set(io2,itprint,"qeebrem",qeebrem[i]);
      tl_prof.set(io2,itprint,"qnpb",qnpb[i]);
      tl_prof.set(io2,itprint,"qplasma",qplasma[i]);
      tl_prof.set(io2,itprint,"qsynch",qsynch[i]);
      tl_prof.set(io2,itprint,"qbubble",qbubble[i]);
      tl_prof.set(io2,itprint,"qpair",qpair[i]);
      tl_prof.set(io2,itprint,"qphoto",qphoto[i]);
      tl_prof.set(io2,itprint,"qbrem_nn",qbrem_nn[i]);
      tl_prof.set(io2,itprint,"qmurca_nucl",qmurca_nucl[i]);
      tl_prof.set(io2,itprint,"qbrem_nucl",qbrem_nucl[i]);
      tl_prof.set(io2,itprint,"qmurca_hyp",qmurca_hyp[i]);
      tl_prof.set(io2,itprint,"qbrem_hyp",qbrem_hyp[i]);
      tl_prof.set(io2,itprint,"qdurca_np",qdurca_np[i]);
      tl_prof.set(io2,itprint,"qdurca_lap",qdurca_lap[i]);
      tl_prof.set(io2,itprint,"qdurca_smn",qdurca_smn[i]);
      tl_prof.set(io2,itprint,"qdurca_smla",qdurca_smla[i]);
      tl_prof.set(io2,itprint,"qdurca_sms0",qdurca_sms0[i]);
      tl_prof.set(io2,itprint,"qfast",qfast[i]);
      tl_prof.set(io2,itprint,"qdurca_q",qdurca_q[i]);
      tl_prof.set(io2,itprint,"qmurca_q",qmurca_q[i]);
      tl_prof.set(io2,itprint,"qpbf_n1s0",qpbf_n1s0[i]);
      tl_prof.set(io2,itprint,"qpbf_p1s0",qpbf_p1s0[i]);
      tl_prof.set(io2,itprint,"qpbf_n3p2",qpbf_n3p2[i]);
      tl_prof.set(io2,itprint,"qpbf_q",qpbf_q[i]);

      {
	std::vector<double> qvec;
	qvec.push_back(fabs(qeebrem[i]));
	qvec.push_back(fabs(qnpb[i]));
	qvec.push_back(fabs(qplasma[i]));
	qvec.push_back(fabs(qsynch[i]));
	qvec.push_back(fabs(qbubble[i]));
	qvec.push_back(fabs(qpair[i]));
	qvec.push_back(fabs(qphoto[i]));
	qvec.push_back(fabs(qbrem_nn[i]));
	qvec.push_back(fabs(qmurca_nucl[i]));
	qvec.push_back(fabs(qbrem_nucl[i]));
	qvec.push_back(fabs(qmurca_hyp[i]));
	qvec.push_back(fabs(qbrem_hyp[i]));
	qvec.push_back(fabs(qdurca_np[i]));
	qvec.push_back(fabs(qdurca_lap[i]));
	qvec.push_back(fabs(qdurca_smn[i]));
	qvec.push_back(fabs(qdurca_smla[i]));
	qvec.push_back(fabs(qdurca_sms0[i]));
	qvec.push_back(fabs(qfast[i]));
	qvec.push_back(fabs(qdurca_q[i]));
	qvec.push_back(fabs(qmurca_q[i]));
	qvec.push_back(fabs(qpbf_n1s0[i]));
	qvec.push_back(fabs(qpbf_p1s0[i]));
	qvec.push_back(fabs(qpbf_n3p2[i]));
	qvec.push_back(fabs(qpbf_q[i]));
	o2scl::permutation order(qvec.size());
	o2scl::vector_sort_index(qvec.size(),qvec,order);
	if (qvec[order[qvec.size()-1]]>0.0) {
	  tl_prof.set(io2,itprint,"qmax",order[qvec.size()-1]+1);
	} else {
	  tl_prof.set(io2,itprint,"qmax",-order[qvec.size()-1]+1);
	}
      }

      {
	double total=qeebrem[i]+qnpb[i]+qplasma[i]+qsynch[i]+
	  qbubble[i]+qpair[i]+qphoto[i]+qbrem_nn[i]+
	  qmurca_nucl[i]+qbrem_nucl[i]+qmurca_hyp[i]+qbrem_hyp[i]+
	  qdurca_np[i]+qdurca_lap[i]+qdurca_smn[i]+qdurca_smla[i]+
	  qdurca_sms0[i]+qfast[i]+qdurca_q[i]+qmurca_q[i]+
	  qpbf_n1s0[i]+qpbf_p1s0[i]+qpbf_n3p2[i]+qpbf_q[i];
	if (fabs(qnu[i]-total)/fabs(qnu[i])>1.0e-6) {
	  std::cout << "Emissivity problem " << qnu[i] << " " << total
		    << std::endl;
	  exit(-1);
	}
      }
      
      if (ptemp>=2.0) {
	std::cout << i << " "
		  << tl_prof.get_grid_x(io2) << " "
		  << rrho[i] << " "
		  << ephi[i] << " "
		  << dvol[i]+dvol[i+1] << " " 
		  << temperature/ephi[i] << " "
		  << lumino/e2phi[i-1] << " "
		  << qnu[i] << " "
		  << qeebrem[i] << " "
		  << qnpb[i] << " "
		  << qplasma[i] << " "
		  << qsynch[i] << " "
		  << qbubble[i] << " "
		  << qpair[i] << " "
		  << qphoto[i] << " "
		  << qbrem_nn[i] << " "
		  << qmurca_nucl[i] << " "
		  << qbrem_nucl[i] << " "
		  << qmurca_hyp[i] << " "
		  << qbrem_hyp[i] << " "
		  << qdurca_np[i] << " "
		  << qdurca_lap[i] << " "
		  << qdurca_smn[i] << " "
		  << qdurca_smla[i] << " "
		  << qdurca_sms0[i] << " "
		  << qfast[i] << " "
		  << qdurca_q[i] << " "
		  << qmurca_q[i] << " "
		  << qpbf_n1s0[i] << " "
		  << qpbf_p1s0[i] << " "
		  << qpbf_n3p2[i] << " "
		  << qpbf_q[i] << std::endl;
      }
    }
    if (ptemp>=2.0) {
      std::cout << std::endl;
    }
    
    return;
  }

  /** \brief Output specific heats
   */
  virtual void print_cv(int itprint, int imax, double *cv,
			double *cv_n, double *cv_p,
			double *cv_e, double *cv_m, double *cv_la, 
			double *cv_sm, double *cv_s0, double *cv_sp,
			double *cv_q) {
    
    if (ptemp>=1.0) {
      for(int i=imax;i>=1;i-=2) {
	int io2=i/2;
	tl_prof.set(io2,itprint,"cv",cv[i]);
	tl_prof.set(io2,itprint,"cv_n",cv_n[i]);
	tl_prof.set(io2,itprint,"cv_p",cv_p[i]);
	tl_prof.set(io2,itprint,"cv_e",cv_e[i]);
	tl_prof.set(io2,itprint,"cv_m",cv_m[i]);
	tl_prof.set(io2,itprint,"cv_la",cv_la[i]);
	tl_prof.set(io2,itprint,"cv_sm",cv_sm[i]);
	tl_prof.set(io2,itprint,"cv_s0",cv_s0[i]);
	tl_prof.set(io2,itprint,"cv_sp",cv_sp[i]);
	tl_prof.set(io2,itprint,"cv_q",cv_q[i]);
	if (ptemp>=2.0) {
	  std::cout << i << " " << cv_n[i] << " " << cv_p[i] << " "
		    << cv_e[i] << " " << cv_m[i] << " "
		    << cv_la[i] << " " << cv_sm[i] << " "
		    << cv_s0[i] << " " << cv_sp[i] << " "
		    << cv_q[i] << std::endl;
	}
      }
    }
    return;
  }
  
  /** \brief Set various numerical parameters

      The parameter \c time0 is the initial time (default 0) \c
      timemax is the maximum time in years (default \f$ 10^{10} \f$).

      \todo More docs here.

      \note If istep is greater than istepmax, the nscool code just
      silently exits without warning.

      This function's wrapper is called in <tt>NSCool.ff</tt>.
  */
  virtual void num_param(double &time0, double &timemax, int &istepmax,
			 int &itrial_max, int &itrial_opt, double &tcut,
			 double &dtime, double &dtlimit, double &scale_dt0,
			 double &scale_dt1, double &repeat, int &istart,
			 double &mratt, double &mratl, double &mrats,
			 double &tvar, double &svar, double &tcon) {
    time0=0.0;
    timemax=2.0e10;
    istepmax=1000000;
    itrial_max=20;
    itrial_opt=12;
    tcut=2.0;
    dtime=1.0e-12;
    dtlimit=3.15e15;
    scale_dt0=1.2;
    scale_dt1=1.5;
    repeat=0.2;
    istart=2;
    mratt=1.0e-12;
    mratl=1.0e-10;
    mrats=1.0e-10;
    tvar=1.20;
    svar=1.05;
    tcon=1.0e12;
    return;
  }

  /** \brief Parameters at the boundary

      The parameter \c ifteff is 1 for the Te-Tb from \ref
      Gudmundsson83 and \c ifteff is 2 for that from \ref Nomoto87.
      These are old, the best is \c ifteff 3 from \ref Potekhin97. If
      \c ifteff is 0, it reads the Te-Tb from a file (probably not
      supported in this version). The value of \c eta determines the
      amount of light elements.

      The value eta is defined by 
      \f[
      \eta = g_{s14}^2 \Delta M / M = P_{\mathrm{light}} / P_0
      \f]
      where \f$ P_0 \equiv 1.193 \times
      10^{34}~\mathrm{dyne}/\mathrm{cm}^2 \f$ where 
      \f$ P_{\mathrm{light}} \f$ is the pressure at the 
      bottom of the light element layer. 

      If \c ifteff is 15, this simulates an accreting neutron
      star where \f$ T_b \f$ is held constant. Then \f$ 
      T_b \f$ is specified in \c tb_acc0 .

      [Dany:] Notice that light elements cannot be present at too high
      densities (e.g., C will burn by pycnonuclear reactions at about
      1010 g cm-3). So there is a maximum value that \f$ \eta \f$ can
      reach. The formula that Potekhin et al. give saturates when \f$
      \eta \f$ grows: thus a value as \f$ \eta=1 \f$ will give the
      maximum possible effect of a light element envelope (even if it
      is physically a wildly unrealistic high value).

      This function's wrapper is called in <tt>NSCool.f</tt>.
  */
  virtual void bound_param(int &ifteff, double &eta_arg,
			   double &mag_coeff, double &tb_acc0) {
    ifteff=3;
    eta_arg=eta;
    mag_coeff=3.0;
    tb_acc0=0.0;
    return;
  }
  
  /** \brief Specify several parameters

      The parameter \c pscreen controls screen output, \c debug is a
      generic debug variable. The parameters \c emnco, \c emncr, and
      \c emp are flags for the automatic computation of the
      effective masses in <tt>precool.f</tt> .

      The parameters \c pteff, \c ptemp, and \c pstar originally
      controlled whether or not iteration information was output to
      files. This file output is now replaced by \ref print_temp()
      and these parameters will be deprecated. 

      This function's wrapper is called in <tt>NSCool.f</tt>.
  */
  virtual void cool_param(int &pscreen, double &debug, int &istep_debug,
			  double &pteff, double &ptemp_arg, double &pstar,
			  int &idump1, int &idump2, int &idump3,
			  double &tempmin, double &tempini,
			  int &icvel_nodeg, double &emnco, double &emncr,
			  double &emp, double &p0, int &itpmax,
			  double *tprint) {

    pscreen=1;
    debug=((double)(nscool_debug));
    //debug=0.5;
    istep_debug=0;
    pteff=0.0;
    ptemp_arg=1.0;
    pstar=0.0;
    idump1=1;
    idump2=111;
    idump3=421;
    tempmin=1.0e4;
    tempini=1.0e10;
    // If this is true, then the function 'cvelec()' is used
    // to handle the electron specific heat rather than
    // the simple expression for degenerate electrons
    icvel_nodeg=0;
    emnco=5.0;
    emncr=5.0;
    emp=3.0;
    p0=0.1;

    // The array 'tprint' the cooling code is actually zero-indexed,
    // but Dany's code ignores the first value.
    itpmax=time_print.size();
    for(int i=1;i<=itpmax;i++) {
      tprint[i]=time_print[i-1];
    }
    return;
  }
  
  /** \brief Specify core composition

      This function provides the core composition by
      copying the information in the \ref nscool_core table
      to the Fortran arrays.

      This function's wrapper is called in <tt>precool.f</tt>.
  */
  virtual void core_comp
    (double *rho_t, double *nbar_t, double *yelect_t, double *ymuon_t,
     double *yneutr_t, double *yprot_t, double *ylambda_t,
     double *ysminus_t, double *yszero_t, double *ysplus_t,
     double *mstp_t, double *mstn_t, double *mstla_t, double *mstsm_t,
     double *msts0_t, double *mstsp_t, int *ix) {

    (*ix)=nscool_core.get_nlines();
    for(size_t i=0;i<nscool_core.get_nlines();i++) {
      // rho is energy density
      rho_t[i]=nscool_core.get("Rho",i);
      nbar_t[i]=nscool_core.get("nbar",i);
      yelect_t[i]=nscool_core.get("Ye",i);
      ymuon_t[i]=nscool_core.get("Ymu",i);
      yneutr_t[i]=nscool_core.get("Yn",i);
      yprot_t[i]=nscool_core.get("Yp",i);
      ylambda_t[i]=nscool_core.get("Yla",i);
      ysminus_t[i]=nscool_core.get("Ysm",i);
      yszero_t[i]=nscool_core.get("Ys0",i);
      ysplus_t[i]=nscool_core.get("Ysp",i);
      
      double diff=yprot_t[i]+ysplus_t[i]-ysminus_t[i]-
	yelect_t[i]-ymuon_t[i];
      if (fabs(diff)>1.0e-5) {
	std::cerr << "Charge problem in nscool_wrap::core_comp()."
		  << std::endl;
	std::cout << i << " " << yprot_t[i] << " " << yelect_t[i] << " "
		  << ymuon_t[i] << " " << ysplus_t[i] << " "
		  << ysminus_t[i] << " " [i] << std::endl;
	(*ix)=0;
	return;
      }
      if (fabs(yprot_t[i]+yneutr_t[i]+ylambda_t[i]+ysminus_t[i]+
	       yszero_t[i]+ysplus_t[i]-1.0)>1.0e-5) {
	std::cerr << "Baryon problem in core_comp." << std::endl;
	std::cout << i << " " << yneutr_t[i] << " "
		  << yprot_t[i] << " " << ylambda_t[i] << " "
		  << ysminus_t[i] << " " << yszero_t[i] << " "
		  << ysplus_t[i] << std::endl;
	(*ix)=0;
	return;
      }
      
      // These are the reduced effective masses, i.e. m^{*}/m
      mstp_t[i]=nscool_core.get("mstp",i);
      mstn_t[i]=nscool_core.get("mstn",i);
      mstla_t[i]=nscool_core.get("mstla",i);
      mstsm_t[i]=nscool_core.get("mstsm",i);
      msts0_t[i]=nscool_core.get("msts0",i);
      mstsp_t[i]=nscool_core.get("mstsp",i);
    }
    
    if (rho_t[1]>rho_t[0]) {
      O2SCL_ERR("Core composition table should be decreasing",
		o2scl::exc_einval);
    }
    
    return;
  };

  /** \brief Fix the settings for direct Urca

      AWS: alpha is the broadening parameter and beta is the
      fractional decrease of the direct Urca threshold
   */
  virtual void urca_settings(double &durca, double &a_durca,
			     double &b_durca) {
    durca=fix_durca;
    a_durca=alpha_durca;
    b_durca=beta_durca;
    return;
  }
  
  /** \brief Specify crust composition
      
      This function provides the crust composition by
      copying the information in the \ref nscool_crust table
      to the Fortran arrays.

      This function's wrapper is called in <tt>precool.f</tt>.
  */
  virtual void crust_comp
    (double *Z_ion_t, double *A_ion_t, double *A_cell_t,
     double *bar_t, double *pres_t, double *rho_t,
     int *jmax) {
    
    *jmax=((int)(nscool_crust.get_nlines()));
    if (nscool_crust.get_nlines()>=500) {
      O2SCL_ERR("Crust table too large in nscool_wrap::crust_comp().",
		o2scl::exc_einval);
    }
    for(size_t i=0;i<nscool_crust.get_nlines();i++) {
      size_t i2=(*jmax)-1-i;
      rho_t[i]=nscool_crust.get("rho",i2);
      pres_t[i]=nscool_crust.get("P",i2);
      bar_t[i]=nscool_crust.get("n",i2);
      A_cell_t[i]=nscool_crust.get("A_cell",i2);
      A_ion_t[i]=nscool_crust.get("A_ion",i2);
      Z_ion_t[i]=nscool_crust.get("Z",i2);
    }
    if (rho_t[1]<rho_t[0]) {
      O2SCL_ERR("Crust composition table should be decreasing",
		o2scl::exc_einval);
    }
    
    return;
  };

  /** \brief Specify crust EOS
      
      This function provides the crust EOS by
      copying the information in the \ref nscool_crust table
      to the Fortran arrays.

      This function's wrapper is called in <tt>precool.f</tt>.
  */
  virtual void crust_eos(double *rho2, double *pres2, int *idata) {

    *idata=((int)(nscool_crust.get_nlines()));
    for(size_t i=0;i<nscool_crust.get_nlines();i++) {
      rho2[i]=nscool_crust.get("rho",(*idata)-1-i);
      pres2[i]=nscool_crust.get("P",(*idata)-1-i);
    }

    return;
  };

  /** \brief Data for superfluid suppression function

      This function's wrapper is called in <tt>neutrino_core.f</tt>.
  */
  virtual void sf_suppress_data(double *lgtau1, double *lgtau2,
				double *lgr) {
    
#include "sf_suppression.h"
    
    return;
  }
  
  /** \brief Data for neutrino pair bremsstrahlung

      This function's wrapper is called in <tt>neutrino_crust.f</tt>.
  */
  virtual void pair_brem_data(double *logt, double *nalpha) {

    nalpha[0]=2.119400e+02;
    nalpha[1]=2.119314e+02;
    nalpha[2]=2.118338e+02;
    nalpha[3]=2.116714e+02;
    nalpha[4]=2.114446e+02;
    nalpha[5]=2.111548e+02;
    nalpha[6]=2.108042e+02;
    nalpha[7]=2.103965e+02;
    nalpha[8]=2.099369e+02;
    nalpha[9]=2.094313e+02;
    nalpha[10]=2.088862e+02;
    nalpha[11]=2.025300e+02;
    nalpha[12]=1.965045e+02;
    nalpha[13]=1.913841e+02;
    nalpha[14]=1.870169e+02;
    nalpha[15]=1.832339e+02;
    nalpha[16]=1.799095e+02;
    nalpha[17]=1.769527e+02;
    nalpha[18]=1.742953e+02;
    nalpha[19]=1.718859e+02;
    nalpha[20]=1.555852e+02;
    nalpha[21]=1.458835e+02;
    nalpha[22]=1.389667e+02;
    nalpha[23]=1.335816e+02;
    nalpha[24]=1.291637e+02;
    nalpha[25]=1.254108e+02;
    nalpha[26]=1.221419e+02;
    nalpha[27]=1.192400e+02;
    nalpha[28]=1.166251e+02;
    nalpha[29]=9.846525e+01;
    nalpha[30]=8.600176e+01;
    nalpha[31]=7.553085e+01;
    nalpha[32]=6.634075e+01;
    nalpha[33]=5.830211e+01;
    nalpha[34]=5.135016e+01;
    nalpha[35]=4.539033e+01;
    nalpha[36]=4.030437e+01;
    nalpha[37]=3.596899e+01;
    nalpha[38]=1.500019e+01;
    nalpha[39]=8.665605e+00;
    nalpha[40]=5.942522e+00;
    nalpha[41]=4.494431e+00;
    nalpha[42]=3.613857e+00;
    nalpha[43]=3.027656e+00;
    nalpha[44]=2.611401e+00;
    nalpha[45]=2.301280e+00;
    nalpha[46]=2.061534e+00;
    nalpha[47]=1.063250e+00;
    nalpha[48]=7.494363e-01;
    nalpha[49]=5.912537e-01;
    nalpha[50]=4.941738e-01;
    nalpha[51]=4.277401e-01;
    nalpha[52]=3.790197e-01;
    nalpha[53]=3.415363e-01;
    nalpha[54]=3.116674e-01;
    nalpha[55]=0.000000e-00;

    logt[0]=log10(1.000000e-10)+9.0;
    logt[1]=log10(1.160093e-04)+9.0;
    logt[2]=log10(2.320186e-04)+9.0;
    logt[3]=log10(3.480278e-04)+9.0;
    logt[4]=log10(4.640371e-04)+9.0;
    logt[5]=log10(5.800464e-04)+9.0;
    logt[6]=log10(6.960557e-04)+9.0;
    logt[7]=log10(8.120650e-04)+9.0;
    logt[8]=log10(9.280742e-04)+9.0;
    logt[9]=log10(1.044084e-03)+9.0;
    logt[10]=log10(1.160093e-03)+9.0;
    logt[11]=log10(2.320186e-03)+9.0;
    logt[12]=log10(3.480278e-03)+9.0;
    logt[13]=log10(4.640371e-03)+9.0;
    logt[14]=log10(5.800464e-03)+9.0;
    logt[15]=log10(6.960557e-03)+9.0;
    logt[16]=log10(8.120650e-03)+9.0;
    logt[17]=log10(9.280742e-03)+9.0;
    logt[18]=log10(1.044084e-02)+9.0;
    logt[19]=log10(1.160093e-02)+9.0;
    logt[20]=log10(2.320186e-02)+9.0;
    logt[21]=log10(3.480278e-02)+9.0;
    logt[22]=log10(4.640371e-02)+9.0;
    logt[23]=log10(5.800464e-02)+9.0;
    logt[24]=log10(6.960557e-02)+9.0;
    logt[25]=log10(8.120650e-02)+9.0;
    logt[26]=log10(9.280742e-02)+9.0;
    logt[27]=log10(1.044084e-01)+9.0;
    logt[28]=log10(1.160093e-01)+9.0;
    logt[29]=log10(2.320186e-01)+9.0;
    logt[30]=log10(3.480278e-01)+9.0;
    logt[31]=log10(4.640371e-01)+9.0;
    logt[32]=log10(5.800464e-01)+9.0;
    logt[33]=log10(6.960557e-01)+9.0;
    logt[34]=log10(8.120650e-01)+9.0;
    logt[35]=log10(9.280742e-01)+9.0;
    logt[36]=log10(1.044084e+00)+9.0;
    logt[37]=log10(1.160093e+00)+9.0;
    logt[38]=log10(2.320186e+00)+9.0;
    logt[39]=log10(3.480278e+00)+9.0;
    logt[40]=log10(4.640371e+00)+9.0;
    logt[41]=log10(5.800464e+00)+9.0;
    logt[42]=log10(6.960557e+00)+9.0;
    logt[43]=log10(8.120650e+00)+9.0;
    logt[44]=log10(9.280742e+00)+9.0;
    logt[45]=log10(1.044084e+01)+9.0;
    logt[46]=log10(1.160093e+01)+9.0;
    logt[47]=log10(2.320186e+01)+9.0;
    logt[48]=log10(3.480278e+01)+9.0;
    logt[49]=log10(4.640371e+01)+9.0;
    logt[50]=log10(5.800464e+01)+9.0;
    logt[51]=log10(6.960557e+01)+9.0;
    logt[52]=log10(8.120650e+01)+9.0;
    logt[53]=log10(9.280742e+01)+9.0;
    logt[54]=log10(1.044084e+02)+9.0;
    logt[55]=log10(1.000000e+10)+9.0;
    
    return;
  };
  
  /** \brief Specify the stellar structure
      
      This function provides the stellar profile by
      copying the information in the \ref nscool_crust table
      to the Fortran arrays. These arrays are unit-indexed 
      in the original Fortran code, but the pointers which
      are sent to this function are pointers to the first
      element, so the C-style arrays in this function are
      zero-indexed.

      See \ref nscool_tov for the proper units and column names for
      the various quantities.

      This function's wrapper is called in <tt>precool.f</tt>.
  */
  virtual void star_struct(int icore, double rhocore,
			   double *rad_t, double *bar_t, double *rho_t,
			   double *pres_t, double *emas_t, double *phi_t,
			   double *rad, int *jmax, int *jcore,
			   double *w1, double *w2) {

    *jmax=nscool_tov.get_nlines();
    if (*jmax>9999) {
      O2SCL_ERR("Table too large.",o2scl::exc_einval);
    }
    *jcore=0;
    for(size_t j=0;j<nscool_tov.get_nlines();j++) {
      rad_t[j]=nscool_tov.get("radius",j)*100.0;
      bar_t[j]=nscool_tov.get("n_baryon",j);
      // This is energy density
      rho_t[j]=nscool_tov.get("density",j);
      pres_t[j]=nscool_tov.get("pressure",j);
      emas_t[j]=nscool_tov.get("emass",j);
      phi_t[j]=nscool_tov.get("phi",j);
      if (rho_t[j]<rhocore && (*jcore)==0) {
	*jcore=j;
      }
    }
    if (rad_t[1]<rad_t[0]) {
      O2SCL_ERR("Structure table should be increasing in radius",
		o2scl::exc_einval);
    }
    if (*jcore==0) {
      O2SCL_ERR2("Variable 'jcore' not set in ",
		 "nscool_wrap::star_struct().",o2scl::exc_einval);
    }

    /*
      double drho=rho_t[(*jcore)-1]-rho_t[(*jcore)];
      *w1=(rhocore-rho_t[(*jcore)])/drho;
      *w2=1.0-(*w1);
      double rad_core=(*w1)*rad_t[(*jcore)-1]+(*w2)*rad_t[(*jcore)];
      for(size_t i=0;i<=icore;i++) {
      rad[i]=cbrt(((float)i)/((float)icore))*rad_core;
      }
    */
    
    return;
  }
  
  /** \brief Specification of nucleon gaps

      This function's wrapper is called in <tt>precool.f</tt>.
  */
  virtual void gaps(int &sfn1s0_arg, double &n1_tc_arg,
		    double &n1_kf_arg, double &n1_dk_arg,
		    int &sfn3p2_arg, double &n3_tc_arg,
		    double &n3_kf_arg, double &n3_dk_arg,
		    int &sfp1s0_arg, double &p1_tc_arg,
		    double &p1_kf_arg, double &p1_dk_arg) {

    // If sfn3p2 and sfp1s0 are 150, then the double parameters
    // specify the Gaussian
    sfn3p2_arg=sfn3p2;
    sfp1s0_arg=sfp1s0;
    sfn1s0_arg=sfn1s0;
    if (sfn3p2_arg==150) {
      n3_tc_arg=n3_tc;
      n3_kf_arg=n3_kf;
      n3_dk_arg=n3_dk;
    }
    if (sfp1s0_arg==150) {
      p1_tc_arg=p1_tc;
      p1_kf_arg=p1_kf;
      p1_dk_arg=p1_dk;
    }
    if (sfn1s0_arg==150) {
      n1_tc_arg=n1_tc;
      n1_kf_arg=n1_kf;
      n1_dk_arg=n1_dk;
    }
    return;
  }

  /** \brief Main output 

      This function stores the main cooling curve output into \ref
      v_time, \ref v_tptr, \ref v_lphot, \ref v_lneut, and \ref
      v_lheat .

      This function's wrapper is called in <tt>NSCool.f</tt>.
      
      This function also flips the stop flag if more than 
      10000 steps are taken by the main cooling loop.
  */
  virtual void main_out(double &time, double &tptr,
			double &lphot, double &lneut, double &lheat,
			int &stop) {

    if (((int)v_time.size())%main_out_it==(main_out_it-1)) {
      std::cout.width(4);
      std::cout << v_time.size() << " "
		<< time << " " << tptr << " " << lphot << " "
		<< lneut << " " << lheat << std::endl;
    }

    v_time.push_back(time);
    v_tptr.push_back(tptr);
    v_lphot.push_back(lphot);
    v_lneut.push_back(lneut);
    v_lheat.push_back(lheat);

    stop=0;
    if (v_time.size()>10000) stop=1;
    
    return;
  }
  //@}
  
 public:

  /** \brief Number of iterations to skip
      for the main output function (default 20)
  */
  int main_out_it;
  
  /** \brief Execute main cooling calculation

      \todo Document in what cases <tt>iret</tt> can be nonzero
  */
  virtual int run(int irank) {
    if (nscool_wrap_ptrs.size()==0) {
      std::cerr << "Object 'nscool_wrap_ptrs' is empty." << std::endl;
      exit(-1);
    }
    if (v_time.size()>0) {
      v_time.clear();
      v_tptr.clear();
      v_lphot.clear();
      v_lneut.clear();
      v_lheat.clear();
    }
    int iret=0;
    // nscool_ is a subroutine and thus has no return value
    nscool_(&irank,&iret,pb_logt,pb_nalpha,pb_n2,
	    sf_lgtau1,sf_lgtau2,sf_lgr,sf_lgr2);
    return iret;
  }

  /** \brief Output the cooling curve to a file
   */
  void write_cool_curve(std::string fname="cool_curve.o2") {
    o2scl_hdf::hdf_file hf;
    hf.open_or_create(fname);
    o2scl::table_units<> t;
    t.line_of_names("t T L_neut L_phot L_heat");
    t.set_unit("t","yr");
    t.set_unit("T","K");
    t.set_unit("L_neut","erg/s");
    t.set_unit("L_phot","erg/s");
    t.set_unit("L_heat","erg/s");
    for(size_t i=0;i<v_time.size();i++) {
      double line[5]={v_time[i],v_tptr[i],v_lneut[i],v_lphot[i],
		      v_lheat[i]};
      t.line_of_data(5,line);
    }
    hdf_output(hf,t,"cool_curve");
    hf.close();
    return;
  }
  
  /** \brief Write temperature and luminosity profiles to a file
      
      This function copies the data into a table with a reorganized
      temperature grid (necessary because the cooling code doesn't
      know at what time the temperature will drop below the minimum).
  */
  void write_tl_prof(std::string fname="tl_prof.o2") {

    if (!tl_prof.is_constant("it_last")) {
      std::cerr << "No table to write in write_tl_prof()."
		<< std::endl;
      return;
    }
    
    o2scl_hdf::hdf_file hf;
    hf.open_or_create(fname);
    
    // Create new table omitting temperatures not stored
    int num_times=1+((int)(tl_prof.get_constant("it_last")+1.0e-6));
    
    // If the time grid hasn't been filled, then restructure
    // the table 
    if (num_times!=((int)time_print.size())) {
      o2scl::table3d tl_prof2;
      const ubvector &r_grid=tl_prof.get_x_data();
      const ubvector &t_grid=tl_prof.get_y_data();
      ubvector t_grid2(num_times);
      o2scl::vector_copy(num_times,t_grid,t_grid2);
      tl_prof2.set_xy("r",r_grid.size(),r_grid,"t",num_times,t_grid2);
      tl_prof2.line_of_names(((std::string)"rho ephi vol Tinf Linf qnu ")+
			     "qeebrem qnpb qplasma qsynch qbubble qpair "+
			     "qphoto qbrem_nn qmurca_nucl qbrem_nucl "+
			     "qmurca_hyp qbrem_hyp qdurca_np qdurca_lap "+
			     "qdurca_smn qdurca_smla qdurca_sms0 qfast "+
			     "qdurca_q qmurca_q qpbf_n1s0 qpbf_p1s0 "+
			     "qpbf_n3p2 qpbf_q qmax cv "+
			     "cv_n cv_p cv_e cv_m cv_la "+
			     "cv_sm cv_s0 cv_sp cv_q");
      for(size_t i=0;i<r_grid.size();i++) {
	for(size_t j=0;j<((size_t)num_times);j++) {
	  tl_prof2.set(i,j,"rho",tl_prof.get(i,j,"rho"));
	  tl_prof2.set(i,j,"ephi",tl_prof.get(i,j,"ephi"));
	  tl_prof2.set(i,j,"vol",tl_prof.get(i,j,"vol"));
	  tl_prof2.set(i,j,"Tinf",tl_prof.get(i,j,"Tinf"));
	  tl_prof2.set(i,j,"Linf",tl_prof.get(i,j,"Linf"));
	  tl_prof2.set(i,j,"qnu",tl_prof.get(i,j,"qnu"));
	  tl_prof2.set(i,j,"qeebrem",tl_prof.get(i,j,"qeebrem"));
	  tl_prof2.set(i,j,"qnpb",tl_prof.get(i,j,"qnpb"));
	  tl_prof2.set(i,j,"qplasma",tl_prof.get(i,j,"qplasma"));
	  tl_prof2.set(i,j,"qsynch",tl_prof.get(i,j,"qsynch"));
	  tl_prof2.set(i,j,"qbubble",tl_prof.get(i,j,"qbubble"));
	  tl_prof2.set(i,j,"qpair",tl_prof.get(i,j,"qpair"));
	  tl_prof2.set(i,j,"qphoto",tl_prof.get(i,j,"qphoto"));
	  tl_prof2.set(i,j,"qbrem_nn",tl_prof.get(i,j,"qbrem_nn"));
	  tl_prof2.set(i,j,"qmurca_nucl",tl_prof.get(i,j,"qmurca_nucl"));
	  tl_prof2.set(i,j,"qbrem_nucl",tl_prof.get(i,j,"qbrem_nucl"));
	  tl_prof2.set(i,j,"qmurca_hyp",tl_prof.get(i,j,"qmurca_hyp"));
	  tl_prof2.set(i,j,"qbrem_hyp",tl_prof.get(i,j,"qbrem_hyp"));
	  tl_prof2.set(i,j,"qdurca_np",tl_prof.get(i,j,"qdurca_np"));
	  tl_prof2.set(i,j,"qdurca_lap",tl_prof.get(i,j,"qdurca_lap"));
	  tl_prof2.set(i,j,"qdurca_smn",tl_prof.get(i,j,"qdurca_smn"));
	  tl_prof2.set(i,j,"qdurca_smla",tl_prof.get(i,j,"qdurca_smla"));
	  tl_prof2.set(i,j,"qdurca_sms0",tl_prof.get(i,j,"qdurca_sms0"));
	  tl_prof2.set(i,j,"qfast",tl_prof.get(i,j,"qfast"));
	  tl_prof2.set(i,j,"qdurca_q",tl_prof.get(i,j,"qdurca_q"));
	  tl_prof2.set(i,j,"qmurca_q",tl_prof.get(i,j,"qmurca_q"));
	  tl_prof2.set(i,j,"qpbf_n1s0",tl_prof.get(i,j,"qpbf_n1s0"));
	  tl_prof2.set(i,j,"qpbf_p1s0",tl_prof.get(i,j,"qpbf_p1s0"));
	  tl_prof2.set(i,j,"qpbf_n3p2",tl_prof.get(i,j,"qpbf_n3p2"));
	  tl_prof2.set(i,j,"qpbf_q",tl_prof.get(i,j,"qpbf_q"));
	  tl_prof2.set(i,j,"qmax",tl_prof.get(i,j,"qmax"));
	  tl_prof2.set(i,j,"cv",tl_prof.get(i,j,"cv"));
	  tl_prof2.set(i,j,"cv_n",tl_prof.get(i,j,"cv_n"));
	  tl_prof2.set(i,j,"cv_p",tl_prof.get(i,j,"cv_p"));
	  tl_prof2.set(i,j,"cv_e",tl_prof.get(i,j,"cv_e"));
	  tl_prof2.set(i,j,"cv_m",tl_prof.get(i,j,"cv_m"));
	  tl_prof2.set(i,j,"cv_la",tl_prof.get(i,j,"cv_la"));
	  tl_prof2.set(i,j,"cv_sm",tl_prof.get(i,j,"cv_sm"));
	  tl_prof2.set(i,j,"cv_s0",tl_prof.get(i,j,"cv_s0"));
	  tl_prof2.set(i,j,"cv_sp",tl_prof.get(i,j,"cv_sp"));
	  tl_prof2.set(i,j,"cv_q",tl_prof.get(i,j,"cv_q"));
	}
      }

      tl_prof2.set_interp_type(o2scl::itp_nearest_neigh);
      std::cout << "Herex." << std::endl;
      o2scl::table3d t3dug=tl_prof2.slice_to_uniform_grid
	("qmax",100,false,100,true);
      std::cout << "Herex2." << std::endl;
      tl_prof2.set_interp_type(o2scl::itp_linear);
      t3dug.set_interp_type(o2scl::itp_linear);
      for(size_t k=0;k<tl_prof2.get_nslices();k++) {
	std::string sl_name=tl_prof2.get_slice_name(k);
	if (sl_name!="qmax") {
	  t3dug.add_slice_from_table(tl_prof2,sl_name,sl_name);
	}
      }
      
      hdf_output(hf,((const o2scl::table3d &)(tl_prof2)),"tl_prof");
      hdf_output(hf,((const o2scl::table3d &)(t3dug)),"tl_prof_ug");

    } else {
      
      tl_prof.set_interp_type(o2scl::itp_nearest_neigh);
      o2scl::table3d t3dug=tl_prof.slice_to_uniform_grid
	("qmax",100,false,100,true);
      tl_prof.set_interp_type(o2scl::itp_linear);
      t3dug.set_interp_type(o2scl::itp_linear);
      for(size_t k=0;k<tl_prof.get_nslices();k++) {
	std::string sl_name=tl_prof.get_slice_name(k);
	if (sl_name!="qmax") {
	  t3dug.add_slice_from_table(tl_prof,sl_name,sl_name);
	}
      }
      
      hdf_output(hf,((const o2scl::table3d &)(tl_prof)),"tl_prof");
      hdf_output(hf,((const o2scl::table3d &)(t3dug)),"tl_prof_ug");
      
    }
    
    hf.close();
    return;
  }

};

/// Set superfluid gaps
extern "C" void nscool_gaps_
(int *irank, double *sfn1s0, double *n1tc, double *n1kf, double *n1dk,
 double *sfn3p2, double *n3tc, double *n3kf, double *n3dk, double *sfp1s0,
 double *p1tc, double *p1kf, double *p1dk) {

  // The superfluid switches are actually double-precision
  // numbers in the cooling code, so we convert from int
  // to double here.
  int isfn3p2, isfp1s0, isfn1s0;
  nscool_wrap_ptrs[*irank]->gaps(isfn1s0,*n1tc,*n1kf,*n1dk,
				 isfn3p2,*n3tc,*n3kf,*n3dk,
				 isfp1s0,*p1tc,*p1kf,*p1dk);
				 
  *sfn3p2=((double)isfn3p2);
  *sfp1s0=((double)isfp1s0);
  *sfn1s0=((double)isfn1s0);
  return;
}

/// Main output function
extern "C" void nscool_main_out_(int *irank, double *time, double *tptr,
				 double *lphot, double *lneut, double *lheat,
				 int *stop) {
  nscool_wrap_ptrs[*irank]->main_out(*time,*tptr,*lphot,*lneut,*lheat,*stop);
  return;
}

/// Obtain stellar structure
extern "C" void nscool_star_struct_
  (int *irank, int *icore,
   double *rhocore, double *rad_t, double *bar_t, double *rho_t,
   double *pres_t, double *emas_t, double *phi_t, double *rad,
   int *jmax, int *jcore, double *w1, double *w2) {
  nscool_wrap_ptrs[*irank]->star_struct(*icore,*rhocore,rad_t,bar_t,rho_t,
					pres_t,emas_t,phi_t,rad,
					jmax,jcore,w1,w2);
  return;
}

/// Obtain core composition
extern "C" void nscool_core_comp_
(int *irank, double *rho_t, double *nbar_t, double *yelect_t, double *ymuon_t,
 double *yneutr_t, double *yprot_t, double *ylambda_t,
 double *ysminus_t, double *yszero_t, double *ysplus_t,
 double *mstp_t, double *mstn_t, double *mstla_t, double *mstsm_t,
 double *msts0_t, double *mstsp_t, int *ix) {
  nscool_wrap_ptrs[*irank]->core_comp(rho_t,nbar_t,yelect_t,ymuon_t,
				      yneutr_t,yprot_t,ylambda_t,ysminus_t,
				      yszero_t,ysplus_t,mstp_t,mstn_t,mstla_t,
				      mstsm_t,msts0_t,mstsp_t,ix);
  return;
}

/// Obtain crust composition
extern "C" void nscool_crust_comp_
(int *irank, double *Z_ion_t, double *A_ion_t, double *A_cell_t,
 double *bar_t, double *pres_t, double *rho_t, int *jmax) {
  nscool_wrap_ptrs[*irank]->crust_comp(Z_ion_t,A_ion_t,A_cell_t,bar_t,
				       pres_t,rho_t,jmax);
  return;
}

/// Obtain crust EOS
extern "C" void nscool_crust_eos_(int *irank, double *rho2, double *pres2,
				  int *idata) {
  nscool_wrap_ptrs[*irank]->crust_eos(rho2,pres2,idata);
  return;
}

/// Set numerical parameters
extern "C" void nscool_num_param_
(int *irank, double *time0, double *timemax, int *istepmax,
 int *itrial_max, int *itrial_opt, double *tcut,
 double *dtime, double *dtlimit, double *scale_dt0,
 double *scale_dt1, double *repeat, int *istart,
 double *mratt, double *mratl, double *mrats,
 double *tvar, double *svar, double *tcon) {

  nscool_wrap_ptrs[*irank]->num_param
    (*time0,*timemax,*istepmax,*itrial_max,*itrial_opt,*tcut,
     *dtime,*dtlimit,*scale_dt0,*scale_dt1,*repeat,*istart,
     *mratt,*mratl,*mrats,*tvar,*svar,*tcon);

  return;
}

/// Set cooling parameters
extern "C" void nscool_cool_param_
(int *irank, int *pscreen, double *debug, int *istep_debug,
 double *pteff, double *ptemp, double *pstar,
 int *idump1, int *idump2, int *idump3,
 double *tempmin, double *tempini,
 int *icvel_nodeg, double *emnco, double *emncr,
 double *emp, double *p0, int *itpmax, double *tprint) {

  nscool_wrap_ptrs[*irank]->cool_param
    (*pscreen,*debug,*istep_debug,*pteff,*ptemp,*pstar,
     *idump1,*idump2,*idump3,*tempmin,*tempini,
     *icvel_nodeg,*emnco,*emncr,*emp,*p0,*itpmax,tprint);
  
  return;
}

/// Set initial temperature profile
extern "C" void nscool_tptr_init_(int *irank, int *ifteff, double *tempini,
				  double *ephi_surf, double *ephi_drip,
				  double *ephi_core, double *tsurface,
				  double *tdrip, double *tcore,
				  double *tb_acc0) {
  nscool_wrap_ptrs[*irank]->tptr_init(*ifteff,*tempini,*ephi_surf,*ephi_drip,
				      *ephi_core,*tsurface,*tdrip,
				      *tcore,*tb_acc0);
  return;
}

/// Set atmosphere parameters
extern "C" double nscool_teff_(int *irank, double *Tb, int *ifteff,
			       double *eta, double *bfield, int *istep,
			       double *time, double *Ts1, double *Ts2,
			       double *Z, double *A, double *Rho,
			       int *debug, double *gs14,
			       double *compactness) {
  return nscool_wrap_ptrs[*irank]->Teff(*Tb,*ifteff,*eta,*bfield,*istep,
					*time,*Ts1,*Ts2,*Z,*A,*Rho,*debug,
					*gs14,*compactness);
}

/// Set boundary parameters
extern "C" void nscool_bound_param_(int *irank, int *ifteff, double *eta,
				    double *mag_coeff, double *tb_acc0) {
  nscool_wrap_ptrs[*irank]->bound_param(*ifteff,*eta,*mag_coeff,*tb_acc0);
  return;
}

/// Print out temperatures and luminosities
extern "C" void nscool_print_temp_
(int *irank, int *istep, int *itprint, double *time,
 double *t_effective, int *imax, double *w1, 
 double *w2, double *otemp, double *temp,
 double *olum, double *lum, double *rad,
 double *rrho, double *ephi, double *dvol,
 double *e2phi, double *qnu, double *qeebrem,
 double *qnpb, double *qplasma,
 double *qsynch, double *qbubble, double *qpair,
 double *qphoto, double *qbrem_nn,
 double *qmurca_nucl, double *qbrem_nucl,
 double *qmurca_hyp, double *qbrem_hyp,
 double *qdurca_np, double *qdurca_lap,
 double *qdurca_smn, double *qdurca_smla,
 double *qdurca_sms0, double *qfast,
 double *qdurca_q, double *qmurca_q, 
 double *qpbf_n1s0, double *qpbf_p1s0,
 double *qpbf_n3p2, double *qpbf_q) {
  nscool_wrap_ptrs[*irank]->print_temp
    (*istep,*itprint,*time,*t_effective,*imax,*w1,*w2,otemp,temp,
     olum,lum,rad,rrho,ephi,dvol,e2phi,qnu,qeebrem,qnpb,qplasma,
     qsynch,qbubble,qpair,qphoto,qbrem_nn,qmurca_nucl,qbrem_nucl,
     qmurca_hyp,qbrem_hyp,qdurca_np,qdurca_lap,qdurca_smn,qdurca_smla,
     qdurca_sms0,qfast,qdurca_q,qmurca_q,qpbf_n1s0,qpbf_p1s0,qpbf_n3p2,
     qpbf_q);
  return;
}

/// Print specific heats
extern "C" void nscool_print_cv_
(int *irank, int *itprint, int *imax, double *cv, double *cv_n, double *cv_p, 
 double *cv_e, double *cv_m, double *cv_la, double *cv_sm, double *cv_s0, 
 double *cv_sp, double *cv_q) {
  nscool_wrap_ptrs[*irank]->print_cv
    (*itprint,*imax,cv,cv_n,cv_p,cv_e,
     cv_m,cv_la,cv_sm,cv_s0,cv_sp,cv_q);
  return;
}

/// Handle Urca settings
extern "C" void nscool_urca_settings_(int *irank, double *fix_durca,
				      double *alpha_durca,
				      double *beta_durca) {
  nscool_wrap_ptrs[*irank]->urca_settings(*fix_durca,*alpha_durca,
					  *beta_durca);
  return;
}

/// Compute critical temperatures
extern "C" void nscool_tc_new_(int *irank, int *imax, int *idrip, int *icore,
			       double *sfn1s0, double *n1tc,
			       double *n1kf, double *n1dk,
			       double *sfp1s0, double *p1tc,
			       double *p1kf, double *p1dk,
			       double *sfn3p2, double *n3tc,
			       double *n3kf, double *n3dk,
			       double *fn1s0, double *fp1s0, double *fn3p2,
			       double *kf_n, double *kf_p, double *kf_la,
			       double *kf_sm, double *kf_s0, double *kf_sp,
			       double *kf_u, double *kf_d, double *kf_s,
			       double *tc_n, double *tc_p, double *tc_la,
			       double *tc_uu, double *tc_dd, double *tc_ss,
			       double *tc_ud, double *tc_us, double *tc_ds,
			       double *tc_u, double *tc_d, double *tc_s,
			       int *isf) {

  /*
    std::cout << "H2: " << *icore << " " << *sfn1s0 << " "
    << *fn1s0 << std::endl;
    std::cout << "H2: " << *icore << " " << *sfp1s0 << " "
    << *fp1s0 << std::endl;
  */
  
  nscool_wrap_ptrs[*irank]->atc.get_tc(*imax,*idrip,*icore,
				       *sfn1s0,*n1tc,*n1kf,*n1dk,
				       *sfp1s0,*p1tc,*p1kf,*p1dk,
				       *sfn3p2,*n3tc,*n3kf,*n3dk,
				       *fn1s0,*fp1s0,*fn3p2,
				       kf_n,kf_p,kf_la,kf_sm,kf_s0,kf_sp,
				       kf_u,kf_d,kf_s,tc_n,tc_p,tc_la,
				       tc_uu,tc_dd,tc_ss,tc_ud,tc_us,
				       tc_ds,tc_u,tc_d,tc_s,isf);

  return;
}

extern "C" void nscool_P_electron_(int *irank, double *ne, double *T,
				   double *pre) {
  nscool_wrap_ptrs[*irank]->P_electron(*ne,*T,pre);
  return;
}

#endif
