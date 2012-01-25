#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include "radioactive.hh"
#include "physical_constants.hh"


// energy (Mev) an probability (# per decay) of gamma rays from a Co56 decay
int n_co56_decays   = 17;
double co56_prob[]  = { 0.360,   1.000,   0.015,   0.137,  0.022,  0.670,
                         0.043,   0.158,   0.031,   0.079,  0.166,  0.058,
                         0.030,   0.074,   0.018,   0.009,  0.000};
double co56_energy[] = { 0.511,   0.847,   0.980,   1.040,  1.180,  1.240,
                         1.360,   1.770,   2.015,   2.030,  2.600,  3.010,
                         3.200,   3.250,   3.270,   3.450,  4.000};
 
int n_ni56_decays    =  7;
double ni56_prob[]   =  {1.000, 0.360, 0.360, 0.500, 0.870, 0.140, 0.000};
double ni56_energy[] =  {0.158, 0.270 ,0.480, 0.750, 0.812, 1.562, 2.000};



int    n_R_proc_fit = 11;
double R_proc_fit[] =  {17.608179, -2.0442059, -0.42565322, 0.39830095,  -0.0059089906,
			-0.054805836, 0.014068697, -0.00086706160, -5.7056758e-05, 2.6401842e-06,
			3.7186979e-07};



//--------------------------------------------------------------------
// returns the total ergs/g from radioactive decay up to time t stop
//---------------------------------------------------------------------
double RADIOACTIVE::Total_Decay_Energy(double t)
{
  // fraction of 56ni left and the end
  double xni = exp(-t/TAU_NI);
  // fraction of 56Co that is left
  double xco = TAU_CO/(TAU_NI-TAU_CO)*(exp(-t/TAU_NI) - exp(-t/TAU_CO));

  // energy per atom in MeV
  double E = AVERAGE_NI_ENERGY*(1-xni) + AVERAGE_CO_ENERGY*(1-xco);
 
  // convert to ergs per gram
  E = E*MEV_TO_ERGS/(56*M_PROTON);
  return E;

}


//--------------------------------------------------------------
// returns the energy from radioactive decay, per time per gram
//--------------------------------------------------------------
double RADIOACTIVE::Decay_Energy_Rate(double time)
{
  // exponential factors to be used
  double e_ni = exp(-time/TAU_NI);
  double e_co = exp(-time/TAU_CO);
  // number divided by decay time)
  double ni56 = e_ni/TAU_NI;
  double co56 = 1.0/(TAU_NI-TAU_CO)*(e_ni - e_co);
  // // get the energy from decays in ergs/s, using unit conversions
  double ni_E = ni56*(AVERAGE_NI_ENERGY*MEV_TO_ERGS);
  double co_E = co56*(AVERAGE_CO_ENERGY*MEV_TO_ERGS);
  
  // energy rate per unit gram
  double dEdt = (ni_E + co_E)/(56*M_PROTON);
  return dEdt;
}



double RADIOACTIVE::Sample_Ni56_Wavelength(double time, gsl_rng *rangen)
{
  // The Ratio Of Energy Coming Out In Nickel
  double E_Ni = exp(-time/TAU_NI);
  double E_Co = exp(-time/TAU_CO);
  double Ni_E = AVERAGE_NI_ENERGY*(E_Ni/TAU_NI);
  double Co_E = AVERAGE_CO_ENERGY/(TAU_NI - TAU_CO)*(E_Ni - E_Co);
  double nico_ratio = Ni_E/(Ni_E + Co_E);
  
  // pick emission wavelength
  double x_val;
  int d;
  double z1 = gsl_rng_uniform(rangen);
  if (z1 < nico_ratio)
  {
    while (true)
    {
      double z2 = gsl_rng_uniform(rangen);
      double z3 = gsl_rng_uniform(rangen);
      d = (int)(n_ni56_decays*z2);
      if (z3 < ni56_prob[d]) break;
    }
    x_val = ni56_energy[d];
  }
  else
  {
    while (true)
    {
      // put 2% of energy into positrons
      double z4 = gsl_rng_uniform(rangen);
      if (z4 < 0.02) return -1;

      double z2 = gsl_rng_uniform(rangen);
      double z3 = gsl_rng_uniform(rangen);
      d = (int)(n_co56_decays*z2);
      if (z3 < co56_prob[d]) break;
    }
    x_val = co56_energy[d];
  }
  
  return x_val;
}




