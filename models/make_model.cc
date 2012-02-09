//--------------------------------------------------------
// simple C++ code to make a initial 3D model file for 
// the smoke monte carlo transport code
//
// header gives:
//  n_x   dv   t_exp
// where
//  n_x  = number of zones in each dimension 
//   dv  = size of each zone (in velocity space, cm/s)
// t_exp = time sinc explosion (in days)
//
// the columns then give (for i,j,k) zones
//   dens  X_Ni   mu_e
// where
// density = density of zone in g/cm^3
//    X_Ni = mass abundance of nickel in the one    
//    mu_e = fraction of electrons per nucleon
//--------------------------------------------------------

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_const_cgsm.h>
#include <gsl/gsl_math.h>

int main(void)
{
  int    nx      = 150;           // number of zones per dimension
  double vmax    = 1e9;           // maximum velocity at edge (cm/s)
  double mass    = 1.4*GSL_CONST_CGSM_SOLAR_MASS;     // total mass of ejecta
  double tstart  = 1.0*60*60*24;  // start time of transport calculation
  double Mni_in  = 0.5;           // mass coordinate of 56Ni
  double Mni_out = 0.75;          // outer mass coordinate of 56Ni

  // output name of model file
  char   fname[] = "lucy_test.mod";

  FILE *out = fopen(fname,"w");
 
  double r0 = vmax*tstart;
  double dx = 2*r0/nx;
  double dv = dx/tstart;
  double rho = mass/(4.0*M_PI/3.0*r0*r0*r0);

  /// print header
  fprintf(out,"%d %e %e\n",nx,dv,tstart/3600./24.);
 
  double tmass = 0;
  for (int i=0;i<nx;i++)
    for (int j=0;j<nx;j++)
      for (int k=0;k<nx;k++)
      {
        double x = (i - nx/2)*dx;
        double y = (j - nx/2)*dx;
        double z = (k - nx/2)*dx;
        double r = sqrt(x*x + y*y + z*z);
        
        double vx = vmax*(x/r0);
        double vy = vmax*(y/r0);
        double vz = vmax*(z/r0);
        double vr = sqrt(vx*vx + vy*vy + vz*vz);

        // set composition
<<<<<<< HEAD
        double m   = rho*4.0*M_PI/3.0*r*r*r/GSL_CONST_CGSM_SOLAR_MASS;
	double ni_comp = 1.0;
	if (m > Mni_in)  ni_comp = 1.0 - (m - Mni_in)/(Mni_out - Mni_in);
        if (m > Mni_out) ni_comp = 0.0;
        
        double mu = 0.5; // 1 electron for every two nucleons
        double d = rho;
        if (vr > vmax) d = rho*1e-10;

        fprintf(out,"%14.4e %14.4e %14.4e\n ",d,ni_comp,mu);
      }
}
