#ifndef _RADIOACTIVE_H
#define _RADIOACTIVE_H

// radioactive decay constonstants
#define TAU_NI  760579         // decay time = 8.803   days
#define TAU_CO  9.81599e+06    // decay time = 113.611 days 

// average energy per decay in MeV
#define AVERAGE_NI_ENERGY 1.728
#define AVERAGE_CO_ENERGY 3.566

class RADIOACTIVE
{
 
public:

  double Total_Decay_Energy(double);
  double Decay_Energy_Rate(double);
  double Sample_Ni56_Wavelength(double time, gsl_rng*);

};


#endif
