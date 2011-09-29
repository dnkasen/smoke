#ifndef _TRANSPORT_H
#define _TRANSPORT_H
#include "particle.hh"
#include "spectrum.hh"
#include "grid.hh"

#define MAX_PARTICLES 1000000

class TRANSPORT
{

private:


  // arrays of particles
  long int   n_particles; 
  PARTICLE  *particle;
  PARTICLE  *pbuffer;
  long int   n_pbuffer;
  long int   n_living_particles;

   // random number generator
  gsl_rng *rangen;

 
  void   Propogate(PARTICLE &p, double tstop);
  void   Emit_Particles(double dt);
  double Klein_Nishina(double);
  void   Compton_Scatter(PARTICLE&);
  void   Rebuffer_Particles();


public:

  // number of total particles
  long int   n_total_particles;
 
  // current time in simulation
  double t_now;
  
  // pointer to grid
  GRID *grid;

  // emmergent spectrum,
  SPECTRUM spectrum;

  int    verbose;
  double n_photons_per;         // number of photons to emit per day
  double grey_opac;
  double emit_min, emit_max;
  double step_size;

 
  TRANSPORT();
  void   Init();
  double Step(double dt);
  int    num_particles()        {return n_particles;}
  int    num_living_particles() {return n_living_particles;}
  

 
};

#endif

