#include <limits>
#include <mpi.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include "physical_constants.hh"
#include "transport.hh"
#include "radioactive.hh"
#include <omp.h>


//--------------------------------------------------------
// constructor
//--------------------------------------------------------
TRANSPORT::TRANSPORT()
{
  step_size = 0.1;
}


//--------------------------------------------------------
// Initialize and allocate
//--------------------------------------------------------
void TRANSPORT::Init()
{


  // default values
  grey_opac = 0.0;
  n_living_particles = 0;
  
  // Setup particle buffers
  n_particles = 0;
  particle = new PARTICLE[MAX_PARTICLES];   // list of all particles
  pbuffer  = new PARTICLE[MAX_PARTICLES];   // used to reshuffle particles
  
  // start at time = 0 with zero particles
  n_particles = 0;
  t_now = 0;
  
  // get mpi rank
  int my_rank; 
  MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
  if (my_rank == 0) verbose = 1; else verbose = 0;

  // setup and seed random number generator
  const gsl_rng_type * TypeR;
  gsl_rng_env_setup();
  gsl_rng_default_seed = (unsigned int)time(NULL) + my_rank;
  TypeR = gsl_rng_default;
  rangen = gsl_rng_alloc (TypeR);
}



//*********************************************************
// --------------------------------------------------------
// Do radiation transport over a single time step
// --------------------------------------------------------
//*********************************************************
double TRANSPORT::Step(double t_step)
{
  int i;

  // set the start timer 
  time_t start_tp,end_tp;
  time(&start_tp);
  
  // emit new particles
  Emit_Particles(t_step);
  
  // set counter to zero
  n_living_particles = 0;

  // --------------------------------------------------------------
  // Propogate the particles
  // --------------------------------------------------------------
  /*
  int my_rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &my_rank ); */
  #pragma omp for private(i)
  /* printf("total number of threads on process %ld: %ld\n", my_rank,
         omp_get_num_threads()); */
  for (i = 0; i < n_particles; i++)
  {
    /* printf("sending particle %ld on process %d on thread %d...\n",
           i, my_rank, omp_get_thread_num()); */
    // check if particles are overflowing buffer
    if (i >= MAX_PARTICLES) {
      printf("particles overflowing buffer: ");
      printf("n_particles = %ld, MAX_PARTICLES = %ld\n",
             n_particles, MAX_PARTICLES);
      printf("stopping...\n");
      exit(1); }

    //    printf("GO %d\n",i);
    // propogate a single particle and count if escaped
    if (particle[i].fate == alive) Propogate(particle[i],t_step);
    
    // count living particles
    if (particle[i].fate == alive) n_living_particles++;
  }
  // --------------------------------------------------------------

  // advance time step
  t_now += t_step;

  // calculate the elapsed time 
  time(&end_tp);
  double time_wasted=difftime(end_tp,start_tp);
  return time_wasted;
  
}




//*********************************************************
//--------------------------------------------------------
// Propogate a single monte carlo particle until
// it  escapes, is absorbed, or the time step ends
//--------------------------------------------------------
//*********************************************************
void TRANSPORT::Propogate(PARTICLE &p, double dt)
{
  enum EVENT {scatter, boundary, tstep};
  EVENT event;

  // time of end of timestep
  double tstop = t_now + dt;
  // local variables
  double tau_r,d_sc,d_tm,this_d;
  double opac;

  // propogate until this flag is set
  int stop = 0;

  // automatically absorb positrons, convert to photons
  if (p.type == positron)  p.type = photon; 

  while (!stop)
  {
    // maximum step size inside zone
    double d_bn = step_size*grid->Get_dx();
    
    // local velocity vector
    double V[3];
    double t_secs = p.t*DAY_TO_SEC;
    V[0] = p.x[0]/t_secs;
    V[1] = p.x[1]/t_secs;
    V[2] = p.x[2]/t_secs;

    // local Lorentz transformation params
    p.beta   = sqrt(V[0]*V[0] + V[1]*V[1] + V[2]*V[2])/C_LIGHT;
    p.gamma  = 1.0/sqrt(1 - p.beta*p.beta); 
    double vdotD  = V[0]*p.D[0] + V[1]*p.D[1] + V[2]*p.D[2];
    double dshift = p.gamma*(1 - vdotD/C_LIGHT);

    // density in this zone
    double rho = grid->Get_Density(p.ind);
    // determine opacity if its a gamma-ray
    if (p.type == gammaray)  
    {
      double ne_dens = rho/M_PROTON*grid->Get_mu_e(p.ind);
      opac = THOMSON_CS*Klein_Nishina(p.E_x)*ne_dens; 
    }
    // determine opacity if it is an optical photon
    else if (p.type == photon) opac = rho*grey_opac;

    // convert opacity from comoving to lab frame
    opac = opac*dshift;

    // random optical depth to next interaction
    tau_r = -1.0*log(1 - gsl_rng_uniform(rangen));
    
    // step size to next interaction event
    d_sc  = tau_r/opac;
    if (opac == 0) d_sc = std::numeric_limits< double >::infinity();
   
    // find distance to end of time step
    d_tm = (tstop - p.t)*DAY_TO_SEC*C_LIGHT;

    // find out what event happens (shortest distance)
    if ((d_sc < d_bn)&&(d_sc < d_tm)) 
      {event = scatter; this_d = d_sc;}
    else if (d_bn < d_tm)
      {event = boundary; this_d = d_bn;}
    else 
      {event = tstep; this_d = d_tm; }

    // move particle the distance
    p.x[0] += this_d*p.D[0];
    p.x[1] += this_d*p.D[1];
    p.x[2] += this_d*p.D[2]; 
    // advance the time
    p.t = p.t + this_d/C_LIGHT/DAY_TO_SEC;
    
    // ---------------------------------
    // Do if scatter
    // ---------------------------------
    if (event == scatter) Compton_Scatter(p);

    // ---------------------------------
    // do if time step end
    // ---------------------------------
    else if (event == tstep) { stop = 1;}

    // find zone index of particle now
    p.ind = grid->Get_Zone_Index(p.x);
  
    // particle escapes if p.ind < 0
    if (p.ind < 0) 
    {
      p.fate = escaped;
      spectrum.Count(p);
      stop = 1;
    }

  }
  
}



//------------------------------------------------------------
// physics of non-isotropic Compton scattering
//------------------------------------------------------------
void TRANSPORT::Compton_Scatter(PARTICLE &p)
{
  // local velocity vector
  double V[3];
  double t_secs = p.t*DAY_TO_SEC;
  V[0] = p.x[0]/t_secs;
  V[1] = p.x[1]/t_secs;
  V[2] = p.x[2]/t_secs;
  double vdotD  = V[0]*p.D[0] + V[1]*p.D[1] + V[2]*p.D[2];
  double dshift_in  = p.gamma*(1 - vdotD/C_LIGHT);

  // transform quantities into comoving frame
  p.energy = p.energy*dshift_in;
  p.E_x    = p.E_x*dshift_in;
  
  // sample new direction by rejection method
  double E_ratio;
  double D_new[3];
  while (true)
  {
    // choose new isotropic direction in comoving frame
    double mu  = 1 - 2.0*gsl_rng_uniform(rangen);
    double phi = 2.0*PI*gsl_rng_uniform(rangen);
    double smu = sqrt(1 - mu*mu);
    D_new[0] = smu*cos(phi);
    D_new[1] = smu*sin(phi);
    D_new[2] = mu;
    
    // use rejection method for differential compton cross-section
    int reject = 0;
    if (p.type == gammaray)
    {
      // angle between old and new directions
      double cost = p.D[0]*D_new[0] + p.D[1]*D_new[1] + p.D[2]*D_new[2];
      // new energy ratio (E_new/E_old) at this angle (assuming lambda in MeV)
      E_ratio = 1/(1 + p.E_x/E_ENERGY_MEV*(1 - cost));
      // klein-nishina differential cross-section
      double diff_cs = 0.5*(E_ratio*E_ratio*(1/E_ratio + E_ratio - 1 + cost*cost));
      // see if this scatter angle OK
      double y = gsl_rng_uniform(rangen);
      if (y > diff_cs) reject = 1;
    }
    if (!reject) break;
  }

  // new gamma-ray energy (i.e., frequency)
  if (p.type == gammaray) {
    grid->add_edep(p.ind, p.E_x - p.E_x*E_ratio);
    p.E_x = p.E_x*E_ratio;
  }
  
  // sample whether we stay alive, if not become a photon
  double y = gsl_rng_uniform(rangen);
  if (y > E_ratio)  {
    p.type = photon; 
     // choose new isotropic direction in comoving frame
    double mu  = 1 - 2.0*gsl_rng_uniform(rangen);
    double phi = 2.0*PI*gsl_rng_uniform(rangen);
    double smu = sqrt(1 - mu*mu);
    D_new[0] = smu*cos(phi);
    D_new[1] = smu*sin(phi);
    D_new[2] = mu;
  }
  
  // outgoing velocity vector
  for (int i=0;i<3;i++) V[i] = -1*V[i];
  
  // doppler shifts outgoing
  double vdp = (D_new[0]*V[0] + D_new[1]*V[1] + D_new[2]*V[2]);
  double vd_out = p.gamma*(1 - vdp/C_LIGHT);
  
  // transformation of direction vector back into lab frame
  p.D[0] = 1.0/vd_out*(D_new[0] - p.gamma*V[0]/C_LIGHT*(1 - p.gamma*vdp/C_LIGHT/(p.gamma+1)));
  p.D[1] = 1.0/vd_out*(D_new[1] - p.gamma*V[1]/C_LIGHT*(1 - p.gamma*vdp/C_LIGHT/(p.gamma+1)));
  p.D[2] = 1.0/vd_out*(D_new[2] - p.gamma*V[2]/C_LIGHT*(1 - p.gamma*vdp/C_LIGHT/(p.gamma+1)));
  double norm = sqrt(p.D[0]*p.D[0] + p.D[1]*p.D[1] + p.D[2]*p.D[2]);
  p.D[0] = p.D[0]/norm;
  p.D[1] = p.D[1]/norm;
  p.D[2] = p.D[2]/norm;
  
  // transformation of energy/wavelength into lab frame
  p.energy = p.energy*vd_out;
  p.E_x    = p.E_x*vd_out;
}
  
  


//-----------------------------------------------------------------
// Klein_Nishina correction to the Compton cross-section
// assumes energy x is in MeV
//-----------------------------------------------------------------
double TRANSPORT::Klein_Nishina(double x)
{
  // divide by m_e c^2 = 0.511 MeV
  x = x/E_ENERGY_MEV;
  double logfac = log(1 + 2*x);
  double term1 = (1+x)/x/x/x*(2*x*(1+x)/(1+2*x) - logfac);
  double term2 = 1.0/2.0/x*logfac;
  double term3 = -1.0*(1 + 3*x)/(1+2*x)/(1+2*x);
  double KN    = .75*(term1 + term2 + term3);
  return KN;
}




//-----------------------------------------------------------------
// Emit particles from radioactive decay
//-----------------------------------------------------------------
void TRANSPORT::Emit_Particles(double dt)
{
  RADIOACTIVE rad;
  int n_x   = grid->Get_n_x();
  double x_cen = grid->Get_x_cen();
  double dx = grid->Get_dx();

  // decay energy rate per unit gram at this time
  double dEdt = rad.Decay_Energy_Rate((t_now + 0.5*dt)*DAY_TO_SEC);
  // total energy emitted over this time step
  double Etot = dEdt*dt*grid->Get_Nickel_Mass();
  // number of photons to add 
  double n_add_step = dt*n_photons_per;
  // energy per photon particle
  double Ep = Etot/n_add_step;
  
  for (int i=0;i<n_x;i++)
  for (int j=0;j<n_x;j++)
  for (int k=0;k<n_x;k++)
  {
    int ind = grid->Get_Index(i,j,k);

    // energy emitted in this zone
    double E = dEdt*dt*grid->Get_Nickel_Mass(ind);
    // number of photons to add
    int n_add = floor(E/Ep);
    // pick up remainder randomly
    if (gsl_rng_uniform(rangen) < E/Ep - n_add) n_add++;
    
    // rebuffer particle list if necessary
    if (n_particles+n_add > MAX_PARTICLES) Rebuffer_Particles();
    // check that we have enough memory to hold particles
    if (n_particles+n_add > MAX_PARTICLES) {
      printf("Ran out of particle space\n");
      return; }

    // setup particles
    for (int q=n_particles;q<n_particles+n_add;q++)
    {
      particle[q].ind  = ind;
      particle[q].fate = alive;

      // randomly sample position in zone
      particle[q].x[0] = i*dx - x_cen + dx*gsl_rng_uniform(rangen);
      particle[q].x[1] = j*dx - x_cen + dx*gsl_rng_uniform(rangen);
      particle[q].x[2] = k*dx - x_cen + dx*gsl_rng_uniform(rangen);

      // emit randomly over time step
      particle[q].t      = t_now + dt*gsl_rng_uniform(rangen);
      
      // emit isotropically
      double mu  = 1 - 2.0*gsl_rng_uniform(rangen);
      double phi = 2.0*PI*gsl_rng_uniform(rangen);
      double smu = sqrt(1 - mu*mu);
      particle[q].D[0] = smu*cos(phi);
      particle[q].D[1] = smu*sin(phi);
      particle[q].D[2] = mu;

      // total energy of particle packet
      particle[q].energy = Ep;

      // automatically count particle into radioactive energy
      particle[q].type = radioactive;
      spectrum.Count(particle[q]);
  
      // energy of particles inside packet (in MeV)
      particle[q].E_x   = rad.Sample_Ni56_Wavelength(particle[q].t,rangen);
      if (particle[q].E_x < 0) particle[q].type = positron;
      else particle[q].type  = gammaray;

      // convert from comoving to lab frame
      double t_secs = t_now*DAY_TO_SEC;
      double V[3];
      V[0] = particle[q].x[0]/t_secs;
      V[1] = particle[q].x[1]/t_secs;
      V[2] = particle[q].x[2]/t_secs;

      // local Lorentz transformation params
      particle[q].beta   = sqrt(V[0]*V[0] + V[1]*V[1] + V[2]*V[2])/C_LIGHT;
      particle[q].gamma  = 1.0/sqrt(1 - particle[q].beta*particle[q].beta); 
      double vdotD  = V[0]*particle[q].D[0] + V[1]*particle[q].D[1] + V[2]*particle[q].D[2];
      double dshift = particle[q].gamma*(1 - vdotD/C_LIGHT);
      
      // transformation of direction vector back into lab frame
      particle[q].D[0] = 1.0/dshift*(particle[q].D[0] - particle[q].gamma*V[0]/C_LIGHT*(1 - particle[q].gamma*vdotD/C_LIGHT/(particle[q].gamma+1)));
      particle[q].D[1] = 1.0/dshift*(particle[q].D[1] - particle[q].gamma*V[1]/C_LIGHT*(1 - particle[q].gamma*vdotD/C_LIGHT/(particle[q].gamma+1)));
      particle[q].D[2] = 1.0/dshift*(particle[q].D[2] - particle[q].gamma*V[2]/C_LIGHT*(1 - particle[q].gamma*vdotD/C_LIGHT/(particle[q].gamma+1)));
      double norm = sqrt(particle[q].D[0]*particle[q].D[0] + particle[q].D[1]*particle[q].D[1] + particle[q].D[2]*particle[q].D[2]);
      particle[q].D[0] = particle[q].D[0]/norm;
      particle[q].D[1] = particle[q].D[1]/norm;
      particle[q].D[2] = particle[q].D[2]/norm;
      
      // transformation of energy/wavelength into lab frame
      particle[q].energy *= dshift;
      particle[q].E_x    *= dshift;

    }
    n_particles += n_add;
  }
}


//*********************************************************
// --------------------------------------------------------
// Function to reshuffle particles and get rid of
// dead ones
// --------------------------------------------------------
//*********************************************************
void TRANSPORT::Rebuffer_Particles()
{
  int i;
  if (verbose) printf("Rebuffering photons ");

  // save particles that are alive
  n_pbuffer = 0;
  for (i=0;i<n_particles;i++)
    if (particle[i].fate == alive) {
      pbuffer[n_pbuffer].copy(particle[i]);
      n_pbuffer++; }
  
  // put particles back into main array
  for (i=0;i<n_pbuffer;i++) particle[i].copy(pbuffer[i]);
  n_particles = n_pbuffer;
  n_living_particles = n_particles;

  if (verbose) printf(" - done %ld \n",n_particles);
}
