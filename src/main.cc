#include <mpi.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fstream>
#include "physical_constants.hh"
#include "grid.hh"
#include "transport.hh"
#include "radioactive.hh"
#include <vector>

#define N_COARSE_VEL_GRID 256
#define COARSE_VEL_MAX 15000

//--------------------------------------------------------
// The main code
//--------------------------------------------------------
int main(int argc, char **argv)
{
  // essential classes to use
  GRID       grid;
  TRANSPORT  transport;
  int verbose = 0;

  // set the global start timer
  time_t start_tp,end_tp;
  time(&start_tp);


  //---------------------------------------------------------------------
  // BEGIN SETTING UP 
  //---------------------------------------------------------------------

  // initialize MPI parallelism
  int my_rank,n_procs;
  MPI_Init( &argc, &argv );
  time_t proc_start_tp, proc_end_tp; // see if workload is balanced
  time(&proc_start_tp);
  MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
  MPI_Comm_size( MPI_COMM_WORLD, &n_procs);
  if (my_rank == 0) verbose = 1;
  if (verbose) printf("\n# Using %d MPI cores\n\n",n_procs);

  // create coarse(r) 1-D spherical velocity grid on which to interpolate
  // gamma-ray deposition
  double *coarse_vel_grid, *coarse_vel_grid_e;
  coarse_vel_grid   = new double[N_COARSE_VEL_GRID];
  coarse_vel_grid_e = new double[N_COARSE_VEL_GRID];

  if (verbose) {
//    printf("%15s %15s", "SHELL INDEX", "VEL (KM/S)\n");
    for (int i = 0; i < N_COARSE_VEL_GRID; i++) {
      coarse_vel_grid[i] = COARSE_VEL_MAX*(double(i)/N_COARSE_VEL_GRID);
      coarse_vel_grid_e[i] = 0.0;
//      printf("%15d %15lf\n", i, coarse_vel_grid[i]);
    }
  }

  // open up the parameter file
  std::string script_file = "param.dat";
  if( argc > 1 ) script_file = std::string( argv[ 1 ] );

  //--------------------------------------
  // default values of input parameters
  double tstep_min   = 0;
  double tstep_max   = 1;
  double tstep_del   = 0.2;
  double opacity     = 0.1;
  double t_delta     = 1;
  double t_stop      = 50;
  int    n_times     = 1000;
  int    n_mu        = 1;
  int    n_phi       = 1;
  long int n_photons = 1e5;
  double step_size   = 0.1;
  //--------------------------------------

  // read paramater file
  std::string line, val,par, mod_file;
  std::string out_file = "light_curve.dat";
  std::ifstream pfile;
  pfile.open(script_file.c_str());
  if (!pfile.is_open()) {
    if (verbose) printf("Can't open param file %s\n",script_file.c_str());
    exit(1); }
  
  while (!pfile.eof())
  {
    std::getline(pfile,line);
    int f1 = line.find_first_of("=");
    int f2 = line.find_first_of("!");
    if (f2 < 0) f2 = line.length();
    if (f1>0) 
    {
      par = line.substr(0,f1);
      val = line.substr(f1+1,f2-f1-1);
      if (!par.find("tstep_min")) tstep_min = atof(val.c_str());
      if (!par.find("tstep_max")) tstep_max = atof(val.c_str());
      if (!par.find("tstep_del")) tstep_del = atof(val.c_str());
      if (!par.find("t_stop"))      t_stop  = atof(val.c_str());
      if (!par.find("t_delta"))     t_delta = atof(val.c_str());
      if (!par.find("n_mu"))          n_mu  = atoi(val.c_str());
      if (!par.find("n_phi"))        n_phi  = atoi(val.c_str());
      if (!par.find("opacity"))     opacity = atof(val.c_str());
      if (!par.find("step_size")) step_size = atof(val.c_str());
      if (!par.find("n_photons")) n_photons = (long int)(atof(val.c_str()));

      if (!par.find("model_file")) 
      {
	int d1 = val.find_first_of("\"");
	int d2 = val.find_last_of("\"");
	mod_file = val.substr(d1+1, d2-d1-1);
      }
      if (!par.find("output_file")) 
      {
	int d1 = val.find_first_of("\"");
	int d2 = val.find_last_of("\"");
	out_file = val.substr(d1+1, d2-d1-1);
      }
    }
  }

  // print out the parameters
  if (verbose)
  {
    printf("------- parameters ---------\n");
    printf("tstep_min = %.4e\n",tstep_min);
    printf("tstep_max = %.4e\n",tstep_max);
    printf("tstep_del = %.4e\n",tstep_del);
    printf("t_stop    = %.4e\n",t_stop);
    printf("t_delta   = %.4e\n",t_delta);
    printf("opacity   = %.4e\n",opacity);
    printf("n_photons = %.4e\n",(float)n_photons);
    printf("step_size = %.4f\n",step_size);
    printf("n_mu      = %d\n",n_mu);
    printf("n_phi     = %d\n",n_phi);
    printf("----------------------------\n\n");
  }

  // Read the model file to setup grid
  grid.Read(mod_file.c_str());

  // printout total energy from radioactive decay
  RADIOACTIVE rad;
  double E_dec = grid.Get_Nickel_Mass()*rad.Total_Decay_Energy(t_stop);
  if (verbose) printf("# total decay energy = %e\n",E_dec);

  
  // setup the transport module
  transport.Init();
  transport.grid = &grid;
  transport.t_now = grid.Get_t_begin();
  transport.grey_opac = opacity;
  transport.step_size = step_size;
  transport.spectrum.Init(grid.Get_t_begin(),t_stop,t_delta,n_mu,n_phi);
  transport.spectrum.Set_Name(out_file.c_str());
  transport.n_photons_per = n_photons/(t_stop)/n_procs;
  

  //---------------------------------------------------------------------
  // SETUP DONE; START CALCULATION
  //---------------------------------------------------------------------

  if (verbose)
    printf("# Sending around %.4e particles per MPI process (%.4e total)\n",
	   1.0*n_photons/n_procs, 1.0*n_photons);

  // start off at beginning time
  double t = grid.Get_t_begin();
  

  if (verbose) 
  {
    printf("#\n");
    printf("#    step      time         dt       n_living   n_total   exec time(sec)"); 
    printf("#\n");
  }
  
  // loop over time steps;   
  for (int it=0;it<n_times;it++)
  {
    // get this time step
    double t_step = tstep_max;
    if (t_step < tstep_min) t_step = tstep_min;
    if (t_step > tstep_max) t_step = tstep_max;
    if ((tstep_del > 0)&&(t > 0)) if (t_step > t*tstep_del) t_step = t*tstep_del;
     
    // printout time step
    if (verbose) printf("%8d %12.4e %12.4e %8d %8d ",it,t,t_step,
			transport.num_living_particles(),
			transport.num_particles());

    // Propogate the particles
    double t_exec = transport.Step(t_step); 
    if (verbose) printf("       %.3e\n",t_exec);
    
    // homologously expand the grid
    double efac = (t + t_step)/t;
    grid.Expand(efac);
 
    // advance time
    t = t + t_step;
    if (t > t_stop) break;
  }
  
  
  //---------------------------------------------------------------------
  // CALCULATION DONE; WRITE OUT AND FINISH
  //---------------------------------------------------------------------
  
  // Finish and write out final spectrum
  transport.spectrum.MPI_Average_All();
  transport.spectrum.Normalize();

  // interpolate from fine grid to coarse grid
  if (verbose) printf("interpolating onto coarse velocity grid...\n");
  // right now this is done with a linear search. trees would be waaaaaaay
  // more efficient...
  for (int i = 0; i < N_COARSE_VEL_GRID; i++) {
    for (int j = 0; j < grid.Get_n_zones(); j++) {
      if ((grid.Get_vel(j)/1.0e5 > coarse_vel_grid[i]) &&
          (grid.Get_vel(j)/1.0e5 <= coarse_vel_grid[i+1]))
        { coarse_vel_grid_e[i] += grid.Get_edep(j); }
    }
  }
  if (verbose) printf("interpolation complete!\n");

  double *tot_coarse_vel_grid_e;
  tot_coarse_vel_grid_e = new double[N_COARSE_VEL_GRID];

  double edep_sum = 0.0;

  for (int i = 0; i < N_COARSE_VEL_GRID; i++) {
    edep_sum += coarse_vel_grid_e[i];
  }
  printf("Before MPI_Reduce, on process %d, edep_sum = %lf\n",
         my_rank, edep_sum);

  int error;
  error = MPI_Reduce(coarse_vel_grid_e, tot_coarse_vel_grid_e,
                     N_COARSE_VEL_GRID, MPI_DOUBLE, MPI_SUM, 0,
		     MPI_COMM_WORLD);
  if (verbose) {
    printf("MPI_Reduce returned error: ");
    switch (error) {
      case MPI_SUCCESS:
        printf("MPI_SUCCESS\n");
        break;
      case MPI_ERR_COMM:
        printf("MPI_ERR_COMM\n");
        break;
      case MPI_ERR_TYPE:
        printf("MPI_ERR_TYPE\n");
        break;
      case MPI_ERR_BUFFER:
        printf("MPI_ERR_BUFFER\n");
        break;
      default:
        printf("Unknown MPI error flag\n");
        break;
    }
  }

  edep_sum = 0.0;

  if (verbose) {
    for (int i = 0; i < N_COARSE_VEL_GRID; i++) {
      edep_sum += tot_coarse_vel_grid_e[i];
    }
    printf("After MPI_Reduce, edep_sum = %lf\n", edep_sum);
  }

  time(&proc_end_tp);
  float proc_time_wasted = difftime(proc_end_tp, proc_start_tp)/60.0;
  printf("process %d took %.2f minutes\n", my_rank, proc_time_wasted);
  if (verbose) transport.spectrum.Print();

  // finish up mpi
  MPI_Finalize();
  

  // calculate the elapsed time 
  time(&end_tp);
  float time_wasted=difftime(end_tp,start_tp)/60.0;
  if (verbose) {
    printf("#\n# CALCULATION took %.3f minutes (%.2f hours)\n",
	   time_wasted,time_wasted/60.0);
    FILE *pfile;
    pfile = fopen("dep.dat", "w");
    grid.show_dep(pfile, coarse_vel_grid, coarse_vel_grid_e);
    fclose(pfile);
  }
  delete coarse_vel_grid, coarse_vel_grid_e, tot_coarse_vel_grid_e;
}
