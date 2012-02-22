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
#include <omp.h>

#define N_COARSE_VEL_GRID 256
#define COARSE_VEL_MAX 15000

/* I think you probably don't have to do this in C++, but I'm still
 * learning C... */
int compare_times(const void *, const void *);

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
  const int rank_root = 0;
  MPI_Init( &argc, &argv );
  // timers to see if workload is balanced
  double proc_time_start, proc_time_end, proc_time;
  proc_time_start = MPI_Wtime();
  MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
  MPI_Comm_size( MPI_COMM_WORLD, &n_procs);
  if (my_rank == rank_root) verbose = 1;
  if (verbose) printf("\n# Using %d MPI cores\n\n",n_procs);

  // create coarse(r) 1-D spherical velocity grid on which to interpolate
  // gamma-ray deposition
  double *coarse_vel_grid, *coarse_vel_grid_e;
  coarse_vel_grid   = new double[N_COARSE_VEL_GRID];
  coarse_vel_grid_e = new double[N_COARSE_VEL_GRID];
  for (int i = 0; i < N_COARSE_VEL_GRID; i++) {
    coarse_vel_grid[i] = COARSE_VEL_MAX*(double(i)/N_COARSE_VEL_GRID);
    coarse_vel_grid_e[i] = 0.0;
  }

  // open up the parameter file
  std::string script_file = "param.dat";
  if( argc > 1 ) script_file = std::string( argv[ 1 ] );

  //--------------------------------------
  // default values of input parameters
  double tstep_min   = 0.0;
  double tstep_max   = 1.0;
  double tstep_del   = 0.2;
  double opacity     = 0.1;
  double t_delta     = 1.0;
  double t_stop      = 50.0;
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
	   (double)n_photons/(double)n_procs, (double)n_photons);

  // start off at beginning time
  double t = grid.Get_t_begin();
  

  if (verbose) 
  {
    printf("#\n");
    printf("%1s %7s %9s %10s %14s %9s %17s\n", "#", "step", "time", "dt",
           "n_living", "n_total", "exec time(sec)#");
  }
  
  // loop over time steps;   
  for (int it = 0; it < n_times; it++)
  {
    // get this time step
    double t_step = tstep_max;
    if (t_step < tstep_min) t_step = tstep_min;
    if (t_step > tstep_max) t_step = tstep_max;
    if ((tstep_del > 0.0) && (t > 0.0))
      if (t_step > t*tstep_del)
        t_step = t*tstep_del;
     
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

  int i, j;

  double *tot_coarse_vel_grid_e;
  tot_coarse_vel_grid_e = new double[N_COARSE_VEL_GRID];

  double *vel_grid_e, *tot_vel_grid_e;
  vel_grid_e     = new double[grid.Get_n_zones()];
  tot_vel_grid_e = new double[grid.Get_n_zones()];

  /* this is stupid to store the deposition data twice, but I don't know
   * how to do an MPI_Reduce on members of a class */
  for (i = 0; i < grid.Get_n_zones(); i++) {
    vel_grid_e[i] = grid.Get_edep(i);
    tot_vel_grid_e[i] = 0.0;
  }

  double edep_sum = 0.0;

  for (i = 0; i < N_COARSE_VEL_GRID; i++) {
    tot_coarse_vel_grid_e[i] = 0.0;
  }

  for (i = 0; i < grid.Get_n_zones(); i++) {
    edep_sum += vel_grid_e[i];
  }

  int error;
  error = MPI_Reduce(vel_grid_e, tot_vel_grid_e,
                     grid.Get_n_zones(), MPI_DOUBLE, MPI_SUM, 0,
		     MPI_COMM_WORLD);
  if (error != MPI_SUCCESS)
    if (verbose)
      printf("ERROR: MPI_Reduce failed!\n");

  edep_sum = 0.0;

  if (verbose) {
    for (i = 0; i < grid.Get_n_zones(); i++) {
      edep_sum += tot_vel_grid_e[i];
    }
  }

  // interpolate from fine grid to coarse grid
  if (verbose) {
    // right now this is done with a linear search. trees would be way
    // faster...
    for (i = 0; i < N_COARSE_VEL_GRID; i++) {
      for (j = 0; j < grid.Get_n_zones(); j++) {
        if ((grid.Get_vel(j)/1.0e5 > coarse_vel_grid[i]) &&
            (grid.Get_vel(j)/1.0e5 <= coarse_vel_grid[i+1]) &&
            (tot_vel_grid_e[j] > 0.0)) {
          tot_coarse_vel_grid_e[i] += tot_vel_grid_e[j];
        }
      }
    }
  }
  if (verbose) transport.spectrum.Print();

  // finish up mpi
  proc_time_end = MPI_Wtime();
  proc_time = proc_time_end - proc_time_start;

  /* Now let's practice with trees. We'll compare the longest and shortest
   * process times using a binary tree sort from libc. */

  /* first we gather all the process times into a single array on the
   * master process */
  int gsize;
  MPI_Comm_size(MPI_COMM_WORLD, &gsize);
  double *all_proc_times;
  if (verbose ) { all_proc_times = new double[gsize]; }
  MPI_Gather(&proc_time, 1, MPI_DOUBLE, all_proc_times, 1, MPI_DOUBLE,
             rank_root, MPI_COMM_WORLD);

  // now sort!
  if (verbose)
    { qsort(all_proc_times, gsize, sizeof(double), compare_times); }

  if (verbose) {
    FILE *pfile;
    pfile = fopen("dep.dat", "w");
    fprintf(pfile, "%15s %15s\n", "VEL (km/s)", "DEP (MeV)");
    for (i = 0; i < N_COARSE_VEL_GRID; i++) {
      fprintf(pfile, "%15e %15e\n", coarse_vel_grid[i],
              tot_coarse_vel_grid_e[i]);
    }
    fclose(pfile);
  }

  // calculate the elapsed time 
  time(&end_tp);
  float time_wasted=difftime(end_tp,start_tp)/60.0;
  if (verbose) {
    printf("#\n# CALCULATION took %.2f minutes (%.3f hours)\n",
	   time_wasted,time_wasted/60.0);
    printf("slowest process time: %.3f min\n", all_proc_times[0] / 60.0);
    printf("fastest process time: %.3f min\n", all_proc_times[gsize-1] / 60.0);
  }
  delete coarse_vel_grid, coarse_vel_grid_e, tot_coarse_vel_grid_e,
         vel_grid_e, tot_vel_grid_e;
  if (verbose) { delete all_proc_times; }

  MPI_Finalize();

  return 0;
}

// function needed by qsort to compare process times
int compare_times(const void *time1, const void *time2) {
  double dtime1 = *(double *)time1;
  double dtime2 = *(double *)time2;
  int result;
  if (dtime2 > dtime1) {
    result = 1;
  } else if (dtime2 < dtime1) {
    result = -1;
  } else if (dtime2 == dtime1) {
    result = 0;
  }
  return result;
}
