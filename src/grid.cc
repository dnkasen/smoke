#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "grid.hh"
#include "physical_constants.hh"

GRID::GRID()
{
}


//------------------------------------------------------------   
// Read in the model file in ascii format
//------------------------------------------------------------ 
void GRID::Read(const char *infile)
{
  // get mpi rank
  int my_rank; 
  MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
  if (my_rank == 0) verbose = 1; else verbose = 0;

  FILE *in = fopen(infile,"r");
  if (in == NULL) {
    printf("ERROR:  Can't open setup file %s\n",infile);
    exit(1); }

  // read in basic grid parameters
  fscanf(in,"%d %lf %lf\n",&n_x,&dv,&t_begin);
  n_zones = n_x*n_x*n_x;
  dx = dv*t_begin*DAY_TO_SEC;
  x_cen   = n_x/2.0*dx;
  t_begin = t_begin;


  // estimate size of grid
  if (verbose)
  {
    double zsize = 1.0*sizeof(ZONE);
    printf("# n_x = %d (%d zones, %.0f bytes each) = %.3f GBs\n",
	   n_x,n_zones,zsize,n_zones*1.0*zsize/1e9);
  }

  // allocate memory
  z = new ZONE[n_zones];

  // calculate sums
  double tmass   = 0;
  double ke      = 0;
  ni_mass        = 0;
  double vol = dx*dx*dx;


  // read size of file
  double x1,x2,x3;
  int ind = 0;
  for (int i=0;i<n_x;i++)
    for (int j=0;j<n_x;j++)
      for (int k=0;k<n_x;k++)
      {
	fscanf(in,"%lf %lf %lf",&x1,&x2,&x3);
	z[ind].rho     = x1;
	z[ind].ni_frac = x2;
	z[ind].mu_e    = x3;
	
	tmass    += z[ind].rho*vol;
	ni_mass  += z[ind].rho*z[ind].ni_frac*vol;
	double vx = (i*dx - x_cen)/(t_begin*DAY_TO_SEC);
	double vy = (j*dx - x_cen)/(t_begin*DAY_TO_SEC);
	double vz = (k*dx - x_cen)/(t_begin*DAY_TO_SEC);
	double vv = vx*vx + vy*vy + vz*vz;
	ke    += 0.5*z[ind].rho*vol*vv;
	ind++;
      }
  fclose(in);

  if (verbose)
  {
    printf("# Model read\n");
    printf("# Total mass = %.3e (%.3e Msun)\n",tmass,tmass/M_SUN);
    printf("# 56ni  mass = %.3e (%.3e Msun)\n",ni_mass,ni_mass/M_SUN);
    printf("# Kinetic E  = %.3e ergs\n",ke);
  }      
} 


//------------------------------------------------------------   
// Get the index of the i,j,k zone
//------------------------------------------------------------ 
int GRID::Get_Index(int i, int j, int k)
{
  return i*n_x*n_x + j*n_x + k;
}

//------------------------------------------------------------   
// Locate the zone index given real x,y,z coordinates
//------------------------------------------------------------ 
int GRID::Get_Zone_Index(double *x)
{
  int ix = round((x[0]+x_cen)/dx); 
  int iy = round((x[1]+x_cen)/dx);
  int iz = round((x[2]+x_cen)/dx);
  
  // check for off grid
  if ((ix < 0)||(ix >= n_x)) return -1;
  if ((iy < 0)||(iy >= n_x)) return -1;
  if ((iz < 0)||(iz >= n_x)) return -1;

  return Get_Index(ix,iy,iz);
}

//------------------------------------------------------------   
// Homologously expand grid by a factor e
//------------------------------------------------------------   
void GRID::Expand(double e)
{
  dx = dx*e;
  x_cen = x_cen*e;
  for (int i=0;i<n_zones;i++)
    z[i].rho = z[i].rho/e/e/e;
}
