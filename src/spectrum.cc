#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <vector>
#include "physical_constants.hh"
#include "spectrum.hh"

//***************************************************************
// Constructors
//***************************************************************

SPECTRUM::SPECTRUM()
{
  strcpy(name,DEFAULT_NAME);
  n_elements = 0;
  n_escaped = 0;
}



//***************************************************************
// Initialization and Allocation
//***************************************************************
void SPECTRUM::Init(double t1, double t2, double dt, int nm, int np)
{
  t_start = t1;
  t_stop  = t2;
  t_delta = dt;
  
  n_times    = (int)((t_stop-t_start)/t_delta);
  n_mu       = nm;
  n_phi      = np;
  n_elements = n_times*n_mu*n_phi;
  
  // allocate
  click  = new double[n_elements];
  count  = new double[n_elements];
  gamma  = new double[n_elements];
  radio  = new double[n_elements];
  
  // clear 
  Wipe();
}



void SPECTRUM::Set_Name(const char *n)
{
  strcpy(name,n);
}

//***************************************************************
// Functional procedure: Wipe
//***************************************************************
void SPECTRUM::Wipe()
{
  n_escaped = 0;
  for (int i=0;i<n_elements;i++) 
  {
    click[i]  = 0;
    count[i]  = 0;
    gamma[i]  = 0;
    radio[i]  = 0;
  }
}



//***************************************************************
//--------------------------------------------------------------
//********************** COUNT FUNCTIONS ************************
//--------------------------------------------------------------
//***************************************************************

int SPECTRUM::Index(int i, int j, int k)
{
  return i*n_mu*n_phi + j*n_phi + k;
}


void SPECTRUM::Count(PARTICLE p)
{
  n_escaped++;

  double mu  = p.D[2];
  double phi = atan2(p.D[1],p.D[0]) + PI;

  // do correction for light crossing time
  double xdot  = p.x[0]*p.D[0] + p.x[1]*p.D[1] + p.x[2]*p.D[2];
  double t_obs = p.t - xdot/C_LIGHT/DAY_TO_SEC;

  // locate bin number in all dimensions
  int t_bin = (int)((t_obs-t_start)/t_delta);
  int m_bin = (int)((mu + 1)/2*n_mu);
  int p_bin = (int)((phi/2/PI*n_phi));

  // if off the grids, just return without counting
  if ((t_bin < 0)||(m_bin < 0)||(p_bin < 0)) return;
  if ((t_bin >= n_times)||(m_bin >= n_mu)||(p_bin >= n_phi)) return;
      
  // add to counters
  int ind      = Index(t_bin,m_bin,p_bin);
  if (p.type == photon)      count[ind]  += p.energy;
  if (p.type == gammaray)    gamma[ind]  += p.energy;
  if (p.type == radioactive) radio[ind]  += p.energy;
  if (p.type == photon)      click[ind]  += 1;
}


//***************************************************************
//--------------------------------------------------------------
//********************** PRINT FUNCTIONS *************************
//--------------------------------------------------------------
//***************************************************************
void SPECTRUM::Print()
{
  FILE *out = fopen(name,"w");
 
  fprintf(out,"# %d %d %d\n",n_times,n_mu,n_phi);

    for (int m=0;m<n_mu;m++)
      for (int p=0;p<n_phi;p++)
	for (int i=0;i<n_times;i++)
	{
	fprintf(out,"%12.5e ",t_start + t_delta*(i + 0.5));
	if (n_mu  > 1) fprintf(out,"%12.5e ",-1.0 + (m+0.5)*2.0/n_mu);
	if (n_phi > 1) fprintf(out,"%12.5e ",2.0*PI*(p+0.5)/n_phi - PI);
	int ind  = Index(i,m,p);
	fprintf(out,"%12.5e %12.5e %12.5e %12.5e %12.5e\n",count[ind], 
		radio[ind]-gamma[ind],gamma[ind],radio[ind],click[ind]);

      }
  fclose(out);

}


//***************************************************************
//--------------------------------------------------------------
//********************** MPI FUNCTIONS **************************
//--------------------------------------------------------------
//***************************************************************

void SPECTRUM::Normalize()
{
  // renormalize flux
  double norm_factor =  1.0*n_mu*n_phi/(1.0*t_delta);
  Rescale(norm_factor);
}


void SPECTRUM::Rescale(double r)
{
  for (int i=0;i<n_elements;i++) count[i] *= r;
}


void SPECTRUM::MPI_Sum_All()
{
  MPI_Allreduce_Array(count);
  MPI_Allreduce_Array(click);
  MPI_Allreduce_Array(gamma);
  MPI_Allreduce_Array(radio);
}

void SPECTRUM::MPI_Average_All()
{
  MPI_Sum_All();

  int mpi_procs;
  MPI_Comm_size( MPI_COMM_WORLD, &mpi_procs );
  for (int i=0;i<n_elements;i++) 
  {
    count[i] /= mpi_procs; 
    gamma[i] /= mpi_procs;
    radio[i] /= mpi_procs;
  }
}


void SPECTRUM::MPI_Allreduce_Array(double *arr)
{
  double *new_ptr;
  int j;          

  // allocate the memory for new pointer
  int chunk = n_elements;
  new_ptr = new double[chunk];
         
  // zero out array
  for (j=0;j<chunk;j++) new_ptr[j] = 0;
  // reduce the stuff
  MPI_Allreduce(arr,new_ptr,chunk,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  // put back into place
  for (j=0;j<chunk;j++) arr[j] = new_ptr[j];

  // free up the memory
  delete new_ptr;
}
