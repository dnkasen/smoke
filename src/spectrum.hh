#ifndef _SPECTRUM_H
#define _SPECTRUM_H 1

#include <string>
#include "particle.hh"
using std::string;

// default values
#define DEFAULT_NAME "spectrum"

class SPECTRUM {
 
private:

  // spectrum name
  char name[1000];

  // number of elements and bin sizes
  int n_elements,n_times,n_mu,n_phi;
  double t_start, t_stop, t_delta;

  // counting arrays
  double *click,*count;
  double *gamma,*radio;
  double n_escaped;

  // Indexing
  int Index(int,int,int);
  
  // Used for MPI
  void MPI_Allreduce_Array(double *arr);
  
    
public:

  // constructors
  SPECTRUM();
  
  // Initialize
  void Init(double, double, double, int, int);
  void Set_Name(const char *n);

  // MPI functions
  void MPI_Sum_All();
  void MPI_Average_All();
  
  // Count and Normalize counted packets
  
  void Count(PARTICLE);
  void Normalize();
  void Rescale(double);
  void Wipe();

  // Print out
  void Print();
  
  double N_Escaped() {return n_escaped;}
};

#endif
