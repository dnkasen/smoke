#ifndef _GRID_H
#define _GRID_H

class ZONE
{

public:

  float rho;      // mass density (g cm^(-3))
  float ni_frac;  // 56ni mass fraction
  float mu_e;     // number of electrons per particle
  double edep;    // energy deposited due to Compton scatters
  double vel;     // velocity of zone

};

//-------------------------------------------------
// Class to hold all the zones
//-------------------------------------------------
class GRID
{

private:

  int verbose;
  ZONE *z;          // array of zones
  int n_zones;      // total number of zones
  int n_x;          // number of zones in x-dimension
  double dx;        // linear size of zone in cm
  double dv;        // size of zone in velocity coords (cm/s)
  double x_cen;     // location of grid center
  double t_begin;   // start time of the model
  double ni_mass;   // total 56ni mass


public:

  GRID();


  void Read(const char *infile);
  int  Get_Zone_Index(double *x);
  void Expand(double);
  
  int    Get_Index(int i, int j, int k);
  double Get_Nickel_Mass(int i) { return dx*dx*dx*z[i].rho*z[i].ni_frac; }
  double Get_Nickel_Mass()      { return ni_mass; }
  int    Get_n_zones() {return n_zones;}
  int    Get_n_x()     {return n_x;}
  double Get_x_cen()   {return x_cen;}
  double Get_t_begin() {return t_begin; }
  double Get_dx() {return dx; }
  double Get_Density(int i) { return z[i].rho; }
  double Get_mu_e(int i) {return z[i].mu_e; }
  double Get_vel(int i) {return z[i].vel; }

  // add energy to the zone during a Compton scatter
  void add_edep(int i, double energy) {
    z[i].edep = z[i].edep + energy;
    return;
  }

  // print gamma-ray energy deposition information
  void show_dep() {
    for (int i = 0; i < n_zones; i++) {
      if (z[i].edep != 0) {
        printf("%d %le %le\n", i, z[i].vel/1.0e5, z[i].edep);
      }
    }
    return;
  }

};


#endif
