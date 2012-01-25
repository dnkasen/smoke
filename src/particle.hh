#ifndef _PARTICLE_H
#define _PARTICLE_H

enum PFATE {alive, escaped};
enum PTYPE {photon, gammaray,positron,radioactive};

// particle class
class PARTICLE 
{

public:
  
  PTYPE type;          // current type (photon or gamma-ray)
  PFATE fate;          // current status (alive or not)

  double x[3];         // x,y,z position
  double D[3];         // direction vector, Dx,Dy,Dz
  long int  ind;       // index of the zone in grid where we are
  double t;            // current time
  double energy;       // total energy in ergs of packet
  double E_x;          // energy of photons in gamma-ray packet in MeV
  double beta;         // local velocity divided by c
  double gamma;        // local lorentz factor (1-v^2/c^2)^(-1/2)

  double r() 
  { return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]); }

   double x_dot_d() 
  {return x[0]*D[0] + x[1]*D[1] + x[2]*D[2]; }

  double mu() 
  {return x_dot_d()/r(); }

  void copy(PARTICLE p)
  {
    for (int i=0;i<3;i++) 
      {x[i] = p.x[i]; D[i] = p.D[i]; }

    ind    = p.ind;
    t      = p.t;
    energy = p.energy;
    fate   = p.fate;
    type   = p.type;
    beta   = p.beta;
    gamma  = p.gamma;
    E_x    = p.E_x;
  }

};

#endif
