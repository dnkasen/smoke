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

    int my_rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );

    verbose = 0;

    // Only the master MPI task reads the model from disk.

    if( my_rank == 0 )
    {

        verbose = 1;

        FILE *in = fopen(infile,"r");
        if( in == NULL )
        {
            printf("ERROR:  Can't open setup file %s\n",infile);
            exit(1); 
        }
    
        // basic grid parameters

        fscanf(in,"%d %lf %lf\n",&n_x,&dv,&t_begin);
        n_zones = n_x*n_x*n_x;
        dx = dv*t_begin*DAY_TO_SEC;
        x_cen   = n_x/2.0*dx;
        t_begin = t_begin;
    
        // estimate size of grid

        double zsize = 1.0 * sizeof( ZONE );
        printf("# n_x = %d (%d zones, %.0f bytes each) = %.3f GBs\n",n_x,n_zones,zsize,n_zones*1.0*zsize/1e9);
    
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
        {
            for (int j=0;j<n_x;j++)
            {
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
                    z[ind].vel = sqrt(vv);
                    ke    += 0.5*z[ind].rho*vol*vv;
                    ind++;
                }
            }
        }

        fclose(in);
    
        printf("# Model read\n");
        printf("# Total mass = %.3e (%.3e Msun)\n",tmass,tmass/M_SUN);
        printf("# 56ni  mass = %.3e (%.3e Msun)\n",ni_mass,ni_mass/M_SUN);
        printf("# Kinetic E  = %.3e ergs\n",ke);

    }

    // Now we exchange the model to all the MPI tasks.

    MPI_Bcast( &n_zones, 1, MPI_INT   , 0, MPI_COMM_WORLD );
    MPI_Bcast( &n_x    , 1, MPI_INT   , 0, MPI_COMM_WORLD );
    MPI_Bcast( &dx     , 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &dv     , 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &x_cen  , 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &t_begin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &ni_mass, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );

    float* rho_buffer     = new float [ n_zones ];
    float* ni_frac_buffer = new float [ n_zones ];
    float* mu_e_buffer    = new float [ n_zones ];

    if( my_rank == 0 )
    {
        for( int i = 0; i < n_zones; ++ i )
        {
            rho_buffer    [ i ] = z[ i ].rho;
            ni_frac_buffer[ i ] = z[ i ].ni_frac;
            mu_e_buffer   [ i ] = z[ i ].mu_e;
        }
    }
    else
    {
        for( int i = 0; i < n_zones; ++ i )
        {
            rho_buffer    [ i ] = 0.0;
            ni_frac_buffer[ i ] = 0.0;
            mu_e_buffer   [ i ] = 0.0;
        }
    }

    MPI_Bcast( rho_buffer    , n_zones, MPI_FLOAT, 0, MPI_COMM_WORLD );
    MPI_Bcast( ni_frac_buffer, n_zones, MPI_FLOAT, 0, MPI_COMM_WORLD );
    MPI_Bcast( mu_e_buffer   , n_zones, MPI_FLOAT, 0, MPI_COMM_WORLD );

    if( my_rank != 0 )
    {
        z = new ZONE[ n_zones ];
        for( int i = 0; i < n_zones; ++ i )
        {
            z[ i ].rho     = rho_buffer    [ i ];
            z[ i ].ni_frac = ni_frac_buffer[ i ];
            z[ i ].mu_e    = mu_e_buffer   [ i ];
        }
    }

    delete [] rho_buffer;
    delete [] ni_frac_buffer;
    delete [] mu_e_buffer;
    
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
