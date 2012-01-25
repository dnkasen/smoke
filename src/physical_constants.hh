#ifndef PHYS_CONST
#define PHYS_CONST 1

// physical constants 
#define PI              3.14159         // just pi
#define C_LIGHT         2.99792458e10   // speed of light (c) in cm/s
#define C_LIGHT_KM      2.99792458e5    // speed of light (c) in km/s
#define FINE_STRUCTURE  7.297352533e-3  // fine structure constant (alpha) dimensionless 
#define K_BOLTZ         1.380658e-16    // boltz const in ergs/K
#define K_BOLTZ_EV      8.615e-5        // boltz const in eV/K
#define H_PLANCK        6.6260755e-27   // plancks constant (h) in ergs-s
#define STEF_BOLTZ      5.67e-5         // stefan-boltz (sigma) erg cm-2 s-1 K-4
#define RAD_CONST       7.5657e-15      // radiation constant (a) erg cm-3 K-4
#define M_ELECTRON      9.10938188e-28  // electron mass in grams
#define M_PROTON        1.67262158e-24  // proton mass in grams
#define E_ELECTRON      4.8e-10         // electron charge 
#define NEWTON_G        6.67e-8         // graviational const gm-1 cm3 s-2
#define THOMSON_CS      0.665e-24       // Thomson cross-section cm-2
#define CLASSICAL_CS    2.6575e-2       // PI*e^2/m_e/c in cm^2
#define SAHA_CONST      4.82940e15      // = (2.0*PI*m_e*k/h^2)^(1.5)
#define RYDBERG         2.1796e-11      // Rydberg constant in ergs
#define RYDBERG_EV      13.6058         // Rydberg constant in eV
#define E_ENERGY_MEV    0.511           // rest energy of electron in Mev

// astronomical constants
#define M_SUN           1.98892e33      // solar mass in grams
#define R_SUN           6.955e10        // solar radius in cm
#define L_SUN           3.827e33        // solar luminosity in  erg s-1
#define M_EARTH         5.98e27         // mass of earth in grams
#define R_EARTH         6.3781e8        // radius of earth in cm
#define R_JUPITER       7.1492e9        // radius of jupiter in cm
#define ABS_MAG_SUN     4.72            // absolute bolometric mag of sun
#define PARSEC          3.08568025e18   // parsec in cm 
#define ASTRO_UNIT      1.49598e13      // astronomical unit

//useful conversions
#define MEV_TO_ERGS    1.60217646e-6            
#define EV_TO_ERGS     1.60217646e-12            
#define ERGS_TO_EV     6.24150974e11
#define DAY_TO_SEC     8.64e4
#define YEAR_TO_SEC  2.2896e7
#define KM_TO_CM        1.0e5
#define ANGS_TO_CM      1.0e-8 
#define CM_TO_ANGS      1.0e8 

// helpful stuff
#define VERY_LARGE_NUMBER 9e99

#endif
