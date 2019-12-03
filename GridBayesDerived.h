#ifndef __GRIDBAYESDERIVED__H
#define __GRIDBAYESDERIVED__H

#include "GridBayes.h"

/*==============================================================================
  DERIVED CLASS DEFINITIONS
==============================================================================*/

// 3 states Lorenz equation ----------------------------------------------------

class Lorenz3D : public GridBayes {

public:
  // Three states Lorenz equation

  // Constructors	and destructors, initialize const values
  Lorenz3D(
    // Base class constant values.
    int n_pmax = 100000
    , int n_upd = 1
    , double t_max = 1
    , double dt = 5e-4
    , double pdf_thresh = 2e-5
    , vector<double> dx = vector<double>(3, 0.4)
    , vector<double> x0 = { -11.5, -10, 9.5 }
    , vector<double> x0_stdev = { 0.5, 0.5, 0.5 }
    , int dn = 5
  )
    : GridBayes(3, n_pmax, n_upd, t_max, dt, pdf_thresh, dx, x0, x0_stdev, dn)
  {}

  // Constant variables (Lorenz model parameters)
  double SIGMA = 4;
  double B = 1;
  double R = 48;

  // Virtual functions
  void model_f(double* f, const double* x, const double* xh)
  {
    // dx/dt = f(x,t) model : 3 states Lorenz equation
    f[0] = SIGMA * (x[1] - (x[0] + xh[0]));
    f[1] = -(x[1] + xh[1]) - x[0] * x[2];
    f[2] = -B * (x[2] + xh[2]) + x[0] * x[1] - B * R;
  }
  double measurement_model_pdf(int id)
  {
    // Measurement model 
    // Returns the measurement pdf for the list member [id]
    
    // Measure z axis only
    double z_meas = xs[2]; // simulated z state value
    double stdev = 0.5;   // measurement standard deviation

    // Current position in the list member [id]
    double z_list = list.pos[id][2] * DX[2];

    // Determine the pdf value based on Gaussian white noise model.
    double pdf = exp( -pow(z_list - z_meas, 2) / 2.0 / stdev / stdev );

    return pdf;
  }
  double flux_limiter(double x)
  {
    // Monotonized Central Difference (MC)
    // min([ (1+x)/2, 2, 2*x])
    double min = (1 + x) / 2.0;
    if (min > 2) min = 2;
    if (min > 2 * x) min = 2 * x;
    if (min > 0) return min;
    else return 0;

    // van Leer
    // return (x + abs(x)) / (1.0 + abs(x)); 
  }

private:

};


// 4 states Lorenz equation ----------------------------------------------------

class Lorenz4D : public GridBayes {

public:
  // Four states Lorenz equation

  // Constructors	and destructors, initialize const values
  Lorenz4D(
    // Base class constant values.
    int n_pmax = 1000000
    , int n_upd = 1
    , double t_max = 0.00075 * 2000
    , double dt = 0.00075
    , double pdf_thresh = 5e-6   // high tol: 5e-5, low tol: 5e-6
    , vector<double> dx = { 0.4, 0.4, 0.4, 0.01 }
    , vector<double> x0 = { -3.5, -3.5, -4, -0.45 }
    , vector<double> x0_stdev = { 0.5, 0.5, 0.5, 0.01 }
    , int dn = 5
  )
    : GridBayes(4, n_pmax, n_upd, t_max, dt, pdf_thresh, dx, x0, x0_stdev, dn)
  {}

  // Constant variables (4D Lorenz model parameters)
  double A = 5;
  double B = 20;
  double C = 1;
  double D = 1;
  double E = 20.6;
  double H = 1;
  double K = 0.1;

  // Virtual functions
  void model_f(double* f, const double* x, const double* xh)
  {
    // dx/dt = f(x,t) model : 4 states Lorenz equation
    // dx/dt = a(y-x) - e w
    // dy/dt = x z - h y
    // dz/dt = b - x y -c z
    // dw/dt = k y - d w

    f[0] = A * (x[1] - (x[0] + xh[0])) - E * x[3];
    f[1] = x[0] * x[2] - H * (x[1] + xh[1]);
    f[2] = B - x[0] * x[1] - C * (x[2] + xh[2]);
    f[3] = K * x[1] - D * (x[3] + xh[3]);
  }
  double measurement_model_pdf(int id)
  {
    // Measurement model 
    // Returns the measurement pdf for the list member [id]
    
    // Measure z axis only
    double z_meas = xs[2]; // simulated z state value
    double stdev = 0.5;   // measurement standard deviation

    // Current position in the list member [id]
    double z_list = list.pos[id][2] * DX[2];

    // Determine the pdf value based on Gaussian white noise model.
    double pdf = exp( -pow(z_list - z_meas, 2) / 2.0 / stdev / stdev );

    return pdf;
  }
  double flux_limiter(double x)
  {
    // Use either model (MC or van Leer)

    // Monotonized Central Difference (MC): min([ (1+x)/2, 2, 2*x])
    double min = (1 + x) / 2.0;
    if (min > 2) min = 2;
    if (min > 2 * x) min = 2 * x;
    if (min > 0) return min;
    else return 0;

    // van Leer
    // return (x + abs(x)) / (1.0 + abs(x));
  }

private:

};

// 6 states planetary model ----------------------------------------------------

class Planetary : public GridBayes {

public:
  // Six states planetary equation:  x = [a, e, i, o, w, t0]
  // a = semi-major axis radius (m)
  // e = eccentricity (rad)
  // i = inclination angle (rad)
  // o = right ascention of ascending node (rad)
  // w = argument of perigee (rad)
  // t0 = time since perigee passage (s)

  // Constructors	and destructors, initialize const values
  Planetary(
    // Base class constant values.
    int n_pmax = 1000000
    , int n_upd = 1
    , double t_max = 1 * 24 * 3600    // 1 day (in seconds)
    , double dt = 60                  // 1 seconds
    , double pdf_thresh = 1e-4
    , vector<double> dx = { 100,   // semi-major axis (m)
                            0.001,   // Eccentricity
                            0.001,   // Inclination angle (rad)
                            0.001,   // RAAN (rad)
                            0.001,   // Argument of perigee (rad) 
                            1 }     // Time of perigee passage (s)
    , vector<double> x0 = { (2000 + 6378) * 1000, // 2000 km altitude (add earth radius)
                            0.1,
                            0.5,
                            0.5,
                            0.5,
                            0 }
    , vector<double> x0_stdev = { 1000,
                                  0.01,
                                  0.01,
                                  0.01,
                                  0.01,
                                  1 }
    , int dn = 2
  )
    : GridBayes(6, n_pmax, n_upd, t_max, dt, pdf_thresh, dx, x0, x0_stdev, dn)
  {
    // Initialize atm data in the constructor
    read_atm_data();

    // Initialize other arrays
    xs_cartesian = new double[6];
    xl_cartesian = new double[6];
    xl_planetary = new double[6];
  }

  ~Planetary()
  {
    // Free dynamic memory (atm data)
    delete[] ATM_DATA_ALTITUDE;
    delete[] ATM_DATA_DENSITY;
    delete[] xs_cartesian;
    delete[] xl_cartesian;
    delete[] xl_planetary;
  }

  // Constant variables (Earth)
  double MU = 3.986004418e14;   // Earth std gravitational parameter (m^3 s^-2)
  double R_E = 6378.13649e3;    // Earth radius (m)
  double J2 = 1.082e-3;         // Earth symmetry about its polar axis

  // Constant variables (space junk)
  double SJ_M = 1;              // mass (kg)
  double SJ_S = 1;              // aerodynamic surface (m^2)
  double SJ_CD = 0.47;          // drag coefficient

  // Constant variables (atmospheric data)
  double* ATM_DATA_ALTITUDE;    // Altitude data (km)
  double* ATM_DATA_DENSITY;     // Density data (kg/m^3)
  int ATM_DATA_SIZE = 0;        // Atm density data size

  // Other variables
  double* xs_cartesian;         // Simulated [pos; vel] vector (m,m/s)
  double* xl_cartesian;         // List member's cartesian [pos;vel] vector
  double* xl_planetary;         // List member's orbital elements

  // Virtual functions
  void model_f(double* f, const double* x, const double* xh);
  double measurement_model_pdf(int id);
  double flux_limiter(double x)
  {
    // Use either model (MC or van Leer)

    // Monotonized Central Difference (MC): min([ (1+x)/2, 2, 2*x])
    double min = (1 + x) / 2.0;
    if (min > 2) min = 2;
    if (min > 2 * x) min = 2 * x;
    if (min > 0) return min;
    else return 0;

    // van Leer
    // return (x + abs(x)) / (1.0 + abs(x));
  }

  // Other functions
  void planetary2cartesian(double* xc, const double* xp, double t);
  void read_atm_data();
  double interpolate_atm_data(double alt);

private:

};

#endif // __GRIDBAYESDERIVED__H