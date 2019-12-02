#include "GridBayesDerived.h"

// 6 states planetary model ----------------------------------------------------

void Planetary::model_f(double* f, const double* x, const double* xh) {

  double n = sqrt(MU / x[0] / x[0] / x[0]);   // Mean motion
  double p = x[0] * (1 - x[1] * x[1]);        // Semi-latus rectum

  // Base planetary equation of motion (no perturbation)
  f[0] = 0;
  f[1] = 0;
  f[2] = 0;
  f[3] = 0;
  f[4] = 0;
  f[5] = n; // dM/dt, temporary

  // J2 perturbation
  f[3] += -3.0 / 2.0 * n * J2 * (R_E / p) * (R_E / p) * cos(x[2]);
  f[4] += -3.0 / 4.0 * n * J2 * (R_E / p) * (R_E / p) *
    (4 - 5 * sin(x[2]) * sin(x[2]));
  f[5] += 3.0 / 4.0 * n * J2 * (R_E / x[0]) * (R_E / x[0]) /
    pow(1 - x[1] * x[1], 1.5) * (2 - 3 * sin(x[2]) * sin(x[2]));

  // Atmospheric drag perturbation

  // M = mean anomaly. Limit M to +- pi
  double M = n * (t_sim - x[5]);
  while (M > PI) M -= 2 * PI;
  while (M < -PI) M += 2 * PI;

  // Solve for eccentric anomaly angle (E) from M using Newton-Raphson
  double tol = 1e-6;
  double E = 0;
  double ff = -M;
  double dff = 1 - x[1];
  while (abs(ff) > tol) {
    E = E - ff / dff;
    ff = E - x[1] * sin(E) - M;
    dff = 1 - x[1] * cos(E);
  }

  // Calculate true anomaly angle (theta)
  double theta = 2.0 * atan(sqrt((1 + x[1]) / (1 - x[1])) * tan(E / 2.0));

  // Satellite's linear distance and velocity
  double r = x[0] * (1 - x[1] * x[1]) / (1 + x[1] * cos(theta));
  double v = sqrt(MU * (2.0 / r - 1.0 / x[0]));

  // Satellite's ballistic coefficient (or space junk)
  double Bc = SJ_CD * SJ_S / 2.0 / SJ_M;

  // Determine the atmospheric density (input: altitude in km)
  double rho = interpolate_atm_data((r - R_E) / 1000.0);

  // Calculate the perturbation
  f[0] += -2 * n * x[0] * x[0] * x[0] * sqrt(1 - x[1] * x[1]) / MU * rho * v * v * Bc;
  f[1] += -n * x[0] * x[0] * sqrt(1 - x[1] * x[1]) / MU
    * rho * v * v * Bc * (cos(theta) + cos(E));
  f[4] += -n * x[0] * x[0] * sqrt(1 - x[1] * x[1]) / MU / x[1]
    * rho * v * v * Bc * (1 + r / p) * sin(theta);

  // Last step: transform dM/dt into dt0/dt
  f[5] = 1 - f[5] / n - 1.5 / x[0] * (t_sim - x[5]) * f[0];

  
}

double Planetary::measurement_model_pdf(int id) {

  // Standard deviation of the measurements
  double stdev_pos = 1000;  // m
  double stdev_vel = 10;    // m/s

  // If the first element in the list, determine the simulated cartesian 
  // position and speed. This will be the measured position and velocity.
  if (id == 0) {
    planetary2cartesian(xs_cartesian, xs, t_sim);
    return 0; // The 1st element of the list is an empty node
  }

  // Determine the cartesian position and velocity if the list at [id].
  for (int it = 0; it < N_DIM; it++) {
    xl_planetary[it] = list.pos[id][it] * DX[it];
  }
  planetary2cartesian(xl_cartesian, xl_planetary, t_sim);

  // Calculate pdf based on the Gaussian noise model in Cartesian coordinates.
// The pdf will be normalized later so that their sum = 1, so no scaling here.
// Assume diagonal covariance matrix.
  double pdf = 0;

  // Position
  for (int it = 0; it < 3; it++) {
    pdf += exp(-pow(xs_cartesian[it] - xl_cartesian[it], 2)
      / 2.0 / stdev_pos / stdev_pos);
  }

  // Velocity
  for (int it = 3; it < 6; it++) {
    pdf += exp(-pow(xs_cartesian[it] - xl_cartesian[it], 2)
      / 2.0 / stdev_vel / stdev_vel);
  }

  return pdf;
}

void Planetary::planetary2cartesian(double* xc, const double* xp, double t) {
  // Outputs xc (Cartesian position and velocity)
  // a  = xp[0]: semi-major axis radius 
  // e  = xp[1]: eccentricity
  // i  = xp[2]: inclination angle
  // o  = xp[3]: right ascention of the ascending node (RAAN)
  // w  = xp[4]: argument of perigee
  // t0 = xp[5]: time of perigee passage

  // Calculate the mean anomaly angle (M), limit M to [-PI,PI]
  double M = (t - xp[5]) * sqrt(MU / pow(xp[0], 3));
  while (M > PI)  M -= 2 * PI;
  while (M < -PI) M += 2 * PI;

  // Calculate the eccentric anomaly angle (E) with Newton-Raphson
  double E = 0;
  double tol = 1e-6;
  double f_nr = E - xp[1] * sin(E) - M;
  double df_nr = 1 - xp[1] * cos(E);
  int counter_nr = 0; // just in case

  // Newton-Raphson iterative method
  while (abs(f_nr) < tol && counter_nr < 100) {
    E -= f_nr / df_nr;
    f_nr = E - xp[1] * sin(E) - M;
    df_nr = 1 - xp[1] * cos(E);
    counter_nr++;
  }

  // Calculate the true anomaly angle (theta)
  double theta = 2.0 * atan(sqrt((1 + xp[1]) / (1 - xp[1])) * tan(E / 2.0));

  // Calculate ECI from Cartesian position and velocity
  double A00 = cos(xp[3]) * cos(xp[4]) - sin(xp[3]) * cos(xp[2]) * sin(xp[4]);
  double A01 = -cos(xp[3]) * sin(xp[4]) - sin(xp[3]) * cos(xp[2]) * cos(xp[4]);
  double A02 = sin(xp[3]) * sin(xp[2]);
  double A10 = sin(xp[3]) * cos(xp[4]) + cos(xp[3]) * cos(xp[2]) * sin(xp[4]);
  double A11 = -sin(xp[3]) * sin(xp[4]) + cos(xp[3]) * cos(xp[2]) * cos(xp[4]);
  double A12 = -cos(xp[3]) * sin(xp[2]);
  double A20 = sin(xp[2]) * sin(xp[4]);
  double A21 = sin(xp[2]) * cos(xp[4]);
  double A22 = cos(xp[2]);

  // Cartesian position
  double temp = xp[0] * (1 - pow(xp[1], 2)) / (1 + xp[1] * cos(theta));
  xc[0] = (A00 * cos(theta) + A01 * sin(theta)) * temp;
  xc[1] = (A10 * cos(theta) + A11 * sin(theta)) * temp;
  xc[2] = (A20 * cos(theta) + A21 * sin(theta)) * temp;

  // Cartesian velocity
  temp = sqrt(MU / xp[0] / (1 - pow(xp[1], 2)));
  xc[3] = (-A00 * sin(theta) + A01 * (xp[1] + cos(theta))) * temp;
  xc[4] = (-A10 * sin(theta) + A11 * (xp[1] + cos(theta))) * temp;
  xc[5] = (-A20 * sin(theta) + A21 * (xp[1] + cos(theta))) * temp;
}

void Planetary::read_atm_data() {
  // Read atmospheric density function from "atm_data.txt"
  // Data format: 
  // 1st row = data size (int)
  // the remaining rows: [altitude (km), atm density (kg/m^3)]

  string line;
  ifstream myfile("atm_data.txt");

  // Read data size and initialize atm data arrays
  myfile >> ATM_DATA_SIZE;
  ATM_DATA_ALTITUDE = new double[ATM_DATA_SIZE];
  ATM_DATA_DENSITY = new double[ATM_DATA_SIZE];

  // Read altitude and density data
  for (int i = 0; i < ATM_DATA_SIZE; i++) {
    myfile >> ATM_DATA_ALTITUDE[i];
    myfile >> ATM_DATA_DENSITY[i];
  }

  myfile.close();
}

double Planetary::interpolate_atm_data(double alt) {

  // Check if the input value is outside the dataset
  if (alt <= ATM_DATA_ALTITUDE[0]) {
    return ATM_DATA_DENSITY[0];
  }
  if (alt >= ATM_DATA_ALTITUDE[ATM_DATA_SIZE - 1]) {
    return ATM_DATA_DENSITY[ATM_DATA_SIZE - 1];
  }

  // If within the dataset, use linear interpolation
  int x0 = floor(alt);
  double y0 = ATM_DATA_DENSITY[x0];
  double y1 = ATM_DATA_DENSITY[x0 + 1];

  return (alt - double(x0)) * (y1 - y0) + y0;
}