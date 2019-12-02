/*==============================================================================

MAIN

==============================================================================*/

#ifndef __GRIDBAYES__H
#define __GRIDBAYES__H

// Standard libraries
#include <vector>
#include <array>
#include <cmath>
#include <iostream>
#include <chrono>
#include <thread>
#include <string>
#include <fstream>

using namespace std;

// Installed libraries


// Gaussian Model 
#define GAUSS_STDEV 0.5
#define PI          3.14159265359

/*==============================================================================
  BASE CLASS DEFINITIONS
==============================================================================*/

class GridBayes {
public:
	// Constructor to initialize the constant member variables
  GridBayes(
    // Default values for the constant variables, usually not used.
    // These values must be specified when decalring the derived class.
    // The vector is a workaround to initialize the array using a constructor.
      int n_dim = 1
    , int n_pmax = 1
    , int n_upd = 1
    , double t_max = 1
    , double dt = 1
    , double pdf_thresh = 1
    , vector<double> dx = { 1 }
		, vector<double> x0 = { 1 }
    , vector<double> x0_stdev = { 1 }
		, int dn = 1
	)
		: N_DIM(n_dim)
		, N_PMAX(n_pmax)
		, N_UPD(n_upd)
		, T_MAX(t_max)
		, DT(dt)
		, PDF_THRESH(pdf_thresh)
		, DN(dn)
	{
    // Assign the initial array values. 
    DX = new double[n_dim];   // State grid size
    X0 = new double[n_dim];   // Initial state
    xs = new double[n_dim];   // Simulated state (actual state)
    X0_STDEV = new double[n_dim];
    for (int i = 0; i < n_dim; i++) 
    {
      DX[i] = dx[i];
      X0[i] = x0[i];
      xs[i] = x0[i];
      X0_STDEV[i] = x0_stdev[i];
    }

    // Initialize list DataStruct (row)
    list.pdf = new double[N_PMAX];  // The only 1D member
    list.pos = new int* [N_PMAX];
    list.nb = new int* [N_PMAX];
    list.ns = new int* [N_PMAX];
    list.v = new double* [N_PMAX];
    list.vp = new double* [N_PMAX];
    list.vn = new double* [N_PMAX];
    list.f = new double* [N_PMAX];

    for (int i = 0; i < N_PMAX; i++) { // column
      list.pos[i] = new int[N_DIM];
      list.nb[i] = new int[N_DIM];
      list.ns[i] = new int[N_DIM];
      list.v[i] = new double[N_DIM];
      list.vp[i] = new double[N_DIM];
      list.vn[i] = new double[N_DIM];
      list.f[i] = new double[N_DIM];
    }

    // Fill in with the default values (0)
    for (int i = 0; i < N_PMAX; i++) {  // column
      reset_list_entry(i);              // Row
    }

    // Dynamic memories for the temporary variables in a function:

    // In RK$

    // Initialize global variable values
    t_sim = 0;
    nn = 0;
    mm = 0;
  }

  ~GridBayes() 
  {
    // Free dynamic memory
    delete[] DX;
    delete[] X0;
    delete[] xs;
    delete[] X0_STDEV;

    for (int i = 0; i < N_PMAX; i++) 
    {
      delete[] list.pos[i];
      delete[] list.nb[i];
      delete[] list.ns[i];
      delete[] list.v[i];
      delete[] list.vp[i];
      delete[] list.vn[i];
      delete[] list.f[i];
    }

    delete[] list.pdf;
    delete[] list.pos;
    delete[] list.nb;
    delete[] list.ns;
    delete[] list.v;
    delete[] list.vp;
    delete[] list.vn;
    delete[] list.f;
  }

  // Variables and structs -----------------------------------------------------

	// Constant member variables, array not defined with const for simplicity. 
	const int N_DIM;			    // State dimension
  const int N_PMAX;			    // Point list initial length
  const int N_UPD;		    	// Measurement update every N_UPD steps
  const double T_MAX;	    	// Simulation time length
  const double DT;		    	// Simulation time step (for RK4)
  const double PDF_THRESH;	// PDF treshold for removal from the list
	double* DX;	              // Grid width, vector of length N_DIM
  double* X0;	              // Starting position
  const int DN;             // Initial PDF up to DN grids away from X0
  double* X0_STDEV;	        // Gaussian stdev for initial grid list pdf values

	// Point list data struct, stored in arrays for maximum speed
	struct DataStruct 
  {
    double* pdf;  // Point PDF value
    int** pos;	  // Position vector (int)
    int** nb;		  // Big neighbor index
    int** ns;		  // Small neighbor index
		double** v;	  // Process model shifted?
    double** vp;	// v+
    double** vn;	// v-
    double** f;	  // Flux
  };
  
  // Member Variables (global variables in the class)
  DataStruct list;  // List struct
  int nn;           // Total elements
  int mm;           // Index to the latest large pdf value in the list
  double* xs;       // Simulated state (true state value)
  double t_sim;     // Simulation time

  // Functions -----------------------------------------------------------------

  // Primary Functions
  void initialize();
  void initialize_vn(int b);
  void modify_list();
  void march_RK4(double* x);
  void march_PDF();
  void measurement_update();

  // Sub Functions
  int neighbor_check(const int* x1, const int* x2);
  void reset_list_entry(int i);
  void record_data(string file_name);
  void swap_list_entries(int k1, int k2);
  void create_new_list_entry(int* pos);
  void check_init();

  // Math functions


  // Virtual functions (make a derived class to use different models).
  // Made to support model of any dimensions.
  virtual void model_f(double* f, const double* x, const double* xh) {  }
  virtual double measurement_model_pdf(int id) { return 0; }
  virtual double flux_limiter(double x) { return 0; }

private:
	// Variables
	


	// Functions

};






#endif // __GRIDBAYES__H