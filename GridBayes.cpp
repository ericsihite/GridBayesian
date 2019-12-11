
/*==============================================================================

  MAIN

==============================================================================*/



#include "GridBayes.h"

/*==============================================================================
	BASE CLASS FUNCTIONS
==============================================================================*/

void GridBayes::initialize() {

  int n_list = 1; // List size
  
  // Determine the initial position grid points.
  // ki = [k_begin, k_end, total_length]
  int** ki = new int*[N_DIM];
  for (int i = 0; i < N_DIM; i++) {

    ki[i] = new int[3];
  }
  
  // Populate the grid around X0 +- DN*DX
  for (int i = 0; i < N_DIM; i++) {

    ki[i][0] = (int)round(X0[i] / DX[i] - DN);
    ki[i][1] = (int)round(X0[i] / DX[i] + DN);
    ki[i][2] = ki[i][1] - ki[i][0] + 1;
  }

  // cout << ki[0][2] << "," << ki[1][2] << "," << ki[2][2] << endl;

  int* kk = new int[N_DIM];
  for (int i = 0; i < N_DIM; i++) {
    kk[i] = 0;
  }

  double temp_d = 0;
  
  // White Gaussian PDF denumerator
  double pdf_den = X0_STDEV[0]*X0_STDEV[0];
  for (int d = 1; d < N_DIM; d++) {
    pdf_den *= X0_STDEV[d]*X0_STDEV[d];
  }
  pdf_den = sqrt(pow(2 * 3.14159265359, N_DIM) * pdf_den);
  
  // Iterate through the entire dimensions.
  // Equivalent to nested for loops for each dimension.
  while (kk[0] != ki[0][2]) {

    // Generate the point values (pdf and position)
    temp_d = 0;
    for (int i = 0; i < N_DIM; i++) {
      // Multivariable Gaussian model (stdev = 0.5)
      temp_d += pow( (ki[i][0] + kk[i]) * DX[i] - X0[i], 2) 
                  / 2.0 / X0_STDEV[i] / X0_STDEV[i];
      
      // Grid position vector (int)
      list.pos[n_list][i] = ki[i][0] + kk[i];
    }
    
    list.pdf[n_list] = exp(-temp_d) / pdf_den;

    n_list++; // Update list counter.

    // Update iterator, starting from the end.
    kk[N_DIM - 1]++;

    if (kk[N_DIM - 1] >= ki[N_DIM - 1][2]) {

      kk[N_DIM - 2]++;
      kk[N_DIM - 1] = 0;

      for (int i = N_DIM - 2; i > 0; i--) {

        if (kk[i] >= ki[i][2]) {
          kk[i - 1]++;
          kk[i] = 0;
        }
      }
    }
  }

  // cout << "List length: " << n_list << endl;

  nn = n_list;    // List entry counter
  mm = n_list-1;  // Latest nontrivial PDF index

  initialize_vn(1);

  // Free the dynamic memory.
  for (int i = 0; i < N_DIM; i++) delete[] ki[i];
  delete[] ki;
  delete[] kk;
}

// Initialize neighbors and process model (nb,ns,v,v+,v-)
void GridBayes::initialize_vn(int b) {

  double* x_pos = new double[N_DIM];

  // Half step forward dx
  double* xh = new double[N_DIM];
  for (int i = 0; i < N_DIM; i++) {
    xh[i] = DX[i] / 2.0;
  }

  // Initialize v, v+, v-
  for (int l = b; l < nn; l++) {

    // v = process model, half step forward
    for (int i = 0; i < N_DIM; i++) {
      x_pos[i] = list.pos[l][i] * DX[i];
    }
    model_f(list.v[l], x_pos, xh);

    // vp, vn = v+ = max(v,0), v- = min(v,0)
    for (int i = 0; i < N_DIM; i++) {
      // max(v,0)
      if (list.v[l][i] > 0) list.vp[l][i] = list.v[l][i];
      else list.vp[l][i] = 0;

      // min(v,0)
      if (list.v[l][i] < 0) list.vn[l][i] = list.v[l][i];
      else list.vn[l][i] = 0;
    }
  }

  // Find neighbor
  int ii = 0;
  for (int l = nn - 1; l >= b; l--) {// From the back
    
    for (int t = 0; t <= l - 1; t++) { // From the front
    
      ii = neighbor_check(list.pos[l], list.pos[t]);  // Neighbor axis
      
      if (ii >= 0) {
        if (  list.pos[l][ii] == list.pos[t][ii] + 1) {
          list.ns[l][ii] = t;
          list.nb[t][ii] = l;
        }
        else if (list.pos[t][ii] == list.pos[l][ii] + 1) {
          list.ns[t][ii] = l;
          list.nb[l][ii] = t;
        }
      }
    }
  }

  // Release memory
  delete[] x_pos;
  delete[] xh;
}


void GridBayes::modify_list() {

  mm = nn-1;                      // Latest large pdf index
  int ll = mm;                    // Temp list index
  int* x_temp = new int[N_DIM];   // temp grid position vector (int)
  int temp_int = 0;

  // Move big elements to the beginning, small elements to the end.
  
  while (list.pdf[ll] < PDF_THRESH) {
    ll--;
    mm--;
  }
  
  ll--; // The next big element. mm = latest big element.

  while (ll > 0) {
    if (list.pdf[ll] < PDF_THRESH) {
      // Found a small element. Swap with data in index mm.
      
      // Fix pointers for the neighbor data list
      for (int d = 0; d < N_DIM; d++) {

        // If switching position with its own neighbor:
        // Small neighbor
        if (list.ns[ll][d] == mm) {
          list.ns[ll][d] = list.ns[mm][d];
          list.ns[mm][d] = ll;
        }
        else if (list.ns[mm][d] == ll) {
          list.ns[mm][d] = list.ns[ll][d];
          list.ns[ll][d] = mm;
        }
        else {

          temp_int = list.ns[ll][d];
          list.ns[ll][d] = list.ns[mm][d];
          list.ns[mm][d] = temp_int;
        }

        // Big neighbor
        if (list.nb[ll][d] == mm) {
          list.nb[ll][d] = list.nb[mm][d];
          list.nb[mm][d] = ll;
        }
        else if (list.nb[mm][d] == ll) {
          list.nb[mm][d] = list.nb[ll][d];
          list.nb[ll][d] = mm;
        }
        else {

          temp_int = list.nb[ll][d];
          list.nb[ll][d] = list.nb[mm][d];
          list.nb[mm][d] = temp_int;
        }

        // Non neighbor switch
        list.ns[list.nb[ll][d]][d] = ll;
        list.ns[list.nb[mm][d]][d] = mm;
        list.nb[list.ns[ll][d]][d] = ll;
        list.nb[list.ns[mm][d]][d] = mm;
      }

      // Swap data other than nb and ns.
      swap_list_entries(ll, mm);
      mm--; 
    }

    ll--;
  }

  
  // Identify the neighbors to the big elements and create entry if neccessary.
  
  // Temporarily borrow first dim of f for neighbor flags
  for (int i = mm + 1; i < nn; i++) {
    list.f[i][0] = 0;   
  }

  // Check the nontrivial pdf for neighbors.
  for (ll = 1; ll <= mm; ll++) {

    for (int d = 0; d < N_DIM; d++) {

      // If no small neighbor, create one.
      if (list.ns[ll][d] == 0) {
        
        // Check the v direction
        if (list.vn[ll][d] < 0) {
          for (int j = 0; j < N_DIM; j++) {
            x_temp[j] = list.pos[ll][j];
          }
          x_temp[d]--;
          create_new_list_entry(x_temp);
        }

      }
      else {
        // 1 indicates the pointer has a big neighbor
        // Check for the upwind flow as well.

        


        list.f[list.ns[ll][d]][0] = 1;
      }

      // If no big neighbor, create one.
      if (list.nb[ll][d] == 0) {
        
        // Check the v direction
        if (list.vp[ll][d] > 0) {
          for (int j = 0; j < N_DIM; j++) {
            x_temp[j] = list.pos[ll][j];
          }
          x_temp[d]++;
          create_new_list_entry(x_temp);
        }

      }
      else {
        // 1 indicates the pointer has a small neighbor
        list.f[list.nb[ll][d]][0] = 1;
      }

      // Compute corner entries and create one if doesn't exist.
      for (int e = d; e < N_DIM; e++) {

        // Small-small corner
        if (  list.ns[list.ns[ll][e]][d] == 0) {

          if (list.vn[list.ns[ll][e]][d] < 0) {
            for (int j = 0; j < N_DIM; j++) {
              x_temp[j] = list.pos[ll][j];
            }
            x_temp[d]--;
            x_temp[e]--;
            create_new_list_entry(x_temp);
          }

        }
        else {
          list.f[ list.ns[ list.ns[ll][e] ][d] ][0] = 1;
        }

        // Small-big corner
        if (  list.ns[list.nb[ll][e]][d] == 0) {
          
          if (list.vn[list.nb[ll][e]][d] < 0) {
            for (int j = 0; j < N_DIM; j++) {
              x_temp[j] = list.pos[ll][j];
            }
            x_temp[d]--;
            x_temp[e]++;
            create_new_list_entry(x_temp);
          }
        }
        else {
          list.f[list.ns[list.nb[ll][e]][d]][0] = 1;
        }

        // Big-small corner
        if (  list.nb[list.ns[ll][e]][d] == 0 ) {

          if (list.vp[list.ns[ll][e]][d] > 0) {
            for (int j = 0; j < N_DIM; j++) {
              x_temp[j] = list.pos[ll][j];
            }
            x_temp[d]++;
            x_temp[e]--;
            create_new_list_entry(x_temp);
          }
        }
        else {
          list.f[list.nb[list.ns[ll][e]][d]][0] = 1;
        }

        // Big-big corner
        if (  list.nb[list.nb[ll][e]][d] == 0) {

          if (list.vp[list.nb[ll][e]][d] > 0) {
            for (int j = 0; j < N_DIM; j++) {
              x_temp[j] = list.pos[ll][j];
            }
            x_temp[d]++;
            x_temp[e]++;
            create_new_list_entry(x_temp);
          }
        }
        else {
          list.f[list.nb[list.nb[ll][e]][d]][0] = 1;
        }

      }
    }
  }

  // Delete small elements not neighboring the big elements
  // and also not in the upwind direction.

  ll = mm + 1;
  bool is_zero = false;

  while (ll < nn) {

    if (list.f[ll][0] < 1) {
      // Not neighboring large element.
      // Swap with the last element and then delete entry.

      // Update pointer
      for (int d = 0; d < N_DIM; d++) {

        list.ns[list.nb[ll][d]][d] = 0;
        if (ll < nn - 1) {
          list.ns[list.nb[nn-1][d]][d] = ll;
        }

        list.nb[list.ns[ll][d]][d] = 0;
        if (ll < nn - 1) {
          list.nb[list.ns[nn-1][d]][d] = ll;
        }

        // Check if the deleted entry has a small neighbor
        for (int j = 0; j < N_DIM; j++) {
          x_temp[j] = list.pos[nn-1][j] - list.pos[ll][j];
        }
        x_temp[d]--;
        is_zero = true;
        for (int j = 0; j < N_DIM; j++) {
          if (x_temp[j] != 0) is_zero = false;
        }
        if (is_zero)  list.ns[ll][d] = 0;         
        else          list.ns[ll][d] = list.ns[nn-1][d];

        // Check if the deleted entry has a big neighbor
        for (int j = 0; j < N_DIM; j++) {
          x_temp[j] = list.pos[nn - 1][j] - list.pos[ll][j];
        }
        x_temp[d]++;
        is_zero = true;
        for (int j = 0; j < N_DIM; j++) {
          if (x_temp[j] != 0) is_zero = false;
        }
        if (is_zero)  list.nb[ll][d] = 0;
        else          list.nb[ll][d] = list.nb[nn-1][d]; 

      }

      // Swap and delete the last entry
      swap_list_entries(ll, nn - 1);
      reset_list_entry(nn - 1);
      nn--;
    }
    else {
      ll++;
    }
  }

  // Normalize PDF such that the sum = 1.
  double pdf_sum = 0;
  for (ll = 1; ll < nn; ll++) {
    if (list.pdf[ll] < 0) list.pdf[ll] = 0;
    pdf_sum += list.pdf[ll];
  }
  for (ll = 1; ll < nn; ll++) {
    list.pdf[ll] = list.pdf[ll] / pdf_sum;
  }

  // Free memory
  delete[] x_temp;

  return;
}


void GridBayes::march_RK4(double* x) {

  // Initialize memory
  double* x_temp = new double[N_DIM];
  double* xh = new double[N_DIM];
  double* f1 = new double[N_DIM];
  double* f2 = new double[N_DIM];
  double* f3 = new double[N_DIM];
  double* f4 = new double[N_DIM];

  // f1
  for (int i = 0; i < N_DIM; i++) {
    x_temp[i] = x[i];
    xh[i] = 0;
  }
  model_f(f1, x_temp, xh);
  
  // f2
  for (int i = 0; i < N_DIM; i++) {
    x_temp[i] = x[i] + f1[i] * DT / 2.0;
  }
  model_f(f2, x_temp, xh);

  // f3
  for (int i = 0; i < N_DIM; i++) {
    x_temp[i] = x[i] + f2[i] * DT / 2.0;
  }
  model_f(f3, x_temp, xh);

  // f4
  for (int i = 0; i < N_DIM; i++) {
    x_temp[i] = x[i] + f3[i] * DT;
  }
  model_f(f4, x_temp, xh);

  // Time march RK4 : x[k+1] = x[k] + (f1/6 + f2/3 + f3/3 + f4/6)*dt
  for (int i = 0; i < N_DIM; i++) {
    x[i] += (f1[i] / 6.0 + f2[i] / 3.0 + f3[i] / 3.0 + f4[i] / 6.0) * DT;
  }

  t_sim += DT;  // Time stamp


  // Release memory
  delete[] x_temp;
  delete[] xh;
  delete[] f1;
  delete[] f2;
  delete[] f3;
  delete[] f4;

  return;
}

void GridBayes::march_PDF() {

  int ll = 0;
  int n1 = 0;
  int n2 = 0;
  int n3 = 0;
  int n4 = 0;
  double flux1 = 0;
  double flux2 = 0;
  double flux3 = 0;
  double temp = 0;
  double theta = 0;
  double* rhs = new double[nn];
  rhs[0] = 0;
  

  /*

  // Calculate the initial flux terms, 1 dim
  for (ll = 1; ll < nn; ll++) {
    for (int d = 0; d < N_DIM; d++) {
      list.f[ll][d] = list.vp[ll][d] * list.pdf[ll] + 
                      list.vn[ll][d] * list.pdf[list.nb[ll][d]];
    }
    rhs[ll] = 0;
  }

  // Corner Transport Upwind (CTU), 2 dim
  for (int d = 0; d < N_DIM; d++) {
    for (ll = 1; ll < nn; ll++) {
      
      i = list.ns[ll][d]; // Small neighbor

      // Big pdf entries only
      if (ll <= mm || (i > 0 && i <= mm)) {

        flux = DT * (list.pdf[ll] - list.pdf[i]) / (2 * DX[d]);
        
        // Compute CTU flux terms
        for (int e = 0; e < N_DIM; e++) {
          if (e != d) {
            list.f[ll][e] -= list.vp[ll][e] * list.vp[i][d] * flux;
            
            j = list.ns[ll][e];
            list.f[j][e] -= list.vn[j][e] * list.vp[i][d] * flux;
            list.f[i][e] -= list.vp[i][e] * list.vn[i][d] * flux;

            j = list.ns[i][e];
            list.f[j][e] -= list.vn[j][e] * list.vn[i][d] * flux;
          }
        }

        // Compute 2nd order correction flux term
        if (list.v[i][d] > 0) {
          theta = (list.pdf[i] - list.pdf[list.ns[i][d]]) / 
                  (list.pdf[ll] - list.pdf[i]);
        }
        else {
          theta = (list.pdf[list.nb[ll][d]] - list.pdf[ll]) /
                  (list.pdf[ll] - list.pdf[i]);
        }

        // Flux limiter function. MC or VL.
        temp = abs(list.v[i][d]);
        list.f[i][d] += temp*(DX[d] / DT - temp)*flux*flux_limiter(theta);
      }

    }
  }

  // Calculate RHS of the PDF equation
  for (ll = 1; ll < nn; ll++) {
    for (int d = 0; d < N_DIM; d++) {
      rhs[ll] -= (list.f[ll][d] - list.f[ list.ns[ll][d] ][d]) / DX[d];
    }
  }

  */

  // Calculate the initial flux terms, 1 dim
  for (ll = 1; ll < nn; ll++) {
    for (int d = 0; d < N_DIM; d++) {
      list.f[ll][d] = list.vp[ll][d] * list.pdf[ll] +
        list.vn[ll][d] * list.pdf[ list.nb[ll][d] ];
    }
    rhs[ll] = 0;
  }

  // Corner Transport Upwind (CTU)
  for (int d1 = 0; d1 < N_DIM; d1++) {
    for (ll = 1; ll < nn; ll++) {

      n1 = list.ns[ll][d1]; // Small neighbor

      // Big pdf entries only
      if (ll <= mm || (n1 > 0 && n1 <= mm)) {

        // Compute CTU flux terms, 2 dim corners--------------------------------

        // Half flux across d1
        flux1 = DT * (list.pdf[ll] - list.pdf[n1]) / (2.0 * DX[d1]);

        for (int d2 = 0; d2 < N_DIM; d2++) {
          if (d2 != d1) {


            n2 = list.ns[ll][d2];

            // pp
            list.f[ll][d2] -= list.vp[ll][d2] * list.vp[n1][d1] * flux1;

            // np
            list.f[n2][d2] -= list.vn[n2][d2] * list.vp[n1][d1] * flux1;

            n2 = list.ns[n1][d2];

            // pn
            list.f[n1][d2] -= list.vp[n1][d2] * list.vn[n1][d1] * flux1;

            // nn
            list.f[n2][d2] -= list.vn[n2][d2] * list.vn[n1][d1] * flux1;
          }
        }

        // Compute 2nd order correction flux term
        if (list.v[n1][d1] > 0) {
          theta = (list.pdf[n1] - list.pdf[list.ns[n1][d1]]) /
            (list.pdf[ll] - list.pdf[n1]);
        }
        else {
          theta = (list.pdf[list.nb[ll][d1]] - list.pdf[ll]) /
            (list.pdf[ll] - list.pdf[n1]);
        }

        // Flux limiter function. MC or VL.
        temp = abs(list.v[n1][d1]);
        list.f[n1][d1] += temp * (DX[d1] / DT - temp)
          * flux1 * flux_limiter(theta);

        // Compute CTU flux terms, 3 dim corners -------------------------------

        if (0) {

          for (int d2 = 0; d2 < N_DIM; d2++) {

            n2 = list.ns[ll][d2];
            n3 = list.ns[n2][d1];

            // Half flux across d1 & d2
            flux1 = list.pdf[ll] - list.pdf[n1];
            flux2 = list.pdf[n2] - list.pdf[n3];
            flux3 = (flux1 - flux2) * DT * DT / (3.0 * DX[d1] * DX[d2]);

            for (int d3 = 0; d3 < N_DIM; d3++) {

              if (d1 != d2 && d1 != d3 && d2 != d3) { // All unique axis

                //  n1 -- ll
                //   |    |
                //  n3 -- n2
                //  n4 = 3rd dimension points

                // ------------------
                n4 = list.ns[ll][d3]; // ll

                // ppp 
                list.f[ll][d3] +=
                  list.vp[ll][d3] * list.vp[n2][d2] * list.vp[n1][d1] * flux3;

                // npp
                list.f[n4][d3] +=
                  list.vn[n4][d3] * list.vp[n2][d2] * list.vp[n1][d1] * flux3;

                // ------------------
                n4 = list.ns[n1][d3]; // n1

                // ppn
                list.f[n1][d3] +=
                  list.vp[n1][d3] * list.vp[n3][d2] * list.vn[n1][d1] * flux3;

                // npn
                list.f[n4][d3] +=
                  list.vn[n4][d3] * list.vp[n3][d2] * list.vn[n1][d1] * flux3;

                // ------------------              
                n4 = list.ns[n2][d3]; // n2

                // pnp
                list.f[n2][d3] +=
                  list.vp[n2][d3] * list.vn[n2][d2] * list.vp[n3][d1] * flux3;

                // nnp
                list.f[n4][d3] +=
                  list.vn[n4][d3] * list.vn[n2][d2] * list.vp[n3][d1] * flux3;


                // ------------------   
                n4 = list.ns[n3][d3]; // n3

                // pnn
                list.f[n3][d3] +=
                  list.vp[n3][d3] * list.vn[n3][d2] * list.vn[n3][d1] * flux3;

                // nnn
                list.f[n4][d3] +=
                  list.vn[n4][d3] * list.vn[n3][d2] * list.vn[n3][d1] * flux3;

              }
            } // end for d3

            // Flux correction? q_xxy, etc?



          } // end for d2

          // Flux correction? q_xxx, etc?


        } // if(1/0)

      } // active cell only
    } // end for ll
  } // end for d1

  // Calculate RHS of the PDF equation
  for (ll = 1; ll < nn; ll++) {
    for (int d = 0; d < N_DIM; d++) {
      rhs[ll] -= (list.f[ll][d] - list.f[ list.ns[ll][d] ][d]) / DX[d];
    }
  }

  // Time march
  double sum = 0; // PDF sum for normalization
  for (ll = 1; ll < nn; ll++) {
    list.pdf[ll] += DT * rhs[ll];
    sum += list.pdf[ll];
  }

  // Normalize PDF
  for (ll = 1; ll < nn; ll++) {
    list.pdf[ll] /= sum;
  }

  // Free memory
  delete[] rhs;

  return;
}

void GridBayes::measurement_update() {

  // Measurements
  //double h[1] = { 0 };
  //model_h(h, xs);

  // Update pdf
  double meas_pdf = 0;
  double sum = 0;

  // Calculate new pdf (measurement model update) and sum the pdf
  for (int i = 0; i < nn; i++) {
    // meas_pdf = exp( -pow(list.pos[i][2] * DX[2] - xs[2], 2) 
    //                / 2.0 / STDEV[2] / STDEV[2]);
    
    //cout << i << endl;

    meas_pdf = measurement_model_pdf(i);
    list.pdf[i] *= meas_pdf;  // Multiply grid pdf with measurement pdf
    sum += list.pdf[i];
  }

  // Normalize PDF (sum of all pdf in the list = 1)
  for (int i = 1; i < nn; i++) {
    list.pdf[i] /= sum;
  }

}

// Sub Functions ---------------------------------------------------------------

int GridBayes::neighbor_check(const int* x1, const int* x2) {

  int d = 0;
  int n = 0;

  for (int i = 0; i < N_DIM; i++) {

    // Potential neighbor (dx = 1)
    if (abs(x1[i] - x2[i]) == 0) {
      n++;
    }
    else {
      d = i;
    }
  }

  if (n == N_DIM - 1) return d;
  else return -1;
}


// Reset the node entry to the empty/default value
void GridBayes::reset_list_entry(int i) {

  list.pdf[i] = 0;
  for (int j = 0; j < N_DIM; j++)   {
    list.pos[i][j] = 0;
    list.nb[i][j] = 0;
    list.ns[i][j] = 0;
    list.v[i][j] = 0;
    list.vp[i][j] = 0;
    list.vn[i][j] = 0;
    list.f[i][j] = 0;
  }
}

void GridBayes::record_data(string file_name) {
	// Record the data into a text file.

  // The first set of data has the following structure:
  // list.pos shows the list.nn instead
  // list.nb shows the list.mm instead.


	ofstream myfile;
	myfile.open(file_name);

	myfile << "GBEES Data." << endl;

	for (int i = 0; i < mm; i++) {
		
		// PDF
    if (i == 0) myfile << t_sim;  // 1st line = time stamp
    else myfile << list.pdf[i];   // pdf
		
		// Position
		for (int j = 0; j < N_DIM; j++) {

      if (i == 0) {
        // 1st line: record misc data
        if (j == 0) myfile << "," << nn;  // total list size
        else myfile << "," << mm;  // active list size
      }
      else 			    myfile << "," << list.pos[i][j];  // state value

		}

		// Neighbors: nb, ns
    /*
		for (int j = 0; j < N_DIM; j++) {

      if (i == 0)   myfile << "," << mm; // 1st line = major pdf list size
      else        	myfile << "," << list.nb[i][j]; // big neighbor

		}

		for (int j = 0; j < N_DIM; j++) {
			myfile << "," << list.ns[i][j]; // small neighbor
		}

		// Process model: v, v+, v-
		for (int j = 0; j < N_DIM; j++) {
			myfile << "," << list.v[i][j]; 
		}
		for (int j = 0; j < N_DIM; j++) {
			myfile << "," << list.vp[i][j];
		}
		for (int j = 0; j < N_DIM; j++) {
			myfile << "," << list.vn[i][j];
		}

		// Fluxes
		for (int j = 0; j < N_DIM; j++) {
			myfile << "," << list.f[i][j];
		}
    */

		myfile << endl;
	}

  myfile.close();
}

// Swap everything but nb and ns.
void GridBayes::swap_list_entries(int k1, int k2) {

  // Temporary values for k1
  double pdf = list.pdf[k1];
  int* pos = new int[N_DIM];
  double* v = new double[N_DIM];
  double* vp = new double[N_DIM];
  double* vn = new double[N_DIM];
  double* f = new double[N_DIM];
  for (int i = 0; i < N_DIM; i++) {
    pos[i] = list.pos[k1][i];
    v[i]   = list.v[k1][i];
    vp[i]  = list.vp[k1][i];
    vn[i]  = list.vn[k1][i];
    f[i]   = list.f[k1][i];
  }

  // Swap values
  list.pdf[k1] = list.pdf[k2];
  list.pdf[k2] = pdf;
  for (int i = 0; i < N_DIM; i++) {
    list.pos[k1][i] = list.pos[k2][i];
    list.v[k1][i]   = list.v[k2][i];
    list.vp[k1][i]  = list.vp[k2][i];
    list.vn[k1][i]  = list.vn[k2][i];
    list.f[k1][i]   = list.f[k2][i];

    list.pos[k2][i] = pos[i];
    list.v[k2][i]   = v[i];
    list.vp[k2][i]  = vp[i];
    list.vn[k2][i]  = vn[i];
    list.f[k2][i]   = f[i];
  }

  // Free memory
  delete[] pos;
  delete[] v;
  delete[] vp;
  delete[] vn;
  delete[] f;
}

void GridBayes::create_new_list_entry(int* pos) {

  nn++;
  list.pdf[nn - 1] = 0;
  
  for (int j = 0; j < N_DIM; j++) {
    list.pos[nn - 1][j] = pos[j];
  }

  initialize_vn(nn-1);  // Initialize the new entry only.
  list.f[nn-1][0] = 1;

}

void GridBayes::check_init() {
  // Check if the class has been initialized properly.

  cout << "Class constant parameters\n\n";
  cout << "Dimension length: " << N_DIM << endl;
  cout << "Maximum list size: " << N_PMAX << endl;        // unused
  cout << "Measurement update steps: " << N_UPD << endl;  // unused
  cout << "Maximum simulation time: " << T_MAX << endl;   // used outside
  cout << "Simulation time step: " << DT << endl;         // used outside
  cout << "PDF value threshold: " << PDF_THRESH << endl;
  cout << "Initial state grid length: " << DN << endl;

  cout << "Grid length: { " << DX[0];
  for (int i = 1; i < N_DIM; i++) {
    cout << " , " << DX[i];
  }
  cout << " }" << endl;

  cout << "Initial state: { " << X0[0];
  for (int i = 1; i < N_DIM; i++) {
    cout << " , " << X0[i];
  }
  cout << " }" << endl;

  cout << "Initial grid standard dev.: { " << X0_STDEV[0];
  for (int i = 1; i < N_DIM; i++) {
    cout << " , " << X0_STDEV[i];
  }
  cout << " }" << endl;

  cout << endl;
}

