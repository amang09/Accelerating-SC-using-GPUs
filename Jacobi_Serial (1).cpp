#include <cmath>
#include <cstdlib>
#include <iostream>

const int NX = 9;
const int NY = 9;

inline int indx(const int i, const int j) { return i + NY * j; }

double linf_norm(const double *u1, const double *u2, const int sz) {
  double max = 0.0;
  for (int i = 0; i < sz; ++i) {
    if (std::abs(u1[i] - u2[i]) > max) {
      max = std::abs(u1[i] - u2[i]);
    }
  }
  return max;
}

double l2_norm(const double *u1, const double *u2, const int sz) {
  double rms = 0.0;
  for (int i=0; i < sz; ++i) {
    rms += std::abs(u1[i] - u2[i])*std::abs(u1[i] - u2[i]);
  }
  rms = sqrt (rms / sz);
  return rms;
}

// Problem: -(d2u/dx2 + d2u/dy2) = f(x,y)
// uexact = sin(pi*x)*cos(pi*y); f(x,y) = 2*pi2*sin(pi*x)*cos(pi*y)

int main(int argc, char *argv[]) {
  double *u = (double *)malloc(NX * NY * sizeof(double));   // solution
  double *ue = (double *)malloc(NX * NY * sizeof(double));  // exact solution
  double *uo =
      (double *)malloc(NX * NY * sizeof(double));  // sol. from prev. iteration
  double *rh = (double *)malloc(
      NX * NY * sizeof(double));  // right hand side of Poisson eq.

  double *xgrid = (double *)malloc(NX * sizeof(double));
  double *ygrid = (double *)malloc(NY * sizeof(double));

  const double xstart = -1.0;
  const double ystart = -1.0;
  const double xend = 1.0;
  const double yend = 1.0;

  const double dx = (xend - xstart) / (NX - 1);
  const double dy = (yend - ystart) / (NY - 1);

  double (*mynorm)(const double *, const double *, const int) = l2_norm;

  // Grid
  for (int i = 0; i < NX; ++i) {
    xgrid[i] = xstart + i * dx;
  }

  for (int j = 0; j < NY; ++j) {
    ygrid[j] = ystart + j * dy;
  }

  // rhs
  for (int j = 0; j < NY; ++j) {
    for (int i = 0; i < NX; ++i) {
      rh[indx(i, j)] =
          2 * (M_PI * M_PI) * sin(M_PI * xgrid[i]) * cos(M_PI * ygrid[j]);
    }
  }

  // Exact solution
  for (int j = 0; j < NY; ++j) {
    for (int i = 0; i < NX; ++i) {
      ue[indx(i, j)] = sin(M_PI * xgrid[i]) * cos(M_PI * ygrid[j]);
    }
  }

  // Initialize
  for (int j = 0; j < NY; ++j) {
    for (int i = 0; i < NX; ++i) {
      if (i == 0 || i == NX - 1 || j == 0 || j == NY - 1) {
        u[indx(i, j)] = ue[indx(i, j)];
      } else {
        u[indx(i, j)] = 0.0;
      }
    }
  }

  // The magnitude of u is of the order of 1, so we will make the convergence
  // criterion as 1.e-14 times 1, i.e., simply 1.e-14
  double tol = 1.e-10;

  // Initialize the difference between old and new solutions
  double diff = 1.0;

  // counter to store number of iterations
  int counter = 0;
  while (diff > tol) { // main while loop
    // copy  u to uold
    for (int i = 0; i < NX * NY; ++i) {
      uo[i] = u[i];
    }

    // Perform the Jacobi iteration (only on interior cells)
    for (int j = 1; j < NY - 1; ++j) {
      for (int i = 1; i < NX - 1; ++i) {
        // here we assume that dx = dy
        u[indx(i, j)] = 0.25 * (uo[indx(i + 1, j)] + uo[indx(i, j + 1)] +
                                uo[indx(i - 1, j)] + uo[indx(i, j - 1)] +
                                (dx * dx) * rh[indx(i, j)]);
      }
    }

    diff = mynorm(u, uo, NX*NY);
    std::cout << "Iteration " << counter << ": diff = " << diff << std::endl;
    counter++;
  }

  double me = mynorm(u, ue, NX*NY); // max. error in sol.

  std::cout << "Finally Iterations = " << counter << "; max error = " << me << std::endl;

  free(u);
  free(ue);
  free(uo);
  free(rh);
  free(xgrid);
  free(ygrid);
}
