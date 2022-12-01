# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>
# include <cmath>

using namespace std;

# define NX 4
# define NY 4

int main ( int argc, char *argv[] );

void rhs ( int nx, int ny, double f[NX][NY] );
void sweep ( int nx, int ny, double dx, double dy, double f[NX][NY], 
double u[NX][NY], double unew[NX][NY] );
double f_exact ( double x, double y );



int main ( int argc, char *argv[] )
{
  bool converged;
  double diff;
  double dx;
  double dy;
  double error;
  double f[NX][NY];
  int i;
  int it;
  int it_max = 1;
  int j;
  int nx = NX;
  int ny = NY;
  double tolerance = 0.1;
  double u[NX][NY];
  double u_norm;
  double udiff[NX][NY];;
  double unew[NX][NY];
  double x;
  double y;
  double uprev[NX][NY];

  dx = 2.0 / ( double ) ( nx - 1 );
  dy = 2.0 / ( double ) ( ny - 1 );
  rhs ( nx, ny, f );

  for ( j = 0; j < ny; j++ )
  {
    for ( i = 0; i < nx; i++ )
    {
      if ( i == 0 || i == nx - 1 || j == 0 || j == ny - 1 )
      {
        unew[i][j] = f[i][j];
      }
      else
      {
        unew[i][j] = 0.0;
      }
    }
  }
  

  converged = false;

  cout << "\n";
  


cout << "\n";

    for(i=0;i<nx;i++){
        for(j=0;j<ny;j++){
        u[i][j] = 0;
        }
    }
    
  for (j = 0; j < ny; j++ )
  {
    for (i = 0; i < nx; i++ )
    {
        cout<<" ";
      for ( it = 1; it <= it_max; it++)
      {
        uprev[i][j] = unew[i][j];

        sweep ( nx, ny, dx, dy, f, u, unew );
        u[i][j] = unew[i][j];
        
        
        udiff[i][j] = u[i][j] - uprev[i][j];
        cout<<u[i][j]<<"        "<<udiff[i][j]<<"\n";
        
        
        
            // if (udiff[i][j] <= tolerance)
            // {
            //     converged = true;
            //     break;
            // }
            
      }
    }
}

  if ( converged )
  {
    cout << "  The iteration has converged.\n";
  }
  else
  {
    cout << "  The iteration has NOT converged.\n";
  }
//
//  Terminate.
//
  cout << "\n";
  cout << "POISSON_SERIAL:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  

  return 0;
}

//****************************************************************************80

void rhs ( int nx, int ny, double f[NX][NY] )
{
  double fnorm;
  int i;
  int j;
  double x;
  double y;

  for ( j = 0; j < ny; j++ )
  {
    y = -1 + 2*( double ) ( j+1 ) / ( double ) ( ny - 1 );
    for ( i = 0; i < nx; i++ )
    {
      x = -1 + 2*( double ) ( i+1) / ( double ) ( nx - 1 );
      f[i][j] = f_exact ( x, y );
    }
  }

  return;
}
//****************************************************************************80

void sweep ( int nx, int ny, double dx, double dy, double f[NX][NY], 
  double u[NX][NY], double unew[NX][NY] )
{
  int i;
  int j;

  for ( j = 0; j < ny; j++ )
  {
    for ( i = 0; i < nx; i++ )
    {
      if ( i == 0 || j == 0 || i == nx - 1 || j == ny - 1 )
      {
        unew[i][j] = f[i][j];
      }
      else
      { 
        unew[i][j] = 0.25 * ( 
          u[i-1][j] + u[i][j+1] + u[i][j-1] + u[i+1][j] + f[i][j] * dx * dy );
      }
    }
  }
  return;
}
//****************************************************************************80

double f_exact ( double x, double y )

{
  double r8_pi = 3.141592653589793;
  double value;

  value = sin ( r8_pi *x )*cos(r8_pi*y);

  return value;
}


# undef NX
# undef NY

