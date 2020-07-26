#! /usr/bin/env python3
#
from dolfin import *

def bvp_02 ( e_num ):

#*****************************************************************************80
#
## bvp_02 solves a boundary value problem.
#
#  Discussion:
#
#    Repeat the calculation in BVP_01, but now concentrate on the following
#    * print a table of the solution values;
#    * compute the L2 norm of the errror;
#    * compute the H1 seminorm of the error;
#    * compute the average value of the solution.
#
#    -u'' + u = x, 0 < x < 1
#    u(0) = 0, u(1) = 0
#
#    Exact solution is u(x) = x - sinh(x)/sinh(1)
#
#  Licensing:
#
#    This code is distributed under the MIT License
#
#  Parameters:
#
#    Input, integer E_NUM, the number of elements.
#
  import numpy as np
#
#  Create a mesh on the unit interval.
#
  mesh = UnitIntervalMesh ( e_num )
#
#  Define the function space to be of Lagrange type
#  using piecewise linear basis functions.
#
  V = FunctionSpace ( mesh, "Lagrange", 1 )
#
#  For convenience, we define the exact solution.
#  In order to compute the H1 seminorm of the error, we
#  also need the derivative of the exact solution.
#
  exact = Expression ( 'x[0] - sinh ( x[0] ) / sinh ( 1.0 )', degree = 10 )
  exact_p = Expression ( '1 - cosh ( x[0] ) / sinh ( 1.0 )', degree = 10 )
#
#  Define U0_BOUNDARY ( ) to indicate boundary points.
#
  def u0_boundary ( x, on_boundary ):
    return on_boundary
#
#  Set the boundary conditions.
#
  bc = DirichletBC ( V, exact, u0_boundary )
#
#  Define the system matrix.
#
  psii = TestFunction ( V )
  psij = TrialFunction ( V )
  A = ( inner ( nabla_grad ( psii ), nabla_grad ( psij ) ) + psii * psij ) * dx
#
#  Set up the right hand side.
#
  f = Expression ( 'x[0]', degree = 10 )
  RHS = psii * f * dx
#
#  Specify that the solution is a linear combination of elements of V.
#
  u = Function ( V )
#
#  Solve the problem for U, with the given boundary conditions.
#
  solve ( A == RHS, u, bc )
#
#  Extract, as the array "X", the coordinates of the mesh.
#
  x = mesh.coordinates ( )
#
#  Extract, as the array "C", the finite element coeficients.
#
  c = u.vector().get_local()
#
#  Set up a 2 point quadrature rule defined on the reference interval [0,1].
#
  xg = np.array ( ( 0.2113248654, 0.7886751346 ) )
  wg = np.array ( ( 0.5, 0.5 ) )
#
#  1) Print a table of the solution at 3*N+1 points.
#
#  Since the U created by the solve command is actually a finite element function,
#  not just a list of coefficients, we can evaluate U anywhere.
#
  print ( '' )
  print ( '    X        U(X)' )
  print ( '' )

  for e in range ( 0, e_num ):

    xl = x[e]
    xr = x[e+1]
#
#  We have to use one extra point in the last interval.
#
    if ( e < e_num - 1 ):
      ihi = 2
    else:
      ihi = 3

    for i in range ( 0, ihi ):
      xi = ( ( 2 - i ) * xl + ( i ) * xr ) / 2
      ux = u ( xi )
      print ( '  %8f  %8f' % ( xi, ux ) )
#
#  2) Compute the L2 norm of the error.
#
#  Square root of the integral of the square of the difference of 
#  the exact and computed solutions.
#
#  If I try to do this explicitly, I get all kinds of nasty errors.
#
  l2_error = ( u - exact ) ** 2 * dx
  l2_error = sqrt ( assemble ( l2_error ) )
  print ( '' )
  print ( '  L2 norm of error is %g' % ( l2_error ) )
#
#  3) Compute the H1 seminorm of the error.
#
#  Square root of the integral of the square of the difference of 
#  the derivatives of the exact and computed solutions.
#
#  Again, I can't write this out myself!
#
  e_V = Function ( V )
  u_e_V = interpolate ( exact, V )
  u_V = interpolate ( u, V )
  e_V.vector()[:] = u_e_V.vector().get_local() - u_V.vector().get_local()
  h1_error = inner ( grad ( e_V ), grad ( e_V ) ) * dx
  h1_error = sqrt ( assemble ( h1_error ) )
  print ( '  H1 seminorm of error is %g' % ( h1_error ) )
#
#  4A) Compute the average value of the solution.
#
#  Integral of u divided by length of interval.
#
  ave = 0.0

  for e in range ( 0, e_num ):

    xl = x[e]
    xr = x[e+1]

    for q in range ( 0, 2 ):

      xq = xl + xg[q] * ( xr - xl )
      phil = ( xr - xq ) / ( xr - xl )
      phir = ( xq - xl ) / ( xr - xl )
      uq = c[e] * phil + c[e+1] * phir
      ave = ave + wg[q] * uq * ( xr - xl )

  ave = ave / 1.0
  print ( '' )
  print ( '  Average value of solution (my way) is %g' % ( ave ) )
#
#  4B) Compute the average value of the solution the FENICS way.
#
  u_V = interpolate ( u, V )
  ave = assemble ( u_V * dx ) / 1.0
  print ( '  Average value of solution (FENICS way) is %g' % ( ave ) )

  return

def bvp_02_test ( ):

#*****************************************************************************80
#
## bvp_02_test tests bvp_02.
#
#  Licensing:
#
#    This code is distributed under the MIT Licesnse
  import dolfin
  import platform
  import time

  print ( time.ctime ( time.time() ) )
#
#  Report level = only warnings or higher.
#
  level = 30
  set_log_level ( level )

  print ( '' )
  print ( 'bvp_02_test:' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  FENICS version %s'% ( dolfin.__version__ ) )
  print ( '  Solve  -u" + u = x, 0 < x < 1' )
  print ( '  u(0) = 0, u(1) = 0' )
  print ( '  Exact solution is u(x) = x - sinh(x)/sinh(1)' )

  e_num = 8
  print ( '  Using %d elements.' % ( e_num ) )

  bvp_02 ( e_num )
#
#  Terminate.
#
  print ( '' )
  print ( 'bvp_02_test:' )
  print ( '  Normal end of execution.' )
  print ( '' )
  print ( time.ctime ( time.time() ) )
  return

if ( __name__ == '__main__' ):

  bvp_02_test ( )
