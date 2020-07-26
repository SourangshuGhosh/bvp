#! /usr/bin/env python3
#
from dolfin import *

def bvp_01 ( e_num ):

#*****************************************************************************80
#
## bvp_01 solves a boundary value problem.
#
#  Discussion:
#
#    -u'' + u = x, 0 < x < 1
#    u(0) = 0, u(1) = 0
#
#    Exact solution is u(x) = x - sinh(x)/sinh(1)
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    23 January 2016
#
#  Author:
#
#    John Burkardt
#
#  Parameters:
#
#    Input, integer E_NUM, the number of elements.
#
  import matplotlib.pyplot as plt
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
#  Define the exact solution.
#
  exact = Expression ( "x[0] - sinh ( x[0] ) / sinh ( 1.0 )", degree = 10 )
#
#  Define U0_BOUNDARY ( ) to indicate boundary points.
#
  def u0_boundary ( x, on_boundary ):
    return on_boundary
#
#  The boundary conditions, stored as "bc", are Dirichlet boundary conditions,
#  defined using the function space V, with value U0, for points for which 
#  U0_BOUNDARY is true.
#
  bc = DirichletBC ( V, exact, u0_boundary )
#
#  The "TEST" function is the one that multiplies the PDE.
#  The "TRIAL" function is the one that is used to build the solution.
#
  psii = TestFunction ( V )
  psij = TrialFunction ( V )
#
#  Define the system matrix, A = integral ( psii' * psij' + psii * psij ) dx
#
  A = ( inner ( nabla_grad ( psii ), nabla_grad ( psij ) ) + psii * psij ) * dx
#
#  Define the function on the right hand side.
#
  f = Expression ( "x[0]", degree = 1 )
#
#  Define the linear system right hand side, RHS = integral ( psii * f ) dx
#
  RHS = psii * f * dx
#
#  Specify that the solution should be a linear combination of elements of V.
#
  u = Function ( V )
#
#  Solve the problem for U, with the given boundary conditions.
#
  solve ( A == RHS, u, bc )
#
#  Plot the solution.
#
  plot ( u, title = 'bvp_01 solution' )
  filename = 'bvp_01_solution.png'
  plt.savefig ( filename )
  print ( '  Graphics saved as "%s"' % ( filename ) )
  plt.close ( )

  return

def bvp_01_test ( ):

#*****************************************************************************80
#
## bvp_01_test tests bvp_01.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    21 October 2018
#
#  Author:
#
#    John Burkardt
#
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
  print ( 'bvp_01_test:' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  FENICS version %s'% ( dolfin.__version__ ) )
  print ( '  Boundary value problem.' )
  print ( '  -u" + u = x, 0 < x < 1' )
  print ( '  u(0) = 0, u(1) = 0' )
  print ( '  Exact solution is u(x) = x - sinh(x)/sinh(1)' )

  e_num = 8
  print ( '  Using %d elements.' % ( e_num ) )

  bvp_01 ( e_num )
#
#  Terminate.
#
  print ( '' )
  print ( 'bvp_01_test:' )
  print ( '  Normal end of execution.' )
  print ( '' )
  print ( time.ctime ( time.time() ) )
  return

if ( __name__ == '__main__' ):

  bvp_01_test ( )
