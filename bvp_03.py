#! /usr/bin/env python3
#
from dolfin import *

def bvp_03 ( header ):

#*****************************************************************************80
#
## bvp_03 solves a boundary value problem on a sequence of meshes.
#
#  Discussion:
#
#    -u'' = 2x/(1+x^2)^2, -8 < x < +8
#    u(0) = u0, u(1) = u1
#
#    Exact solution is u_exact(x) = arctangent ( x )
#
#    The point of this problem is to have the mesh read in from a file.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    02 October 2014
#
#  Author:
#
#    John Burkardt
#
#  Parameters:
#
#    Input, string HEADER, specifies file names.
#    HEADER + '_mesh.xml' is the mesh file to be read.
#    HEADER + '_solution.png' is the solution plot file to be created.
#
  import matplotlib.pyplot as plt
#
#  Read the mesh from an XML file.
#  The number of elements and their extents are determined from the file.
#
  filename = header + '_mesh.xml'
  print ( '' )
  print ( '  Reading mesh file "%s"' % ( filename ) )
  mesh = Mesh ( filename )
#
#  Define the exact solution.
#
  u_exact = Expression ( "atan ( x[0] )", degree = 10 )
#
#  Define the right hand side function.
#
  f = Expression ( "2 * x[0] / pow ( 1 + x[0] * x[0], 2 )", degree = 10 )
#
#  Define the function X_DIRICHLET ( ) to be TRUE for points X where
#  a Dirichlet condition is to be applied.
#
  def x_dirichlet ( x, on_boundary ):
    return on_boundary
#
#  Define the function space to be of type 'Lagrange'
#  using piecewise linear basis functions (1).
#
  V = FunctionSpace ( mesh, "Lagrange", 1 )
#
#  The "TEST" function is the one that multiplies the PDE.
#  The "TRIAL" function is the one that is used to build the solution.
#  The solution should be a function from the function space.
#
  psii = TestFunction ( V )
  psij = TrialFunction ( V )
#
#  Define the system matrix, Aij = integral ( d/dx psii * d/dx psij ) dx
#
  A = inner ( grad ( psii ), grad ( psij ) ) * dx
#
#  The system right hand side.
#
  RHS = psii * f * dx
#
#  The solution should be a function from the function space.
#
  u = Function ( V )
#
#  The boundary conditions, stored as "bc", are Dirichlet boundary conditions,
#  defined using the function space V, with value U_EXACT, for points for which 
#  the X_DIRICHLET function is true.
#
  bc = DirichletBC ( V, u_exact, x_dirichlet )
#
#  Solve the problem for U, with the given boundary conditions.
#
  solve ( A == RHS, u, bc )
#
#  Plot the solution.
#
  label = header 
  plot ( u, title = header )
  filename = header + '_solution.png'
  plt.savefig ( filename )
  print ( '  Graphics saved as "%s"' % ( filename ) )
  plt.close ( )
#
#  Terminate.
#
  return

def bvp_03_test ( ):

#*****************************************************************************80
#
## bvp_03_test tests bvp_03.
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
  print ( 'bvp_03_test:' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  FENICS version %s'% ( dolfin.__version__ ) )
  print ( '  Solve  -u'' = 2x/(1+x^2)^2, -8 < x < +8' )
  print ( '  u(0) = u0, u(1) = u1 ' )
  print ( '  Exact solution is u(x) = arctangent(x)' )
  print ( '  Solve on a sequence of meshes.' )

  header = 'bvp_03_04'
  bvp_03 ( header )

  header = 'bvp_03_08'
  bvp_03 ( header )

  header = 'bvp_03_16'
  bvp_03 ( header )

  header = 'bvp_03_32'
  bvp_03 ( header )

  print ( '' )
  print ( 'bvp_03_test:' )
  print ( '  Normal end of execution.' )
  print ( '' )
  print ( time.ctime ( time.time() ) )
  return

if ( __name__ == '__main__' ):

  bvp_03_test ( )
