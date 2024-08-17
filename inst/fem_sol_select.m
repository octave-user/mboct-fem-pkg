## Copyright (C) 2020(-2021) Reinhard <octave-user@a1.net>
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; If not, see <http://www.gnu.org/licenses/>.

##-*- texinfo -*-
## @deftypefn{Function File} @var{solver} = fem_sol_select(@var{blambda})
## @deftypefnx{} @var{solver} = fem_sol_select(@var{blambda}, @var{solver})
## Return the name of a solver which can be used to solve a system of linear equations.
## If @var{blambda} is true, return a solver which can do also pivoting.
## The preferred solver can be provided as input argument.
## @end deftypefn

function solver = fem_sol_select(blambda, solver)
  if (nargin < 1 || nargin > 2 || nargout > 1)
    print_usage();
  endif

  if (nargin < 2 || ~fem_sol_check_func(solver))
    solvers = {"pardiso", "pastix", "mumps", "umfpack", "chol", "lu"};
    
    for i=1:numel(solvers)
      if (fem_sol_check_func(solvers{i}))
        solver = solvers{i};
        break;
      endif
    endfor
  endif
  
  switch (solver)
    case "chol"
      if (blambda)
        solver = "lu";
      endif
  endswitch
endfunction
