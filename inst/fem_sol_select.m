## Copyright (C) 2020(-2024) Reinhard <octave-user@a1.net>
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

  persistent solvers = [];

  if (isempty(solvers))
    ## FIXME: work around suspected regression in PaStiX 6.4.0 (1151c30a25e2014ff4b39bf8c7ac4b381913f193)
    solvers = struct("name",   {"pastix", "pardiso", "mumps", "umfpack", "chol",  "lu", "mldivide"}, ...
                     "lambda", {    true,      true,    true,      true,  false,  true,       true});

    idx_func = cellfun(@fem_sol_check_func, {solvers.name}, "UniformOutput", true);

    solvers = solvers(idx_func);
  endif

  solsel = solvers;

  if (blambda)
    solsel = solsel([solvers.lambda]);
  endif

  if (nargin < 2)
    solver = solsel(1).name;
  endif

  valid_solver = any(cellfun(@(name) strcmp(solver, name), {solsel.name}, "UniformOutput", true));

  if (~valid_solver)
    solver = solsel(1).name;
  endif
endfunction
