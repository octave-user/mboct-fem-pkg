## Copyright (C) 2011(-2020) Reinhard <octave-user@a1.net>
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

## -*- texinfo -*-
## @deftypefn {Function File} @var{U} = fem_sol_linsolve(@var{K}, @var{R})
## @deftypefnx {} @dots{} = fem_sol_linsolve(@dots{}, @var{options})
##
## Choose an appropriate linear solver and solve  @var{K} * @var{U} = @var{R}
##
## @var{K} @dots{} stiffness matrix
##
## @var{R} @dots{} matrix of load vectors
##
## @var{U} @dots{} matrix of displacement vectors
##
## @var{options}.number_of_threads @dots{} set the number of threads for parallel solution
##
## @var{options}.refine_max_iter @dots{} set the maximum number of iterations for refinement of the solution
##
## @end deftypefn

function U = fem_sol_linsolve(K, R, options)
  if (nargin < 2 || nargin > 3 || nargout > 1)
    print_usage();
  endif

  if (nargin < 3)
    options = struct();
  endif

  if (~isfield(options, "solver"))
    options.solver = struct();
  endif

  if (~isfield(options, "refine_max_iter"))
    options.refine_max_iter = int32(3);
  endif

  if (~isfield(options.solver, "number_of_threads"))
    options.solver.number_of_threads = int32(1);
  endif
  
  blambda = full(any(diag(K) <= 0));
  
  if (fem_sol_check_func("pastix"))
    opt_pastix.matrix_type = PASTIX_API_SYM_YES;      
    opt_pastix.factorization = PASTIX_API_FACT_LDLT;
    opt_pastix.verbose = PASTIX_API_VERBOSE_NOT;
    opt_pastix.bind_thread_mode = PASTIX_API_BIND_NO;
    opt_pastix.number_of_threads = options.solver.number_of_threads;
    opt_pastix.refine_max_iter = options.refine_max_iter; ## Do not use pastix without refinement (see fem_pre_mesh_import:test13)!

    Kfact = fem_fact_pastix(K, opt_pastix);
  elseif (fem_sol_check_func("mumps"))
    if (blambda)
      opt_mumps.matrix_type = MUMPS_MAT_SYM;
    else
      opt_mumps.matrix_type = MUMPS_MAT_DEF;
    endif
    opt_mumps.refine_max_iter = options.refine_max_iter;
    opt_mumps.verbose = MUMPS_VER_ERR;
    Kfact = fem_fact_mumps(K, opt_mumps);    
  elseif (~blambda) ## If we are using Lagrange multipliers, there will be always zeros at the main diagonal
    opt_chol = options;
    Kfact = fem_fact_chol(K, opt_chol);
  elseif (fem_sol_check_func("umfpack"))
    opt_umfpack.refine_max_iter = options.refine_max_iter;
    Kfact = fem_fact_umfpack(K, opt_umfpack);
  else
    Kfact = fem_fact(K);
  endif

  U = Kfact \ R;
endfunction
