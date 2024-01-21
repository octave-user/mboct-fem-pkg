## Copyright (C) 2020(-2023) Reinhard <octave-user@a1.net>
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
## @deftypefn{Function File} @var{Afact} = fem_sol_factor(@var{A})
## @deftypefnx{} @var{Afact} = fem_sol_factor(@var{A}, @var{options})
## Build a factor object @var{Afact} for matrix @var{A} which can be used to solve
## a system of linear equations @var{A} * @var{x} = @var{b}
## by means of the expression @var{x} = @var{Afact} \ @var{b}.
##
## @var{A} @dots{} Quadratic full rank matrix.
##
## @var{options}.solver @dots{} Character string name of the linear solver to use.
##
## Use one of @{"pastix" | "pardiso" | "mumps" | "umfpack" | "lu" | "chol" | "mldivide"@}.
##
## @var{options}.refine_max_iter @dots{} Maximum number of refinement iterations.
##
## @var{options}.number_of_threads @dots{} Number of threads to use of a SMP solution.
##
## @var{options}.symmetric @dots{} Define if @var{A} is a symmetric matrix even if only the
##                                 upper or lower triangular part of @var{A} is passed
##
## @var{options}.verbose @dots{} A value greater than zero enables verbose output.
##
## @var{options}.bind_thread_mode @dots{} Define if thread affinity should used
##
## @var{options}.workspace_inc @dots{} Memory allocation flag used only if solver == "mumps".
##
## @end deftypefn

function Afact = fem_sol_factor(A, options)
  if (nargin < 1 || nargin > 2 || nargout > 1)
    print_usage();
  endif

  if (nargin < 2)
    options = struct();
  endif

  if (~isfield(options, "solver"))
    options.solver = fem_sol_select(false);
  endif

  if (~isfield(options, "refine_max_iter"))
    options.refine_max_iter = int32(10);
  endif

  if (~isfield(options, "epsilon_refinement"))
    options.epsilon_refinement = -1;
  endif

  if (~isfield(options, "number_of_threads"))
    options.number_of_threads = int32(1);
  endif

  if (~isfield(options, "pre_scaling"))
    options.pre_scaling = false;
  endif

  if (~isfield(options, "scale_tol"))
    options.scale_tol = sqrt(eps);
  endif

  if (~isfield(options, "scale_max_iter"))
    options.scale_max_iter = int32(100);
  endif

  if (~isfield(options, "symmetric"))
    options.symmetric = true;
  endif

  if (~isfield(options, "bind_thread_mode"))
    options.bind_thread_mode = false;
  endif

  if (~isfield(options, "verbose"))
    options.verbose = false;
  endif

  if (options.pre_scaling)
    Afact = fem_fact_scale(A, options);
    return;
  endif

  switch (options.solver)
    case {"chol", "lu", "mldivide"}
      if (options.refine_max_iter > 0)
        Afact = fem_fact_refine(A, options);
        return;
      endif
  endswitch

  blambda = full(any(diag(A) <= 0));

  options.solver = fem_sol_select(blambda, options.solver);

  if (~issparse(A))
    switch (options.solver)
      case {"chol", "lu", "mldivide"}
        ## solvers which can handle also dense matrices
      otherwise
        options.solver = "lu";
    endswitch
  endif

  switch (options.solver)
    case "pastix"
      if (options.symmetric)
        options.matrix_type = PASTIX_API_SYM_YES;
        options.factorization = PASTIX_API_FACT_LDLT;
      else
        options.matrix_type = PASTIX_API_SYM_NO;
        options.factorization = PASTIX_API_FACT_LU;
      endif

      if (options.verbose)
        options.verbose = PASTIX_API_VERBOSE_YES;
      else
        options.verbose = PASTIX_API_VERBOSE_NOT;
      endif

      if (options.bind_thread_mode)
        options.bind_thread_mode = PASTIX_API_BIND_AUTO;
      else
        options.bind_thread_mode = PASTIX_API_BIND_NO;
      endif

      linear_solver = @fem_fact_pastix;
    case "pardiso"
      if (~isfield(options, "scaling"))
        options.scaling = blambda;
      endif

      if (~isfield(options, "weighted_matching"))
        options.weighted_matching = blambda;
      endif

      linear_solver = @fem_fact_pardiso;
    case "mumps"
      if (options.symmetric)
        options.matrix_type = MUMPS_MAT_SYM;
      else
        options.matrix_type = MUMPS_MAT_GEN;
      endif

      if (~isfield(options, "workspace_inc"))
        options.workspace_inc = int32(50);
      endif

      if (options.verbose)
        options.verbose = MUMPS_VER_ALL;
      else
        options.verbose = MUMPS_VER_ERR;
      endif

      linear_solver = @fem_fact_mumps;
    case "umfpack"
      linear_solver = @fem_fact_umfpack;
    case "chol"
      ## If we are using Lagrange multipliers, there will be always zeros at the main diagonal
      linear_solver = @fem_fact_chol;
    case "lu"
      linear_solver = @fem_fact_lu;
    case "mldivide"
      linear_solver = @fem_fact;
    otherwise
      error("unknown linear solver: \"%s\"", options.solver);
  endswitch

  Afact = linear_solver(A, options);
endfunction

