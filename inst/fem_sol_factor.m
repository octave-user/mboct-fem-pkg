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
## @var{options}.factorization @dots{} Kind of factorization to use.
##
## If solver == "pastix" use one of @{PASTIX_API_FACT_LDLT | PASTIX_API_FACT_LLT | PASTIX_API_FACT_LU@}
##
## @var{options}.matrix_type @dots{} Type of input matrix.
##
## If solver == "pastix" use one of @{PASTIX_API_SYM_YES | PASTIX_API_SYM_NO@}
##
## If solver == "mumps" use one of @{MUMPS_MAT_DEF | MUMPS_MAT_GEN | MUMPS_MAT_SYM@}
##
## @var{options}.verbose @dots{} A value greater than zero enables verbose output.
##
## @var{options}.bind_thread_mode @dots{} If solver == "pastix" use one of @{PASTIX_API_BIND_NO, PASTIX_API_BIND_AUTO@}
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
    options.solver = "pastix";
  endif

  if (~isfield(options, "refine_max_iter"))
    options.refine_max_iter = int32(3);
  endif

  if (~isfield(options, "number_of_threads"))
    options.number_of_threads = int32(1);
  endif

  if (~isfield(options, "pre_scaling"))
    options.pre_scaling = true;
  endif

  if (~isfield(options, "scale_tol"))
    options.scale_tol = sqrt(eps);
  endif

  if (~isfield(options, "scale_max_iter"))
    options.scale_max_iter = int32(100);
  endif

  if (options.pre_scaling)
    Afact = fem_fact_scale(A, options);
  else
    blambda = full(any(diag(A) <= 0));

    options.solver = fem_sol_select(blambda, options.solver);

    switch (options.solver)
      case {"pastix", "pastix_ref"}
        if (~isfield(options, "matrix_type"))
          if (isfield(options, "symmetric"))
            if (options.symmetric)
              options.matrix_type = PASTIX_API_SYM_YES;
            else
              options.matrix_type = PASTIX_API_SYM_NO;
            endif
          else
            options.matrix_type = PASTIX_API_SYM_YES;
          endif
        endif

        if (~isfield(options, "factorization"))
          options.factorization = PASTIX_API_FACT_LDLT;
        endif

        if (~isfield(options, "verbose"))
          options.verbose = PASTIX_API_VERBOSE_NOT;
        endif

        if (~isfield(options, "bind_thread_mode"))
          options.bind_thread_mode = PASTIX_API_BIND_NO;
        endif

        switch (options.solver)
          case "pastix"
            linear_solver = @fem_fact_pastix;
          case "pastix_ref"
            linear_solver = @fem_fact_pastix_ref;
        endswitch
      case "pardiso"
        if (~isfield(options, "symmetric"))
          options.symmetric = true;
        endif

        if (~isfield(options, "verbose"))
          options.verbose = false;
        endif

        if (~isfield(options, "scaling"))
          options.scaling = blambda;
        endif

        if (~isfield(options, "weighted_matching"))
          options.weighted_matching = blambda;
        endif

        linear_solver = @fem_fact_pardiso;
      case "mumps"
        if (blambda)
          options.matrix_type = MUMPS_MAT_SYM;
        else
          options.matrix_type = MUMPS_MAT_DEF;
        endif

        if (~isfield(options, "workspace_inc"))
          options.workspace_inc = int32(50);
        endif

        if (~isfield(options, "verbose"))
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
  endif
endfunction

%!test
%! state = rand("state");
%! unwind_protect
%!   rand("seed", 0);
%!   solvers = {"pastix", "pardiso", "mumps", "umfpack", "lu", "chol", "mldivide"};
%!   M = 20;
%!   options.refine_max_iter = 10;
%!   options.number_of_threads = int32(4);
%!   tol = sqrt(eps);
%!   for s=[false, true]
%!     options.pre_scaling = s;
%!     for N=2.^(1:8)
%!       A = gallery("poisson", N);
%!       A += A.';
%!       b = rand(columns(A), M);
%!       for i=1:numel(solvers)
%!         options.solver = solvers{i};
%!         Afact = fem_sol_factor(A, options);
%!         x = Afact \ b;
%!         f = max(norm(A * x - b, "cols") ./ norm(A * x + b, "cols"));
%!         assert(f < tol);
%!       endfor
%!     endfor
%!   endfor
%! unwind_protect_cleanup
%!   rand("state", state);
%! end_unwind_protect
