## Copyright (C) 2019(-2024) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{U}] = fem_sol_harmonic_gmres(@var{omega}, @var{M}, @var{D}, @var{K}, @var{R}, @var{idx_output}, @var{opt_sol})
## @deftypefnx {} [@dots{}] = fem_sol_harmonic_gmres(@var{omega}, @var{callback})
## Solve -@var{omega}^2 * @var{M} + 1j * @var{omega} * @var{D} + @var{K} = @var{R} using a preconditioned gmres solver
## @end deftypefn

function U = fem_sol_harmonic_gmres(omega, varargin)
  if (numel(varargin) == 1 && is_function_handle(varargin{1}))
    func = @fem_sol_harmonic_gmres_callback;
  elseif (numel(varargin) == 6 && ismatrix(varargin{1}) && ismatrix(varargin{2}) && ismatrix(varargin{3}) && ismatrix(varargin{4}) && isvector(varargin{5}) && isstruct(varargin{6}))
    func = @fem_sol_harmonic_gmres_matrix;
  else
    print_usage();
  endif

  U = func(omega, varargin{:});
endfunction

function U = fem_sol_harmonic_gmres_matrix(omega, M, D, K, R, idx_output, opt_sol)
  cb = @(stage, opt_dummy, idx, omega, Ui, status) fem_sol_harmonic_gmres_matrix_callback(stage, opt_sol, idx, omega, Ui, status, M, D, K, R, idx_output);
  U = fem_sol_harmonic_gmres_callback(omega, cb);
endfunction

function varargout = fem_sol_harmonic_gmres_matrix_callback(stage, opt_sol, idx, omega, Ui, status, M, D, K, R, idx_output)
  switch (stage)
    case "init"
      U = complex(zeros(numel(idx_output), numel(omega)));
      varargout = {U, opt_sol};
    case "pre"
      varargout = {M, D, K, R};
    case "post"
      varargout = {Ui(idx_output, :)};
  endswitch
endfunction

function U = fem_sol_harmonic_gmres_callback(omega, callback)
  opt_sol = struct();

  status.gmres.iter = repmat(intmax(), 1, 2);
  status.gmres.flag = -1;
  status.cpu.factor = 0;
  status.cpu.solve = 0;

  [U, opt_sol] = callback("init", opt_sol, [], omega, [], status);

  for i=1:numel(omega)
    [M, D, K, R] = callback("pre", opt_sol, i, omega(i), [], status);

    opt_sol.gmres.restart = min([opt_sol.gmres.restart, columns(K)]);

    do
      status.do_fact = (status.gmres.flag ~= 0 || status.gmres.iter(1) > 1 || status.gmres.iter(2) >= opt_sol.gmres.maxiter || status.cpu.solve > opt_sol.gmres.alpha * status.cpu.factor);

      if (status.do_fact)
        start = cputime();
        A = -omega(i)^2 * M + 1j * omega(i) * D + K;
        [AS, D1, D2] = fem_sol_matrix_scale(A, opt_sol.scale.tol, opt_sol.scale.maxiter);
        clear KSfact;
        ASfact = fem_sol_factor(AS, opt_sol.factor);
        status.cpu.factor = cputime() - start;
      endif

      Z = ASfact \ (diag(D1) * R);

      start = cputime();

      status.gmres.relres = -realmax();

      for j=1:columns(Z)
        [Z(:, j), status.gmres.flag, relres, status.gmres.iter] = gmres(@(Z) mat_product(M, D, K, D1, D2, ASfact, omega(i), Z), Z(:, j), opt_sol.gmres.restart, opt_sol.gmres.tol, opt_sol.gmres.maxiter, [], [], Z(:, j));

        status.gmres.relres = max([relres, status.gmres.relres]);

        if (status.gmres.flag ~= 0)
          break;
        endif
      endfor

      if (status.gmres.flag ~= 0)
        continue;
      endif

      Ui = zeros(size(Z));
      Ui = diag(D2) * Z;
      status.cpu.solve = cputime() - start;
    until (status.gmres.flag == 0 || status.do_fact)

    U(:, i) = callback("post", opt_sol, i, omega(i), Ui, status);
  endfor
endfunction

function Y = mat_product(M, D, K, D1, D2, ASfact, omega, Z)
  Y = diag(D2) * Z;
  Y = ASfact \ (diag(D1) * ((-omega^2) * (M * Y) + (1j * omega) * (D * Y) + K * Y));
endfunction
