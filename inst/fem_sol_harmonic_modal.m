## Copyright (C) 2023(-2023) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} @var{U} = fem_sol_harmonic_modal(@var{h}, @var{lambda}, @var{Phi}, @var{omega})
## Compute the harmonic solution based on eigenvalues @var{lambda}, normalized mode shapes @var{Phi} and modal excitation @var{h}
##
## @var{h} @dots{} Modal excitation returned from fem_sol_modes_scale
##
## @var{lambda} @dots{} Complex eigenvalues returned from fem_sol_modal_damped
##
## @var{Phi} @dots{} Complex mode shapes returned from fem_sol_modal_damped
##
## @var{omega} @dots{} Angular velocity of harmonic excitation @var{h}
##
## @seealso{fem_sol_modal_damped}
## @end deftypefn

function U = fem_sol_harmonic_modal(h, lambda, Phi, omega)
  if (nargout > 1 || nargin ~= 4)
    print_usage();
  endif

  if (~(rows(h) == columns(Phi) && columns(Phi) == numel(lambda)))
    error("size of h, lambda and Phi is not consistent");
  endif

  if (rows(omega) ~= 1)
    error("omega must be a row vector");
  endif

  U = zeros(rows(Phi), numel(omega), columns(h));

  for j=1:columns(h)
    for i=1:columns(Phi)
      qi = h(i, j) ./ (lambda(i) - 1j * omega);
      U(:, :, j) += Phi(:, i) * qi;
    endfor
  endfor
endfunction

%!test
%! m1 = 0.5;
%! m2 = 30;
%! m3 = 0;
%! m4 = 0;
%! m5 = 0;
%! m6 = 0;
%! k1 = 1e-3;
%! k2 = 10000;
%! k3 = 1e10;
%! k4 = 1e-10;
%! d1 = 50;
%! d2 = 1e2;
%! d3 = 0;
%! d4 = 0;
%! r1 = 1e6;
%! r2 = -0.5;
%! r3 = 1e15;
%! r4 = -1e10;
%! r5 = 1e10;
%! r6 = 1e7;
%! c5 = c6 = norm([k1, k2, k3, k4]);
%! omega = linspace(0, 2 * pi * 10, 10000);
%! num_modes = 4;
%! mat_ass.M = [m1,  0,  0,  0,  0,  0;
%!              0,  m2,  0,  0,  0,  0;
%!              0,   0, m3,  0,  0,  0;
%!              0,   0,  0, m4,  0,  0;
%!              0,   0,  0,  0, m5,  0;
%!              0,   0,  0,  0,  0, m6];
%! mat_ass.K = [k1, 0,  0,  0,  0,  0;
%!              0, k2,  0,  0,  0,  0;
%!              0,  0, k3,  0, c5,  0;
%!              0,  0,  0, k4,  0, c6;
%!              0,  0, c5,  0,  0,  0;
%!              0,  0,  0, c6,  0,  0];
%! mat_ass.D = [d1,   0,   0,  0, 0, 0;
%!              0,  d2,   0,  0, 0, 0;
%!              0,   0,  d3,  0, 0, 0;
%!              0,   0,   0, d4, 0, 0;
%!              0,   0,   0,  0, 0, 0;
%!              0,   0,   0,  0, 0, 0];
%! mat_ass.R = [r1;
%!              r2;
%!              r3;
%!              r4;
%!              r5;
%!              r6];
%! mat_ass.M = sparse(mat_ass.M);
%! mat_ass.D = sparse(mat_ass.D);
%! mat_ass.K = sparse(mat_ass.K);
%! mat_ass.R = sparse(mat_ass.R);
%! mesh.nodes = zeros(2, 6);
%! dof_map.domain = FEM_DO_STRUCTURAL;
%! dof_map.totdof = rows(mat_ass.M);
%! dof_map.ndof = zeros(4, 6, "int32");
%! dof_map.ndof(1, 1:2) = [1:2];
%! dof_map.ndof(2, 1:2) = [3:4];
%! solvers = {"pastix", "pardiso", "umfpack", "lu", "mldivide"};
%! for i=1:numel(solvers)
%!   opt_sol.solver = solvers{i};
%!   opt_sol.refine_max_iter = int32(100);
%!   opt_sol.pre_scaling = true;
%!   opt_sol.verbose = int32(0);
%!   opt_sol.scaling = false;
%!   opt_sol.weighted_matching = false; ## FIXME: Pardiso is not able to solve it with weighted matching enabled
%!   opt_eig.disp = int32(0);
%!   opt_eig.maxit = int32(10);
%!   opt_sol.symmetric = true; 
%!   [sol, Phi] = fem_sol_modal_damped(mesh, dof_map, mat_ass, num_modes, opt_sol, opt_eig);
%!   [Phi, h] = fem_sol_modes_scale(mat_ass.M, mat_ass.K, sol.lambda, Phi, mat_ass.R);
%!   U = fem_sol_harmonic_modal(h, sol.lambda, Phi(dof_map.ndof(1, 1:2), :), omega);
%!   Uref = [r1 ./ (-omega.^2 * m1 + 1j * omega * d1 + k1);
%!           r2 ./ (-omega.^2 * m2 + 1j * omega * d2 + k2)];
%!   lambda1ref = -d1/(2*m1) + [1, -1] * sqrt((d1/(2*m1))^2 - k1/m1);
%!   lambda2ref = -d2/(2*m2) + [1, -1] * sqrt((d2/(2*m2))^2 - k2/m2);
%!   omega01 = sqrt(k1 / m1);
%!   omega02 = sqrt(k2 / m2);
%!   delta1 = d1 / (2 * m1);
%!   delta2 = d2 / (2 * m2);
%!   D1 = delta1 / omega01;
%!   D2 = delta2 / omega02;
%!   Dref = [D1, D1, D2, D2];
%!   lambdaref = [lambda1ref, lambda2ref];
%!   [~, idx] = sortrows([imag(lambdaref)(:), real(lambdaref)(:)]);
%!   lambdaref = lambdaref(idx);
%!   Dref = min(1, Dref(idx));
%!   tol = eps^0.9;
%!   assert(U, Uref, tol * norm(Uref));
%!   assert(sol.lambda, lambdaref, tol * norm(lambdaref));
%!   assert(sol.D, Dref, tol * norm(Dref));
%! endfor
