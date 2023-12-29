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
## @deftypefn {Function File} @var{U} = fem_sol_harmonic_modal2(@var{dgen}, @var{kgen}, @var{rgen}, @var{Phi}, @var{omega})
## Compute the harmonic solution based on modal damping @var{dgen}, modal stiffness @var{kgen} modal excitation @var{rgen} and mode shapes @var{Phi}.
## It is assumed, that the modal mass @var{mgen} is equal to ones(size(kgen)).
##
## @var{dgen} @dots{} Modal damping returned from fem_sol_modes_scale2
##
## @var{kgen} @dots{} Modal stiffness returned from fem_sol_modes_scale2
##
## @var{rgen} @dots{} Modal excitation returned from fem_sol_modes_scale2
##
## @var{Phi} @dots{} Real mode shapes returned from fem_sol_modal
##
## @var{omega} @dots{} Angular velocity of harmonic excitation @var{rgen}
##
## @seealso{fem_sol_modal, fem_sol_modes_scale2}
## @end deftypefn

function U = fem_sol_harmonic_modal2(dgen, kgen, rgen, Phi, omega)
  U = zeros(rows(Phi), numel(omega), columns(rgen));

  for j=1:columns(rgen)
    for i=1:columns(Phi)
      qi = rgen(i, j) ./ (-omega.^2 + 1j * omega * dgen(i) + kgen(i)); ## assume mgen == 1
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
%! omega = linspace(0, 2 * pi * 10, 100000);
%! num_modes = 2;
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
%! solvers = {"pastix", "pardiso", "mumps", "chol", "umfpack", "lu", "mldivide"};
%! for i=1:numel(solvers)
%!   fprintf(stderr, "linear solver: \"%s\"\n", solvers{i});
%!   if (~fem_sol_check_func(solvers{i}))
%!     fprintf(stderr, "solver %s is not available\n", solvers{i});
%!     continue;
%!   endif
%!   opt_sol.solver = solvers{i};
%!   opt_sol.verbose = int32(0);
%!   opt_sol.refine_max_iter = int32(100);
%!   opt_sol.pre_scaling = true;
%!   opt_sol.scaling = true;
%!   opt_eig.disp = int32(0);
%!   opt_eig.maxit = int32(10);
%!   switch (solvers{i})
%!   case "pardiso"
%!     opt_sol.symmetric = false; ## FIXME: Pardiso is not able to solve it in symmetric mode with weighted matching enabled
%!     opt_sol.weighted_matching = false;
%!   otherwise
%!     opt_sol.symmetric = true;
%!     opt_sol.weighted_matching = true;
%!   endswitch
%!   opt_sol.algorithm = "generic";
%!   [sol, Phi] = fem_sol_modal(mesh, dof_map, mat_ass, num_modes, opt_sol, opt_eig);
%!   [dgen, kgen, rgen] = fem_sol_modes_scale2(mat_ass.M, mat_ass.D, mat_ass.K, Phi, mat_ass.R);
%!   U = fem_sol_harmonic_modal2(dgen, kgen, rgen, Phi(dof_map.ndof(1, 1:2), :), omega);
%!   Uref = [r1 ./ (-omega.^2 * m1 + 1j * omega * d1 + k1);
%!           r2 ./ (-omega.^2 * m2 + 1j * omega * d2 + k2)];
%!   omega01 = sqrt(k1 / m1);
%!   omega02 = sqrt(k2 / m2);
%!   delta1 = d1 / (2 * m1);
%!   delta2 = d2 / (2 * m2);
%!   D1 = delta1 / omega01;
%!   D2 = delta2 / omega02;
%!   Dref = [D1, D1, D2, D2];
%!   lambdaref = 1j * [omega01, omega02];
%!   [~, idx] = sortrows([imag(lambdaref)(:), real(lambdaref)(:)]);
%!   lambdaref = lambdaref(idx);
%!   Dref = min(1, Dref(idx));
%!   tol = eps^0.9;
%!   assert_simple(U, Uref, tol * norm(Uref));
%!   assert_simple(imag(sol.lambda), imag(lambdaref), tol * norm(imag(lambdaref)));
%! endfor
