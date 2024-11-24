## fem_sol_harmonic_modal.m:01
%!test
%! try
%! ## lumped model
%! m1 = 27;
%! k1 = 1700;
%! k2 = 2900;
%! k3 = 1800;
%! d1 = 87;
%! r1 = 10;
%! ## condensed model
%! mt = m1;
%! dt = d1;
%! kt = k1 + 1/(1/k2 + 1/k3);
%! rt = r1;
%! omega = linspace(0, 2 * pi * 10, 10000);
%! num_modes = 2;
%! mat_ass.M = [m1,  0,  0,  0;
%!              0,   0,  0,  0;
%!              0,   0,  0,  0;
%!              0,   0,  0,  0];
%! mat_ass.K = [k1 + k2, -k2,   0,  0;
%!                  -k2,  k2,   0,  1;
%!                    0,   0,  k3, -1;
%!                    0,   1,  -1,  0];
%! mat_ass.D = [d1,   0,   0,  0;
%!               0,   0,   0,  0;
%!               0,   0,   0,  0;
%!               0,   0,   0,  0];
%! mat_ass.R = [r1;
%!               0;
%!               0;
%!               0];
%! mat_ass.M = sparse(mat_ass.M);
%! mat_ass.D = sparse(mat_ass.D);
%! mat_ass.K = sparse(mat_ass.K);
%! mat_ass.R = sparse(mat_ass.R);
%! mesh.nodes = zeros(1, 3);
%! dof_map.domain = FEM_DO_STRUCTURAL;
%! dof_map.totdof = rows(mat_ass.M);
%! dof_map.ndof = zeros(3, 6, "int32");
%! dof_map.ndof(1:3, 1) = 1:3;
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
%!   [~, idxlambda] = sortrows([imag(sol.lambda(:)), real(sol.lambda(:))], [-2, -1]);
%!   Phi = Phi(:, idxlambda);
%!   sol.def = sol.def(:,:, idxlambda);
%!   sol.lambda = sol.lambda(idxlambda);
%!   [Phi, h] = fem_sol_modes_scale(mat_ass.M, mat_ass.K, sol.lambda, Phi, mat_ass.R);
%!   U = fem_sol_harmonic_modal(h, sol.lambda, Phi(dof_map.ndof(1, 1), :), omega);
%!   Uref = [rt ./ (-omega.^2 * mt + 1j * omega * dt + kt)];
%!   lambdaref = -dt/(2*mt) + [1, -1] * sqrt((dt/(2*mt))^2 - kt/mt);
%!   omega0ref = sqrt(kt / mt);
%!   deltaref = dt / (2 * mt);
%!   Dref = deltaref / omega0ref;
%!   tol = eps^0.8;
%!   assert_simple(U, Uref, tol * norm(Uref));
%!   assert_simple(sol.lambda, lambdaref, tol * norm(lambdaref));
%!   assert_simple(sol.D, repmat(Dref, 1, 2), tol * norm(Dref));
%! endfor
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
