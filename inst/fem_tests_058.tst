## fem_tests.m:58
%!test
%! try
%! ## TEST 58
%! close all;
%! number_of_modes = 5;
%! material.E = 210000e6;
%! material.nu = 0.3;
%! material.rho = 7850;
%! Fx = -1;
%! h = 10e-3/8;
%! geometry.l = 250e-3;
%! geometry.w = 10e-3;
%! geometry.h = 10e-3;
%! Iz = geometry.w * geometry.h^3 / 12;
%! Iy = geometry.h * geometry.w^3 / 12;
%! Imin = min(Iy, Iz);
%! ## Elastic Euler buckling (W. Beitz, K.H.Grote, 1997 Dubbel C42)
%! lK = 2 * geometry.l;
%! FK = pi^2 * material.E * Imin / lK^2;
%! mesh_size.num_elem_l = ceil(geometry.l / h);
%! mesh_size.num_elem_w = ceil(geometry.w / h);
%! mesh_size.num_elem_h = ceil(geometry.h / h);
%! e1 = [0.7; 0.2; 0.3];
%! e2 = [0.5; 0.6; 0.2];
%! e3 = cross(e1, e2);
%! e2 = cross(e3, e1);
%! R1 = [e1, e2, e3];
%! R1 = R1 * diag(1 ./ norm(R1, "cols"));
%! f = R1 * [ Fx; 0; 0 ];
%! [mesh, load_case] = fem_pre_mesh_cube_create(geometry, mesh_size, material, f);
%! mesh.nodes = [mesh.nodes(:, 1:3) * R1.', zeros(rows(mesh.nodes), 3)];
%! [dof_map] = fem_ass_dof_map(mesh, load_case);
%! dof_map.parallel.threads_ass = mbdyn_solver_num_threads_default();
%! [mat_ass.K, ...
%!  mat_ass.R] = fem_ass_matrix(mesh, ...
%!                              dof_map, ...
%!                              [FEM_MAT_STIFFNESS, ...
%!                               FEM_VEC_LOAD_CONSISTENT], ...
%!                              load_case);
%!
%! [sol_stat] = fem_sol_static(mesh, dof_map, mat_ass);
%!
%! sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                  dof_map, ...
%!                                  [FEM_VEC_STRESS_CAUCH], ...
%!                                  load_case, ...
%!                                  sol_stat);
%! load_case.tau0 = sol_stat.stress.tau;
%! mat_ass.KTAU0 = fem_ass_matrix(mesh, ...
%!                              dof_map, ...
%!                              [FEM_MAT_STIFFNESS_TAU0], ...
%!                              load_case);
%! opts.number_of_threads = dof_map.parallel.threads_ass;
%! opts.problem = "buckling";
%! [U, sol_buck.lambda] = fem_sol_eigs(mat_ass.K, mat_ass.KTAU0, number_of_modes, opts);
%! sol_buck.def = fem_post_def_nodal(mesh, dof_map, U);
%! tol = 1e-2;
%! assert_simple(sol_buck.lambda(1), FK, tol * abs(FK));
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
