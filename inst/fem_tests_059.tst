## fem_tests.m:59
%!test
%! try
%! ## TEST 59
%! close all;
%! material.E = 210000e6;
%! material.nu = 0.3;
%! material.rho = 7850;
%! num_modes = int32(10);
%! h = 10e-3;
%! geometry.l = 1000e-3;
%! geometry.w = 20e-3;
%! geometry.h = 20e-3;
%! r = linspace(0.05 * geometry.l, 0.97 * geometry.l, 100);
%! omegay = 1000;
%! sigmaxx_ref = material.rho * omegay^2 * (geometry.l^2 - r.^2) / 2;
%! mesh_size.num_elem_l = ceil(geometry.l / h);
%! mesh_size.num_elem_w = ceil(geometry.w / h);
%! mesh_size.num_elem_h = ceil(geometry.h / h);
%! f = zeros(3, 1);
%! [mesh, load_case_dof] = fem_pre_mesh_cube_create(geometry, mesh_size, material, f);
%! mesh.nodes(:, 2) -= 0.5 * geometry.w;
%! mesh.nodes(:, 3) -= 0.5 * geometry.h;
%! xi = mesh.nodes(:, 1)(mesh.elements.iso8);
%! yi = mesh.nodes(:, 2)(mesh.elements.iso8);
%! zi = mesh.nodes(:, 3)(mesh.elements.iso8);
%! e1 = [1; 0.7; 0.9];
%! e2 = [0.8; -0.5; 0.6];
%! e3 = cross(e1, e2);
%! e2 = cross(e3, e1);
%! R1 = [e1, e2, e3];
%! R1 *= diag(1 ./ norm(R1, "cols"));
%! dof_map = fem_ass_dof_map(mesh, load_case_dof);
%! load_case(1).omega = [0; omegay; 0];
%! load_case(1).omega = R1 * load_case.omega;
%! load_case(2).omega = -load_case(1).omega;
%! mesh.nodes = [mesh.nodes(:, 1:3) * R1.', zeros(rows(mesh.nodes), 3)];
%! [mat_ass.M, ...
%!  mat_ass.K, ...
%!  mat_ass.R] = fem_ass_matrix(mesh, ...
%!                              dof_map, ...
%!                              [FEM_MAT_MASS, ...
%!                               FEM_MAT_STIFFNESS, ...
%!                               FEM_VEC_LOAD_CONSISTENT], ...
%!                              load_case);
%! sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%! sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                  dof_map, ...
%!                                  [FEM_VEC_STRESS_CAUCH], ...
%!                                  load_case, ...
%!                                  sol_stat);
%! for i=1:2
%!   load_case(i).tau0.iso8 = sol_stat.stress.tau.iso8(:, :, :, i);
%! endfor
%! [mat_ass.KTAU0] = fem_ass_matrix(mesh, ...
%!                                dof_map, ...
%!                                [FEM_MAT_STIFFNESS_TAU0], ...
%!                                load_case);
%! sol_eig = fem_sol_modal(mesh, dof_map, mat_ass, num_modes);
%! mat_ass.K += mat_ass.KTAU0;
%! sol_eig_preload = fem_sol_modal(mesh, dof_map, mat_ass, num_modes);
%! TAU = zeros(3, 3, rows(sol_stat.stress.taum.iso8), columns(sol_stat.stress.taum.iso8));
%! for k=1:size(TAU, 4)
%!   for j=1:size(TAU, 3)
%!     for i=1:3
%!       TAU(i, i, j, k) = sol_stat.stress.taum.iso8(j, k, i);
%!     endfor
%!     TAU(1, 2, j, k) = TAU(2, 1, j, k) = sol_stat.stress.taum.iso8(j, k, 4);
%!     TAU(2, 3, j, k) = TAU(3, 2, j, k) = sol_stat.stress.taum.iso8(j, k, 5);
%!     TAU(1, 3, j, k) = TAU(3, 1, j, k) = sol_stat.stress.taum.iso8(j, k, 6);
%!     TAU(:, :, j, k) = R1.' * TAU(:, :, j, k) * R1;
%!   endfor
%! endfor
%! sigmaxx = griddata3(xi, yi, zi, reshape(TAU(1, 1, :, :), size(xi)), r, zeros(size(r)), zeros(size(r)));
%! tol = 1e-3;
%! assert_simple(sigmaxx(:), sigmaxx_ref(:), tol * max(abs(sigmaxx_ref)));
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
