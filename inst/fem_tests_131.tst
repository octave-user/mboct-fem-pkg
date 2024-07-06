## fem_tests.m:19
%!test
%! try
%! ##########################################################################################
%! ## TEST 19: Test case for sfncon4s
%! ##########################################################################################
%! close all;
%! p = 1.25e6;
%! a = 25e-3;
%! b = 15e-3;
%! c = 15e-3;
%! A = 35e-3;
%! B = 20e-3;
%! C = 20e-3;
%! function K = stiffness(R, h0, dA, k, sigma_delta, sigma0)
%!   c1 = 4.4086e-5;
%!   c2 = 6.804;
%!   Hmax = 4.;
%!   H0 = h0 / sigma_delta;
%!   F5_2 = (H0 <= Hmax) * c1 * (Hmax - H0).^c2;
%!   dF5_2_dH0 = -(Hmax - H0).^(c2 - 1.) * c1 * c2;
%!   p = real(k) * F5_2;
%!   kh = -k / sigma_delta * dF5_2_dH0;
%!   K = diag([p * sigma0, p * sigma0, kh]);
%! endfunction
%! scale = 10e-3;
%! num_modes = int32(6);
%! do_plot = false;
%! X = [-0.5 * a, -0.5 * b, c + C;
%!      0.5 * a, -0.5 * b, c + C;
%!      0.5 * a,  0.5 * b, c + C;
%!      -0.5 * a,  0.5 * b, c + C;
%!      -0.5 * a, -0.5 * b, C;
%!      0.5 * a, -0.5 * b, C;
%!      0.5 * a,  0.5 * b, C;
%!      -0.5 * a,  0.5 * b, C;
%!      -0.5 * A, -0.5 * B, C;
%!      0.5 * A, -0.5 * B, C;
%!      0.5 * A,  0.5 * B, C;
%!      -0.5 * A,  0.5 * B, C;
%!      -0.5 * A, -0.5 * B, 0;
%!      0.5 * A, -0.5 * B, 0;
%!      0.5 * A,  0.5 * B, 0;
%!      -0.5 * A,  0.5 * B, 0];
%! e1 = [1; 0.5; 0.3];
%! e2 = [-0.5; 0.9; 1];
%! e3 = cross(e1, e2);
%! e2 = cross(e3, e1);
%! R = [e1, e2, e3];
%! R *= diag(1 ./ norm(R, "cols"));
%! R = eye(3);
%! X *= R.';
%! mesh.nodes = [X, zeros(rows(X), 3)];
%! mesh.elements.iso8 = int32([1:8;
%!                             9:16]);
%! mesh.elements.sfncon4s.master = int32(9:12);
%! mesh.elements.sfncon4s.slave = int32(5:8).';
%! mesh.elements.sfncon4s.maxdist = sqrt(eps) * max(norm(X, "rows"));
%! sigma0 = 1e5;
%! sigma_delta = 1e-6;
%! k = 100000e6;
%! mesh.elements.sfncon4s.k = @(R, h0, dA) stiffness(R, h0, dA, k, sigma_delta, sigma0);
%! E(1) = 210000e6;
%! nu(1) = 0.3;
%! mesh.material_data(1).rho = 7850;
%! E(2) = 70000e6;
%! nu(2) = 0.3;
%! mesh.material_data(2).rho = 2700;
%! for i=1:numel(mesh.material_data)
%!   mesh.material_data(i).C = fem_pre_mat_isotropic(E(i), nu(i));
%! endfor
%! mesh.materials.iso8 = int32([1; 2]);
%! load_case_dof.locked_dof = false(size(mesh.nodes));
%! load_case_dof.locked_dof(13:16, 1:3) = true;
%! load_case(1).pressure.iso4.elements = int32([1,2,3,4;
%!                                              9,10,11,12]);
%! load_case(1).pressure.iso4.p = repmat(p, 2, 4);
%! load_case(1).g = [0; 0; -9.81];
%! load_case(2).pressure.iso4.elements = int32([5,6,7,8;
%!                                              13,14,15,16]);
%! load_case(2).pressure.iso4.p = repmat(p, 2, 4);
%! load_case(2).g = [0; 0; -9.81];
%! dof_map = fem_ass_dof_map(mesh, load_case_dof);
%! [mat_ass.K, ...
%!  mat_ass.M, ...
%!  mat_ass.R, ...
%!  mat_ass.dm, ...
%!  mat_ass.S, ...
%!  mat_ass.J, ...
%!  mat_ass.surface, ...
%!  mat_ass.mat_info, ...
%!  mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_MAT_STIFFNESS, ...
%!                                       FEM_MAT_MASS, ...
%!                                       FEM_VEC_LOAD_CONSISTENT, ...
%!                                       FEM_SCA_TOT_MASS, ...
%!                                       FEM_VEC_INERTIA_M1, ...
%!                                       FEM_MAT_INERTIA_J, ...
%!                                       FEM_VEC_SURFACE_AREA], ...
%!                                      load_case);
%! opt_sol.solver = "umfpack";
%! opt_sol.pre_scaling = true;
%! opt_sol.refine_max_iter = 100;
%! sol_stat = fem_sol_static(mesh, dof_map, mat_ass, opt_sol);
%! sol_eig = fem_sol_modal(mesh, dof_map, mat_ass, num_modes, opt_sol);
%! Aref = 2 * (a * b + A * B);
%! assert_simple(sum(sum(mat_ass.surface.iso4)), Aref, eps^0.8 * Aref);
%! if (do_plot)
%!   figure("visible", "off");
%!   fem_post_sol_plot(mesh, sol_stat, scale / max(norm(sol_stat.def(:, 1:3), "rows")), 1);
%!   xlabel("x [m]");
%!   ylabel("y [m]");
%!   zlabel("z [m]");
%!   grid on;
%!   grid minor on;
%!   title("static deflection");
%!   for i=1:numel(sol_eig.f)
%!     figure("visible", "off");
%!     fem_post_sol_plot(mesh, sol_eig, scale / max(norm(sol_eig.def(:, 1:3, i), "rows")), i);
%!     xlabel("x [m]");
%!     ylabel("y [m]");
%!     zlabel("z [m]");
%!     grid on;
%!     grid minor on;
%!     title(sprintf("mode %d: %.1fHz", i, sol_eig.f(i)));
%!   endfor
%!   figure_list();
%! endif
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
