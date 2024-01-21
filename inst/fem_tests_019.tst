## fem_tests.m:19
%!test
%! ##########################################################################################
%! ## TEST 19: Test case for sfncon4
%! ##########################################################################################
%! close all;
%! p = 1.25e6;
%! a = 15e-3;
%! b = 25e-3;
%! c = 20e-3;
%! A = 20e-3;
%! B = 30e-3;
%! C = 40e-3;
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
%! mesh.nodes = [X, zeros(rows(X), 3)];
%! mesh.elements.iso8 = int32([1:8;
%!                             9:16]);
%! mesh.elements.sfncon4.master = int32(9:12);
%! mesh.elements.sfncon4.slave = int32(5:8).';
%! mesh.elements.sfncon4.maxdist = sqrt(eps) * max(norm(X, "rows"));
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
%! sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%! sol_eig = fem_sol_modal(mesh, dof_map, mat_ass, num_modes);
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
