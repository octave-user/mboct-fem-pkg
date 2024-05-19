## fem_pre_load_case_merge.m:03
%!test
%! try
%! close all;
%! a = 70e-3;
%! b = 20e-3;
%! c = 10e-3;
%! d = 140e-3;
%! scale_def = 20e-3;
%! Fx = 100;
%! Fy = 150;
%! Fz = -230;
%! X = [ a,  0.5 * b,  0.5 * c;
%!       0,  0.5 * b,  0.5 * c;
%!       0, -0.5 * b,  0.5 * c;
%!       a, -0.5 * b,  0.5 * c;
%!       a,  0.5 * b, -0.5 * c;
%!       0,  0.5 * b, -0.5 * c;
%!       0, -0.5 * b, -0.5 * c;
%!       a, -0.5 * b, -0.5 * c];
%! mesh.nodes = [X, zeros(rows(X), 3)];
%! mesh.elements.iso8 = int32(1:8);
%! mesh.materials.iso8 = int32(1);
%! E = 210000e6;
%! nu = 0.3;
%! mesh.material_data.rho = 7850;
%! mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%! load_case1.locked_dof = false(rows(mesh.nodes), 6);
%! load_case1.locked_dof([2; 3; 6; 7], 1:3) = true;
%! load_case1.loaded_nodes = int32([1; 4; 5; 8]);
%! load_case1.loads = Fx * [0.25, 0, 0, 0, 0, 0;
%!                          0.25, 0, 0, 0, 0, 0;
%!                          0.25, 0, 0, 0, 0, 0;
%!                          0.25, 0, 0, 0, 0, 0];
%! load_case2.loaded_nodes = int32([1; 4; 5; 8]);
%! load_case2.loads = Fy * [0, 0.25, 0, 0, 0, 0;
%!                          0, 0.25, 0, 0, 0, 0;
%!                          0, 0.25, 0, 0, 0, 0;
%!                          0, 0.25, 0, 0, 0, 0];
%! load_case3.loaded_nodes = int32([1; 4; 5; 8]);
%! load_case3.loads = Fz * [0, 0, 0.25, 0, 0, 0;
%!                          0, 0, 0.25, 0, 0, 0;
%!                          0, 0, 0.25, 0, 0, 0;
%!                          0, 0, 0.25, 0, 0, 0];
%! load_case = fem_pre_load_case_merge(load_case1, load_case2, load_case3);
%! dof_map = fem_ass_dof_map(mesh, load_case(1));  
%! [mat_ass.K, ...
%!  mat_ass.R] = fem_ass_matrix(mesh, ...
%!                              dof_map, ...
%!                              [FEM_MAT_STIFFNESS, ...
%!                               FEM_VEC_LOAD_CONSISTENT], ...
%!                              load_case);
%! sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%! for i=1:size(sol_stat.def, 3)
%!   figure("visible", "off");
%!   fem_post_sol_plot(mesh, sol_stat, scale_def / max(max(abs(sol_stat.def(:, :, i)))), i);
%!   xlabel("x [m]");
%!   ylabel("y [m]");
%!   zlabel("z [m]");
%!   grid on;
%!   grid minor on;
%!   title(sprintf("load case %d", i));
%!   view(30, 30);
%! endfor
%! figure_list();
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
