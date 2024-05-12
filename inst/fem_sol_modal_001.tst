## fem_sol_modal.m:01
%!test
%! try
%! close all;
%! material.E = 210000e6;
%! material.nu = 0.3;
%! material.rho = 7850;
%! h = 10e-3 / 2;
%! geometry.l = 1000e-3;
%! geometry.w = 10e-3;
%! geometry.h = 50e-3;
%! mesh_size.num_elem_l = ceil(geometry.l / h);
%! mesh_size.num_elem_w = ceil(geometry.w / h);
%! mesh_size.num_elem_h = ceil(geometry.h / h);
%! number_of_modes = 10;
%! f = zeros(3, 1);
%! [mesh, load_case] = fem_pre_mesh_cube_create(geometry, mesh_size, material, f);
%! [dof_map] = fem_ass_dof_map(mesh, load_case);
%! [mat_ass.M, ...
%!  mat_ass.K] = fem_ass_matrix(mesh, ...
%!                              dof_map, ...
%!                              [FEM_MAT_MASS, ...
%!                               FEM_MAT_STIFFNESS], ...
%!                              load_case);
%! [sol_eig] = fem_sol_modal(mesh, dof_map, mat_ass, number_of_modes);
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
