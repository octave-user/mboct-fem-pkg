%!test
%! try
%! ## TEST 1
%! close all;
%! material.E = 210000e6;
%! material.nu = 0.3;
%! material.rho = 7850;
%! Fx = 0;
%! Fy = 0;
%! Fz = 30000;
%! h = 10e-3;
%! a = 500e-3;
%! geometry.l = 1000e-3;
%! geometry.w = 10e-3;
%! geometry.h = 50e-3;
%! mesh_size.num_elem_l = ceil(geometry.l / h);
%! mesh_size.num_elem_w = ceil(geometry.w / h);
%! mesh_size.num_elem_h = ceil(geometry.h / h);
%! f = [ Fx; Fy; Fz ];
%! [mesh, load_case] = fem_pre_mesh_cube_create(geometry, mesh_size, material, f);
%! idx_node = find(mesh.nodes(:, 1) == geometry.l & mesh.nodes(:, 2) == 0 & mesh.nodes(:, 3) == 0);
%! [dof_map] = fem_ass_dof_map(mesh, load_case);
%! opt_sol.pre_scaling = true;
%! t = 0.01:0.01:1;
%! sol_incr = fem_sol_static_incr(mesh, dof_map, load_case, t, opt_sol);
%! [mat_ass.K, ...
%!  mat_ass.R] = fem_ass_matrix(mesh, ...
%!                              dof_map, ...
%!                              [FEM_MAT_STIFFNESS, ...
%!                               FEM_VEC_LOAD_CONSISTENT], ...
%!                              load_case);
%! sol_lin = fem_sol_static(mesh, dof_map, mat_ass, load_case);

%! for i=1:3
%!   figure("visible", "off");
%!   hold on;
%!   plot(t, sol_incr.def(idx_node, i, :)(:), sprintf("-;incr %s;", {"ux", "uy", "uz"}{i}));
%!   plot(t, t * sol_lin.def(idx_node, i), sprintf("-x;lin %s;", {"ux", "uy", "uz"}{i}));
%! endfor
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
