## fem_tests.m:16
%!test
%! try
%! ## TEST 16
%! material.E = 210000e6;
%! material.nu = 0.3;
%! material.rho = 7850;
%! geometry.l = 1000e-3;
%! geometry.w = 10e-3;
%! geometry.h = 50e-3;
%! elem_size = 5e-3;
%! mesh_size.num_elem_l = ceil(geometry.l / elem_size);
%! mesh_size.num_elem_w = ceil(geometry.w / elem_size);
%! mesh_size.num_elem_h = ceil(geometry.h / elem_size);
%! number_of_modes = 10;
%! number_of_modes_disp = 3;
%! plot_def = false;
%! f = [ 0; 0; 15000];
%! [mesh, load_case] = fem_pre_mesh_cube_create(geometry, mesh_size, material, f);
%! [dof_map] = fem_ass_dof_map(mesh, load_case);
%! [mat_ass_sym.K, ...
%!  mat_ass_sym.R, ...
%!  mat_ass.K, ...
%!  mat_ass.R] = fem_ass_matrix(mesh, ...
%!                              dof_map, ...
%!                              [FEM_MAT_STIFFNESS, ...
%!                               FEM_VEC_LOAD_CONSISTENT, ...
%!                               FEM_MAT_STIFFNESS, ...
%!                               FEM_VEC_LOAD_CONSISTENT], ...
%!                              load_case);
%! opt_sol.solver = "chol";
%! [sol_stat_sym] = fem_sol_static(mesh, dof_map, mat_ass_sym, opt_sol);
%! [sol_stat] = fem_sol_static(mesh, dof_map, mat_ass);
%! assert_simple(sol_stat_sym.def, sol_stat.def, sqrt(eps) * max(norm(sol_stat.def, "rows")));
%! z = linspace(0,geometry.l,100);
%! I = [ geometry.w * geometry.h, geometry.h * geometry.w^3 / 12, geometry.w * geometry.h^3 / 12 ];
%! y(1,:) = f(1) * geometry.l / ( material.E * I(1) ) * ( 1 - z / geometry.l );
%! for i=2:3
%!   y(i,:) = f(i) * geometry.l^3 / ( 6 * material.E * I(i) ) * ( 2 - 3 * z / geometry.l + ( z / geometry.l ).^3 );
%! endfor
%! uz = griddata3(mesh.nodes(:, 1), ...
%!                mesh.nodes(:, 2), ...
%!                mesh.nodes(:, 3), ...
%!                sol_stat.def(:, 3),  ...
%!                geometry.l - z, ...
%!                zeros(size(z)), ...
%!                zeros(size(z)), ...
%!                "linear");
%! assert_simple(uz, y(3, :).', 1e-2 * max(abs(y(3, :))));
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
