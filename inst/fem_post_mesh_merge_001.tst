## fem_post_mesh_merge.m:01
%!test
%! try
%! close all;
%! c = 10e-3;
%!
%! material1.E = 210000e6;
%! material1.nu = 0.3;
%! material1.rho = 7850;
%!
%! material2.E = material1.E;
%! material2.nu = material1.nu;
%! material2.rho = material1.rho;
%!
%! h1 = 4.5e-3;
%! h2 = 5.5e-3;
%!
%! geometry1.l = 20e-3;
%! geometry1.w = 20e-3;
%! geometry1.h = 40e-3;
%! geometry2.l = 60e-3;
%! geometry2.w = 30e-3;
%! geometry2.h = 50e-3;
%!
%! mesh_size1.num_elem_l = ceil(geometry1.l / h1);
%! mesh_size1.num_elem_w = ceil(geometry1.w / h1);
%! mesh_size1.num_elem_h = ceil(geometry1.h / h1);
%! mesh_size2.num_elem_l = ceil(geometry2.l / h2);
%! mesh_size2.num_elem_w = ceil(geometry2.w / h2);
%! mesh_size2.num_elem_h = ceil(geometry2.h / h2);
%!
%! [data(1).mesh] = fem_pre_mesh_cube_create(geometry1, mesh_size1, material1, zeros(3, 1));
%!
%! data(1).mesh.nodes(:, 2) -= 0.5 * geometry1.w;
%! data(1).mesh.nodes(:, 3) -= 0.5 * geometry1.h;
%!
%! [data(2).mesh] = fem_pre_mesh_cube_create(geometry2, mesh_size2, material2, zeros(3, 1));
%!
%! data(2).mesh.nodes(:, 1) += geometry1.l + c;
%! data(2).mesh.nodes(:, 2) -= 0.5 * geometry2.w;
%! data(2).mesh.nodes(:, 3) -= 0.5 * geometry2.h;
%!
%! for i=1:numel(data)
%!  data(i).load_case.locked_dof = false(size(data(i).mesh.nodes));
%!  data(i).dof_map = fem_ass_dof_map(data(i).mesh, data(i).load_case);
%! endfor
%!
%! [mesh, dof_map] = fem_post_mesh_merge(data);
%! figure("visible", "off");
%! fem_post_sol_plot(mesh);
%! xlabel("x [m]");
%! ylabel("y [m]");
%! zlabel("z [m]");
%! grid on;
%! grid minor on;
%! title("merged mesh");
%! view(30, 30);
%! figure_list();
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
