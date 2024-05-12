## fem_post_sol_plot.m:01
%!test
%! try
%! ## TEST 1
%! close all;
%! material.E = 210000e6;
%! material.nu = 0.3;
%! material.rho = 7850;
%! h = 10e-3;
%! geometry.l = 60e-3;
%! geometry.w = 20e-3;
%! geometry.h = 50e-3;
%! mesh_size.num_elem_l = ceil(geometry.l / h);
%! mesh_size.num_elem_w = ceil(geometry.w / h);
%! mesh_size.num_elem_h = ceil(geometry.h / h);
%! f = [1; 0; 0];
%! [mesh, load_case] = fem_pre_mesh_cube_create(geometry, mesh_size, material, f);
%! figure("visible", "off");
%! fem_post_sol_plot(mesh);
%! view(30, 30);
%! xlabel("x [m]");
%! ylabel("y [m]");
%! zlabel("z [m]");
%! view(15, 15);
%! grid on;
%! grid minor on;
%! title("undeformed cube mesh");
%! figure_list();
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
