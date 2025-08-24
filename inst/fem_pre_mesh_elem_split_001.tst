## fem_pre_mesh_elem_split.m:01
%!test
%! try
%! pkg load mbdyn_util_oct;
%! do_plot = false;
%! if (do_plot)
%!   close all;
%! endif
%! a = 70e-3;
%! b = 20e-3;
%! c = 10e-3;
%! d = 120e-3;
%! rho = 7850;
%! m = 3 / 8 * a * b * c * rho;
%! tol = eps^0.9;
%!
%! X = [0.5 * a, 0.5 * b, c;
%!            0, 0.5 * b, c;
%!            0,       0, c;
%!            a,       0, c;
%!      0.5 * a, 0.5 * b, 0;
%!            0, 0.5 * b, 0;
%!            0,       0, 0;
%!            a,       0, 0];
%! Phi1 = [0, 20, 45, 330] * pi / 180;
%! Phi2 = [0, 60, 270, 285] * pi / 180;
%! Phi3 = [0, -15, 30, -195] * pi / 180;
%! options.iso8.max_angle = 95 * pi / 180;
%! rand("seed", 0);
%! idx_elem = int32([1:8;
%!                   5, 1, 4, 8, 6, 2, 3, 7;
%!                   4, 3, 7, 8, 1, 2, 6, 5;
%!                   3, 2, 6, 7, 4, 1, 5, 8;
%!                   1, 4, 8, 5, 2, 3, 7, 6;
%!                   1, 5, 6, 2, 4, 8, 7, 3;
%!                   3, 2, 6, 7, 4, 1, 5, 8;
%!                   8, 7, 6, 5, 4, 3, 2, 1;
%!                   7, 6, 5, 8, 3, 2, 1, 4;
%!                   5, 8, 7, 6, 1, 4, 3, 2]);
%! idx_mesh = int32(0);
%! for i=1:rows(idx_elem)
%! for j=1:length(Phi1)
%!   R1 = euler123_to_rotation_matrix([Phi1(j); Phi2(j); Phi3(j)]);
%!   mesh(++idx_mesh).nodes = [X * R1.', zeros(rows(X), 3)];
%!   mesh(idx_mesh).elements.iso8 = idx_elem(i, :);
%!   mesh(idx_mesh).materials.iso8 = int32([1]);
%!   load_case(idx_mesh).locked_dof = false(rows(mesh(idx_mesh).nodes), 6);
%!   E = 210000e6;
%!   nu = 0.3;
%!   mesh(idx_mesh).material_data.rho = rho;
%!   mesh(idx_mesh).material_data.C = fem_pre_mat_isotropic(E, nu);
%!   dof_map(idx_mesh) = fem_ass_dof_map(mesh(idx_mesh), load_case(idx_mesh));
%!   [rndval, idx_node] = sort(rand(rows(mesh(idx_mesh).nodes), 1));
%!   idx_node_inv(idx_node) = 1:rows(mesh(idx_mesh).nodes);
%!   for k=1:rows(mesh(j).elements.iso8)
%!     mesh(idx_mesh).elements.iso8(k, :) = idx_node_inv(mesh(idx_mesh).elements.iso8(k, :));
%!   endfor
%!   mesh(idx_mesh).nodes = mesh(idx_mesh).nodes(idx_node, :);
%!   mesh_split(idx_mesh) = fem_pre_mesh_elem_split(mesh(idx_mesh), options);
%!   mesh_split(idx_mesh).nodes(:, 3) += d;
%!   mat_ass(idx_mesh).m = fem_ass_matrix(mesh(idx_mesh), ...
%!                                        dof_map(idx_mesh), ...
%!                                        [FEM_SCA_TOT_MASS], ...
%!                                        load_case(idx_mesh));
%!   mat_ass_split(idx_mesh).m = fem_ass_matrix(mesh_split(idx_mesh), ...
%!                                              dof_map(idx_mesh), ...
%!                                              [FEM_SCA_TOT_MASS], ...
%!                                              load_case(idx_mesh));
%!   if (do_plot)
%!   figure("visible","off");
%!   fem_post_sol_plot(mesh(idx_mesh));
%!   fem_post_sol_plot(mesh_split(idx_mesh));
%!   title(sprintf("mesh %d: elements before %d, elements after %d", j, rows(mesh(idx_mesh).elements.iso8), rows(mesh_split(idx_mesh).elements.iso8)));
%!   xlabel("x [m]");
%!   ylabel("y [m]");
%!   zlabel("z [m]");
%!   grid on;
%!   grid minor on;
%!   endif
%! endfor
%! endfor
%! assert_simple(all(abs([mat_ass.m] - m) < tol * m));
%! assert_simple(all(abs([mat_ass_split.m] - m) < tol * m));
%! if (do_plot)
%! figure_list();
%! endif
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
