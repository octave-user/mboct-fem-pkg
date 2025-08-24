## fem_tests.m:18
%!test
%! try
%! pkg load mbdyn_util_oct;
%! ##########################################################################################
%! ## TEST 18: Test case for pressure load
%! ##########################################################################################
%! a = 0.5;
%! b = 0.3;
%! c = 0.7;
%! p = 2.5e9;
%! do_plot = false;
%! scale = -0.5 * a;
%! X = [      0,       0,       0;
%!            a,       0,       0;
%!            0,       b,       0;
%!            0,       0,       c;
%!            0.5 * a,       0,       0;
%!            0.5 * a, 0.5 * b,       0;
%!            0, 0.5 * b,       0;
%!            0,       0, 0.5 * c;
%!            0.5 * a,       0, 0.5 * c;
%!            0, 0.5 * b, 0.5 * c].';
%! Phi1 = [0, 30, 120] * pi / 180;
%! Phi2 = [0, -45, 270] * pi / 180;
%! Phi3 = [0, 170, 310] * pi / 180;
%! if (do_plot)
%!   close all;
%! endif
%! for i=1:numel(Phi1)
%!   R1 = euler123_to_rotation_matrix([Phi1(i); Phi2(i); Phi3(i)]);
%!   data(i).mesh.nodes = [(R1 * X).', zeros(columns(X), 3)];
%!   data(i).mesh.elements.tet10 = int32(1:10);
%!   data(i).mesh.elements.joints = struct("C", cell(1,0), "nodes", cell(1,0));
%!   for j=[1,2,4,5,8,9]
%!     data(i).mesh.elements.joints(end + 1).nodes = j;
%!     data(i).mesh.elements.joints(end).C = [eye(3), zeros(3, 3)];
%!   endfor
%!   data(i).mesh.materials.tet10 = int32(1);
%!   E = 210000e6;
%!   nu = 0.3;
%!   data(i).mesh.material_data.rho = 7850;
%!   data(i).mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%!   data(i).load_case.pressure.tria6.elements = int32([1,2,3,5,6,7]);
%!   data(i).load_case.pressure.tria6.p = [p, p, 0, p, 0.5 * p, 0.5 * p];
%!   data(i).load_case.locked_dof = false(rows(data(i).mesh.nodes), 6);
%!   [data(i).dof_map] = fem_ass_dof_map(data(i).mesh, data(i).load_case);
%!   [data(i).mat_ass.K, ...
%!    data(i).mat_ass.R, ...
%!    data(i).mat_ass.Rlumped] = fem_ass_matrix(data(i).mesh, ...
%!                                              data(i).dof_map, ...
%!                                              [FEM_MAT_STIFFNESS, ...
%!                                               FEM_VEC_LOAD_CONSISTENT, ...
%!                                               FEM_VEC_LOAD_LUMPED], ...
%!                                              data(i).load_case);
%!   data(i).Flumped = fem_post_def_nodal(data(i).mesh, data(i).dof_map, data(i).mat_ass.Rlumped);
%!   data(i).Fcon = fem_post_def_nodal(data(i).mesh, data(i).dof_map, data(i).mat_ass.R);
%!   assert_simple(R1.' * sum(data(i).Fcon(:, 1:3), 1).', a * b * p / 3 * [0; 0; -1], eps * p);
%!   assert_simple(R1.' * sum(data(i).Flumped(:, 1:3), 1).', a * b * p / 3 * [0; 0; -1], eps * p);
%!   data(i).sol_stat = fem_sol_static(data(i).mesh, data(i).dof_map, data(i).mat_ass);
%!   data(i).sol_stat_lumped = fem_sol_static(data(i).mesh, data(i).dof_map, setfield(data(i).mat_ass, "R", data(i).mat_ass.Rlumped));
%!   if (do_plot)
%!     figure("visible", "off");
%!     fem_post_sol_plot(data(i).mesh, data(i).sol_stat, scale / max(norm(data(i).sol_stat.def, "rows")), 1);
%!     xlabel("x [m]");
%!     ylabel("y [m]");
%!     zlabel("z [m]");
%!     grid on;
%!     grid minor on;
%!     title(sprintf("consistent load case %d", i));
%!     figure("visible", "off");
%!     fem_post_sol_plot(data(i).mesh, data(i).sol_stat_lumped, scale / max(norm(data(i).sol_stat_lumped.def, "rows")), 1);
%!     xlabel("x [m]");
%!     ylabel("y [m]");
%!     zlabel("z [m]");
%!     grid on;
%!     grid minor on;
%!     title(sprintf("lumped load case %d", i));
%!   endif
%! endfor
%! if (do_plot)
%!   figure_list();
%! endif
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
