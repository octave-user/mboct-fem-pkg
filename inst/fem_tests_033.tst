## fem_tests.m:33
%!test
%! try
%! ##########################################################################################
%! ## TEST 33: Test case for sfncon4 using joints
%! ##########################################################################################
%! close all;
%! Fx = 1250;
%! c = 1e-6;
%! scale = 20e-3;
%! num_modes = 6;
%! material1.E = 210000e6;
%! material1.nu = 0.3;
%! material1.rho = 7850;
%! contour_plot = do_plot = false;
%! h1 = 5e-3;
%! geometry1.l = 20e-3;
%! geometry1.w = 30e-3;
%! geometry1.h = 50e-3;
%! mesh_size1.num_elem_l = ceil(geometry1.l / h1);
%! mesh_size1.num_elem_w = ceil(geometry1.w / h1);
%! mesh_size1.num_elem_h = ceil(geometry1.h / h1);
%! [data(1).mesh] = fem_pre_mesh_cube_create(geometry1, mesh_size1, material1, zeros(3, 1));
%! data(1).mesh.nodes(:, 2) -= 0.5 * geometry1.w;
%! data(1).mesh.nodes(:, 3) -= 0.5 * geometry1.h;
%! material2.E = material1.E;
%! material2.nu = material1.nu;
%! material2.rho = material1.rho;
%! h2 = 2.5e-3;
%! geometry2.l = 20e-3;
%! geometry2.w = 30e-3;
%! geometry2.h = 50e-3;
%! mesh_size2.num_elem_l = ceil(geometry2.l / h2);
%! mesh_size2.num_elem_w = ceil(geometry2.w / h2);
%! mesh_size2.num_elem_h = ceil(geometry2.h / h2);
%! [data(2).mesh] = fem_pre_mesh_cube_create(geometry2, mesh_size2, material2, zeros(3, 1));
%! data(2).mesh.nodes(:, 1) += geometry1.l + c;
%! data(2).mesh.nodes(:, 2) -= 0.5 * geometry2.w;
%! data(2).mesh.nodes(:, 3) -= 0.5 * geometry2.h;
%! for i=1:numel(data)
%!   data(i).load_case.locked_dof = false(size(data(i).mesh.nodes));
%!   data(i).dof_map = fem_ass_dof_map(data(i).mesh, data(i).load_case);
%! endfor
%! [data(3).mesh, data(3).dof_map] = fem_post_mesh_merge(data);
%! idx_clamp = find(data(3).mesh.nodes(:, 1) == 0);
%! idx_force = find(data(3).mesh.nodes(:, 1) == geometry1.l + geometry2.l + c);
%! idx_master = find(data(3).mesh.nodes(:, 1) == geometry1.l);
%! idx_slave = find(data(3).mesh.nodes(:, 1) == geometry1.l + c);
%! data(3).load_case.locked_dof = false(size(data(3).mesh.nodes));
%! data(3).load_case.locked_dof(idx_clamp, 1:3) = true;
%! data(3).load_case.loaded_nodes = int32(idx_force);
%! data(3).load_case.loads = [repmat(Fx / numel(idx_force), numel(idx_force), 1), zeros(numel(idx_force), 5)];
%! elements.sfncon4.slave = int32(idx_slave);
%! elements.sfncon4.master = zeros(0, 4, "int32");
%! for i=1:numel(idx_master)
%!   [ielem, inode] = find(data(3).mesh.elements.iso8 == idx_master(i));
%!   elements.sfncon4.master(end + (1:numel(ielem)), :) = data(3).mesh.elements.iso8(ielem, [4, 1, 5, 8]);
%! endfor
%! elements.sfncon4.maxdist = c * (1 + sqrt(eps));
%! data(3).mesh.elements.joints = fem_pre_mesh_constr_surf_to_node(data(3).mesh.nodes, elements);
%! data(3).dof_map = fem_ass_dof_map(data(3).mesh, data(3).load_case);
%! [data(3).mat_ass.K, ...
%!  data(3).mat_ass.M, ...
%!  data(3).mat_ass.R, ...
%!  data(3).mat_ass.dm, ...
%!  data(3).mat_ass.S, ...
%!  data(3).mat_ass.J, ...
%!  data(3).mat_ass.mat_info, ...
%!  data(3).mat_ass.mesh_info] = fem_ass_matrix(data(3).mesh, ...
%!                                              data(3).dof_map, ...
%!                                              [FEM_MAT_STIFFNESS, ...
%!                                               FEM_MAT_MASS, ...
%!                                               FEM_VEC_LOAD_CONSISTENT, ...
%!                                               FEM_SCA_TOT_MASS, ...
%!                                               FEM_VEC_INERTIA_M1, ...
%!                                               FEM_MAT_INERTIA_J], ...
%!                                              data(3).load_case);
%! data(3).sol_stat = fem_sol_static(data(3).mesh, data(3).dof_map, data(3).mat_ass);
%! data(3).sol_eig = fem_sol_modal(data(3).mesh, data(3).dof_map, data(3).mat_ass, num_modes);
%! data(3).sol_eig.stress = fem_ass_matrix(data(3).mesh, ...
%!                                         data(3).dof_map, ...
%!                                         [FEM_SCA_STRESS_VMIS], ...
%!                                         data(3).load_case, ...
%!                                         data(3).sol_eig);
%! if (do_plot)
%!   figure("visible", "off");
%!   fem_post_sol_plot(data(3).mesh, data(3).sol_stat, scale / max(norm(data(3).sol_stat.def(:, 1:3), "rows")), 1);
%!   view(30,30);
%!   xlabel("x [m]");
%!   ylabel("y [m]");
%!   zlabel("z [m]");
%!   grid on;
%!   grid minor on;
%!   title("static deflection");
%!   for i=1:numel(data(3).sol_eig.f)
%!     figure("visible", "off");
%!     fem_post_sol_plot(data(3).mesh, data(3).sol_eig, scale / max(norm(data(3).sol_eig.def(:, 1:3, i), "rows")), i);
%!     view(30,30);
%!     xlabel("x [m]");
%!     ylabel("y [m]");
%!     zlabel("z [m]");
%!     grid on;
%!     grid minor on;
%!     title(sprintf("%d: patched model eigenmode %.1fHz", i, data(3).sol_eig.f(i)));
%!   endfor
%! endif
%! material.E = 210000e6;
%! material.nu = 0.3;
%! material.rho = 7850;
%! h = min([h1,h2]);
%! geometry.l = geometry1.l + geometry2.l;
%! geometry.w = mean([geometry1.w, geometry2.w]);
%! geometry.h = mean([geometry1.h, geometry2.h]);
%! mesh_size.num_elem_l = ceil(geometry.l / h);
%! mesh_size.num_elem_w = ceil(geometry.w / h);
%! mesh_size.num_elem_h = ceil(geometry.h / h);
%! [mesh, load_case] = fem_pre_mesh_cube_create(geometry, mesh_size, material, zeros(3,1));
%! mesh.nodes(:, 2) -= 0.5 * geometry.w;
%! mesh.nodes(:, 3) -= 0.5 * geometry.h;
%! [dof_map] = fem_ass_dof_map(mesh, load_case);
%! [mat_ass.M, ...
%!  mat_ass.K] = fem_ass_matrix(mesh, ...
%!                              dof_map, ...
%!                              [FEM_MAT_MASS, ...
%!                               FEM_MAT_STIFFNESS], ...
%!                              load_case);
%! [sol_eig] = fem_sol_modal(mesh, dof_map, mat_ass, num_modes);
%! sol_eig.stress = fem_ass_matrix(mesh, ...
%!                                 dof_map, ...
%!                                 [FEM_SCA_STRESS_VMIS], ...
%!                                 load_case, ...
%!                                 sol_eig);
%! if (do_plot)
%!   for i=1:length(sol_eig.f)
%!     figure("visible", "off");
%!     hold on;
%!     fem_post_sol_plot(mesh, sol_eig, scale / max(norm(sol_eig.def(:, 1:3, i), "rows")), i);
%!     view(30,30);
%!     xlabel('x [m]');
%!     ylabel('y [m]');
%!     zlabel('z [m]');
%!     grid on;
%!     grid minor on;
%!     title(sprintf("%d. uniform model: eigenmode: %gHz", i, sol_eig.f(i)));
%!   endfor
%! endif
%! if (contour_plot)
%!   for i=1:size(sol_eig.stress.vmis.iso8,3)
%!     figure("visible", "off");
%!     subplot(2, 1, 1);
%!     xi = mesh.nodes(mesh.elements.iso8,1);
%!     yi = mesh.nodes(mesh.elements.iso8,2);
%!     zi = mesh.nodes(mesh.elements.iso8,3);
%!     taui = sol_eig.stress.vmis.iso8(:,:,i)(:);
%!     x = linspace (min(xi), max(xi), 100);
%!     y = linspace(min(yi),max(yi),100);
%!     z = max(zi);
%!     tau = griddata3(xi, yi, zi, taui, x, y, z);
%!     contourf(x, y, tau);
%!     colormap("jet");
%!     grid minor on;
%!     colorbar;
%!     title(sprintf("continuous model mode %d", i));
%!     xi = data(3).mesh.nodes(data(3).mesh.elements.iso8, 1);
%!     yi = data(3).mesh.nodes(data(3).mesh.elements.iso8, 2);
%!     zi = data(3).mesh.nodes(data(3).mesh.elements.iso8, 3);
%!     x = linspace(min(xi), max(xi), 100);
%!     y = linspace(min(yi), max(yi), 100);
%!     z = max(zi);
%!     taui = data(3).sol_eig.stress.vmis.iso8(:,:,i)(:);
%!     tau = griddata3(xi,yi,zi,taui, x, y, z);
%!     subplot(2, 1, 2);
%!     contourf(x, y, tau);
%!     colormap("jet");
%!     grid minor on;
%!     colorbar;
%!     title(sprintf("patched model mode %d", i));
%!   endfor
%! endif
%! if (do_plot || contour_plot)
%!   figure_list();
%! endif
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
