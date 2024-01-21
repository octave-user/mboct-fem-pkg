## fem_tests.m:120
%!test
%! ## DEMO 1
%! ## Cantilever beam with rectangular cross section and lateral load
%! ## W.Beitz, K.-H.Grothe, 1997, Dubbel, section 2.4.6, page C17, figure 23
%! close all;
%! do_post_pro = false;
%! num_iso = 10;
%! material.E = 210000e6;
%! material.nu = 0.3;
%! material.rho = 7850;
%! Fz = -15000;
%! h = 2e-3;
%! geometry.l = 150e-3;
%! geometry.w = h;
%! geometry.h = 60e-3;
%! mesh_size.num_elem_l = ceil(geometry.l / h);
%! mesh_size.num_elem_w = ceil(geometry.w / h);
%! mesh_size.num_elem_h = ceil(geometry.h / h);
%! f = [ 0; 0; Fz ];
%! [mesh, load_case] = fem_pre_mesh_cube_create(geometry, mesh_size, material, f);
%! mesh.nodes(:, 3) -= 0.5 * geometry.h;
%! [dof_map] = fem_ass_dof_map(mesh, load_case);
%! [mat_ass.K, ...
%!  mat_ass.R] = fem_ass_matrix(mesh, ...
%!                              dof_map, ...
%!                              [FEM_MAT_STIFFNESS, ...
%!                               FEM_VEC_LOAD_CONSISTENT], ...
%!                              load_case);
%! [sol_stat] = fem_sol_static(mesh, dof_map, mat_ass);
%! sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                  dof_map, ...
%!                                  [FEM_VEC_STRESS_CAUCH], ...
%!                                  load_case, ...
%!                                  sol_stat);
%! x = linspace(0, geometry.l, 100);
%! z = linspace(-0.5 * geometry.h, 0.5 * geometry.h, 50);
%! [xx, zz] = meshgrid(x, z);
%! xtauel = mesh.nodes(:, 1)(mesh.elements.iso8);
%! ytauel = mesh.nodes(:, 1)(mesh.elements.iso8);
%! ztauel = mesh.nodes(:, 3)(mesh.elements.iso8);
%! tauxxel = sol_stat.stress.taum.iso8(:, :, 1);
%! tauxzel = sol_stat.stress.taum.iso8(:, :, 6);
%! tauxx = griddata(xtauel(:), ztauel(:), tauxxel(:), xx, zz);
%! tauxz = griddata(xtauel(:), ztauel(:), tauxzel(:), xx, zz);
%! Iy = geometry.w * geometry.h^3 / 12;
%! tauxx_a = -Fz / Iy * (geometry.l - xx) .* zz;
%! tauxz_a = 3 / 2 * Fz / (geometry.h * geometry.w) * (1 - (zz / (0.5 * geometry.h)).^2);
%! scale_tauxx = linspace(min(min(tauxx_a)), max(max(tauxx_a)), num_iso + 1);
%! scale_tauxz = linspace(min(min(tauxz_a)), max(max(tauxz_a)), num_iso + 1);
%! figure("visible", "off");
%! subplot(2, 1, 1);
%! contourf(xx, zz, tauxx_a, scale_tauxx);
%! daspect([1,1,1]);
%! colormap jet;
%! colorbar;
%! xlabel("x [m]");
%! ylabel("z [m]");
%! title("stress component tauxx [Pa]");
%! grid on;
%! grid minor on;
%! subplot(2, 1, 2);
%! contourf(xx, zz, tauxx, scale_tauxx);
%! daspect([1,1,1]);
%! colormap jet;
%! colorbar;
%! xlabel("x [m]");
%! ylabel("z [m]");
%! title("stress component tauxx [Pa]");
%! grid on;
%! grid minor on;
%! figure("visible", "off");
%! subplot(2, 1, 1);
%! contourf(xx, zz, tauxz_a, scale_tauxz);
%! daspect([1,1,1]);
%! colormap jet;
%! colorbar;
%! xlabel("x [m]");
%! ylabel("z [m]");
%! title("stress component tauxz [Pa]");
%! grid on;
%! grid minor on;
%! subplot(2, 1, 2);
%! contourf(xx, zz, tauxz, scale_tauxz);
%! daspect([1,1,1]);
%! colormap jet;
%! colorbar;
%! xlabel("x [m]");
%! ylabel("z [m]");
%! title("stress component tauxz [Pa]");
%! grid on;
%! grid minor on;
%! if (do_post_pro)
%!   opts.scale_def = 0.3 * geometry.l / max(max(abs(sol_stat.def)));
%!   opts.print_and_exit = true;
%!   opts.print_to_file = "";
%!   opts.skin_only = true;
%!   opts.show_element = true;
%!   opts.print_to_file = "";
%!   unwind_protect
%!     opts.print_to_file = tempname();
%!     opts.rotation_angle = [pi/2, 0, 0];
%!     fem_post_sol_external(mesh, sol_stat, opts);
%!     [img, map, alpha] = imread([opts.print_to_file, "_001.jpg"]);
%!     figure("visible", "off");
%!     imshow(img, map);
%!     title("Gmsh - deformed mesh / continuous stress tensor");
%!   unwind_protect_cleanup
%!     if (numel(opts.print_to_file))
%!       [~] = unlink([opts.print_to_file, "_001.jpg"]);
%!     endif
%!   end_unwind_protect
%! endif
%! figure_list();
%! idx_x = find((xx(:) > 0.4 * geometry.l) & (xx(:) < 0.6 * geometry.l));
%! assert_simple(tauxx(:)(idx_x), tauxx_a(:)(idx_x), 1e-2 * max(tauxx_a(:)(idx_x)));
%! assert_simple(tauxz(:)(idx_x), tauxz_a(:)(idx_x), 6e-2 * max(abs(tauxz_a(:)(idx_x))));
