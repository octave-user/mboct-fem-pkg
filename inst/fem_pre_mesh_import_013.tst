## fem_pre_mesh_import.m:13
%!test
%! try
%! ### TEST13
%! do_plot = false;
%! if (do_plot)
%!   close all;
%! endif
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   unwind_protect
%!     [fd, msg] = fopen([filename, ".geo"], "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s.geo\"", filename);
%!     endif
%!     a = 60e-3;
%!     c = 10e-3;
%!     d = 0e-3;
%!     h = 2e-3;
%!     b = h;
%!     Fx = 0;
%!     Fz = -1000;
%!     My = 0;
%!     num_iso = 10;
%!     do_post_pro = false;
%!     fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "a=%g;\n", a);
%!     fprintf(fd, "b=%g;\n", b);
%!     fprintf(fd, "c=%g;\n", c);
%!     fprintf(fd, "h = %g;\n", h);
%!     fputs(fd, "Point(1) = {0.0, -0.5 * b, -0.5 * c, h};\n");
%!     fputs(fd, "Point(2) = {  a, -0.5 * b, -0.5 * c, h};\n");
%!     fputs(fd, "Point(3) = {  a,  0.5 * b, -0.5 * c, h};\n");
%!     fputs(fd, "Point(4) = {0.0,  0.5 * b, -0.5 * c, h};\n");
%!     fputs(fd, "Line(1) = {4,3};\n");
%!     fputs(fd, "Line(2) = {3,2};\n");
%!     fputs(fd, "Line(3) = {2,1};\n");
%!     fputs(fd, "Line(4) = {1,4};\n");
%!     fputs(fd, "Line Loop(5) = {2,3,4,1};\n");
%!     fputs(fd, "Plane Surface(6) = {5};\n");
%!     fputs(fd, "tmp[] = Extrude {0,0.0,c} {\n");
%!     fputs(fd, "  Surface{6};\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%!     fputs(fd, "Physical Surface(\"clamp\",1) = {tmp[4]};\n");
%!     fputs(fd, "Physical Surface(\"load\",2) = {tmp[2]};\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   [~] = unlink([filename, ".geo"]);
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!   master_node_idx = int32(rows(mesh.nodes) + 1);
%!   mesh.nodes(master_node_idx, 1:3) = [a + d, 0, 0];
%!   mesh.elements.rbe3 = fem_pre_mesh_rbe3_from_surf(mesh, 2, master_node_idx);
%!   [~] = unlink([filename, ".msh"]);
%!   load_case.locked_dof = false(rows(mesh.nodes), 6);
%!   load_case.locked_dof(mesh.groups.tria6(find([[mesh.groups.tria6].id] == 1)).nodes, :) = true;
%!   load_case.loaded_nodes = [master_node_idx];
%!   load_case.loads = [Fx, 0, Fz, 0, My, 0];
%!   mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%!   E = 210000e6;
%!   nu = 0.3;
%!   mesh.material_data.rho = 7850;
%!   mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.K, ...
%!    mat_ass.R, ...
%!    mat_ass.dm] = fem_ass_matrix(mesh, ...
%!                                 dof_map, ...
%!                                 [FEM_MAT_STIFFNESS, ...
%!                                  FEM_VEC_LOAD_CONSISTENT, ...
%!                                  FEM_SCA_TOT_MASS], ...
%!                                 load_case);
%!   sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%!   sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_VEC_STRESS_CAUCH], ...
%!                                    load_case, ...
%!                                    sol_stat);
%!   x = linspace(0, a, 100);
%!   z = linspace(-0.5 * c, 0.5 * c, 50);
%!   [xx, zz] = meshgrid(x, z);
%!   xtauel = mesh.nodes(:, 1)(mesh.elements.tet10);
%!   ytauel = mesh.nodes(:, 1)(mesh.elements.tet10);
%!   ztauel = mesh.nodes(:, 3)(mesh.elements.tet10);
%!   tauxxel = sol_stat.stress.taum.tet10(:, :, 1);
%!   tauxzel = sol_stat.stress.taum.tet10(:, :, 6);
%!   tauxx = griddata(xtauel(:), ztauel(:), tauxxel(:), xx, zz);
%!   tauxz = griddata(xtauel(:), ztauel(:), tauxzel(:), xx, zz);
%!   Iy = b * c^3 / 12;
%!   tauxx_a = -Fz / Iy * (a - xx) .* zz + Fx / (b * c) + My / Iy * zz;
%!   tauxz_a = 3 / 2 * Fz / (b * c) * (1 - (zz / (0.5 * c)).^2);
%!   scale_tauxx = linspace(min(min(tauxx_a)), max(max(tauxx_a)), num_iso + 1);
%!   if (do_plot)
%!     figure("visible", "off");
%!     subplot(2, 1, 1);
%!     contourf(xx, zz, tauxx_a, scale_tauxx);
%!     daspect([1,1,1]);
%!     colormap jet;
%!     colorbar;
%!     xlabel("x [m]");
%!     ylabel("z [m]");
%!     title("stress component tauxx [Pa]");
%!     grid on;
%!     grid minor on;
%!     subplot(2, 1, 2);
%!     contourf(xx, zz, tauxx, scale_tauxx);
%!     daspect([1,1,1]);
%!     colormap jet;
%!     colorbar;
%!     xlabel("x [m]");
%!     ylabel("z [m]");
%!     title("stress component tauxx [Pa]");
%!     grid on;
%!     grid minor on;
%!     scale_tauxz = linspace(min(min(tauxz_a)), max(max(tauxz_a)), num_iso + 1);
%!     figure("visible", "off");
%!     subplot(2, 1, 1);
%!     contourf(xx, zz, tauxz_a, scale_tauxz);
%!     daspect([1,1,1]);
%!     colormap jet;
%!     colorbar;
%!     xlabel("x [m]");
%!     ylabel("z [m]");
%!     title("stress component tauxz [Pa]");
%!     grid on;
%!     grid minor on;
%!     subplot(2, 1, 2);
%!     contourf(xx, zz, tauxz, scale_tauxz);
%!     daspect([1,1,1]);
%!     colormap jet;
%!     colorbar;
%!     xlabel("x [m]");
%!     ylabel("z [m]");
%!     title("stress component tauxz [Pa]");
%!     grid on;
%!     grid minor on;
%!     figure_list();
%!   endif
%!   if (do_post_pro)
%!     opts.scale_def = 0.3 * a / max(max(abs(sol_stat.def)));
%!     fem_post_sol_external(mesh, sol_stat, opts);
%!   endif
%!   idx_x = find((xx(:) > 0.1 * a) & (xx(:) < 0.9 * a));
%!   assert_simple(tauxx(:)(idx_x), tauxx_a(:)(idx_x), 0.3e-2 * max(max(max(abs(sol_stat.stress.taum.tet10)))));
%!   assert_simple(tauxz(:)(idx_x), tauxz_a(:)(idx_x), 0.3e-2 * max(max(max(abs(sol_stat.stress.taum.tet10)))));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       [~] = unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
