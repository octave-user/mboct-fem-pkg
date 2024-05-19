## fem_pre_mesh_import.m:266
%!test
%! try
%! ## TEST266
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
%!   num_modes = 6;
%!   shift = 0;
%!   scale_eig = 1.5e-3;
%!   E = [210000e6; 70000e6];
%!   nu = [0.3, 0.3];
%!   rho = [7850; 2700];
%!   ri = [4e-3; 2.5e-3];
%!   ro = [5e-3; ri(1)];
%!   h = 12e-3;
%!   scale_def = 5e-3;
%!   toldist = 1e-3;
%!   mesh_size = 2e-3;
%!   num_nodes = zeros(1, numel(mesh_size));
%!   assert_simple(numel(ri), numel(ro));
%!   f = nan(numel(mesh_size), num_modes);
%!   for m=1:numel(mesh_size)
%!     clear data;
%!     clear mesh;
%!     clear mat_ass;
%!     clear sol_eig;
%!     for i=1:numel(ri)
%!       fd = -1;
%!       unwind_protect
%!      [fd, msg] = fopen([filename, ".geo"], "w");
%!      if (fd == -1)
%!           error("failed to open file \"%s.geo\"", filename);
%!      endif
%!      fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!      fprintf(fd, "ri = %g;\n", ri(i));
%!      fprintf(fd, "ro = %g;\n", ro(i));
%!      fprintf(fd, "h = %g;\n", h);
%!      fputs(fd, "Point(1) = {ri,0.0,0.0};\n");
%!      fputs(fd, "Point(2) = {ro,0.0,0.0};\n");
%!      fputs(fd, "Point(5) = {ro,0.0,h};\n");
%!      fputs(fd, "Point(6) = {ri,0.0,h};\n");
%!      fputs(fd, "Line(1) = {1,2};\n");
%!      fputs(fd, "Line(4) = {2,5};\n");
%!      fputs(fd, "Line(5) = {5,6};\n");
%!      fputs(fd, "Line(8) = {6,1};\n");
%!      fputs(fd, "Line Loop(5) = {1,4,5,8};\n");
%!      fputs(fd, "Plane Surface(6) = {5};\n");
%!      fputs(fd, "tmp[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} { Surface{6}; };\n");
%!      fprintf(fd, "Physical Volume(\"volume\",%d) = {tmp[1]};\n", (i - 1) * 100 + 1);
%!      fprintf(fd, "Physical Surface(\"bottom\",%d) = {tmp[2]};\n", (i - 1) * 100 + 1);
%!      fprintf(fd, "Physical Surface(\"outside\",%d) = {tmp[3]};\n", (i - 1) * 100 + 2);
%!      fprintf(fd, "Physical Surface(\"inside\",%d) = {tmp[5]};\n", (i - 1) * 100 + 3);
%!      fprintf(fd, "Physical Surface(\"top\",%d) = {tmp[4]};\n", (i - 1) * 100 + 4);
%!      fprintf(fd, "Physical Surface(\"left\",%d) = {tmp[0]};\n", (i - 1) * 100 + 5); ##x
%!      fprintf(fd, "Physical Surface(\"right\",%d) = {6};\n", (i - 1) * 100 + 6);  ##y
%!      fprintf(fd, "Mesh.HighOrderOptimize=2;\n");
%!       unwind_protect_cleanup
%!         if (fd ~= -1)
%!           fclose(fd);
%!         endif
%!       end_unwind_protect
%!       [~] = unlink([filename, ".msh"]);
%!       pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "3", "-clmin", sprintf("%g", 0.75 * mesh_size(m)), "-clmax", sprintf("%g", 1.25 *mesh_size(m)), [filename, ".geo"]});
%!       status = spawn_wait(pid);
%!       if (status ~= 0)
%!         warning("gmsh failed with status %d", status);
%!       endif
%!       [~] = unlink([filename, ".geo"]);
%!       data(i).mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!       [~] = unlink([filename, ".msh"]);
%!     endfor
%!     for i=1:numel(data)
%!       data(i).mesh.materials.tet20 = ones(rows(data(i).mesh.elements.tet20), 1, "int32");
%!       data(i).mesh.material_data.rho = rho(i);
%!       data(i).mesh.material_data.C = fem_pre_mat_isotropic(E(i), nu(i));
%!     endfor
%!     mesh = fem_post_mesh_merge(data);
%!     mesh.elements.sfncon10.master = mesh.elements.tria10(mesh.groups.tria10(find([mesh.groups.tria10.id]==3)).elements, :);
%!     mesh.elements.sfncon10.slave = mesh.groups.tria10(find([mesh.groups.tria10.id]==102)).nodes(:);
%!     mesh.elements.sfncon10.maxdist = toldist * ri(1);
%!     mesh.elements.sfncon10.constraint = FEM_CT_SLIDING;
%!     group_id = [[mesh.groups.tria10].id];
%!     node_constr1 = [mesh.groups.tria10(find((group_id == 1)|(group_id == 101))).nodes];
%!     node_constr5 = [mesh.groups.tria10(find(((group_id == 5)|(group_id == 105)))).nodes];
%!     node_constr6 = [mesh.groups.tria10(find(((group_id == 6)|(group_id == 106)))).nodes];
%!     node_constr = [node_constr1, node_constr5, node_constr6];
%!     mesh.elements.joints = repmat(struct("nodes",[],"C",[]), 1, numel(node_constr));
%!     idx_joint = int32(0);
%!     for i=1:numel(node_constr1)
%!       mesh.elements.joints(++idx_joint).C = [0,0,1, zeros(1, 3)];
%!       mesh.elements.joints(idx_joint).nodes = node_constr1(i);
%!     endfor
%!     for i=1:numel(node_constr5)
%!       mesh.elements.joints(++idx_joint).C = [1,0,0, zeros(1, 3)];
%!       mesh.elements.joints(idx_joint).nodes = node_constr5(i);
%!     endfor
%!     for i=1:numel(node_constr6)
%!       mesh.elements.joints(++idx_joint).C = [0,1,0, zeros(1, 3)];
%!       mesh.elements.joints(idx_joint).nodes = node_constr6(i);
%!     endfor
%!     load_case.locked_dof = false(rows(mesh.nodes), 6);
%!     dof_map = fem_ass_dof_map(mesh, load_case);
%!     num_nodes(m) = rows(mesh.nodes);
%!     [mat_ass.K, ...
%!      mat_ass.M, ...
%!      mat_ass.mtot] = fem_ass_matrix(mesh, ...
%!                                     dof_map, ...
%!                                     [FEM_MAT_STIFFNESS, ...
%!                                      FEM_MAT_MASS, ...
%!                                      FEM_SCA_TOT_MASS], ...
%!                                     load_case);
%!     sol_eig = fem_sol_modal(mesh, dof_map, mat_ass, num_modes, shift);
%!     f(m, :) = sol_eig.f;
%!     if (do_plot)
%!       figure("visible", "off");
%!       for i=1:numel(data)
%!      fem_post_sol_plot(data(i).mesh);
%!       endfor
%!       xlabel("x [m]");
%!       ylabel("y [m]");
%!       zlabel("z [m]");
%!       grid on;
%!       grid minor on;
%!       title("undeformed mesh");
%!     endif
%!     fref = [52248
%!             65901
%!             105250
%!             113070
%!             132630
%!             201780];
%!     assert_simple(sol_eig.f(:), fref, 0.1e-2 * max(fref));
%!     if (do_plot)
%!      for i=1:numel(sol_eig.f)
%!      figure("visible", "off");
%!      hold on;
%!      fem_post_sol_plot(mesh, sol_eig, scale_eig/max(norm(sol_eig.def(:, 1:3, i), "rows")), i);
%!      view(30,30);
%!      xlabel('x [m]');
%!      ylabel('y [m]');
%!      zlabel('z [m]');
%!      grid on;
%!      grid minor on;
%!      title(sprintf("mode %d: %.1fHz", i, sol_eig.f(i)));
%!       endfor
%!     endif
%!   endfor
%!   if (do_plot)
%!     figure("visible", "off");
%!     hold on;
%!     for i=1:columns(f)
%!       plot(num_nodes.', f(:, i), sprintf("-;mode %d;%d", i, i));
%!     endfor
%!     xlabel("nodes [1]");
%!     ylabel("f [Hz]");
%!     grid on;
%!     grid minor on;
%!     title("natural frequencies versus mesh size");
%!     figure_list();
%!   endif
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
