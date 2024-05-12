## fem_pre_mesh_import.m:11
%!test
%! try
%! ## TEST11
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
%!   ri = [8e-3; 5e-3];
%!   ro = [10e-3; ri(1)];
%!   h = 12e-3;
%!   scale_def = 5e-3;
%!   toldist = 1e-3;
%!   tolf = 8e-2;
%!   mesh_size = linspace(0.6e-3, 0.3e-3, 5)(1);
%!   num_nodes = zeros(1, numel(mesh_size));
%!   assert_simple(numel(ri), numel(ro));
%!   orange = [2,1];
%!   rrange = [false,true];
%!   f = nan(numel(mesh_size), num_modes, numel(rrange), numel(orange));
%!   for o=1:numel(orange)
%!     for r=1:numel(rrange)
%!       if (orange(o) == 2 && rrange(r))
%!         continue;
%!       endif
%!       for m=1:numel(mesh_size)
%!         clear data;
%!         clear mesh;
%!         clear mat_ass;
%!         clear sol_eig;
%!         for i=1:numel(ri)
%!           fd = -1;
%!           unwind_protect
%!             [fd, msg] = fopen([filename, ".geo"], "w");
%!             if fd == -1
%!               error("failed to open file \"%s.geo\"", filename);
%!             endif
%!             fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!             fprintf(fd, "ri = %g;\n", ri(i));
%!             fprintf(fd, "ro = %g;\n", ro(i));
%!             fprintf(fd, "h = %g;\n", h);
%!             fputs(fd, "Point(1) = {ri,0.0,0.0};\n");
%!             fputs(fd, "Point(2) = {ro,0.0,0.0};\n");
%!             fputs(fd, "Point(5) = {ro,0.0,h};\n");
%!             fputs(fd, "Point(6) = {ri,0.0,h};\n");
%!             fputs(fd, "Line(1) = {1,2};\n");
%!             fputs(fd, "Line(4) = {2,5};\n");
%!             fputs(fd, "Line(5) = {5,6};\n");
%!             fputs(fd, "Line(8) = {6,1};\n");
%!             fputs(fd, "Line Loop(5) = {1,4,5,8};\n");
%!             fputs(fd, "Plane Surface(6) = {5};\n");
%!             if (orange(o) == 1 || i == 1)
%!               if (rrange(r) || i == 1)
%!              fprintf(fd, "tmp[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} { Surface{6}; Layers{%d}; Recombine; };\n", ceil(0.5 * pi * ri(1) / mesh_size(m)));
%!              fprintf(fd, "Recombine Surface{6,tmp[0]};\n");
%!               else
%!              fprintf(fd, "tmp[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} { Surface{6}; Layers{%d};};\n", ceil(0.5 * pi * ri(1) / mesh_size(m)));
%!               endif
%!             else
%!               fputs(fd, "tmp[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} { Surface{6};};\n");
%!             endif
%!             fprintf(fd, "Physical Volume(\"volume\",%d) = {tmp[1]};\n", (i - 1) * 100 + 1);
%!             fprintf(fd, "Physical Surface(\"bottom\",%d) = {tmp[2]};\n", (i - 1) * 100 + 1);
%!             fprintf(fd, "Physical Surface(\"outside\",%d) = {tmp[3]};\n", (i - 1) * 100 + 2);
%!             fprintf(fd, "Physical Surface(\"inside\",%d) = {tmp[5]};\n", (i - 1) * 100 + 3);
%!             fprintf(fd, "Physical Surface(\"top\",%d) = {tmp[4]};\n", (i - 1) * 100 + 4);
%!             fprintf(fd, "Physical Surface(\"left\",%d) = {tmp[0]};\n", (i - 1) * 100 + 5); ##x
%!             fprintf(fd, "Physical Surface(\"right\",%d) = {6};\n", (i - 1) * 100 + 6);  ##y
%!           unwind_protect_cleanup
%!             if (fd ~= -1)
%!               fclose(fd);
%!             endif
%!           end_unwind_protect
%!           [~] = unlink([filename, ".msh"]);
%!           if (orange(o) == 1 || i == 1)
%!             optargs = {"-order", "1"};
%!           else
%!             optargs = {"-order", "2"};
%!           endif
%!           pid = spawn("gmsh", {"-format", "msh2", "-3", optargs{:}, "-clmin", sprintf("%g", 0.75 * mesh_size(m) * orange(o)), "-clmax", sprintf("%g", 1.25 *mesh_size(m) * orange(o)), [filename, ".geo"]});
%!           status = spawn_wait(pid);
%!           if status ~= 0
%!             warning("gmsh failed with status %d", status);
%!           endif
%!           [~] = unlink([filename, ".geo"]);
%!           data(i).mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!           [~] = unlink([filename, ".msh"]);
%!         endfor
%!         for i=1:numel(data)
%!           if (isfield(data(i).mesh.elements, "iso8"))
%!             data(i).mesh.materials.iso8 = ones(rows(data(i).mesh.elements.iso8), 1, "int32");
%!           else
%!             data(i).mesh.materials.tet10 = ones(rows(data(i).mesh.elements.tet10), 1, "int32");
%!           endif
%!           data(i).mesh.material_data.rho = rho(i);
%!           data(i).mesh.material_data.C = fem_pre_mat_isotropic(E(i), nu(i));
%!         endfor
%!         mesh = fem_post_mesh_merge(data);
%!         mesh.elements.sfncon4.master = mesh.elements.iso4(mesh.groups.iso4(find([mesh.groups.iso4.id]==3)).elements, :);
%!         if (orange(o) == 1)
%!           mesh.elements.sfncon4.slave = mesh.groups.iso4(find([mesh.groups.iso4.id]==102)).nodes(:);
%!         else
%!           mesh.elements.sfncon4.slave = mesh.groups.tria6(find([mesh.groups.tria6.id]==102)).nodes(:);
%!         endif
%!         mesh.elements.sfncon4.maxdist = toldist * ri(1);
%!         mesh.elements.sfncon4.constraint = FEM_CT_SLIDING;
%!         if (orange(o) == 1)
%!           group_id = [[mesh.groups.iso4].id];
%!           node_constr1 = [mesh.groups.iso4(find((group_id == 1)|(group_id == 101))).nodes];
%!           node_constr5 = [mesh.groups.iso4(find(((group_id == 5)|(group_id == 105)))).nodes];
%!           node_constr6 = [mesh.groups.iso4(find(((group_id == 6)|(group_id == 106)))).nodes];
%!         else
%!           group_id4 = [[mesh.groups.iso4].id];
%!           group_id6 = [[mesh.groups.tria6].id];
%!           node_constr1 = [mesh.groups.iso4(find((group_id4 == 1))).nodes,   mesh.groups.tria6(find((group_id6 == 101))).nodes];
%!           node_constr5 = [mesh.groups.iso4(find(((group_id4 == 5)))).nodes, mesh.groups.tria6(find(((group_id6 == 105)))).nodes];
%!           node_constr6 = [mesh.groups.iso4(find(((group_id4 == 6)))).nodes, mesh.groups.tria6(find(((group_id6 == 106)))).nodes];
%!         endif
%!         node_constr = [node_constr1, node_constr5, node_constr6];
%!         mesh.elements.joints = repmat(struct("nodes",[],"C",[]), 1, numel(node_constr));
%!         idx_joint = int32(0);
%!         for i=1:numel(node_constr1)
%!           mesh.elements.joints(++idx_joint).C = [0,0,1, zeros(1, 3)];
%!           mesh.elements.joints(idx_joint).nodes = node_constr1(i);
%!         endfor
%!         for i=1:numel(node_constr5)
%!           mesh.elements.joints(++idx_joint).C = [1,0,0, zeros(1, 3)];
%!           mesh.elements.joints(idx_joint).nodes = node_constr5(i);
%!         endfor
%!         for i=1:numel(node_constr6)
%!           mesh.elements.joints(++idx_joint).C = [0,1,0, zeros(1, 3)];
%!           mesh.elements.joints(idx_joint).nodes = node_constr6(i);
%!         endfor
%!         load_case.locked_dof = false(rows(mesh.nodes), 6);
%!         dof_map = fem_ass_dof_map(mesh, load_case);
%!         num_nodes(m) = rows(mesh.nodes);
%!         [mat_ass.K, ...
%!          mat_ass.M, ...
%!          mat_ass.mtot] = fem_ass_matrix(mesh, ...
%!                                         dof_map, ...
%!                                         [FEM_MAT_STIFFNESS, ...
%!                                          FEM_MAT_MASS, ...
%!                                          FEM_SCA_TOT_MASS], ...
%!                                         load_case);
%!         sol_eig = fem_sol_modal(mesh, dof_map, mat_ass, num_modes, shift);
%!         fem_post_mesh_export([filename, ".msh"], mesh);
%!         f(m, :, r, o) = sol_eig.f;
%!         if (do_plot)
%!           figure("visible", "off");
%!           for i=1:numel(data)
%!             fem_post_sol_plot(data(i).mesh);
%!           endfor
%!           xlabel("x [m]");
%!           ylabel("y [m]");
%!           zlabel("z [m]");
%!           grid on;
%!           grid minor on;
%!           title("undeformed mesh");
%!           for i=1:numel(sol_eig.f)
%!             figure("visible", "off");
%!             hold on;
%!             fem_post_sol_plot(mesh, sol_eig, scale_eig/max(norm(sol_eig.def(:, 1:3, i), "rows")), i);
%!             view(30,30);
%!             xlabel('x [m]');
%!             ylabel('y [m]');
%!             zlabel('z [m]');
%!             grid on;
%!             grid minor on;
%!             title(sprintf("mode %d: %.1fHz", i, sol_eig.f(i)));
%!           endfor
%!         endif
%!       endfor
%!       if (do_plot)
%!      figure("visible", "off");
%!      hold on;
%!      for i=1:columns(f)
%!           plot(num_nodes.', f(:, i, r, o), sprintf("-;mode %d;%d", i, i));
%!      endfor
%!      xlabel("nodes [1]");
%!      ylabel("f [Hz]");
%!      grid on;
%!      grid minor on;
%!      title("natural frequencies versus mesh size");
%!       endif
%!     endfor
%!   endfor
%!   fref = [26023  59514  91469  1.0372e+05  1.1294e+05  1.154e+05];
%!   for o=1:numel(orange)
%!     for r=1:numel(rrange)
%!       for i=1:rows(f)
%!         if (~all(isnan(f(i, :, r, o))))
%!           assert_simple(all(abs(f(i, :, r, o) ./ fref - 1) < tolf));
%!         endif
%!       endfor
%!     endfor
%!   endfor
%!   figure_list();
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
