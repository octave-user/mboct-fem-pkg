## fem_pre_mesh_import.m:14
%!test
%! ### TEST14
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
%!     N = 3;
%!     a = 10e-3;
%!     b = 5e-3;
%!     c = 5e-3;
%!     h = 5e-3;
%!     fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "a=%g;\n", a);
%!     fprintf(fd, "b=%g;\n", b);
%!     fprintf(fd, "c=%g;\n", c);
%!     fprintf(fd, "h = %g;\n", h);
%!     fputs(fd, "Point(1) = {0.0,0.0,0.0,h};\n");
%!     fputs(fd, "Point(2) = {a,0.0,0.0,h};\n");
%!     fputs(fd, "Point(3) = {a,b,0.0,h};\n");
%!     fputs(fd, "Point(4) = {0,b,0.0,h};\n");
%!     fputs(fd, "Line(1) = {4,3};\n");
%!     fputs(fd, "Line(2) = {3,2};\n");
%!     fputs(fd, "Line(3) = {2,1};\n");
%!     fputs(fd, "Line(4) = {1,4};\n");
%!     fputs(fd, "Line Loop(5) = {2,3,4,1};\n");
%!     fputs(fd, "Plane Surface(6) = {5};\n");
%!     fputs(fd, "tmp[] = Extrude {0,0.0,c} {\n");
%!     fputs(fd, "  Surface{6};\n");
%!     fprintf(fd, "  Layers{%d}; Recombine;\n", ceil(c / h));
%!     fputs(fd, "};\n");
%!     fprintf(fd, "Recombine Surface{6,tmp[0]};\n");
%!     fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%!     fputs(fd, "Physical Surface(\"clamp\",1) = {tmp[4]};\n");
%!     fputs(fd, "Physical Surface(\"load\",2) = {tmp[2]};\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "1", [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   [~] = unlink([filename, ".geo"]);
%!   mesh1 = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!   [~] = unlink([filename, ".msh"]);
%!   mesh1.materials.iso8 = ones(rows(mesh1.elements.iso8), 1, "int32");
%!   E = 210000e6;
%!   nu = 0.3;
%!   mesh1.material_data.rho = 7850;
%!   mesh1.material_data.C = fem_pre_mat_isotropic(E, nu);
%!   data.mesh = mesh1;
%!   data = repmat(data, 1, 3);
%!   for i=2:3
%!     data(i).mesh.nodes(:, 1) += a * (i - 1);
%!     for j=1:numel(data(i).mesh.groups.iso4)
%!       data(i).mesh.groups.iso4(j).id += 100 * (i - 1);
%!       data(i).mesh.groups.iso4(j).name = sprintf("%s[%d]", data(i).mesh.groups.iso4(j).name, data(i).mesh.groups.iso4(j).id);
%!     endfor
%!     for j=1:numel(data(i).mesh.groups.iso8)
%!       data(i).mesh.groups.iso8(j).id += 100 * (i - 1);
%!       data(i).mesh.groups.iso8(j).name = sprintf("%s[%d]", data(i).mesh.groups.iso8(j).name, data(i).mesh.groups.iso8(j).id);
%!     endfor
%!   endfor
%!   [mesh] = fem_post_mesh_merge(data);
%!   mesh.nodes(:, 2) -= 0.5 * b;
%!   mesh.nodes(:, 3) -= 0.5 * c;
%!   for i=1:N
%!     if (i > 1)
%!       unwind_protect
%!         fem_post_mesh_export([filename, "_in.msh"], data(i - 1).mesh);
%!         pid = spawn("gmsh", {"-refine", "-format", "msh2", "-o", [filename, "_out.msh"], [filename, "_in.msh"]});
%!         status = spawn_wait(pid);
%!         if (status ~= 0)
%!           error("gmsh failed with status %d", status);
%!         endif
%!         data(i).mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, "_out.msh"]));
%!         data(i).mesh.material_data = mesh.material_data;
%!         data(i).mesh.materials.iso8 = zeros(rows(data(i).mesh.elements.iso8), 1, "int32");
%!         for j=1:numel(data(i).mesh.groups.iso8)
%!           data(i).mesh.materials.iso8(data(i).mesh.groups.iso8(j).elements) = j;
%!         endfor
%!       unwind_protect_cleanup
%!         [~] = unlink([filename, "_in.msh"]);
%!         [~] = unlink([filename, "_out.msh"]);
%!       end_unwind_protect
%!     else
%!       data(i).mesh = mesh;
%!     endif
%!     data(i).h = h / i;
%!     node_idx_rbe3 = int32(rows(data(i).mesh.nodes) + 1);
%!     data(i).mesh.nodes(node_idx_rbe3, 1:3) = [3 * a, 0, 0];
%!     grp_idx = find([data(i).mesh.groups.iso4.id] == 202);
%!     data(i).mesh.elements.rbe3.nodes = [node_idx_rbe3, data(i).mesh.groups.iso4(grp_idx).nodes];
%!     data(i).mesh.elements.rbe3.weight = ones(numel(data(i).mesh.elements.rbe3.nodes) - 1, 1);
%!     for j=1:2
%!       grp_idx_slave = find([[data(i).mesh.groups.iso4].id] == (j - 1) * 100 + 2);
%!       grp_idx_master = data(i).mesh.groups.iso4(find([[data(i).mesh.groups.iso4].id] == j * 100 + 1)).elements;
%!       data(i).mesh.elements.sfncon4(j).slave = data(i).mesh.groups.iso4(grp_idx_slave).nodes(:);
%!       data(i).mesh.elements.sfncon4(j).master = data(i).mesh.elements.iso4(grp_idx_master, :);
%!       data(i).mesh.elements.sfncon4(j).maxdist = sqrt(eps) * max(abs([a,b,c]));
%!     endfor
%!     data(i).load_case.locked_dof = false(rows(data(i).mesh.nodes), 6);
%!     data(i).load_case.locked_dof(data(i).mesh.groups.iso4(find([[data(i).mesh.groups.iso4].id] == 1)).nodes, :) = true;
%!     data(i).load_case.loaded_nodes = node_idx_rbe3;
%!     data(i).load_case.loads = [0,-1000, 0, 0, 0, 0];
%!     data(i).dof_map = fem_ass_dof_map(data(i).mesh, data(i).load_case);
%!     [data(i).mat_ass.K, ...
%!      data(i).mat_ass.R, ...
%!      data(i).mat_ass.mtot] = fem_ass_matrix(data(i).mesh, ...
%!                                             data(i).dof_map, ...
%!                                             [FEM_MAT_STIFFNESS, ...
%!                                              FEM_VEC_LOAD_CONSISTENT, ...
%!                                              FEM_SCA_TOT_MASS], ...
%!                                             data(i).load_case);
%!     assert_simple(data(i).mat_ass.mtot, ...
%!            a * b * c * sum([data(i).mesh.material_data.rho]), ...
%!            sqrt(eps) * a * b * c * sum([data(i).mesh.material_data.rho]));
%!     [data(i).sol_stat, data(i).sol_stat.U] = fem_sol_static(data(i).mesh, data(i).dof_map, data(i).mat_ass);
%!     data(i).sol_stat.stress = fem_ass_matrix(data(i).mesh, ...
%!                                              data(i).dof_map, ...
%!                                              [FEM_VEC_STRESS_CAUCH], ...
%!                                              data(i).load_case, ...
%!                                              data(i).sol_stat);
%!     data(i).W = data(i).sol_stat.U.' * data(i).mat_ass.K * data(i).sol_stat.U;
%!   endfor
%!   if (do_plot)
%!     figure("visible", "off");
%!     loglog([data.h], [data.W], "-x;W(h) [J];1");
%!     xlabel("h [m]");
%!     ylabel("W [J]");
%!     grid on;
%!     grid minor on;
%!     title("strain energy");
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
