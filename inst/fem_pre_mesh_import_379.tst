## fem_pre_mesh_import.m:379
%!test
%! ## TEST 379
%! close all;
%! a = 5e-3;
%! b = 5e-3;
%! c = 10e-3;
%! h = 0.5e-3;
%! E = 210000e6;
%! nu = 0.3;
%! rho = 7850;
%! gamma = 10 * pi / 180;
%! tolstat = 4e-2;
%! scale = 5e-3;
%! maxdist = 1e-2 * max([a,b,c]);
%! eliminate = false;
%! animate = false;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   for i=1:4
%!     fd = -1;
%!     unwind_protect
%!       [fd, msg] = fopen([filename, ".geo"], "w");
%!       if (fd == -1)
%!         error("failed to open file \"%s.geo\"", filename);
%!       endif
%!       switch (i)
%!         case {1,2}
%!           order = 1;
%!         otherwise
%!           order = 2;
%!       endswitch
%!       fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!       fprintf(fd, "a=%g;\n", a);
%!       fprintf(fd, "b=%g;\n", b);
%!       fprintf(fd, "c=%g;\n", c);
%!       fprintf(fd, "h=%g;\n", h * order);
%!       fputs(fd, "Point(1) = {0.0,0.0,0.0,h};\n");
%!       fputs(fd, "Point(2) = {a,0.0,0.0,h};\n");
%!       fputs(fd, "Point(3) = {a,b,0.0,h};\n");
%!       fputs(fd, "Point(4) = {0,b,0.0,h};\n");
%!       fputs(fd, "Line(1) = {1,2};\n");
%!       fputs(fd, "Line(2) = {2,3};\n");
%!       fputs(fd, "Line(3) = {3,4};\n");
%!       fputs(fd, "Line(4) = {4,1};\n");
%!       fputs(fd, "Line Loop(5) = {1,2,3,4};\n");
%!       fputs(fd, "Plane Surface(6) = {5};\n");
%!       fputs(fd, "tmp[] = Extrude {0.0,0.0,c} {\n");
%!       switch (i)
%!         case 1
%!           fprintf(fd, "  Surface{6}; Layers{%d}; Recombine;\n", ceil(c/h));
%!         otherwise
%!           fputs(fd, "  Surface{6};\n");
%!       endswitch
%!       fputs(fd, "};\n");
%!       switch (i)
%!         case 1
%!           fprintf(fd, "Recombine Surface{6,tmp[0]};\n");
%!       endswitch
%!       fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%!       fprintf(fd, "Physical Surface(\"right\",%d) = {tmp[2]};\n", 100 * i + 1);
%!       fprintf(fd, "Physical Surface(\"rear\",%d) = {tmp[3]};\n", 100 * i + 2);
%!       fprintf(fd, "Physical Surface(\"left\",%d) = {tmp[4]};\n", 100 * i + 3);
%!       fprintf(fd, "Physical Surface(\"front\",%d) = {tmp[5]};\n", 100 * i + 4);
%!       fprintf(fd, "Physical Surface(\"top\",%d) = {tmp[0]};\n", 100 * i + 5);
%!       fprintf(fd, "Physical Surface(\"bottom\",%d) = {6};\n", 100 * i + 6);
%!     unwind_protect_cleanup
%!       if (fd ~= -1)
%!         fclose(fd);
%!       endif
%!     end_unwind_protect
%!     pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", sprintf("%d", order), [filename, ".geo"]});
%!     status = spawn_wait(pid);
%!     if (status ~= 0)
%!       warning("gmsh failed with status %d", status);
%!     endif
%!     unlink([filename, ".geo"]);
%!     data(i).mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!     unlink([filename, ".msh"]);
%!   endfor
%!   data(2).mesh.nodes(:, 1) += a;
%!   data(3).mesh.nodes(:, 1) += a;
%!   data(3).mesh.nodes(:, 2) += b;
%!   data(4).mesh.nodes(:, 2) += b;
%!   for i=1:numel(data)
%!     if (isfield(data(i).mesh.elements, "iso8"))
%!       data(i).mesh.materials.iso8 = ones(rows(data(i).mesh.elements.iso8), 1, "int32");
%!     else
%!       data(i).mesh.materials.tet10 = ones(rows(data(i).mesh.elements.tet10), 1, "int32");
%!     endif
%!     data(i).mesh.material_data.rho = rho;
%!     data(i).mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%!   endfor
%!   mesh = fem_post_mesh_merge(data);
%!   mesh.nodes(:, 3) -= 0.5 * c;
%!   constr = FEM_CT_FIXED;
%!   mesh.elements.sfncon4(1).master = mesh.elements.iso4(mesh.groups.iso4(find([[mesh.groups.iso4].id] == 103)).elements, :);
%!   mesh.elements.sfncon4(1).slave = mesh.groups.tria6(find([[mesh.groups.tria6].id] == 401)).nodes(:);
%!   mesh.elements.sfncon4(1).maxdist = maxdist;
%!   mesh.elements.sfncon4(1).constraint = constr;
%!   mesh.elements.sfncon4(2).master = mesh.elements.iso4(mesh.groups.iso4(find([[mesh.groups.iso4].id] == 102)).elements, :);
%!   mesh.elements.sfncon4(2).slave = mesh.groups.iso4(find([[mesh.groups.iso4].id] == 204)).nodes(:);
%!   mesh.elements.sfncon4(2).maxdist = maxdist;
%!   mesh.elements.sfncon4(2).constraint = constr;
%!   mesh.elements.sfncon6(1).master = mesh.elements.tria6(mesh.groups.tria6(find([[mesh.groups.tria6].id] == 402)).elements, :);
%!   mesh.elements.sfncon6(1).slave = mesh.groups.tria6(find([[mesh.groups.tria6].id] == 304)).nodes(:);
%!   mesh.elements.sfncon6(1).maxdist = maxdist;
%!   mesh.elements.sfncon6(1).constraint = constr;
%!   mesh.elements.sfncon6(2).master = mesh.elements.tria6(mesh.groups.tria6(find([[mesh.groups.tria6].id] == 301)).elements, :);
%!   mesh.elements.sfncon6(2).slave = mesh.groups.iso4(find([[mesh.groups.iso4].id] == 203)).nodes(:);
%!   mesh.elements.sfncon6(2).maxdist = maxdist;
%!   mesh.elements.sfncon6(2).constraint = constr;
%!   slave_nodes = [];
%!   for i=1:numel(mesh.elements.sfncon4)
%!     idx_slave = [];
%!     for j=1:numel(slave_nodes)
%!       idx_slave = [idx_slave; find(mesh.elements.sfncon4(i).slave == slave_nodes(j))];
%!     endfor
%!     if (numel(idx_slave))
%!       mesh.elements.sfncon4(i).slave(idx_slave) = 0;
%!       mesh.elements.sfncon4(i).slave = mesh.elements.sfncon4(i).slave(find(mesh.elements.sfncon4(i).slave));
%!     endif
%!     slave_nodes = [slave_nodes; mesh.elements.sfncon4(i).slave];
%!   endfor
%!   for i=1:numel(mesh.elements.sfncon6)
%!     idx_slave = [];
%!     for j=1:numel(slave_nodes)
%!       idx_slave = [idx_slave; find(mesh.elements.sfncon6(i).slave == slave_nodes(j))];
%!     endfor
%!     if (numel(idx_slave))
%!       mesh.elements.sfncon6(i).slave(idx_slave) = 0;
%!       mesh.elements.sfncon6(i).slave = mesh.elements.sfncon6(i).slave(find(mesh.elements.sfncon6(i).slave));
%!     endif
%!     slave_nodes = [slave_nodes; mesh.elements.sfncon6(i).slave];
%!   endfor
%!   group_id4 = [[mesh.groups.iso4].id];
%!   group_id6 = [[mesh.groups.tria6].id];
%!   node_front = [[mesh.groups.iso4(find(group_id4 == 104))].nodes, [mesh.groups.tria6(find(group_id6 == 404))].nodes];
%!   node_rear = [[mesh.groups.iso4(find(group_id4 == 202))].nodes, [mesh.groups.tria6(find(group_id6 == 302))].nodes];
%!   node_right = [mesh.groups.iso4(find((group_id4 == 101) | (group_id4 == 201))).nodes, ...
%!                 mesh.groups.tria6(find((group_id6==303)|(group_id6==403))).nodes];
%!   node_bottom = [[mesh.groups.iso4(find(mod(group_id4, 100) == 6))].nodes, ...
%!                  [mesh.groups.tria6(find(mod(group_id6,100) == 6))].nodes];
%!   node_top = [[mesh.groups.iso4(find(mod(group_id4, 100) == 5))].nodes, ...
%!               [mesh.groups.tria6(find(mod(group_id6,100) == 5))].nodes];
%!   node_constr = [node_front, node_rear, node_right, node_bottom, node_top];
%!   mesh.elements.joints = repmat(struct("nodes",[],"C",[]), 1, numel(node_constr));
%!   load_case.joints = repmat(struct("U",[]), 1, numel(node_constr));
%!   idx_joint = int32(0);
%!   for i=1:numel(node_front)
%!     if (~numel(find(slave_nodes == node_front(i))))
%!       mesh.elements.joints(++idx_joint).C = [[1, 0, 0; 0, 0, 1],  zeros(2, 3)];
%!       mesh.elements.joints(idx_joint).nodes = node_front(i);
%!       load_case.joints(idx_joint).U = [gamma * mesh.nodes(node_front(i), 3); 0];
%!     endif
%!   endfor
%!   for i=1:numel(node_rear)
%!     if (~numel(find(slave_nodes == node_rear(i))))
%!       mesh.elements.joints(++idx_joint).C = [[1, 0, 0; 0, 0, 1], zeros(2, 3)];
%!       mesh.elements.joints(idx_joint).nodes = node_rear(i);
%!       load_case.joints(idx_joint).U = [gamma * mesh.nodes(node_rear(i), 3); 0];
%!     endif
%!   endfor
%!   for i=1:numel(node_right)
%!     if (~numel(find(slave_nodes == node_right(i))))
%!       mesh.elements.joints(++idx_joint).C = [0,1,0, zeros(1, 3)];
%!       mesh.elements.joints(idx_joint).nodes = node_right(i);
%!       load_case.joints(idx_joint).U = 0;
%!     endif
%!   endfor
%!   for i=1:numel(node_bottom)
%!     xi = mesh.nodes(node_bottom(i), 1);
%!     if (~(numel(find(slave_nodes == node_bottom(i))) || xi == 0 || xi == 2 * a))
%!       mesh.elements.joints(++idx_joint).C = [0,0,1, zeros(1, 3)];
%!       mesh.elements.joints(idx_joint).nodes = node_bottom(i);
%!       load_case.joints(idx_joint).U = 0;
%!     endif
%!   endfor
%!   for i=1:numel(node_top)
%!     xi = mesh.nodes(node_top(i), 1);
%!     if (~(numel(find(slave_nodes == node_top(i))) || xi == 0 || xi == 2 * a))
%!       mesh.elements.joints(++idx_joint).C = [0,0,1, zeros(1, 3)];
%!       mesh.elements.joints(idx_joint).nodes = node_top(i);
%!       load_case.joints(idx_joint).U = 0;
%!     endif
%!   endfor
%!   mesh.elements.joints = mesh.elements.joints(1:idx_joint);
%!   load_case.joints = load_case.joints(1:idx_joint);
%!   load_case.locked_dof = false(size(mesh.nodes));
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   if (eliminate)
%!     K_mat_type = FEM_MAT_STIFFNESS;
%!   else
%!     K_mat_type = FEM_MAT_STIFFNESS;
%!   endif
%!   [mat_ass.K, ...
%!    mat_ass.R, ...
%!    mat_ass.mat_info] = fem_ass_matrix(mesh, ...
%!                                       dof_map, ...
%!                                       [K_mat_type, ...
%!                                        FEM_VEC_LOAD_CONSISTENT], ...
%!                                       load_case);
%!   if (eliminate)
%!     [Tred, Kred, Mred, Rred] = fem_cms_constr_elim(mesh, dof_map, mat_ass);
%!     opt_ls.refine_max_iter = int32(100);
%!     Kfact = fem_sol_factor(Kred, opt_ls);
%!     Ured = Kfact \ Rred;
%!     sol_stat.def = fem_post_def_nodal(mesh, dof_map, Tred * Ured);
%!   else
%!     sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%!   endif
%!   [sol_stat.stress] = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_VEC_STRESS_CAUCH], ...
%!                                      load_case, ...
%!                                      sol_stat);
%!   if (animate)
%!     opt_anim.scale_def = 1;
%!     opt_anim.animation_delay = 1;
%!     opt_anim.print_and_exit = true;
%!     opt_anim.print_to_file = filename;
%!     opt_anim.rotation_angle = [286, 2, 205] * pi / 180;
%!     opt_anim.skin_only = true;
%!     opt_anim.show_element = true;
%!     unwind_protect
%!       fem_post_sol_external(mesh, sol_stat, opt_anim);
%!       [img, map, alpha] = imread([opt_anim.print_to_file, "_001.jpg"]);
%!       figure("visible", "off");
%!       imshow(img, map);
%!       title("Gmsh - deformed mesh / continuous stress tensor");
%!     unwind_protect_cleanup
%!       unlink([opt_anim.print_to_file, "_001.jpg"]);
%!     end_unwind_protect
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
