## fem_pre_mesh_import.m:12
%!test
%! ## TEST12
%! do_plot = false;
%! if (do_plot)
%!   close all;
%! endif
%! a = 10e-3;
%! b = 20e-3;
%! c = 10e-3;
%! h = 2e-3;
%! E = 210000e6;
%! nu = 0.3;
%! rho = 7850;
%! px = -10e6;
%! py = -2e6;
%! pz = -3e6;
%! tolstat = 4e-2;
%! scale = 30e-3;
%! maxdist = 1e-2 * max([a,b,c]);
%! eliminate = false;
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
%!     [~] = unlink([filename, ".geo"]);
%!     data(i).mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!     [~] = unlink([filename, ".msh"]);
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
%!   node_bottom = [[mesh.groups.iso4(find(mod(group_id4, 100) == 6))].nodes,   [mesh.groups.tria6(find(mod(group_id6,100)==6))].nodes];
%!   node_front = [[mesh.groups.iso4(find(group_id4 == 104))].nodes, [mesh.groups.tria6(find(group_id6 == 404))].nodes];
%!   node_right = [mesh.groups.iso4(find((group_id4 == 101) | (group_id4 == 201))).nodes];
%!   node_constr = [node_bottom, node_front, node_right];
%!   mesh.elements.joints = repmat(struct("nodes",[],"C",[]), 1, numel(node_constr));
%!   idx_joint = int32(0);
%!   for i=1:numel(node_bottom)
%!     if (~numel(find(slave_nodes == node_bottom(i))))
%!       mesh.elements.joints(++idx_joint).C = [0,0,1, zeros(1, 3)];
%!       mesh.elements.joints(idx_joint).nodes = node_bottom(i);
%!     endif
%!   endfor
%!   for i=1:numel(node_front)
%!     if (~numel(find(slave_nodes == node_front(i))))
%!       mesh.elements.joints(++idx_joint).C = [1,0,0, zeros(1, 3)];
%!       mesh.elements.joints(idx_joint).nodes = node_front(i);
%!     endif
%!   endfor
%!   for i=1:numel(node_right)
%!     if (~numel(find(slave_nodes == node_right(i))))
%!       mesh.elements.joints(++idx_joint).C = [0,1,0, zeros(1, 3)];
%!       mesh.elements.joints(idx_joint).nodes = node_right(i);
%!     endif
%!   endfor
%!   mesh.elements.joints = mesh.elements.joints(1:idx_joint);
%!   iso4_top = [[mesh.groups.iso4(find(mod(group_id4, 100) == 5))].elements];
%!   tria6_top = [[mesh.groups.tria6(find(mod(group_id6, 100) == 5))].elements];
%!   iso4_rear = mesh.groups.iso4(find(group_id4 == 202)).elements;
%!   tria6_rear = mesh.groups.tria6(find(group_id6 == 302)).elements;
%!   tria6_left = [[mesh.groups.tria6(find((group_id6 == 303) | (group_id6 == 403)))].elements];
%!   load_case.locked_dof = false(size(mesh.nodes));
%!   load_case.pressure.tria6.elements = mesh.elements.tria6([tria6_left, tria6_rear, tria6_top], :);
%!   load_case.pressure.tria6.p = [repmat(py, numel(tria6_left), 6); repmat(px, numel(tria6_rear), 6); repmat(pz, numel(tria6_top), 6)];
%!   load_case.pressure.iso4.elements = mesh.elements.iso4([iso4_rear, iso4_top], :);
%!   load_case.pressure.iso4.p = [repmat(px, numel(iso4_rear), 4); repmat(pz, numel(iso4_top), 4)];
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.K, ...
%!    mat_ass.R, ...
%!    mat_ass.mat_info] = fem_ass_matrix(mesh, ...
%!                                       dof_map, ...
%!                                       [FEM_MAT_STIFFNESS, ...
%!                                        FEM_VEC_LOAD_CONSISTENT], ...
%!                                       load_case);
%!   if (eliminate)
%!     [Tred, Kred, Mred, Rred] = fem_cms_constr_elim(mesh, dof_map, mat_ass);
%!     opt_ls.refine_max_iter = int32(100);
%!     Kfact = fem_sol_factor(Kred, opt_ls);
%!     Ured = Kfact \ Rred;
%!     sol_stat.def = fem_post_def_nodal(mesh, dof_map, Tred * Ured);
%!   else
%!     opt_sol.verbose = int32(0);
%!     [sol_stat, U] = fem_sol_static(mesh, dof_map, mat_ass, opt_sol);
%!   endif
%!   sol_stat.F = fem_post_def_nodal(mesh, dof_map, mat_ass.R);
%!   Ftot = sum(sol_stat.F, 1);
%!   if (do_plot)
%!     figure("visible", "off");
%!     fem_post_sol_plot(mesh);
%!     xlabel("x [m]");
%!     ylabel("y [m]");
%!     zlabel("z [m]");
%!     grid on;
%!     grid minor on;
%!     title("undeformed mesh");
%!     figure("visible", "off");
%!     hold on;
%!     fem_post_sol_plot(mesh, sol_stat, scale/max(norm(sol_stat.def, "rows")), 1);
%!     view(30,30);
%!     xlabel('x [m]');
%!     ylabel('y [m]');
%!     zlabel('z [m]');
%!     grid on;
%!     grid minor on;
%!     title("deformed mesh");
%!   endif
%!   sigma_a = [-px; -py; -pz; zeros(3, 1)];
%!   epsilon_a = mesh.material_data(1).C \ sigma_a;
%!   U_a = zeros(rows(mesh.nodes), 3);
%!   for i=1:3
%!     U_a(:, i) = epsilon_a(i) * mesh.nodes(:, i);
%!   endfor
%!   Ftot_a = [-2 * b * c * px, -2 * a * c * py, -4 * a * b * pz, zeros(1, 3)];
%!   fprintf(stderr, "max(err)=%g\n", max(max(abs(sol_stat.def(:, 1:3) - U_a))) / max(max(abs(U_a))));
%!   fprintf(stderr, "mean(err)=%g\n", mean(mean(abs(sol_stat.def(:, 1:3) - U_a))) / max(max(abs(U_a))));
%!   assert_simple(sol_stat.def(:, 1:3), U_a, tolstat * max(max(abs(U_a))));
%!   assert_simple(Ftot, Ftot_a, sqrt(eps) * norm(Ftot_a));
%!   if (do_plot)
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
