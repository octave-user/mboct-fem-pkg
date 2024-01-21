## fem_ehd_pre_comp_mat_unstruct.m:01
%!test
%! close all;
%! fd = -1;
%! filename = "";
%! unwind_protect
%! filename = tempname();
%! if (ispc())
%!   filename(filename == "\\") = "/";
%! endif
%! unwind_protect
%! [fd, msg] = fopen([filename, ".geo"], "w");
%! if (fd == -1)
%!   error("failed to open file \"%s.geo\"", filename);
%! endif
%! ri = 8e-3;
%! ro = 10e-3;
%! h = 12e-3;
%! c = 2e-3;
%! b = h - 2 * c;
%! mesh_size = 1.5e-3;
%! fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%! fprintf(fd, "ri = %g;\n", ri);
%! fprintf(fd, "ro = %g;\n", ro);
%! fprintf(fd, "h = %g;\n", h);
%! fprintf(fd, "c = %g;\n", c);
%! fputs(fd, "Point(1) = {ri,0.0,0.0};\n");
%! fputs(fd, "Point(2) = {ro,0.0,0.0};\n");
%! fputs(fd, "Point(3) = {ro,0.0,c};\n");
%! fputs(fd, "Point(4) = {ro,0.0,h - c};\n");
%! fputs(fd, "Point(5) = {ro,0.0,h};\n");
%! fputs(fd, "Point(6) = {ri,0.0,h};\n");
%! fputs(fd, "Point(7) = {ri,0.0,h - c};\n");
%! fputs(fd, "Point(8) = {ri,0.0,c};\n");
%! fputs(fd, "Line(1) = {1,2};\n");
%! fputs(fd, "Line(2) = {2,3};\n");
%! fputs(fd, "Line(3) = {3,4};\n");
%! fputs(fd, "Line(4) = {4,5};\n");
%! fputs(fd, "Line(5) = {5,6};\n");
%! fputs(fd, "Line(6) = {6,7};\n");
%! fputs(fd, "Line(7) = {7,8};\n");
%! fputs(fd, "Line(8) = {8,1};\n");
%! fputs(fd, "Line Loop(5) = {1,2,3,4,5,6,7,8};\n");
%! fputs(fd, "Plane Surface(6) = {5};\n");
%! fputs(fd, "tmp[] = Extrude {{0, 0, 1}, {0, 0, 0}, 2*Pi} { Surface{6}; };\n");
%! fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%! fputs(fd, "Physical Surface(\"clamp\",1) = {tmp[2]};\n");
%! fputs(fd, "Physical Surface(\"load1\",2) = {tmp[4]};\n");
%! fputs(fd, "Physical Surface(\"load2\",3) = {tmp[8]};\n");
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     fclose(fd);
%!   endif
%! end_unwind_protect
%! pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-clmin", sprintf("%g", 0.75 * mesh_size), "-clmax", sprintf("%g", 1.25 *mesh_size), [filename, ".geo"]});
%! status = spawn_wait(pid);
%! if (status ~= 0)
%!  error("gmsh failed with status %d", status);
%! endif
%! mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%! fprintf(stderr, "%d nodes\n", rows(mesh.nodes));
%! grp_id_clamp = find([[mesh.groups.tria6].id] == 1);
%! grp_id_p1 = find([[mesh.groups.tria6].id] == 3);
%! grp_id_p2 = find([[mesh.groups.tria6].id] == 2);
%! node_idx_modal = rows(mesh.nodes) + 1;
%! node_idx_itf1 = rows(mesh.nodes) + 2;
%! node_idx_itf2 = rows(mesh.nodes) + 3;
%! mesh.nodes(node_idx_itf2, 1:3) = [0, 0, c + 0.5 * b];
%! mesh.nodes(node_idx_itf1, 1:3) = [0, 0, c + 0.5 * b];
%! mesh.nodes(node_idx_modal, :) = 0;
%! mesh.elements.rbe3 = fem_pre_mesh_rbe3_from_surf(mesh, [grp_id_p1, grp_id_p2], [node_idx_itf1, node_idx_itf2], "tria6");
%! mat_type = {"nodal", "nodal substruct", "nodal substruct total", "modal substruct total"};
%! for imat=1:numel(mat_type)
%! data(imat).mesh = mesh;
%! data(imat).bearing_surf(1).group_idx = grp_id_p1;
%! data(imat).bearing_surf(1).options.reference_pressure = 1e9;
%! data(imat).bearing_surf(1).options.mesh_size = 1.5 * mesh_size;
%! data(imat).bearing_surf(1).options.bearing_type = "shell";
%! data(imat).bearing_surf(1).options.matrix_type = mat_type{imat};
%! data(imat).bearing_surf(1).master_node_no = node_idx_itf1;
%! data(imat).bearing_surf(1).nodes = mesh.groups.tria6(grp_id_p1).nodes;
%! data(imat).bearing_surf(1).r = ri;
%! data(imat).bearing_surf(1).w = b;
%! data(imat).bearing_surf(1).X0 = [0; 0; b/2 + c];
%! data(imat).bearing_surf(1).R = eye(3);
%! data(imat).bearing_surf(1).relative_tolerance = 0;
%! data(imat).bearing_surf(1).absolute_tolerance = sqrt(eps) * ri;
%! data(imat).bearing_surf(2).group_idx = grp_id_p2;
%! data(imat).bearing_surf(2).options.reference_pressure = 1e9;
%! data(imat).bearing_surf(2).options.mesh_size = 1.5 * mesh_size;
%! data(imat).bearing_surf(2).options.bearing_type = "journal";
%! data(imat).bearing_surf(2).options.matrix_type = mat_type{imat};
%! data(imat).bearing_surf(2).master_node_no = node_idx_itf2;
%! data(imat).bearing_surf(2).nodes = mesh.groups.tria6(grp_id_p2).nodes;
%! data(imat).bearing_surf(2).r = ro;
%! data(imat).bearing_surf(2).w = b;
%! data(imat).bearing_surf(2).X0 = [0; 0; b/2 + c];
%! data(imat).bearing_surf(2).R = eye(3);
%! data(imat).bearing_surf(2).relative_tolerance = 0;
%! data(imat).bearing_surf(2).absolute_tolerance = sqrt(eps) * ri;
%! [data(imat).load_case_press, data(imat).bearing_surf] = fem_ehd_pre_comp_mat_load_case(data(imat).mesh, data(imat).bearing_surf);
%! data(imat).load_case.locked_dof = false(rows(data(imat).mesh.nodes), 6);
%! data(imat).load_case.locked_dof(data(imat).mesh.groups.tria6(grp_id_clamp).nodes, :) = true;
%! data(imat).load_case = fem_pre_load_case_merge(data(imat).load_case, data(imat).load_case_press);
%! data(imat).mesh.materials.tet10 = ones(rows(data(imat).mesh.elements.tet10), 1, "int32");
%! E = 210000e6;
%! nu = 0.3;
%! data(imat).mesh.material_data.rho = 7850;
%! data(imat).mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%! data(imat).cms_opt.nodes.modal.number = node_idx_modal;
%! data(imat).cms_opt.nodes.modal.name = "node_id_modal";
%! data(imat).cms_opt.nodes.interfaces(1).number = node_idx_itf1;
%! data(imat).cms_opt.nodes.interfaces(1).name = "node_id_itf1";
%! data(imat).cms_opt.nodes.interfaces(2).number = node_idx_itf2;
%! data(imat).cms_opt.nodes.interfaces(2).name = "node_id_itf2";
%! data(imat).cms_opt.number_of_threads = mbdyn_solver_num_threads_default();
%! [data(imat).mesh, data(imat).mat_ass, data(imat).dof_map, data(imat).sol_eig, data(imat).cms_opt] = fem_cms_create(data(imat).mesh, data(imat).load_case, data(imat).cms_opt);
%! data(imat).comp_mat = fem_ehd_pre_comp_mat_unstruct(data(imat).mesh, data(imat).mat_ass, data(imat).dof_map, data(imat).cms_opt, data(imat).bearing_surf);
%! data(imat).opt_plot.plot_nodal = true;
%! data(imat).opt_plot.contour_levels = int32(15);
%! data(imat).opt_plot.plot_step.x = int32(2);
%! data(imat).opt_plot.plot_step.z = int32(2);
%! fem_ehd_pre_comp_mat_plot(data(imat).comp_mat, data(imat).opt_plot);
%! for i=1:numel(data(imat).comp_mat)
%!   fem_ehd_pre_comp_mat_export(data(imat).comp_mat(i), data(imat).bearing_surf(i).options, sprintf("%s_%03d.dat", filename, i));
%! endfor
%! endfor
%! figure_list();
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       status = unlink(fullfile(fn(i).folder, fn(i).name));
%!       if (status ~= 0)
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect
