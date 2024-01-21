## fem_pre_mesh_rbe3_from_surf.m:02
%!test
%! ## TEST2
%! close all;
%! filename = "";
%! unwind_protect
%! filename = tempname();
%! if (ispc())
%!  filename(filename == "\\") = "/";
%! endif
%! fd = -1;
%! unwind_protect
%! [fd, msg] = fopen([filename, ".geo"], "w");
%! if (fd == -1)
%!   error("failed to open file \"%s.geo\"", filename);
%! endif
%! a = 30e-3;
%! b = 20e-3;
%! c = 10e-3;
%! d = -5e-3;
%! e = 35e-3;
%! h = 4e-3;
%! fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%! fprintf(fd, "a=%g;\n", a);
%! fprintf(fd, "b=%g;\n", b);
%! fprintf(fd, "c=%g;\n", c);
%! fprintf(fd, "h = %g;\n", h);
%! fputs(fd, "Point(1) = {0.0,0.0,0.0,h};\n");
%! fputs(fd, "Point(2) = {a,0.0,0.0,h};\n");
%! fputs(fd, "Point(3) = {a,b,0.0,h};\n");
%! fputs(fd, "Point(4) = {0,b,0.0,h};\n");
%! fputs(fd, "Line(1) = {4,3};\n");
%! fputs(fd, "Line(2) = {3,2};\n");
%! fputs(fd, "Line(3) = {2,1};\n");
%! fputs(fd, "Line(4) = {1,4};\n");
%! fputs(fd, "Line Loop(5) = {2,3,4,1};\n");
%! fputs(fd, "Plane Surface(6) = {5};\n");
%! fputs(fd, "tmp[] = Extrude {0,0.0,c} {\n");
%! fputs(fd, "  Surface{6};\n");
%! fputs(fd, "};\n");
%! fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%! fputs(fd, "Physical Surface(\"modal\",1) = {tmp[4]};\n");
%! fputs(fd, "Physical Surface(\"interfaces\",2) = {tmp[2]};\n");
%! unwind_protect_cleanup
%! if (fd ~= -1)
%! fclose(fd);
%! endif
%! end_unwind_protect
%! pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "1", [filename, ".geo"]});
%! status = spawn_wait(pid);
%! if (status ~= 0)
%!  warning("gmsh failed with status %d", status);
%! endif
%! unlink([filename, ".geo"]);
%! mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%! node_id_ground = rows(mesh.nodes) + 1;
%! node_id_load = rows(mesh.nodes) + 2;
%! mesh.nodes(node_id_ground, :) = [d, 0.5 * b, 0.5 * c, 0, 0, 0];
%! mesh.nodes(node_id_load, :) = [e, 0.5 * b, 0.5 * c, 0, 0, 0];
%! load_case.locked_dof = false(rows(mesh.nodes), 6);
%! mesh.materials.iso8 = ones(rows(mesh.elements.iso8), 1, "int32");
%! mesh.material_data.E = 210000e6;
%! mesh.material_data.nu = 0.3;
%! mesh.material_data.rho = 7850;
%! mesh.elements.rbe3 = fem_pre_mesh_rbe3_from_surf(mesh, ...
%!                                                [1, 2], ...
%!                                                [node_id_ground, ...
%!                                                 node_id_load], "iso4");
%! mesh.elements.joints(1).nodes = node_id_ground;
%! mesh.elements.joints(1).C = eye(6);
%! mesh.elements.joints(2).nodes = node_id_load;
%! mesh.elements.joints(2).C = [0, 0, 0, 0, 1, 0;];
%! load_case.joints(1).U = zeros(6, 1);
%! load_case.joints(2).U = [1];
%! dof_map = fem_ass_dof_map(mesh, load_case);
%! [mat_ass.K, ...
%!  mat_ass.R] = fem_ass_matrix(mesh, ...
%!                              dof_map, ...
%!                              [FEM_MAT_STIFFNESS, ...
%!                               FEM_VEC_LOAD_CONSISTENT], ...
%!                              load_case);
%! sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%! figure("visible", "off");
%! opts.elem_types = {"iso8"};
%! fem_post_sol_plot(mesh, sol_stat, 20e-3 / max(norm(sol_stat.def(:, 1:3), "rows")), 1, opts);
%! xlabel("x [m]");
%! ylabel("y [m]");
%! zlabel("z [m]");
%! grid on;
%! grid minor on;
%! title("deformed mesh");
%! figure_list();
%! unwind_protect_cleanup
%! if (numel(filename))
%!   fn = dir([filename, "*"]);
%!   for i=1:numel(fn)
%!     unlink(fullfile(fn(i).folder, fn(i).name));
%!   endfor
%! endif
%! end_unwind_protect
