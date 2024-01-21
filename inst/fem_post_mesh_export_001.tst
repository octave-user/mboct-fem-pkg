## fem_post_mesh_export.m:01
%!test
%! ## TEST 1
%! ## K.J.Bathe 2002, page 328 4.20a
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
%! mesh_size = 20;
%! p1 = 0.006;
%! E = 55;
%! nu = 0.3;
%! rho = 1000e-12;
%! fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%! fputs(fd, "Point(1) = { 0.0, 0.0,-20.0};\n");
%! fputs(fd, "Point(2) = {10.0, 0.0,-20.0};\n");
%! fputs(fd, "Point(3) = {10.0, 0.0,-10.0};\n");
%! fputs(fd, "Point(4) = {15.0, 0.0, -5.0};\n");
%! fputs(fd, "Point(5) = {65.0, 0.0, -5.0};\n");
%! fputs(fd, "Point(6) = {65.0, 0.0,  5.0};\n");
%! fputs(fd, "Point(7) = {15.0, 0.0,  5.0};\n");
%! fputs(fd, "Point(8) = {10.0, 0.0, 10.0};\n");
%! fputs(fd, "Point(9) = {10.0, 0.0, 20.0};\n");
%! fputs(fd, "Point(10)= { 0.0, 0.0, 20.0};\n");
%! fputs(fd, "Point(11)= {15.0, 0.0,-10.0};\n");
%! fputs(fd, "Point(12)= {15.0, 0.0, 10.0};\n");
%! fputs(fd, "Line(1) = {1,2};\n");
%! fputs(fd, "Line(2) = {2,3};\n");
%! fputs(fd, "Circle(3) = {3,11,4};\n");
%! fputs(fd, "Line(4) = {4,5};\n");
%! fputs(fd, "Line(5) = {5,6};\n");
%! fputs(fd, "Line(6) = {6,7};\n");
%! fputs(fd, "Circle(7) = {7,12,8};\n");
%! fputs(fd, "Line(8) = {8,9};\n");
%! fputs(fd, "Line(9) = {9,10};\n");
%! fputs(fd, "Line(10) = {10,1};\n");
%! fputs(fd, "Line Loop(11) = {1,2,3,4,5,6,7,8,9,10};\n");
%! fputs(fd, "Plane Surface(14) = {11};\n");
%! fprintf(fd, "tmp[] = Extrude {0, %g, 0}{ Surface{14}; };\n", mesh_size);
%! fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%! fputs(fd, "Physical Surface(\"clamp\",2) = {tmp[11]};\n");
%! fputs(fd, "Physical Surface(\"pressure\",3) = {tmp[7],tmp[8],tmp[9]};\n");
%! fputs(fd, "Physical Surface(\"displacement\",4) = {tmp[6]};\n");
%! fputs(fd, "Physical Surface(\"stress\",5) = {tmp[8]};\n");
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     fclose(fd);
%!   endif
%! end_unwind_protect
%! fprintf(stderr, "meshing ...\n");
%! [~] = unlink([filename, ".msh"]);
%! pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-clmin", sprintf("%g", 0.75 * mesh_size), "-clmax", sprintf("%g", 1.25 * mesh_size), [filename, ".geo"]});
%! status = spawn_wait(pid);
%! if (status ~= 0)
%!  warning("gmsh failed with status %d", status);
%! endif
%! [~] = unlink([filename, ".geo"]);
%! fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%! mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%! fprintf(stderr, "%d nodes\n", rows(mesh.nodes));
%! [~] = unlink([filename, ".msh"]);
%! load_case.locked_dof = false(rows(mesh.nodes), 6);
%! grp_id_clamp = find([[mesh.groups.tria6].id] == 2);
%! load_case.locked_dof(mesh.groups.tria6(grp_id_clamp).nodes, 1:3) = true;
%! grp_id_p1 = find([[mesh.groups.tria6].id] == 3);
%! elem_id_p1 = mesh.groups.tria6(grp_id_p1).elements;
%! elno_p1 = mesh.elements.tria6(elem_id_p1, :);
%! load_case.pressure.tria6.elements = elno_p1;
%! load_case.pressure.tria6.p = [repmat(p1, rows(elno_p1), columns(elno_p1))];
%! mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%! mesh.material_data.rho = rho;
%! mesh.material_data.E = E;
%! mesh.material_data.nu = nu;
%! dof_map = fem_ass_dof_map(mesh, load_case);
%! opts.format = "apdl";
%! fem_post_mesh_export([filename, ".dat"], mesh, opts, load_case, dof_map);
%! spawn_wait(spawn("cat", {[filename, ".dat"]}));
%! unwind_protect_cleanup
%! if (numel(filename))
%!   fn = dir([filename, "*"]);
%!   for i=1:numel(fn)
%!     [~] = unlink(fullfile(fn(i).folder, fn(i).name));
%!   endfor
%! endif
%! end_unwind_protect
