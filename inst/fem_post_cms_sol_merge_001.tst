## fem_post_cms_sol_merge.m:01
%!test
%! try
%! ## TEST
%! close all;
%! number_of_modes = 10;
%! a = [30e-3, 100e-3];
%! b = [20e-3, 20e-3];
%! c = [10e-3, 5e-3];
%! d = [0, 40e-3];
%! h = 3.5e-3;
%! p = [25e6, 25e6];
%! f_run_post_proc = false;
%! filename = "";
%! unwind_protect
%! filename = tempname();
%! if (ispc())
%!   filename(filename == "\\") = "/";
%! endif
%! for i=1:numel(a)
%!   unwind_protect
%!     fd = -1;
%!     unwind_protect
%!       [fd, msg] = fopen([filename, ".geo"], "w");
%!       if (fd == -1)
%!         error("failed to open file \"%s.geo\"", filename);
%!       endif
%!       unwind_protect
%!         fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!         fprintf(fd, "a=%g;\n", a(i));
%!         fprintf(fd, "b=%g;\n", b(i));
%!         fprintf(fd, "c=%g;\n", c(i));
%!         fprintf(fd, "d=%g;\n", d(i));
%!         fprintf(fd, "h = %g;\n", h);
%!         fputs(fd, "Point(1) = {d,0.0,0.0,h};\n");
%!         fputs(fd, "Point(2) = {a+d,0.0,0.0,h};\n");
%!         fputs(fd, "Point(3) = {a+d,b,0.0,h};\n");
%!         fputs(fd, "Point(4) = {d,b,0.0,h};\n");
%!         fputs(fd, "Line(1) = {4,3};\n");
%!         fputs(fd, "Line(2) = {3,2};\n");
%!         fputs(fd, "Line(3) = {2,1};\n");
%!         fputs(fd, "Line(4) = {1,4};\n");
%!         fputs(fd, "Line Loop(5) = {2,3,4,1};\n");
%!         fputs(fd, "Plane Surface(6) = {5};\n");
%!         fputs(fd, "tmp[] = Extrude {0,0.0,c} {\n");
%!         fputs(fd, "  Surface{6};\n");
%!         fputs(fd, "};\n");
%!         fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%!         fputs(fd, "Physical Surface(\"clamp\",1) = {tmp[4]};\n");
%!         fputs(fd, "Physical Surface(\"load\",2) = {tmp[3]};\n");
%!       unwind_protect_cleanup
%!         if (fd ~= -1)
%!           fclose(fd);
%!         endif
%!       end_unwind_protect
%!       fprintf(stderr, "meshing ...\n");
%!       pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", [filename, ".geo"]});
%!       status = spawn_wait(pid);
%!       if (status ~= 0)
%!         warning("gmsh failed with status %d", status);
%!       endif
%!     unwind_protect_cleanup
%!       [~] = unlink([filename, ".geo"]);
%!     end_unwind_protect
%!     fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%!     mesh_data(i).mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!     fprintf(stderr, "%d nodes\n", rows(mesh_data(1).mesh.nodes));
%!   unwind_protect_cleanup
%!     [~] = unlink([filename, ".msh"]);
%!   end_unwind_protect
%! endfor
%! for i=1:numel(mesh_data)
%!   fprintf(stderr, "assembling matrices ...\n");
%!   mesh_data(i).load_case.locked_dof = false(rows(mesh_data(i).mesh.nodes), 6);
%!   mesh_data(i).load_case.locked_dof(mesh_data(i).mesh.groups.tria6(find([[mesh_data(i).mesh.groups.tria6].id] == 1)).nodes, :) = true;
%!   mesh_data(i).load_case.pressure.tria6.elements = mesh_data(i).mesh.elements.tria6(mesh_data(i).mesh.groups.tria6(find([mesh_data(i).mesh.groups.tria6.id] == 2)).elements, :);
%!   mesh_data(i).load_case.pressure.tria6.p = repmat(p(i),size(mesh_data(i).load_case.pressure.tria6.elements));
%!   mesh_data(i).mesh.materials.tet10 = ones(rows(mesh_data(i).mesh.elements.tet10), 1, "int32");
%!   E = 210000e6;
%!   nu = 0.3;
%!   mesh_data(i).mesh.material_data.rho = 7850;
%!   mesh_data(i).mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%!   mesh_data(i).dof_map = fem_ass_dof_map(mesh_data(i).mesh, mesh_data(i).load_case);
%!   [mesh_data(i).mat_ass.K, ...
%!    mesh_data(i).mat_ass.R] = fem_ass_matrix(mesh_data(i).mesh, ...
%!                                             mesh_data(i).dof_map, ...
%!                                             [FEM_MAT_STIFFNESS, ...
%!                                              FEM_VEC_LOAD_CONSISTENT], mesh_data(i).load_case);
%!   sol.bodies(i).def = fem_sol_static(mesh_data(i).mesh, mesh_data(i).dof_map, mesh_data(i).mat_ass).def;
%!   sol.bodies(i).stress = fem_ass_matrix(mesh_data(i).mesh, ...
%!                                         mesh_data(i).dof_map, ...
%!                                         [FEM_VEC_STRESS_CAUCH], ...
%!                                         mesh_data(i).load_case, ...
%!                                         sol.bodies(i));
%! endfor
%! [mesh_comb, dof_map_comb] = fem_post_mesh_merge(mesh_data);
%! [sol_comb] = fem_post_cms_sol_merge(mesh_comb, dof_map_comb, sol);
%! opts.rotation_angle = [0, 0, 0];
%! opts.print_and_exit = 1;
%! opts.print_to_file = filename;
%! opts.skin_only = true;
%! opts.show_element = true;
%! if (f_run_post_proc)
%! unwind_protect
%!   fem_post_sol_external(mesh_comb, sol_comb, opts);
%!   [img, map, alpha] = imread([opts.print_to_file, "_001.jpg"]);
%!   figure("visible", "off");
%!   imshow(img, map);
%!   title("Gmsh - stress and deformation");
%! unwind_protect_cleanup
%!   [~] = unlink([opts.print_to_file, "_001.jpg"]);
%! end_unwind_protect
%! figure_list();
%! endif
%! unwind_protect_cleanup
%! if (numel(filename))
%!   fn = dir([filename, "*"]);
%!   for i=1:numel(fn)
%!     [~] = unlink(fullfile(fn(i).folder, fn(i).name));
%!   endfor
%! endif
%! end_unwind_protect
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
