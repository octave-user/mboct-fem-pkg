## fem_pre_mesh_groups_create.m:01
%!test
%! try
%! close all;
%! filename = "";
%! unwind_protect
%! filename = tempname();
%! if (ispc())
%!   filename(filename == "\\") = "/";
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
%! h = 3.5e-3;
%! p = 25e6;
%! fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
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
%! fputs(fd, "Physical Surface(\"clamp\",1) = {tmp[4]};\n");
%! fputs(fd, "Physical Surface(\"load\",2) = {tmp[2]};\n");
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     fclose(fd);
%!   endif
%! end_unwind_protect
%! pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", [filename, ".geo"]});
%! status = spawn_wait(pid);
%! if (status ~= 0)
%!  warning("gmsh failed with status %d", status);
%! endif
%! [~] = unlink([filename, ".geo"]);
%! mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%! [~] = unlink([filename, ".msh"]);
%! group_defs(1).id = 1;
%! group_defs(1).name = "box1";
%! group_defs(1).R = eye(3);
%! group_defs(1).X0 = zeros(3, 1);
%! group_defs(1).type = "box";
%! group_defs(1).geometry.xmin = 0;
%! group_defs(1).geometry.xmax = 0;
%! group_defs(1).geometry.ymin = 0;
%! group_defs(1).geometry.ymax = b;
%! group_defs(1).geometry.zmin = 0;
%! group_defs(1).geometry.zmax = c;
%! group_defs(1).elem_type = "tria6";
%! group_defs(2).id = 2;
%! group_defs(2).name = "cylinder1";
%! group_defs(2).R = [-1, 0, 0;
%!                     0, 0, 1;
%!                     0, 1, 0];
%! group_defs(2).X0 = [0; 0.5 * b; 0.5 * c];
%! group_defs(2).type = "cylinder";
%! group_defs(2).geometry.rmin = 0;
%! group_defs(2).geometry.rmax = 0.5 * c;
%! group_defs(2).geometry.zmin = -0.5 * b;
%! group_defs(2).geometry.zmax = 0.5 * b;
%! group_defs(2).elem_type = "tria6";
%! group_defs(3).id = 3;
%! group_defs(3).name = "cylinder2";
%! group_defs(3).R = [-1, 0, 0;
%!                     0, 0, 1;
%!                     0, 1, 0];
%! group_defs(3).X0 = [a; 0.5 * b; 0.5 * c];
%! group_defs(3).type = "cylinder";
%! group_defs(3).geometry.rmin = 0;
%! group_defs(3).geometry.rmax = 0.5 * c;
%! group_defs(3).geometry.zmin = -0.5 * b;
%! group_defs(3).geometry.zmax = 0.5 * b;
%! group_defs(3).elem_type = "tria6";
%! group_defs(4).id = 4;
%! group_defs(4).name = "box2";
%! group_defs(4).R = eye(3);
%! group_defs(4).X0 = [a; 0.5 * b; 0.5 * c];
%! group_defs(4).type = "box";
%! group_defs(4).geometry.xmin = 0;
%! group_defs(4).geometry.xmax = 0;
%! group_defs(4).geometry.ymin = -0.5 * b;
%! group_defs(4).geometry.ymax = 0.5 * b;
%! group_defs(4).geometry.zmin = -0.5 * c;
%! group_defs(4).geometry.zmax = 0.5 * c;
%! group_defs(4).elem_type = "tria6";
%! tol_rel = sqrt(eps);
%! mesh.groups = fem_pre_mesh_groups_create(mesh, group_defs, tol_rel);
%! figure("visible", "off");
%! fem_post_sol_plot(mesh);
%! xlabel("x [m]");
%! ylabel("y [m]");
%! zlabel("z [m]");
%! grid on;
%! grid minor on;
%! title("undeformed mesh");
%! view(30, 30);
%! figure_list();
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
