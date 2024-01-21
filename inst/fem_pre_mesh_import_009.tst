## fem_pre_mesh_import.m:09
%!test
%! ### TEST9
%! do_plot = false;
%! if (do_plot)
%!   close all;
%! endif
%! number_of_modes = 10;
%! scale_eig = 1e-3;
%! tol = 1e-2;
%! do_rotate = false;
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
%!     a = 15e-3;
%!     b = 20e-3;
%!     c = 10e-3;
%!     h = 3.5e-3;
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
%!     fputs(fd, "};\n");
%!     if (do_rotate)
%!       fputs(fd, "Rotate{{1, 0, 0},{0, 0, 0}, 30 * Pi / 180}{ Volume{tmp[1]}; }\n");
%!       fputs(fd, "Rotate{{0, 1, 0},{0, 0, 0}, 7 * Pi / 180}{ Volume{tmp[1]}; }\n");
%!       fputs(fd, "Rotate{{0, 0, 1},{0, 0, 0}, 15 * Pi / 180}{ Volume{tmp[1]}; }\n");
%!     endif
%!     fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%!     fputs(fd, "Physical Surface(\"clamp\",1) = {tmp[4]};\n");
%!     fputs(fd, "Physical Surface(\"load\",2) = {tmp[2]};\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   [~] = unlink([filename, ".geo"]);
%!   mesh1 = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!   [~] = unlink([filename, ".msh"]);
%!   if (~do_rotate)
%!     group_defs(1).id = 1;
%!     group_defs(1).name = "box1";
%!     group_defs(1).R = eye(3);
%!     group_defs(1).X0 = zeros(3, 1);
%!     group_defs(1).type = "box";
%!     group_defs(1).geometry.xmin = 0;
%!     group_defs(1).geometry.xmax = 0;
%!     group_defs(1).geometry.ymin = 0;
%!     group_defs(1).geometry.ymax = b;
%!     group_defs(1).geometry.zmin = 0;
%!     group_defs(1).geometry.zmax = c;
%!     group_defs(2).id = 2;
%!     group_defs(2).name = "cylinder1";
%!     group_defs(2).R = [-1, 0, 0;
%!                        0, 0, 1;
%!                        0, 1, 0];
%!     group_defs(2).X0 = [0; 0.5 * b; 0.5 * c];
%!     group_defs(2).type = "cylinder";
%!     group_defs(2).geometry.rmin = 0;
%!     group_defs(2).geometry.rmax = 0.5 * c;
%!     group_defs(2).geometry.zmin = -0.5 * b;
%!     group_defs(2).geometry.zmax = 0.5 * b;
%!     group_defs(3).id = 3;
%!     group_defs(3).name = "cylinder2";
%!     group_defs(3).R = [-1, 0, 0;
%!                        0, 0, 1;
%!                        0, 1, 0];
%!     group_defs(3).X0 = [a; 0.5 * b; 0.5 * c];
%!     group_defs(3).type = "cylinder";
%!     group_defs(3).geometry.rmin = 0;
%!     group_defs(3).geometry.rmax = 0.5 * c;
%!     group_defs(3).geometry.zmin = -0.5 * b;
%!     group_defs(3).geometry.zmax = 0.5 * b;
%!     group_defs(4).id = 4;
%!     group_defs(4).name = "box2";
%!     group_defs(4).R = eye(3);
%!     group_defs(4).X0 = [a; 0.5 * b; 0.5 * c];
%!     group_defs(4).type = "box";
%!     group_defs(4).geometry.xmin = 0;
%!     group_defs(4).geometry.xmax = 0;
%!     group_defs(4).geometry.ymin = -0.5 * b;
%!     group_defs(4).geometry.ymax = 0.5 * b;
%!     group_defs(4).geometry.zmin = -0.5 * c;
%!     group_defs(4).geometry.zmax = 0.5 * c;
%!   endif
%!   mesh1.materials.tet10 = ones(rows(mesh1.elements.tet10), 1, "int32");
%!   E = 210000e6;
%!   nu = 0.3;
%!   mesh1.material_data.rho = 7850;
%!   mesh1.material_data.C = fem_pre_mat_isotropic(E, nu);
%!   mesh2 = mesh1;
%!   mesh2.nodes(:, 1) += a;
%!   for i=1:numel(mesh2.groups.tria6)
%!     mesh2.groups.tria6(i).id += 100;
%!   endfor
%!   data(1).mesh = mesh1;
%!   data(2).mesh = mesh2;
%!   [mesh] = fem_post_mesh_merge(data);
%!   mesh.elements.sfncon6.slave = mesh.groups.tria6(find([[mesh.groups.tria6].id] == 2)).nodes(:);
%!   mesh.elements.sfncon6.master = mesh.elements.tria6(mesh.groups.tria6(find([[mesh.groups.tria6].id] == 101)).elements, :);
%!   mesh.elements.sfncon6.maxdist = eps^0.4 * max(abs([a,b,c]));
%!   mesh.elements.sfncon6.constraint = FEM_CT_SLIDING;
%!   load_case.locked_dof = false(rows(mesh.nodes), 6);
%!   load_case.locked_dof(mesh.groups.tria6(find([[mesh.groups.tria6].id] == 1)).nodes, :) = true;
%!   load_case.locked_dof(mesh.groups.tria6(find([[mesh.groups.tria6].id] == 102)).nodes, :) = true;
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.K, ...
%!    mat_ass.M, ...
%!    mtot] = fem_ass_matrix(mesh, ...
%!                           dof_map, ...
%!                           [FEM_MAT_STIFFNESS, ...
%!                            FEM_MAT_MASS, ...
%!                            FEM_SCA_TOT_MASS], ...
%!                           load_case);
%!   assert_simple(mtot, a * b * c * sum([mesh.material_data.rho]), sqrt(eps) * a * b * c * sum([mesh.material_data.rho]));
%!   [sol_eig] = fem_sol_modal(mesh, dof_map, mat_ass, number_of_modes);
%!   if (do_plot)
%!     for i=1:numel(sol_eig.f)
%!       figure("visible", "off");
%!       hold on;
%!       fem_post_sol_plot(mesh, sol_eig, scale_eig/max(norm(sol_eig.def(:, 1:3, i), "rows")), i);
%!       view(30,30);
%!       xlabel('x [m]');
%!       ylabel('y [m]');
%!       zlabel('z [m]');
%!       grid on;
%!       grid minor on;
%!       title(sprintf("mode %d: %.1fHz", i, sol_eig.f(i)));
%!     endfor
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
