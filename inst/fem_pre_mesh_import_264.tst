## fem_pre_mesh_import.m:264
%!test
%! try
%! ### TEST 264
%! do_plot = false;
%! if (do_plot)
%!   close all;
%! endif
%! number_of_modes = 10;
%! scale_eig = 10e-3;
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
%!     a = 10e-3;
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
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "3", [filename, ".geo"]});
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
%!   mesh1.materials.tet20 = ones(rows(mesh1.elements.tet20), 1, "int32");
%!   E = 210000e6;
%!   nu = 0.3;
%!   mesh1.material_data.rho = 7850;
%!   mesh1.material_data.C = fem_pre_mat_isotropic(E, nu);
%!   data.mesh = mesh1;
%!   data = repmat(data, 1, 3);
%!   for i=2:3
%!     data(i).mesh.nodes(:, 1) += a * (i - 1);
%!     for j=1:numel(data(i).mesh.groups.tria10)
%!       data(i).mesh.groups.tria10(j).id += 100 * (i - 1);
%!     endfor
%!   endfor
%!   [mesh] = fem_post_mesh_merge(data);
%!   for i=1:2
%!     grp_idx_slave = find([[mesh.groups.tria10].id] == (i - 1) * 100 + 2);
%!     grp_idx_master = mesh.groups.tria10(find([[mesh.groups.tria10].id] == i * 100 + 1)).elements;
%!     mesh.elements.sfncon10(i).slave = mesh.groups.tria10(grp_idx_slave).nodes(:);
%!     mesh.elements.sfncon10(i).master = mesh.elements.tria10(grp_idx_master, :);
%!     mesh.elements.sfncon10(i).maxdist = 1e-6;
%!   endfor
%!   load_case.locked_dof = false(rows(mesh.nodes), 6);
%!   load_case.locked_dof(mesh.groups.tria10(find([[mesh.groups.tria10].id] == 1)).nodes, :) = true;
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
%!   sol_eig.stress = fem_ass_matrix(mesh, ...
%!                                   dof_map, ...
%!                                   [FEM_VEC_STRESS_CAUCH], ...
%!                                   load_case, ...
%!                                   sol_eig);
%!   if (do_plot)
%!     for i=1:numel(sol_eig.f)
%!       figure("visible", "off");
%!       hold on;
%!       fem_post_sol_plot(mesh, sol_eig, scale_eig/max(norm(sol_eig.def(:, 1:3,i), "rows")), i);
%!       view(30,30);
%!       xlabel('x [m]');
%!       ylabel('y [m]');
%!       zlabel('z [m]');
%!       grid on;
%!       grid minor on;
%!       title(sprintf("mode %d: %.1fHz", i, sol_eig.f(i)));
%!     endfor
%!   endif
%!   f_ref = [8768.74;
%!            14636.1;
%!            21145.7;
%!            39712.8;
%!            43555.5;
%!            47909;
%!            62270.4;
%!            84324.4;
%!            92665.1;
%!            94563];
%!   for i=1:length(sol_eig.f)
%!     fprintf(stderr, "mode %d f=%.0f\n", i, sol_eig.f(i));
%!   endfor
%!   assert_simple(sol_eig.f(:), f_ref, tol * max(f_ref));
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
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
