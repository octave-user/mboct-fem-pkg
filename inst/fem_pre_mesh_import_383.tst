## fem_pre_mesh_import.m:383
%!test
%! try
%! ### TEST383
%! do_plot = false;
%! if (do_plot)
%!   close all;
%! endif
%! number_of_modes = 10;
%! scale_eig = 10e-3;
%! tol = 1e-2;
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
%!     a = 30e-3;
%!     b = 20e-3;
%!     c = 10e-3;
%!     h = 2e-3;
%!     p = 25e6;
%!     fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "a=%g;\n", a);
%!     fprintf(fd, "b=%g;\n", b);
%!     fprintf(fd, "c=%g;\n", c);
%!     fprintf(fd, "h = %g;\n", h);
%!     fputs(fd, "Mesh.SecondOrderIncomplete=1;\n");
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
%!     fputs(fd, "  Surface{6}; Layers{Ceil(c/h)}; Recombine;\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "Recombine Surface{6, tmp[0]};\n");
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
%!   opt_msh.elem_type = {"iso20r", "quad8"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_msh));
%!   [~] = unlink([filename, ".msh"]);
%!   group_defs(1).id = 1;
%!   group_defs(1).name = "box1";
%!   group_defs(1).R = eye(3);
%!   group_defs(1).X0 = zeros(3, 1);
%!   group_defs(1).type = "box";
%!   group_defs(1).geometry.xmin = 0;
%!   group_defs(1).geometry.xmax = 0;
%!   group_defs(1).geometry.ymin = 0;
%!   group_defs(1).geometry.ymax = b;
%!   group_defs(1).geometry.zmin = 0;
%!   group_defs(1).geometry.zmax = c;
%!   group_defs(1).elem_type = "quad8";
%!   group_defs(2).id = 2;
%!   group_defs(2).name = "cylinder1";
%!   group_defs(2).R = [-1, 0, 0;
%!                      0, 0, 1;
%!                      0, 1, 0];
%!   group_defs(2).X0 = [0; 0.5 * b; 0.5 * c];
%!   group_defs(2).type = "cylinder";
%!   group_defs(2).geometry.rmin = 0;
%!   group_defs(2).geometry.rmax = 0.5 * c;
%!   group_defs(2).geometry.zmin = -0.5 * b;
%!   group_defs(2).geometry.zmax = 0.5 * b;
%!   group_defs(2).elem_type = "quad8";
%!   group_defs(3).id = 3;
%!   group_defs(3).name = "cylinder2";
%!   group_defs(3).R = [-1, 0, 0;
%!                      0, 0, 1;
%!                      0, 1, 0];
%!   group_defs(3).X0 = [a; 0.5 * b; 0.5 * c];
%!   group_defs(3).type = "cylinder";
%!   group_defs(3).geometry.rmin = 0;
%!   group_defs(3).geometry.rmax = 0.5 * c;
%!   group_defs(3).geometry.zmin = -0.5 * b;
%!   group_defs(3).geometry.zmax = 0.5 * b;
%!   group_defs(3).elem_type = "quad8";
%!   group_defs(4).id = 4;
%!   group_defs(4).name = "box2";
%!   group_defs(4).R = eye(3);
%!   group_defs(4).X0 = [a; 0.5 * b; 0.5 * c];
%!   group_defs(4).type = "box";
%!   group_defs(4).geometry.xmin = 0;
%!   group_defs(4).geometry.xmax = 0;
%!   group_defs(4).geometry.ymin = -0.5 * b;
%!   group_defs(4).geometry.ymax = 0.5 * b;
%!   group_defs(4).geometry.zmin = -0.5 * c;
%!   group_defs(4).geometry.zmax = 0.5 * c;
%!   group_defs(4).elem_type = "quad8";
%!   groups = fem_pre_mesh_groups_create(mesh, group_defs, sqrt(eps));
%!   assert_simple(numel(groups.quad8), 4);
%!   assert_simple([groups.quad8.id], [group_defs.id]);
%!   assert_simple(sort(groups.quad8(1).nodes), sort(mesh.groups.quad8(1).nodes));
%!   assert_simple(sort(groups.quad8(2).nodes), sort(mesh.groups.quad8(1).nodes));
%!   assert_simple(sort(groups.quad8(3).nodes), sort(mesh.groups.quad8(2).nodes));
%!   assert_simple(sort(groups.quad8(4).nodes), sort(mesh.groups.quad8(2).nodes));
%!   assert_simple(groups.quad8(1).elements, mesh.groups.quad8(1).elements);
%!   assert_simple(groups.quad8(2).elements, mesh.groups.quad8(1).elements);
%!   assert_simple(groups.quad8(3).elements, mesh.groups.quad8(2).elements);
%!   assert_simple(groups.quad8(4).elements, mesh.groups.quad8(2).elements);
%!   load_case.locked_dof = false(rows(mesh.nodes), 6);
%!   load_case.locked_dof(mesh.groups.quad8(find([[mesh.groups.quad8].id] == 1)).nodes, :) = true;
%!   load_case.pressure.quad8.elements = mesh.elements.quad8(mesh.groups.quad8(find([mesh.groups.quad8.id] == 2)).elements, :);
%!   mesh.materials.iso20r = ones(rows(mesh.elements.iso20r), 1, "int32");
%!   E = 210000e6;
%!   nu = 0.3;
%!   mesh.material_data.rho = 7850;
%!   mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.K, ...
%!    mat_ass.M, ...
%!    mtot] = fem_ass_matrix(mesh, ...
%!                           dof_map, ...
%!                           [FEM_MAT_STIFFNESS, ...
%!                            FEM_MAT_MASS, ...
%!                            FEM_SCA_TOT_MASS], load_case);
%!   assert_simple(mtot, a * b * c * mesh.material_data.rho, sqrt(eps) * a * b * c * mesh.material_data.rho);
%!   [sol_eig] = fem_sol_modal(mesh, dof_map, mat_ass, number_of_modes);
%!   sol_eig.stress = fem_ass_matrix(mesh, ...
%!                                   dof_map, ...
%!                                   [FEM_VEC_STRESS_CAUCH], ...
%!                                   load_case, ...
%!                                   sol_eig);
%!   X = mesh.nodes(unique(load_case.pressure.quad8.elements), 1:3).';
%!   dof_idx = dof_map.ndof(unique(load_case.pressure.quad8.elements), 1:3);
%!   f = sol_eig.f(:);
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
%!   for i=1:length(f)
%!     fprintf(stderr, "mode %d f=%.0f\n", i, f(i));
%!   endfor
%!   assert_simple(f, f_ref, tol * max(f_ref));
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
