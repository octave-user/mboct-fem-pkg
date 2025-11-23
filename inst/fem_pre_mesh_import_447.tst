## fem_pre_mesh_import.m:447
%!test
%! try
%! ### TEST1
%! do_plot = false;
%! if (do_plot)
%!   close all;
%! endif
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   F = 1000;
%!   M = 25000*0;
%!   R = 762e-3;
%!   da = 50.8e-3;
%!   s = 20e-3;
%!   rm = (da - s) / 2;
%!   K = R * s / (1.65 * rm^2);
%!   h = 20e-3;
%!   E = 210000e6;
%!   nu = 0.3;
%!   rho = 7850;
%!   N = 10;
%!   di = da - 2 * s;
%!   Iy = (da^4 - di^4) * pi / 64;
%!   unwind_protect
%!     [fd, msg] = fopen([filename, ".geo"], "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s.geo\"", filename);
%!     endif
%!     fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "R = %g;\n", R);
%!     fprintf(fd, "da = %g;\n", da);
%!     fprintf(fd, "di = %g;\n", di);
%!     fprintf(fd, "s = %g;\n", s);
%!     fprintf(fd, "h = %g;\n", h);
%!     fputs(fd, "Point(1) = {0.0, 0.0, 0.0, h};\n");
%!     fputs(fd, "Point(2) = {0.0, 0.5 * da, 0.0, h};\n");
%!     fputs(fd, "Point(3) = {-0.5 * da, 0.0, 0.0, h};\n");
%!     fputs(fd, "Point(4) = {0.0, -0.5 * da, 0.0, h};\n");
%!     fputs(fd, "Point(5) = {0.5 * da, 0.0, 0.0, h};\n");
%!     fputs(fd, "Point(6) = {0.0, 0.5 * da - s, 0.0, h};\n");
%!     fputs(fd, "Point(7) = {-(0.5 * da - s), 0.0, 0.0, h};\n");
%!     fputs(fd, "Point(8) = {0.0, -(0.5 * da - s), 0.0, h};\n");
%!     fputs(fd, "Point(9) = {0.5 * da - s, 0.0, 0.0, h};\n");
%!     fputs(fd, "Point(10) = {0,-R,R};\n");
%!     fputs(fd, "Circle(1) = {2, 1, 3};\n");
%!     fputs(fd, "Circle(2) = {3, 1, 4};\n");
%!     fputs(fd, "Circle(3) = {4, 1, 5};\n");
%!     fputs(fd, "Circle(4) = {5, 1, 2};\n");
%!     fputs(fd, "Circle(5) = {6, 1, 7};\n");
%!     fputs(fd, "Circle(6) = {7, 1, 8};\n");
%!     fputs(fd, "Circle(7) = {8, 1, 9};\n");
%!     fputs(fd, "Circle(8) = {9, 1, 6};\n");
%!     fputs(fd, "Line(9) = {2, 6};\n");
%!     fputs(fd, "Line(10) = {5, 9};\n");
%!     fputs(fd, "Line(11) = {4, 8};\n");
%!     fputs(fd, "Line(12) = {3, 7};\n");
%!     fputs(fd, "Transfinite Curve(1) = Round(0.25 * da * Pi / h) + 1;\n");
%!     fputs(fd, "Transfinite Curve(2) = Round(0.25 * da * Pi / h) + 1;\n");
%!     fputs(fd, "Transfinite Curve(3) = Round(0.25 * da * Pi / h) + 1;\n");
%!     fputs(fd, "Transfinite Curve(4) = Round(0.25 * da * Pi / h) + 1;\n");
%!     fputs(fd, "Transfinite Curve(5) = Round(0.25 * da * Pi / h) + 1;\n");
%!     fputs(fd, "Transfinite Curve(6) = Round(0.25 * da * Pi / h) + 1;\n");
%!     fputs(fd, "Transfinite Curve(7) = Round(0.25 * da * Pi / h) + 1;\n");
%!     fputs(fd, "Transfinite Curve(8) = Round(0.25 * da * Pi / h) + 1;\n");
%!     fputs(fd, "Transfinite Curve(9) = Round(0.5 * (da - di)  / h) + 1;\n");
%!     fputs(fd, "Transfinite Curve(10) = Round(0.5 * (da - di)  / h) + 1;\n");
%!     fputs(fd, "Transfinite Curve(11) = Round(0.5 * (da - di)  / h) + 1;\n");
%!     fputs(fd, "Transfinite Curve(12) = Round(0.5 * (da - di)  / h) + 1;\n");
%!     fputs(fd, "Curve Loop(100) = {4,10,-8,-9};\n");
%!     fputs(fd, "Curve Loop(101) = {3,11,7,10};\n");
%!     fputs(fd, "Curve Loop(102) = {2,12,-6,-11};\n");
%!     fputs(fd, "Curve Loop(103) = {1,9,-5,-12};\n");
%!     fputs(fd, "Plane Surface(200) = {100};\n");
%!     fputs(fd, "Plane Surface(201) = {101};\n");
%!     fputs(fd, "Plane Surface(202) = {102};\n");
%!     fputs(fd, "Plane Surface(203) = {103};\n");
%!     fputs(fd, "Transfinite Surface(13) = {2, 5, 9, 6};\n");
%!     fputs(fd, "Transfinite Surface(14) = {5, 4, 8, 9};\n");
%!     fputs(fd, "Transfinite Surface(15) = {4, 3, 7, 8};\n");
%!     fputs(fd, "Transfinite Surface(16) = {2, 3, 7, 6};\n");
%!     fputs(fd, "tmp[] = Extrude {{1, 0, 0}, {0, -R, 0}, Pi / 2} { Surface{200,201,202,203}; Layers{Round(0.5 * R * Pi / h) + 1}; Recombine; };\n");
%!     fputs(fd, "Recombine Surface{200,tmp[0]};\n");
%!     fputs(fd, "Recombine Surface{201,tmp[6]};\n");
%!     fputs(fd, "Recombine Surface{202,tmp[12]};\n");
%!     fputs(fd, "Recombine Surface{203,tmp[18]};\n");
%!     fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1],tmp[7],tmp[13],tmp[19]};\n");
%!     fputs(fd, "Physical Surface(\"clamp\",1) = {200,201,202,203};\n");
%!     fputs(fd, "Physical Surface(\"load\",2) = {tmp[0],tmp[6],tmp[12],tmp[18]};\n");
%!     fputs(fd, "Physical Point(\"deformation\", 3) = {10};\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete=1;\n");
%!     fputs(fd, "Mesh.Algorithm = 3;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   # spawn_wait(spawn("gmsh", {[filename, ".geo"]}));return;
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   [~] = unlink([filename, ".geo"]);
%!   opt_msh.elem_type = {"iso20r", "quad8r", "point1"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_msh));
%!   [~] = unlink([filename, ".msh"]);
%!   grp_idx_volume = find([mesh.groups.iso20r.id] == 1);
%!   grp_idx_clamp = find([mesh.groups.quad8r.id] == 1);
%!   grp_idx_load = find([mesh.groups.quad8r.id] == 2);
%!   grp_idx_rbe3 = find([mesh.groups.point1.id] == 3);
%!   node_idx_rbe3 = mesh.groups.point1(grp_idx_rbe3).nodes;
%!   mesh.material_data.E = E;
%!   mesh.material_data.nu = nu;
%!   mesh.material_data.rho = rho;
%!   mesh.materials.iso20r = zeros(rows(mesh.elements.iso20r), 1, "int32");
%!   mesh.materials.iso20r(mesh.groups.iso20r(grp_idx_volume).elements) = 1;
%!   mesh.elements.rbe3 = fem_pre_mesh_rbe3_from_surf(mesh, 2, node_idx_rbe3, "quad8r");
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   load_case_dof.locked_dof(mesh.groups.quad8r(grp_idx_clamp).nodes, 1:3) = true;
%!   load_case.loads = [0, 0, -F, M, 0, 0];
%!   load_case.loaded_nodes = int32(node_idx_rbe3);
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!   [mat_ass.M, ...
%!    mat_ass.K, ...
%!    mat_ass.R, ...
%!    mat_ass.mat_info, ...
%!    mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                        dof_map, ...
%!                                        [FEM_MAT_MASS, ...
%!                                         FEM_MAT_STIFFNESS, ...
%!                                         FEM_VEC_LOAD_CONSISTENT], ...
%!                                        load_case);
%!   sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%!   sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_VEC_STRESS_CAUCH], ...
%!                                    load_case, ...
%!                                    sol_stat);
%!   sol_eig = fem_sol_modal(mesh, dof_map, mat_ass, N);
%!   sol_eig.stress = fem_ass_matrix(mesh, ...
%!                                   dof_map, ...
%!                                   [FEM_SCA_STRESS_VMIS], ...
%!                                   load_case, ...
%!                                   sol_eig);
%!   uref = F * R^3 / (2 * E * Iy); ## W. Beitz, K.-H. Grote, 1997, Dubbel page C25, figure 34b
%!   tauref = F * R / Iy * da / 2;
%!   u = -sol_stat.def(node_idx_rbe3, 2);
%!   tau = zeros(rows(mesh.nodes), 6);
%!   for i=1:6
%!     tau(mesh.elements.iso20r(:), i) = sol_stat.stress.taum.iso20r(:, :, i)(:);
%!   endfor
%!   taumax = max(tau(mesh.groups.quad8r(grp_idx_clamp).nodes, 3));
%!   tol = 0.5e-2;
%!   assert_simple(u, uref, tol * abs(uref));
%!   assert_simple(taumax, tauref, 0.2 * abs(tauref));
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
