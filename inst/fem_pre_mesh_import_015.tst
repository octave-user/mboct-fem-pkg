## fem_pre_mesh_import.m:15
%!test
%! try
%! ### TEST 15
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
%!   unwind_protect
%!     [fd, msg] = fopen([filename, ".geo"], "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s.geo\"", filename);
%!     endif
%!     h1 = 10e-3;
%!     h2 = 2e-3;
%!     w = 4e-3;
%!     l = 30e-3;
%!     h = 2e-3;
%!     do_post_pro = false;
%!     scale_def = 20e-3;
%!     fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "h1=%g;\n", h1);
%!     fprintf(fd, "h2=%g;\n", h2)
%!     fprintf(fd, "w=%g;\n", w);
%!     fprintf(fd, "l=%g;\n", l);
%!     fprintf(fd, "h=%g;\n", h);
%!     fputs(fd, "Point(1) = {0.0, -0.5 * w, -0.5 * h1 - h2, h};\n");
%!     fputs(fd, "Point(2) = {0.0, -0.5 * w, -0.5 * h1, h};\n");
%!     fputs(fd, "Point(3) = {0.0,  0.5 * w, -0.5 * h1, h};\n");
%!     fputs(fd, "Point(4) = {0.0,  0.5 * w, -0.5 * h1 - h2, h};\n");
%!     fputs(fd, "Point(5) = {0.0, -0.5 * w, 0.5 * h1, h};\n");
%!     fputs(fd, "Point(6) = {0.0, -0.5 * w, 0.5 * h1 + h2, h};\n");
%!     fputs(fd, "Point(7) = {0.0,  0.5 * w, 0.5 * h1 + h2, h};\n");
%!     fputs(fd, "Point(8) = {0.0,  0.5 * w, 0.5 * h1, h};\n");
%!     fputs(fd, "Line(1) = {1,2};\n");
%!     fputs(fd, "Line(2) = {2,3};\n");
%!     fputs(fd, "Line(3) = {3,4};\n");
%!     fputs(fd, "Line(4) = {4,1};\n");
%!     fputs(fd, "Line(5) = {2,5};\n");
%!     fputs(fd, "Line(6) = {5,8};\n");
%!     fputs(fd, "Line(7) = {8,3};\n");
%!     fputs(fd, "Line(8) = {3,2};\n");
%!     fputs(fd, "Line(9) = {5,6};\n");
%!     fputs(fd, "Line(10) = {6,7};\n");
%!     fputs(fd, "Line(11) = {7,8};\n");
%!     fputs(fd, "Line(12) = {8,5};\n");
%!     fputs(fd, "Line Loop(13) = {1,2,3,4};\n");
%!     fputs(fd, "Line Loop(14) = {5,6,7,8};\n");
%!     fputs(fd, "Line Loop(15) = {9,10,11,12};\n");
%!     fputs(fd, "Plane Surface(16) = {13};\n");
%!     fputs(fd, "Plane Surface(17) = {14};\n");
%!     fputs(fd, "Plane Surface(18) = {15};\n");
%!     fputs(fd, "v1[] = Extrude {l,0.0,0.0} {\n");
%!     fputs(fd, "  Surface{16};\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "v2[] = Extrude {l,0.0,0.0} {\n");
%!     fputs(fd, "  Surface{17};\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "v3[] = Extrude {l,0.0,0.0} {\n");
%!     fputs(fd, "  Surface{18};\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "v = newv;\n");
%!     fputs(fd, "BooleanUnion(v) = {Volume{v1[1]};}{ Volume{v2[1],v3[1]};};\n");
%!     fputs(fd, "Coherence;\n");
%!     fputs(fd, "Physical Volume(\"volume1\",1) = {v1[1]};\n");
%!     fputs(fd, "Physical Volume(\"volume2\",2) = {v2[1]};\n");
%!     fputs(fd, "Physical Volume(\"volume3\",3) = {v3[1]};\n");
%!     fputs(fd, "Physical Surface(\"surface1\",1) = {v1[0]};\n");
%!     fputs(fd, "Physical Surface(\"surface2\",2) = {v2[0]};\n");
%!     fputs(fd, "Physical Surface(\"surface3\",3) = {v3[0]};\n");
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
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!   E(1) = 80000e6;
%!   nu(1) = 0.4;
%!   mesh.material_data(1).rho = 1700;
%!   mesh.material_data(1).alpha = 1e-8;
%!   mesh.material_data(1).beta = 1e-6;
%!   E(2) = 10000e6;
%!   nu(2) = 0.4;
%!   mesh.material_data(2).rho = 500;
%!   mesh.material_data(2).alpha = 0;
%!   mesh.material_data(2).beta = 0;
%!   for i=1:numel(mesh.material_data)
%!     mesh.material_data(i).C = fem_pre_mat_isotropic(E(i), nu(i));
%!   endfor
%!   mesh.materials.tet10 = zeros(rows(mesh.elements.tet10),1);
%!   mesh.materials.tet10(mesh.groups.tet10(1).elements) = 1;
%!   mesh.materials.tet10(mesh.groups.tet10(2).elements) = 2;
%!   mesh.materials.tet10(mesh.groups.tet10(3).elements) = 1;
%!   load_case.locked_dof = false(size(mesh.nodes));
%!   load_case.locked_dof(find(mesh.nodes(:, 1)==0), 1:3) = true;
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.M, ...
%!    mat_ass.D, ...
%!    mat_ass.K, ...
%!    mat_ass.dm, ...
%!    mat_ass.mat_info] = fem_ass_matrix(mesh, ...
%!                                       dof_map, ...
%!                                       [FEM_MAT_MASS, ...
%!                                        FEM_MAT_DAMPING, ...
%!                                        FEM_MAT_STIFFNESS, ...
%!                                        FEM_SCA_TOT_MASS], ...
%!                                       load_case);
%!   [sol_eig, ~, err] = fem_sol_modal(mesh, dof_map, mat_ass, 10);
%!   [sol_eig.stress] = fem_ass_matrix(mesh, ...
%!                                     dof_map, ...
%!                                     [FEM_VEC_STRESS_CAUCH], ...
%!                                     load_case, ...
%!                                     sol_eig);
%!   opts.scale_def = scale_def / max(max(max(abs(sol_eig.def))));
%!   rho1 = mesh.material_data(1).rho;
%!   rho2 = mesh.material_data(2).rho;
%!   m1 = rho1 * h2 * w * l;
%!   m2 = rho2 * h1 * w * l;
%!   m = 2 * m1 + m2;
%!   assert_simple(mat_ass.dm, m, sqrt(eps) * m);
%!   opts.skin_only = true;
%!   if (do_post_pro)
%!     fem_post_sol_external(mesh, sol_eig, opts);
%!   endif
%! unwind_protect_cleanup
%!   [~] = unlink([filename, ".msh"]);
%!   [~] = unlink([filename, ".geo"]);
%! end_unwind_protect
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
