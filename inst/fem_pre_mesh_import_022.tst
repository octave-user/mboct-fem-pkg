## fem_pre_mesh_import.m:22
%!test
%! try
%! ## TEST 22: static patch test of iso20
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
%!     a = 8e-3;
%!     b = 15e-3;
%!     c = 12e-3;
%!     px = 25.79e6;
%!     py = 7.83e6;
%!     pz = 1.3758e6;
%!     mesh_size = 7e-3;
%!     scale_def = 10e-3;
%!     fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "a = %g;\n", a);
%!     fprintf(fd, "b = %g;\n", b);
%!     fprintf(fd, "c = %g;\n", c);
%!     fprintf(fd, "h = %g;\n", mesh_size);
%!     fputs(fd, "Point(1) = {0,0,0,h};\n");
%!     fputs(fd, "Point(2) = {a,0,0,h};\n");
%!     fputs(fd, "Point(3) = {a,b,0,h};\n");
%!     fputs(fd, "Point(4) = {0,b,0,h};\n");
%!     fputs(fd, "Line(1) = {1,2};\n");
%!     fputs(fd, "Line(2) = {2,3};\n");
%!     fputs(fd, "Line(3) = {3,4};\n");
%!     fputs(fd, "Line(4) = {4,1};\n");
%!     fputs(fd, "Line Loop(5) = {1,2,3,4};\n");
%!     fputs(fd, "Plane Surface(6) = {5};\n");
%!     fputs(fd, "tmp[] = Extrude {0, 0, c}{ Surface{6}; Layers{Ceil(c / h)}; Recombine; };\n");
%!     fputs(fd, "Recombine Surface{6, tmp[0]};\n");
%!     fputs(fd, "Physical Volume(\"volume\", 1) = {tmp[1]};\n");
%!     fputs(fd, "Physical Surface(\"load-x\", 2) = {tmp[3],tmp[5]};\n");
%!     fputs(fd, "Physical Surface(\"load-y\", 3) = {tmp[2],tmp[4]};\n");
%!     fputs(fd, "Physical Surface(\"load-z\", 4) = {6,tmp[0]};\n");
%!     fputs(fd, "Physical Point(\"A\",1) = {1};\n");
%!     fputs(fd, "Physical Point(\"B\",2) = {2};\n");
%!     fputs(fd, "Physical Point(\"C\",3) = {4};\n");
%!     fputs(fd, "ReorientMesh Volume{tmp[1]};\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete = 1;\n");
%!     fputs(fd, "Mesh.ElementOrder = 2;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   [~] = unlink([filename, ".msh"]);
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   [~] = unlink([filename, ".geo"]);
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!   [~] = unlink([filename, ".msh"]);
%!   load_case.locked_dof = false(rows(mesh.nodes), 6);
%!   grp_id_A = find([[mesh.groups.point1].id] == 1);
%!   grp_id_B = find([[mesh.groups.point1].id] == 2);
%!   grp_id_C = find([[mesh.groups.point1].id] == 3);
%!   mesh.elements.joints(1).nodes = mesh.groups.point1(grp_id_A).nodes;
%!   mesh.elements.joints(1).C = eye(6)(1:3, :);
%!   mesh.elements.joints(2).nodes = mesh.groups.point1(grp_id_B).nodes;
%!   mesh.elements.joints(2).C = eye(6)([2, 3], :);
%!   mesh.elements.joints(3).nodes = mesh.groups.point1(grp_id_C).nodes;
%!   mesh.elements.joints(3).C = eye(6)(3, :);
%!   grp_id_px = find([[mesh.groups.quad8].id] == 2);
%!   grp_id_py = find([[mesh.groups.quad8].id] == 3);
%!   grp_id_pz = find([[mesh.groups.quad8].id] == 4);
%!   elem_id_px = mesh.groups.quad8(grp_id_px).elements;
%!   elem_id_py = mesh.groups.quad8(grp_id_py).elements;
%!   elem_id_pz = mesh.groups.quad8(grp_id_pz).elements;
%!   elno_px = mesh.elements.quad8(elem_id_px, :);
%!   elno_py = mesh.elements.quad8(elem_id_py, :);
%!   elno_pz = mesh.elements.quad8(elem_id_pz, :);
%!   load_case.pressure.quad8.elements = [elno_px; elno_py; elno_pz];
%!   load_case.pressure.quad8.p = [repmat(px, rows(elno_px), columns(elno_px));
%!                                 repmat(py, rows(elno_py), columns(elno_py));
%!                                 repmat(pz, rows(elno_pz), columns(elno_pz))];
%!   mesh.materials.iso20 = ones(rows(mesh.elements.iso20), 1, "int32");
%!   E = 210000e6;
%!   nu = 0.3;
%!   mesh.material_data.rho = 7850;
%!   mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.K, ...
%!    mat_ass.R, ...
%!    mat_ass.dm, ...
%!    mat_ass.S, ...
%!    mat_ass.J, ...
%!    mat_ass.mat_info, ...
%!    mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                       dof_map, ...
%!                                       [FEM_MAT_STIFFNESS, ...
%!                                        FEM_VEC_LOAD_CONSISTENT, ...
%!                                       FEM_SCA_TOT_MASS, ...
%!                                       FEM_VEC_INERTIA_M1, ...
%!                                       FEM_MAT_INERTIA_J], ...
%!                                       load_case);
%!   sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%!   F = zeros(rows(mesh.nodes), 3);
%!   for i=1:3
%!     idof = dof_map.ndof(:, i);
%!     idx = idof > 0;
%!     F(idx, i) = mat_ass.R(idof(idx));
%!   endfor
%!   [sol_stat.stress, ...
%!    sol_stat.strain] = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_VEC_STRESS_CAUCH, ...
%!                                       FEM_VEC_STRAIN_TOTAL], ...
%!                                      load_case, ...
%!                                      sol_stat);
%!   if (do_plot)
%!     figure("visible", "off");
%!     fem_post_sol_plot(mesh, sol_stat, scale_def / max(norm(sol_stat.def(:, 1:3), "rows")));
%!     xlabel("x [m]");
%!     ylabel("y [m]");
%!     zlabel("z [m]");
%!     grid on;
%!     grid minor on;
%!     title("static deflection - consistent pressure load");
%!     figure_list();
%!   endif
%!   tol = eps^0.7;
%!   Fref = [-px * b * c;
%!           -py * a * c;
%!           -pz * a * b];
%!   p = [px, py, pz];
%!   Xcg = mat_ass.S / mat_ass.dm;
%!   Jcg = mat_ass.J + skew(Xcg) * skew(Xcg) * mat_ass.dm;
%!   dmref = mesh.material_data.rho * a * b * c;
%!   Xcgref = 0.5 * [a; b; c];
%!   Jxxref = dmref * (b^2 + c^2) / 12;
%!   Jyyref = dmref * (a^2 + c^2) / 12;
%!   Jzzref = dmref * (a^2 + b^2) / 12;
%!   Jcgref = diag([Jxxref, Jyyref, Jzzref]);
%!   assert_simple(mat_ass.dm, dmref, tol * dmref);
%!   assert_simple(Xcg, Xcgref, tol * norm(Xcgref));
%!   assert_simple(Jcg, Jcgref, tol * norm(Jcgref));
%!   epsilon = mesh.material_data.C \ [-px; -py; -pz; zeros(3, 1)];
%!   for i=1:3
%!     assert_simple(sum(F(:, i)), 0, tol * abs(Fref(i)));
%!   endfor
%!   for i=1:3
%!     assert(sol_stat.def(:, i), mesh.nodes(:, i) * epsilon(i), tol * a * abs(epsilon(i)));
%!   endfor
%!   for i=1:3
%!     assert_simple(max(max(max(abs(sol_stat.stress.tau.iso20(:, :, i) / -p(i) - 1)))) < tol);
%!   endfor
%!   assert_simple(max(max(max(abs(sol_stat.stress.tau.iso20(:, :, 4:6) / max(abs(p)))))) < tol);
%!   for i=1:3
%!     assert_simple(max(max(max(abs(sol_stat.strain.epsilon.iso20(:, :, i) / epsilon(i) - 1)))) < tol);
%!   endfor
%!   assert_simple(max(max(max(abs(sol_stat.strain.epsilon.iso20(:, :, 4:6) / max(abs(epsilon)))))) < tol);
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
