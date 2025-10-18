## fem_pre_mesh_import.m:100
%!test
%! try
%! l = 2500e-3;
%! w = 20e-3;
%! h = 150e-3;
%! zf = 0e-3;
%! sx = h / 4;
%! sy = w / 4;
%! sz = h / 4;
%! E = 210000e6;
%! nu = 0.3;
%! G = E / (2 * (1 + nu));
%! rho = 7850;
%! elem_type = "rbe2";
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
%!     fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "l = %g;\n", l);
%!     fprintf(fd, "w = %g;\n", w);
%!     fprintf(fd, "h = %g;\n", h);
%!     fprintf(fd, "zf = %g;\n", zf);
%!     fprintf(fd, "sx = %g;\n", sx);
%!     fprintf(fd, "sy = %g;\n", sy);
%!     fprintf(fd, "sz = %g;\n", sz);
%!     fputs(fd, "Mesh.SecondOrderIncomplete=1;\n");
%!     fputs(fd, "Point(1) = {0,  0.5 * w,  0.5 * h};\n");
%!     fputs(fd, "Point(2) = {0, -0.5 * w,  0.5 * h};\n");
%!     fputs(fd, "Point(3) = {0, -0.5 * w, -0.5 * h};\n");
%!     fputs(fd, "Point(4) = {0,  0.5 * w, -0.5 * h};\n");
%!     fputs(fd, "Point(5) = {l,        0,       zf};\n");
%!     fputs(fd, "Line(1) = {1,2};\n");
%!     fputs(fd, "Line(2) = {2,3};\n");
%!     fputs(fd, "Line(3) = {3,4};\n");
%!     fputs(fd, "Line(4) = {4,1};\n");
%!     fputs(fd, "Transfinite Curve(1) = Round(w / sy) + 1;\n");
%!     fputs(fd, "Transfinite Curve(2) = Round(h / sz) + 1;\n");
%!     fputs(fd, "Transfinite Curve(3) = Round(w / sy) + 1;\n");
%!     fputs(fd, "Transfinite Curve(4) = Round(h / sz) + 1;\n");
%!     fputs(fd, "Line Loop(5) = {1,2,3,4};\n");
%!     fputs(fd, "Plane Surface(6) = {5};\n");
%!     fputs(fd, "Transfinite Surface(6) = {1, 2, 3, 4};\n");
%!     fputs(fd, "tmp[] = Extrude {l, 0, 0}{ Surface{6}; Layers{Ceil(l / sx)}; Recombine; };\n");
%!     fputs(fd, "Recombine Surface{6, tmp[0]};\n");
%!     fputs(fd, "Physical Volume(\"volume\", 1) = {tmp[1]};\n");
%!     fputs(fd, "Physical Surface(\"clamp\", 2) = {6};\n");
%!     fputs(fd, "Physical Surface(\"rbe2\", 3) = {tmp[0]};\n");
%!     fputs(fd, "Physical Point(\"load\", 4) = 5;\n");
%!     fputs(fd, "Mesh.ElementOrder = 2;\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete=1;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   [~] = unlink([filename, ".msh"]);
%!   #spawn_wait(spawn("gmsh", {[filename, ".geo"]}));
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "1", [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   [~] = unlink([filename, ".geo"]);
%!   opt_msh.elem_type = {"iso20r", "quad8r", "point1"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_msh));
%!   [~] = unlink([filename, ".msh"]);
%!   grp_idx_volume = find([mesh.groups.iso20r.id] == 1);
%!   grp_idx_clamp = find([mesh.groups.quad8r.id] == 2);
%!   grp_idx_rbe2 = find([mesh.groups.quad8r.id] == 3);
%!   grp_idx_load = find([mesh.groups.point1.id] == 4);
%!   mesh.materials.iso20r = zeros(rows(mesh.elements.iso20r), 1, "int32");
%!   mesh.materials.iso20r(mesh.groups.iso20r(grp_idx_volume).elements) = 1;
%!   mesh.material_data.E = E;
%!   mesh.material_data.nu = nu;
%!   mesh.material_data.rho = rho;
%!   load_case_dof.locked_dof = zeros(size(mesh.nodes));
%!   load_case_dof.locked_dof(mesh.groups.quad8r(grp_idx_clamp).nodes, :) = true;
%!   node_idx_rbe2 = mesh.groups.quad8r(grp_idx_rbe2).nodes;
%!   node_idx_load = mesh.groups.point1(grp_idx_load).nodes;
%!   switch (elem_type)
%!   case "rbe2"
%!     mesh.elements.rbe2 = struct("nodes", cell(1, numel(node_idx_rbe2)));
%!     for i=1:numel(node_idx_rbe2)
%!       mesh.elements.rbe2(i).nodes = [node_idx_load, node_idx_rbe2(i)];
%!     endfor
%!   case "rbe3"
%!     mesh.elements.rbe3 = fem_pre_mesh_rbe3_from_surf(mesh, 3, node_idx_load, "quad8r");
%!   endswitch
%!   load_case.loads = [0, 0, -1, 0, 0, 0];
%!   load_case.loaded_nodes = node_idx_load;
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!   [mat_ass.K, ...
%!    mat_ass.R, ...
%!    mat_ass.mat_info, ...
%!    mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                        dof_map, ...
%!                                        [FEM_MAT_STIFFNESS, ...
%!                                         FEM_VEC_LOAD_CONSISTENT], ...
%!                                        load_case);
%!   opt_sol.solver = "pardiso";
%!   opt_sol.refine_max_iter = 100;
%!   sol_stat = fem_sol_static(mesh, dof_map, mat_ass, opt_sol);
%!   sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_VEC_STRESS_CAUCH], ...
%!                                    load_case, ...
%!                                    sol_stat);
%!   load_case.lambda = sol_stat.lambda;
%!   load_case.tau0 = sol_stat.stress.tau;
%!   [mat_ass.KTAU0] = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_MAT_STIFFNESS_TAU0], ...
%!                                    load_case);
%!   opt_sol.problem = "buckling";
%!   opt_sol.p = 20;
%!   number_of_modes = 2;
%!   [U, sol_eig.lambda, err] = fem_sol_eigs(mat_ass.K, mat_ass.KTAU0, number_of_modes, opt_sol);
%!   sol_eig.def = fem_post_def_nodal(mesh, dof_map, U);
%!   sol_eig.stress = fem_ass_matrix(mesh, ...
%!                                   dof_map, ...
%!                                   [FEM_VEC_STRESS_CAUCH], ...
%!                                   load_case, ...
%!                                   sol_eig);
%!   Iy = w * h^3 / 12;
%!   Iz = w^3 * h / 12;
%!   c1 = interp1([1, 1.5, 2, 3, 4, 6, 8, 10, realmax], [0.141, 0.196, 0.229, 0.263, 0.281, 0.298, 0.307, 0.312, 0.333], h / w);
%!   It = c1 * max([h, w]) * min([h, w])^3;
%!   K = sqrt(E * Iz * G * It * (Iy - Iz) / Iy);
%!   FKref = 4.013 / l^2 * K * (1 - zf / l * sqrt(E * Iz / (G * It)));
%!   tol = 4.5e-2;
%!   assert(sol_eig.lambda(2), FKref, tol * FKref);
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
