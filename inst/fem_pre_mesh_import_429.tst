## fem_pre_mesh_import.m:214
%!test
%! try
%! ## TEST 214
%! ## M. Maeder, R. D'Auria, E. Grasso, G. Petrone b, S. De Rosa, M. Klaerner, L. Kroll, S. Marburg
%! ## Numerical analysis of sound radiation from rotating discs
%! ## Journal of Sound and Vibration 468 (2020) 115085
%! close all;
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
%!     unit_meters = 1e-3;
%!     unit_second = 1e4;
%!     unit_kilograms = 1e-3;
%!     unit_newton = unit_kilograms * unit_meters / unit_second^2;
%!     unit_pascal = unit_newton / unit_meters^2;
%!     unit_watt = unit_newton * unit_meters / unit_second;
%!     P0 = 1e-12 / unit_watt;
%!     d1 = 700e-3 / unit_meters;
%!     d2 = 120e-3 / unit_meters;
%!     d3 = 60e-3 / unit_meters;
%!     b = 200e-3 / unit_meters;
%!     t = 3.5e-3 / unit_meters;
%!     r = 1500e-3 / unit_meters;
%!     h1 = 30 * t;
%!     h2 = 400e-3 / unit_meters;
%!     h3 = 400e-3 / unit_meters;
%!     E1 = 210000e6 / unit_pascal;
%!     rho1 = 7800 / (unit_kilograms / unit_meters^3);
%!     nu1 = 0.3;
%!     alpha1 = 0.1826 / (unit_second^-1);
%!     beta1 = 5.0125e-6 / unit_second;
%!     c2 = 340 / (unit_meters / unit_second);
%!     rho2 = 1.225 / (unit_kilograms / unit_meters^3);
%!     Fz = (1 + 0j) / unit_newton;
%!     f = (100) / (unit_second^-1);
%!     km = 2 * pi * mean(f) / c2;
%!     deltaPML = h3;
%!     alphaPML = 1;
%!     nPML = 1;
%!     solver = "precond";
%!     f_enable_PML = true;
%!     f_enable_plot = false;
%!     smoothPML = 0;
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "d1 = %g;\n", d1);
%!     fprintf(fd, "d2 = %g;\n", d2);
%!     fprintf(fd, "d3 = %g;\n", d3);
%!     fprintf(fd, "b = %g;\n", b);
%!     fprintf(fd, "t = %g;\n", t);
%!     fprintf(fd, "r = %g;\n", r);
%!     fprintf(fd, "deltaPML = %g;\n", deltaPML);
%!     fprintf(fd, "h1 = %g;\n", h1);
%!     fprintf(fd, "h2 = %g;\n", h2);
%!     fprintf(fd, "h3 = %g;\n", h3);
%!     fprintf(fd, "nPML = %d;\n", ceil(nPML));
%!     fputs(fd, "r1 = r + b;\n");
%!     fputs(fd, "v1 = newv;\n");
%!     fputs(fd, "Cylinder(v1) = {0, 0, -0.5 * t, 0, 0, t, 0.5 * d1};\n");
%!     fputs(fd, "v2 = newv;\n");
%!     fputs(fd, "Cylinder(v2) = {0, 0, -0.5 * t, 0, 0, t, 0.5 * d2};\n");
%!     fputs(fd, "v3 = newv;\n");
%!     fputs(fd, "Cylinder(v3) = {0, 0, -0.5 * t, 0, 0, t, 0.5 * d3};\n");
%!     fputs(fd, "v4 = newv;\n");
%!     fputs(fd, "Sphere(v4) = {0, 0, 0, r};\n");
%!     fputs(fd, "v5 = newv;\n");
%!     fputs(fd, "Box(v5) = {-r1, -r1, -r1, 2 * r1, 2 * r1, 2 * r1};\n");
%!     fputs(fd, "Extrude {deltaPML, 0, 0} {\n");
%!     fputs(fd, "  Surface{12}; Layers {nPML}; Recombine;\n");
%!     fputs(fd, "}\n");
%!     fputs(fd, "Extrude {0, deltaPML, 0} {\n");
%!     fputs(fd, "  Surface{14}; Layers {nPML}; Recombine;\n");
%!     fputs(fd, "}\n");
%!     fputs(fd, "Extrude {0, 0, deltaPML} {\n");
%!     fputs(fd, "  Surface{16}; Layers {nPML}; Recombine;\n");
%!     fputs(fd, "}\n");
%!     fputs(fd, "Extrude {-deltaPML, 0, 0} {\n");
%!     fputs(fd, "  Surface{11}; Layers {nPML}; Recombine;\n");
%!     fputs(fd, "}\n");
%!     fputs(fd, "Extrude {0, -deltaPML, 0} {\n");
%!     fputs(fd, "  Surface{13}; Layers {nPML}; Recombine;\n");
%!     fputs(fd, "}\n");
%!     fputs(fd, "Extrude {0, 0, -deltaPML} {\n");
%!     fputs(fd, "  Surface{15}; Layers {nPML}; Recombine;\n");
%!     fputs(fd, "}\n");
%!     fputs(fd, "Extrude {0, 0, deltaPML} {\n");
%!     fputs(fd, "  Surface{18}; Layers {nPML}; Recombine;\n");
%!     fputs(fd, "}\n");
%!     fputs(fd, "Extrude {0, 0, -deltaPML} {\n");
%!     fputs(fd, "  Surface{20}; Layers {nPML}; Recombine;\n");
%!     fputs(fd, "}\n");
%!     fputs(fd, "Extrude {0, deltaPML, 0} {\n");
%!     fputs(fd, "  Surface{19}; Layers {nPML}; Recombine;\n");
%!     fputs(fd, "}\n");
%!     fputs(fd, "Extrude {0, -deltaPML, 0} {\n");
%!     fputs(fd, "  Surface{17}; Layers {nPML}; Recombine;\n");
%!     fputs(fd, "}\n");
%!     fputs(fd, "Extrude {0, 0, deltaPML} {\n");
%!     fputs(fd, "  Surface{33}; Layers {nPML}; Recombine;\n");
%!     fputs(fd, "}\n");
%!     fputs(fd, "Extrude {0, 0, -deltaPML} {\n");
%!     fputs(fd, "  Surface{35}; Layers {nPML}; Recombine;\n");
%!     fputs(fd, "}\n");
%!     fputs(fd, "Extrude {0, deltaPML, 0} {\n");
%!     fputs(fd, "  Surface{34}; Layers {nPML}; Recombine;\n");
%!     fputs(fd, "}\n");
%!     fputs(fd, "Extrude {0, -deltaPML, 0} {\n");
%!     fputs(fd, "  Surface{32}; Layers {nPML}; Recombine;\n");
%!     fputs(fd, "}\n");
%!     fputs(fd, "Extrude {0, 0, deltaPML} {\n");
%!     fputs(fd, "  Surface{57}; Surface{24}; Surface{77}; Surface{82}; Surface{39}; Surface{62}; Layers {nPML}; Recombine;\n");
%!     fputs(fd, "}\n");
%!     fputs(fd, "Extrude {0, 0, -deltaPML} {\n");
%!     fputs(fd, "  Surface{58}; Surface{22}; Surface{78}; Surface{83}; Surface{37}; Surface{63}; Layers {nPML}; Recombine;\n");
%!     fputs(fd, "}\n");
%!     fputs(fd, "Coherence;\n");
%!     fputs(fd, "Physical Volume(\"solid\", 1) = {40, 41};\n");
%!     fputs(fd, "Physical Volume(\"fluid\", 2) = {7, 42, 43};\n");
%!     fputs(fd, "Physical Volume(\"PML-x\", 3) = {14, 28, 20, 33, 22, 23, 34, 21, 39, 30, 24, 31, 26, 17, 27, 36, 25, 37};\n");
%!     fputs(fd, "Physical Volume(\"PML-y\", 4) = {30, 29, 28, 26, 15, 22, 36, 35, 34, 37, 38, 39, 27, 18, 23, 31, 32, 33};\n");
%!     fputs(fd, "Physical Volume(\"PML-z\", 5) = {30, 29, 28, 24, 16, 20, 31, 32, 33, 36, 35, 34, 25, 19, 21, 37, 38, 39};\n");
%!     fputs(fd, "Physical Surface(\"fluid-struct\", 1) = {38, 39, 41, 42, 7, 37};\n");
%!     fputs(fd, "Physical Surface(\"clamp\", 2) = {41, 42};\n");
%!     fputs(fd, "Physical Surface(\"fluid-boundary\", 3) = {43};\n");
%!     fputs(fd, "Physical Surface(\"PML-boundary\", 4) = {113, 87, 109, 106, 57, 103, 71, 119, 116, 108, 105, 102, 95, 79, 52, 125, 122, 128, 117, 101, 70, 82, 21, 78, 137, 74, 121, 99, 132, 112, 135, 138, 62, 115, 118, 83, 133, 136, 139, 75, 123, 126, 129, 91, 67, 127, 90, 131, 94, 107, 86, 111, 36, 98};\n");
%!     fputs(fd, "ReverseMesh Surface{7};\n");
%!     fputs(fd, "Mesh.HighOrderOptimize = 1;\n");
%!     fputs(fd, "MeshSize{PointsOf{Volume{14, 28, 20, 33, 22, 23, 34, 21, 39, 30, 24, 31, 26, 17, 27, 36, 25, 37, 30, 29, 28, 26, 15, 22, 36, 35, 34, 37, 38, 39, 27, 18, 23, 31, 32, 33, 30, 29, 28, 24, 16, 20, 31, 32, 33, 36, 35, 34, 25, 19, 21, 37, 38, 39}; } } = h3;\n");
%!     fputs(fd, "MeshSize{PointsOf{Volume{7, 42, 43}; } } = h2;\n");
%!     fputs(fd, "MeshSize{PointsOf{Volume{40, 41}; } } = h1;\n");
%!     fputs(fd, "Mesh.HighOrderIterMax = 100;\n");
%!     fputs(fd, "Mesh.HighOrderPassMax = 100;\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete = 1;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-ho_min", "0.5", "-ho_max", "1.5",  [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   [~] = unlink([filename, ".geo"]);
%!   opt_mesh.elem_type = {"iso20r", "penta15", "tet10", "tria6h", "quad8r"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_mesh));
%!   [~] = unlink([filename, ".msh"]);
%!   grp_idx_v_solid = find([mesh.groups.tet10.id] == 1);
%!   grp_idx_v_fluid = find([mesh.groups.tet10.id] == 2);
%!   grp_idx_v_p_PML_x = find([mesh.groups.penta15.id] == 3);
%!   grp_idx_v_p_PML_y = find([mesh.groups.penta15.id] == 4);
%!   grp_idx_v_p_PML_z = find([mesh.groups.penta15.id] == 5);
%!   grp_idx_v_h_PML_x = find([mesh.groups.iso20r.id] == 3);
%!   grp_idx_v_h_PML_y = find([mesh.groups.iso20r.id] == 4);
%!   grp_idx_v_h_PML_z = find([mesh.groups.iso20r.id] == 5);
%!   grp_idx_s_fsi = find([mesh.groups.tria6h.id] == 1);
%!   grp_idx_s_clamp = find([mesh.groups.tria6h.id] == 2);
%!   grp_idx_s_boundary = find([mesh.groups.tria6h.id] == 3);
%!   grp_idx_s_q_boundaryPML = find([mesh.groups.quad8r.id] == 4);
%!   grp_idx_s_t_boundaryPML = find([mesh.groups.tria6h.id] == 4);
%!   elem_idx_p_PML_x = mesh.groups.penta15(grp_idx_v_p_PML_x).elements;
%!   elem_idx_p_PML_y = mesh.groups.penta15(grp_idx_v_p_PML_y).elements;
%!   elem_idx_p_PML_z = mesh.groups.penta15(grp_idx_v_p_PML_z).elements;
%!   elem_idx_h_PML_x = mesh.groups.iso20r(grp_idx_v_h_PML_x).elements;
%!   elem_idx_h_PML_y = mesh.groups.iso20r(grp_idx_v_h_PML_y).elements;
%!   elem_idx_h_PML_z = mesh.groups.iso20r(grp_idx_v_h_PML_z).elements;
%!   node_idx_constr = [mesh.groups.quad8r(grp_idx_s_q_boundaryPML).nodes, ...
%!                      mesh.groups.tria6h(grp_idx_s_t_boundaryPML).nodes];
%!   empty_cell = cell(1, 2);
%!   mesh.material_data = struct("E", empty_cell, ...
%!                               "nu", empty_cell, ...
%!                               "rho", empty_cell, ...
%!                               "c", empty_cell, ...
%!                               "alpha", empty_cell, ...
%!                               "beta", empty_cell);
%!   mesh.material_data(1).E = E1;
%!   mesh.material_data(1).rho = rho1;
%!   mesh.material_data(1).nu = nu1;
%!   mesh.material_data(1).alpha = alpha1;
%!   mesh.material_data(1).beta = beta1;
%!   mesh.material_data(2).c = c2;
%!   mesh.material_data(2).rho = rho2;
%!   mesh.materials.tet10 = zeros(rows(mesh.elements.tet10), 1, "int32");
%!   mesh.materials.penta15 = zeros(rows(mesh.elements.penta15), 1, "int32");
%!   mesh.materials.iso20r = zeros(rows(mesh.elements.iso20r), 1, "int32");
%!   mesh.materials.tet10(mesh.groups.tet10(grp_idx_v_solid).elements) = 1;
%!   mesh.materials.tet10(mesh.groups.tet10(grp_idx_v_fluid).elements) = 2;
%!   mesh.materials.penta15(mesh.groups.penta15(grp_idx_v_p_PML_x).elements) = 2;
%!   mesh.materials.penta15(mesh.groups.penta15(grp_idx_v_p_PML_y).elements) = 2;
%!   mesh.materials.penta15(mesh.groups.penta15(grp_idx_v_p_PML_z).elements) = 2;
%!   mesh.materials.iso20r(mesh.groups.iso20r(grp_idx_v_h_PML_x).elements) = 2;
%!   mesh.materials.iso20r(mesh.groups.iso20r(grp_idx_v_h_PML_y).elements) = 2;
%!   mesh.materials.iso20r(mesh.groups.iso20r(grp_idx_v_h_PML_z).elements) = 2;
%!   mesh.elements.fluid_struct_interface.tria6h = mesh.elements.tria6h(mesh.groups.tria6h(grp_idx_s_fsi).elements, :);
%!   mesh.elements.acoustic_boundary.tria6h = mesh.elements.tria6h(mesh.groups.tria6h(grp_idx_s_boundary).elements, :);
%!   mesh.materials.acoustic_boundary.tria6h = repmat(int32(2), rows(mesh.elements.acoustic_boundary.tria6h), 1);
%!   load_case_dof.locked_dof = false(rows(mesh.nodes), 7);
%!   load_case_dof.locked_dof(mesh.groups.tria6h(grp_idx_s_clamp).nodes, 1:3) = true;
%!   if (f_enable_PML)
%!     load_case_dof.locked_dof(node_idx_constr, 7) = true;
%!   endif
%!   load_case_dof.domain = FEM_DO_FLUID_STRUCT;
%!   empty_cell = cell(1, 2);
%!   load_case = struct("loads", empty_cell, "loaded_nodes", empty_cell);
%!   tol = sqrt(eps) * d1;
%!   node_id_load = find((abs(mesh.nodes(:, 1) - 0.5 * d1) < tol) & ...
%!                       (abs(mesh.nodes(:, 2)) < tol) & ...
%!                       (abs(mesh.nodes(:, 3) - 0.5 * t) < tol));
%!   node_id_load = node_id_load(1);
%!   load_case(1).loaded_nodes = node_id_load;
%!   load_case(1).loads = [zeros(1, 2), real(Fz), zeros(1, 3)];
%!   load_case(2).loaded_nodes = node_id_load;
%!   load_case(2).loads = [zeros(1, 2), imag(Fz), zeros(1, 3)];
%!
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!
%!   if (f_enable_PML)
%!     mesh.elements.perfectly_matched_layers.penta15.f = zeros(3, columns(mesh.elements.penta15), rows(mesh.elements.penta15));
%!     mesh.elements.perfectly_matched_layers.penta15.e1 = zeros(3, columns(mesh.elements.penta15), rows(mesh.elements.penta15));
%!     mesh.elements.perfectly_matched_layers.penta15.e2 = zeros(3, columns(mesh.elements.penta15), rows(mesh.elements.penta15));
%!     mesh.elements.perfectly_matched_layers.penta15.e1(1, :, :) = 1;
%!     mesh.elements.perfectly_matched_layers.penta15.e2(2, :, :) = 1;
%!     mesh.elements.perfectly_matched_layers.iso20r.f = zeros(3, columns(mesh.elements.iso20r), rows(mesh.elements.iso20r));
%!     mesh.elements.perfectly_matched_layers.iso20r.e1 = zeros(3, columns(mesh.elements.iso20r), rows(mesh.elements.iso20r));
%!     mesh.elements.perfectly_matched_layers.iso20r.e2 = zeros(3, columns(mesh.elements.iso20r), rows(mesh.elements.iso20r));
%!     mesh.elements.perfectly_matched_layers.iso20r.e1(1, :, :) = 1;
%!     mesh.elements.perfectly_matched_layers.iso20r.e2(2, :, :) = 1;
%!     for i=1:3
%!       mesh.elements.perfectly_matched_layers.penta15.f(i, :, :) = 1;
%!       mesh.elements.perfectly_matched_layers.iso20r.f(i, :, :) = 1;
%!     endfor
%!     xi_p = abs(mesh.nodes(:, 1)(mesh.elements.penta15(elem_idx_p_PML_x, :).')) - (r + b);
%!     yi_p = abs(mesh.nodes(:, 2)(mesh.elements.penta15(elem_idx_p_PML_y, :).')) - (r + b);
%!     zi_p = abs(mesh.nodes(:, 3)(mesh.elements.penta15(elem_idx_p_PML_z, :).')) - (r + b);
%!     xi_h = abs(mesh.nodes(:, 1)(mesh.elements.iso20r(elem_idx_h_PML_x, :).')) - (r + b);
%!     yi_h = abs(mesh.nodes(:, 2)(mesh.elements.iso20r(elem_idx_h_PML_y, :).')) - (r + b);
%!     zi_h = abs(mesh.nodes(:, 3)(mesh.elements.iso20r(elem_idx_h_PML_z, :).')) - (r + b);
%!     sigmax_p_k = alphaPML ./ (deltaPML - xi_p) - smoothPML * alphaPML / deltaPML;
%!     sigmay_p_k = alphaPML ./ (deltaPML - yi_p) - smoothPML * alphaPML / deltaPML;
%!     sigmaz_p_k = alphaPML ./ (deltaPML - zi_p) - smoothPML * alphaPML / deltaPML;
%!     sigmax_h_k = alphaPML ./ (deltaPML - xi_h) - smoothPML * alphaPML / deltaPML;
%!     sigmay_h_k = alphaPML ./ (deltaPML - yi_h) - smoothPML * alphaPML / deltaPML;
%!     sigmaz_h_k = alphaPML ./ (deltaPML - zi_h) - smoothPML * alphaPML / deltaPML;
%!   endif
%!   P = T = zeros(1, numel(f));
%!   ITER = repmat(intmax(), 1, 2);
%!   MAXITER = 20;
%!   TOL = eps^0.5;
%!   RESTART = 10;
%!   MAXITERS = int32(100);
%!   TOLS = eps^0.5;
%!   perm = [];
%!   Phi = zeros(dof_map.totdof, 1);
%!   cpu_factor = 0;
%!   cpu_solve = 0;
%!   FLAG = -1;
%!   alpha = 0.75;
%!   unwind_protect
%!     for i=1:numel(f)
%!       omega = 2 * pi * f(i);
%!       k = omega / c2;
%!       if (f_enable_PML)
%!         mesh.elements.perfectly_matched_layers.penta15.f(1, :, elem_idx_p_PML_x) = 1 ./ (1 - 1j * sigmax_p_k / k);
%!         mesh.elements.perfectly_matched_layers.penta15.f(2, :, elem_idx_p_PML_y) = 1 ./ (1 - 1j * sigmay_p_k / k);
%!         mesh.elements.perfectly_matched_layers.penta15.f(3, :, elem_idx_p_PML_z) = 1 ./ (1 - 1j * sigmaz_p_k / k);
%!         mesh.elements.perfectly_matched_layers.iso20r.f(1, :, elem_idx_h_PML_x) = 1 ./ (1 - 1j * sigmax_h_k / k);
%!         mesh.elements.perfectly_matched_layers.iso20r.f(2, :, elem_idx_h_PML_y) = 1 ./ (1 - 1j * sigmay_h_k / k);
%!         mesh.elements.perfectly_matched_layers.iso20r.f(3, :, elem_idx_h_PML_z) = 1 ./ (1 - 1j * sigmaz_h_k / k);
%!       endif
%!       [mat_ass.Mfs_re, ...
%!        mat_ass.Mfs_im, ...
%!        mat_ass.Kfs_re, ...
%!        mat_ass.Kfs_im, ...
%!        mat_ass.Dfs_re, ...
%!        mat_ass.Dfs_im, ...
%!        mat_ass.Rfs, ...
%!        mat_ass.mat_info, ...
%!        mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                         dof_map, ...
%!                                         [FEM_MAT_MASS_FLUID_STRUCT_RE, ...
%!                                          FEM_MAT_MASS_FLUID_STRUCT_IM, ...
%!                                          FEM_MAT_STIFFNESS_FLUID_STRUCT_RE, ...
%!                                          FEM_MAT_STIFFNESS_FLUID_STRUCT_IM, ...
%!                                          FEM_MAT_DAMPING_FLUID_STRUCT_RE, ...
%!                                          FEM_MAT_DAMPING_FLUID_STRUCT_IM, ...
%!                                          FEM_VEC_LOAD_FLUID_STRUCT], ...
%!                                         load_case);
%!       Keff = -omega^2 * complex(mat_ass.Mfs_re, mat_ass.Mfs_im) ...
%!            + 1j * omega * complex(mat_ass.Dfs_re, mat_ass.Dfs_im) ...
%!            + complex(mat_ass.Kfs_re, mat_ass.Kfs_im);
%!       Reff = complex(mat_ass.Rfs(:, 1), mat_ass.Rfs(:, 2));
%!       opt_sol.number_of_threads = mbdyn_solver_num_threads_default();
%!       opt_sol.solver = "pastix";
%!       opt_sol.compress_when = int32(0);
%!       opt_sol.compress_min_ratio = 1;
%!       opt_sol.compress_tolerance = 1e-2;
%!       opt_sol.verbose = int32(0);
%!       switch (solver)
%!         case "direct"
%!           opt_sol.refine_max_iter = MAXITER;
%!           opt_sol.epsilon_refinement = TOL;
%!           opt_sol.pre_scaling = true;
%!           Phi = fem_sol_factor(Keff, opt_sol) \ Reff;
%!         case "precond"
%!           opt_sol.refine_max_iter = int32(0);
%!           opt_sol.pre_scaling = false;
%!           if (isempty(perm))
%!             perm = csymamd(Keff);
%!           endif
%!           do
%!             do_fact = (FLAG ~= 0 || cpu_solve > alpha * cpu_factor);
%!             if (do_fact)
%!               start = cputime();
%!               [KS, D1, D2] = fem_sol_matrix_scale(Keff(perm, perm), TOLS, MAXITERS);
%!               clear KSfact;
%!               KSfact = fem_sol_factor(KS, opt_sol);
%!               cpu_factor = cputime() - start;
%!             else
%!               KS = diag(D1) * Keff(perm, perm) * diag(D2);
%!             endif
%!             start = cputime();
%!             Z0 = KSfact \ (diag(D1) * Reff(perm));
%!             [Z, FLAG, RELRES, ITER] = gmres(@(Z) KSfact \ (KS * Z), Z0, RESTART, TOL, MAXITER, [], [], Z0);
%!             Phi(perm) = diag(D2) * Z;
%!             cpu_solve = cputime() - start;
%!           until (FLAG == 0 || do_fact)
%!         otherwise
%!           error("unknown solver: \"%s\"", solver);
%!       endswitch
%!       KeffPhi = Keff * Phi;
%!       f1 = norm(KeffPhi - Reff);
%!       f2 = norm(KeffPhi + Reff);
%!       sol.Phi = sol.PhiP = zeros(rows(mesh.nodes), 1);
%!       idx = dof_map.ndof(:, 7);
%!       iact = find(idx > 0);
%!       sol.Phi(iact) = Phi(idx(iact));
%!       sol.PhiP(iact) = 1j * omega * Phi(idx(iact));
%!       sol.p = -sol.PhiP;
%!       [sol.particle_velocity, ...
%!        sol.acoustic_intensity] = fem_ass_matrix(mesh, ...
%!                                                 dof_map, ...
%!                                                 [FEM_VEC_PARTICLE_VELOCITY_C, ...
%!                                                  FEM_SCA_ACOUSTIC_INTENSITY_C], ...
%!                                                 load_case, ...
%!                                                 sol);
%!       k = omega / c2;
%!       z = rho2 * c2 * k * r .* exp(1j * acot(k * r)) ./ sqrt(1 + (k * r).^2); ## according Jont Allen equation 5.11.9
%!       p = sol.p(mesh.elements.acoustic_boundary.tria6h);
%!       vx = sol.particle_velocity.vn.tria6h;
%!       R = abs(mean(mean((p - z * vx) ./ (p + z * vx))));
%!       T(i) = 1 ./ (1 - R.^2);
%!       P(i) = sum(sol.acoustic_intensity.P.tria6h);
%!       fprintf(stderr, "%d/%d %.2fHz P=%.2fdB TL=%.2fdB res=%.2e iter=%d/%d\n", i, numel(f), f(i) * unit_second^-1, 10 * log10(P(i) / P0), 10 * log10(T(i)), f1 / f2, ITER(1), ITER(2));
%!       solt.t = linspace(0, 2 * pi / omega, 36);
%!       solt.p = real(-sol.PhiP .* exp(1j * omega * solt.t));
%!       solt.def = zeros(rows(mesh.nodes), 6, numel(solt.t));
%!       solt.v = zeros(rows(mesh.nodes), 3, numel(solt.t));
%!       for j=1:6
%!         idx = dof_map.ndof(:, j);
%!         iact = find(idx > 0);
%!         for k=1:numel(solt.t)
%!           solt.def(iact, j, k) = real(Phi(idx(iact)) * exp(1j * omega * solt.t(k)));
%!         endfor
%!       endfor
%!       for j=1:columns(mesh.elements.tet10)
%!         for k=1:3
%!           for l=1:numel(solt.t)
%!             solt.v(mesh.elements.tet10(:, j), k, l) = real(sol.particle_velocity.v.tet10(:, j, k) .* exp(1j * omega * solt.t(l)));
%!           endfor
%!         endfor
%!       endfor
%!       if (f_enable_plot)
%!         Zc = linspace(0, r + b + deltaPML , 300);
%!         Xc = Yc = zeros(size(Zc));
%!         delta0 = 0.1 / unit_meters;
%!         idxc = find((mesh.nodes(:, 3) >= 0) & (sqrt(mesh.nodes(:, 1).^2 + mesh.nodes(:, 2).^2) < delta0));
%!         pc = zeros(numel(Zc), columns(solt.p));
%!         for j=1:columns(pc)
%!           pc(:, j) = griddata3(mesh.nodes(idxc, 1), mesh.nodes(idxc, 2), mesh.nodes(idxc, 3), solt.p(idxc, j), Xc, Yc, Zc);
%!         endfor
%!         for j=1:columns(pc)
%!           figure("visible","off");hold on
%!           plot(Zc * unit_meters, pc(:, j) * unit_pascal, "-;p(X=0,Y=0,Z);r");
%!           xlabel("r [m]");
%!           ylabel("p [Pa]");
%!           grid on;
%!           grid minor on;
%!           title(sprintf("sound pressure Phi=%.1fdeg f=%.1fHz", 180 / pi * (solt.t(i) * omega), f(i) * unit_second^-1));
%!           ylim([min(min(pc)),max(max(pc))] * unit_pascal);
%!         endfor
%!       endif
%!     endfor
%!   unwind_protect_cleanup
%!     clear KSfact;
%!   end_unwind_protect
%!   figure("visible", "off");
%!   hold on;
%!   plot(f * unit_second^-1, 10 * log10(P / P0), "-;P(f);r");
%!   xlabel("f [Hz]");
%!   ylabel("Lw [dB]");
%!   grid on;
%!   grid minor on;
%!   title("sound power level");
%!   figure("visible", "off");
%!   hold on;
%!   plot(f * unit_second^-1, 10 * log10(T), "-;TL(f);r");
%!   xlabel("f [Hz]");
%!   ylabel("TL [dB]");
%!   title("transmission loss");
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
