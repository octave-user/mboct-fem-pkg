## fem_pre_mesh_import.m:253
%!test
%! try
%! ## TEST 253
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
%!     unit_meters = 1e-1;
%!     unit_second = 1e-2;
%!     unit_kilograms = 1;
%!     unit_newton = unit_kilograms * unit_meters / unit_second^2;
%!     unit_pascal = unit_newton / unit_meters^2;
%!     unit_watt = unit_newton * unit_meters / unit_second;
%!     e = 10e-3 / unit_meters;
%!     D = 200e-3 / unit_meters;
%!     h = 50e-3 / unit_meters;
%!     c1 = 1400 / (unit_meters / unit_second);
%!     rho1 = 1000 / (unit_kilograms / unit_meters^3);
%!     E2 = 1e6 / unit_pascal;
%!     nu2 = 0.3;
%!     rho2 = 0 / (unit_kilograms / unit_meters^3)
%!     l = 5000e-3 / unit_meters;
%!     vn = 50 / (unit_meters / unit_second);
%!     K = rho1 * c1^2;
%!     ## Twyman, J. (2016). Wave Speed Calculation For Water Hammer Analysis.
%!     Psi = 1 / (1 + e / D) * (5 / 4 - nu2 + 2 * e / D * (1 + nu2) * (1 + e / D)); ## Table 1
%!     a = sqrt(K / (rho1 * (1 + D / e * K / E2 * Psi))); ## Equation (3)
%!     l = h * ceil(l / h);
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "e = %.16g;\n", e);
%!     fprintf(fd, "D = %.16g;\n", D);
%!     fprintf(fd, "h = %.16g;\n", h);
%!     fprintf(fd, "l = %.16g;\n", l);
%!     fputs(fd, "Point(1) = {0, 0, 0};\n");
%!     fputs(fd, "Point(2) = {0, 0.5 * D, 0};\n");
%!     fputs(fd, "Point(3) = {0, 0.5 * D + e, 0};\n");
%!     fputs(fd, "Line(1) = {1, 2};\n");
%!     fputs(fd, "Line(2) = {2, 3};\n");
%!     fputs(fd, "tmp1 = Extrude{l, 0, 0}{Line{1,2}; Layers{Ceil(l/h)}; Recombine; };\n");
%!     fputs(fd, "tmp2 = Extrude{{1,0,0},{0, 0, 0},Pi/2}{Surface{tmp1[1],tmp1[5]}; Layers{Ceil(0.5 * D * Pi / 2/ h)}; Recombine; };\n");
%!     fputs(fd, "Physical Volume(\"fluid\", 1) = {1};\n");
%!     fputs(fd, "Physical Volume(\"solid\", 2) = {2};\n");
%!     fputs(fd, "Physical Surface(\"fsi-boundary\", 3) = {3};\n");
%!     fputs(fd, "Physical Surface(\"solid-sym-xz\", 4) = {10};\n");
%!     fputs(fd, "Physical Surface(\"solid-sym-xy\", 5) = {2};\n");
%!     fputs(fd, "Physical Surface(\"solid-clamp-yz\", 6) = {8};\n");
%!     fputs(fd, "Physical Surface(\"fluid-sym-xz\", 7) = {6};\n");
%!     fputs(fd, "Physical Surface(\"fluid-sym-xy\", 8) = {1};\n");
%!     fputs(fd, "Physical Surface(\"fluid-inlet\", 9) = {4};\n");
%!     fputs(fd, "Physical Surface(\"fluid-outlet\", 10) = {5};\n");
%!     fputs(fd, "ReorientMesh Volume{2};\n");
%!     fputs(fd, "MeshSize{PointsOf{Volume{1};}} = h;\n");
%!     fputs(fd, "MeshSize{PointsOf{Volume{2};}} = h;\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete=1;\n");
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
%!   opt_mesh.elem_type = {"iso20", "quad8", "penta15", "tria6h"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_mesh));
%!   [~] = unlink([filename, ".msh"]);
%!   grp_idx_fluid_h = find([mesh.groups.iso20.id] == 1);
%!   grp_idx_solid_h = find([mesh.groups.iso20.id] == 2);
%!   if (isfield(mesh.elements, "penta15"))
%!     grp_idx_fluid_p = find([mesh.groups.penta15.id] == 1);
%!     grp_idx_solid_p = find([mesh.groups.penta15.id] == 2);
%!   endif
%!   grp_idx_fsi = find([mesh.groups.quad8.id] == 3);
%!   grp_idx_solid_sym_xz = find([mesh.groups.quad8.id] == 4);
%!   grp_idx_solid_sym_xy = find([mesh.groups.quad8.id] == 5);
%!   grp_idx_solid_clamp_yz = find([mesh.groups.quad8.id] == 6);
%!   grp_idx_fluid_sym_xz = find([mesh.groups.quad8.id] == 7);
%!   grp_idx_fluid_sym_xy = find([mesh.groups.quad8.id] == 8);
%!   grp_idx_fluid_inlet_q = find([mesh.groups.quad8.id] == 9);
%!   grp_idx_fluid_outlet_q = find([mesh.groups.quad8.id] == 10);
%!   if (isfield(mesh.groups, "tria6h"))
%!     grp_idx_fluid_inlet_t = find([mesh.groups.tria6h.id] == 9);
%!     grp_idx_fluid_outlet_t = find([mesh.groups.tria6h.id] == 10);
%!   endif
%!   empty_cell = cell(1, 2);
%!   mesh.material_data = struct("E", empty_cell, "nu", empty_cell, "rho", empty_cell, "c", empty_cell);
%!   mesh.material_data(1).c = c1;
%!   mesh.material_data(1).rho = rho1;
%!   mesh.material_data(2).E = E2;
%!   mesh.material_data(2).nu = nu2;
%!   mesh.material_data(2).rho = rho2;
%!   mesh.materials.iso20 = zeros(rows(mesh.elements.iso20), 1, "int32");
%!   mesh.materials.iso20(mesh.groups.iso20(grp_idx_fluid_h).elements) = 1;
%!   mesh.materials.iso20(mesh.groups.iso20(grp_idx_solid_h).elements) = 2;
%!   if (isfield(mesh.elements, "penta15"))
%!     mesh.materials.penta15 = zeros(rows(mesh.elements.penta15), 1, "int32");
%!     if (~isempty(grp_idx_fluid_p))
%!       mesh.materials.penta15(mesh.groups.penta15(grp_idx_fluid_p).elements) = 1;
%!     endif
%!     if (~isempty(grp_idx_solid_p))
%!       mesh.materials.penta15(mesh.groups.penta15(grp_idx_solid_p).elements) = 2;
%!     endif
%!   endif
%!   if (isfield(mesh.groups, "tria6h"))
%!     node_idx_inlet = [mesh.groups.quad8(grp_idx_fluid_inlet_q).nodes, ...
%!                       mesh.groups.tria6h(grp_idx_fluid_inlet_t).nodes];
%!     node_idx_outlet = [mesh.groups.quad8(grp_idx_fluid_outlet_q).nodes, ...
%!                        mesh.groups.tria6h(grp_idx_fluid_outlet_t).nodes];
%!   else
%!     node_idx_inlet = mesh.groups.quad8(grp_idx_fluid_inlet_q).nodes;
%!     node_idx_outlet = mesh.groups.quad8(grp_idx_fluid_outlet_q).nodes;
%!   endif
%!   mesh.elements.fluid_struct_interface.quad8 = mesh.elements.quad8(mesh.groups.quad8(grp_idx_fsi).elements, :);
%!   mesh.elements.particle_velocity.quad8.nodes = mesh.elements.quad8(mesh.groups.quad8(grp_idx_fluid_inlet_q).elements, :);
%!   mesh.materials.particle_velocity.quad8 = ones(rows(mesh.elements.particle_velocity.quad8.nodes), 1, "int32");
%!   load_case.particle_velocity.quad8.vn = repmat(-vn, size(mesh.elements.particle_velocity.quad8.nodes));
%!   if (isfield(mesh.groups, "tria6h"))
%!     mesh.elements.particle_velocity.tria6h.nodes = mesh.elements.tria6h(mesh.groups.tria6h(grp_idx_fluid_inlet_t).elements, :);
%!     mesh.materials.particle_velocity.tria6h = ones(rows(mesh.elements.particle_velocity.tria6h.nodes), 1, "int32");
%!     load_case.particle_velocity.tria6h.vn = repmat(-vn, size(mesh.elements.particle_velocity.tria6h.nodes));
%!   endif
%!   mesh.elements.acoustic_impedance.quad8.nodes = mesh.elements.quad8(mesh.groups.quad8(grp_idx_fluid_outlet_q).elements, :);
%!   mesh.elements.acoustic_impedance.quad8.z = repmat(rho1 * a, size(mesh.elements.acoustic_impedance.quad8.nodes));
%!   mesh.materials.acoustic_impedance.quad8 = ones(rows(mesh.elements.acoustic_impedance.quad8.nodes), 1, "int32");
%!   if (isfield(mesh.groups, "tria6h"))
%!     mesh.elements.acoustic_impedance.tria6h.nodes = mesh.elements.tria6h(mesh.groups.tria6h(grp_idx_fluid_outlet_t).elements, :);
%!     mesh.elements.acoustic_impedance.tria6h.z = repmat(rho1 * a, size(mesh.elements.acoustic_impedance.tria6h.nodes));
%!     mesh.materials.acoustic_impedance.tria6h = ones(rows(mesh.elements.acoustic_impedance.tria6h.nodes), 1, "int32");
%!   endif
%!   load_case_dof.domain = FEM_DO_FLUID_STRUCT;
%!   load_case_dof.locked_dof = false(rows(mesh.nodes), 7);
%!   load_case_dof.locked_dof(mesh.groups.quad8(grp_idx_solid_sym_xz).nodes, 2) = true;
%!   load_case_dof.locked_dof(mesh.groups.quad8(grp_idx_solid_sym_xy).nodes, 3) = true;
%!   load_case_dof.locked_dof(mesh.groups.quad8(grp_idx_solid_clamp_yz).nodes, 1) = true;
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!   [mat_ass.Mfs_re, ...
%!    mat_ass.Dfs_re, ...
%!    mat_ass.Kfs_re, ...
%!    mat_ass.Rfs, ...
%!    mat_ass.mat_info, ...
%!    mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                        dof_map, ...
%!                                        [FEM_MAT_MASS_FLUID_STRUCT_RE, ...
%!                                         FEM_MAT_DAMPING_FLUID_STRUCT_RE, ...
%!                                         FEM_MAT_STIFFNESS_FLUID_STRUCT_RE, ...
%!                                         FEM_VEC_LOAD_FLUID_STRUCT], ...
%!                                        load_case);
%!   T = 1.5 * l / a;
%!   dt = h / a;
%!   T = ceil(T / dt) * dt;
%!   sol.t = 0:dt:T;
%!   tau = l / a;
%!   f = sin(pi * sol.t / tau).^2 .* (sol.t <= tau / 2) + (sol.t > tau / 2);
%!   Z = ZP = ZPP = zeros(dof_map.totdof, 1);
%!   opt_sol.solver = "pardiso";
%!   opt_sol.number_of_threads = mbdyn_solver_num_threads_default();
%!   opt_sol.delta = 0.5;
%!   opt_sol.pre_scaling = true;
%!   opt_sol.refine_max_iter = int32(1000);
%!   solver = fem_sol_transient_init(mat_ass.Mfs_re, mat_ass.Dfs_re, mat_ass.Kfs_re, dt, opt_sol);
%!   sol.p = zeros(rows(mesh.nodes), numel(sol.t));
%!   sol.def = zeros(rows(mesh.nodes), 6, numel(sol.t));
%!   idx_p = find(dof_map.ndof(:, 7) > 0);
%!   for i=2:numel(sol.t)
%!     [Z, ZP, ZPP, res] = fem_sol_transient_step(Z, ...
%!                                                ZP, ...
%!                                                ZPP, ...
%!                                                mat_ass.Rfs * f(i), ...
%!                                                solver);
%!     fprintf(stderr, "Step %d/%d: %.2e\n", i - 1, numel(sol.t) - 1, res);
%!     sol.p(idx_p, i) = -ZP(dof_map.ndof(idx_p, 7));
%!     for j=1:6
%!       idx_def = find(dof_map.ndof(:, j) > 0);
%!       sol.def(idx_def, j, i) = Z(dof_map.ndof(idx_def, j));
%!     endfor
%!   endfor
%!   pin = mean(sol.p(node_idx_inlet, :), 1);
%!   pout = mean(sol.p(node_idx_outlet, :), 1);
%!   figure("visible", "off");
%!   hold on;
%!   plot(sol.t * unit_second, pin * unit_pascal, "-;inlet;r");
%!   plot((sol.t - l / a) * unit_second, pout * unit_pascal, "-;outlet;b");
%!   plot(sol.t * unit_second, rho1 * a * vn * f * unit_pascal, "-;z * vn;k");
%!   xlabel("t [s]");
%!   ylabel("p [Pa]");
%!   grid on;
%!   grid minor on;
%!   title("inlet/outlet pressure versus time");
%!   pref = 0.5 * pin(end);
%!   tref0 = sol.t(find(pin >= pref)(1));
%!   tref1 = sol.t(find(pout >= pref)(1));
%!   tol = 3e-2;
%!   assert_simple(tref1 - tref0, l / a, tol * l / a);
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
