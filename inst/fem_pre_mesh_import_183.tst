## fem_pre_mesh_import.m:183
%!test
%! ### TEST 183
%! ### rigid body motion of a solid domain between two fluid domains with prescribed pressure within the fluid
%! close all;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   unwind_protect
%!   [fd, msg] = fopen([filename, ".geo"], "w");
%!   if (fd == -1)
%!     error("failed to open file \"%s.geo\"", filename);
%!   endif
%!   unit_meters = 1e-3;
%!   unit_second = 1e2;
%!   unit_kilograms = 1e3;
%!   unit_newton = unit_kilograms * unit_meters / unit_second^2;
%!   unit_pascal = unit_newton / unit_meters^2;
%!   l1 = 100e-3 / unit_meters;
%!   l2 = 10e-3 / unit_meters;
%!   l3 = 100e-3 / unit_meters;
%!   c1 = 1400 / (unit_meters / unit_second);
%!   rho1 = 1000 / (unit_kilograms / unit_meters^3);
%!   eta1 = 1e-3 / (unit_pascal * unit_second);
%!   zeta1 = 0 / (unit_pascal * unit_second);
%!   E2 = 210000e6 / unit_pascal;
%!   rho2 = 7850 / (unit_kilograms / unit_meters^3);
%!   nu2 = 0.3;
%!   alpha2 = 1e-6 / (unit_second^-1);
%!   beta2 = 1e-8 / unit_second;
%!   c3 = 1400 / (unit_meters / unit_second);
%!   rho3 = 1000 / (unit_kilograms / unit_meters^3);
%!   eta3 = 1e-3 / (unit_pascal * unit_second);
%!   zeta3 = 0 / (unit_pascal * unit_second);
%!   dx = 1e-3 / unit_meters;
%!   w = dx;
%!   h = dx;
%!   omega = 100 / (1 / unit_second);
%!   pinput = (2 + 1j) / unit_pascal;
%!   poutput = (2 - 1j) / unit_pascal;
%!   Uref = -(pinput - poutput) / (omega^2 * rho2 * l2);
%!   fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!   fprintf(fd, "l1=%g;\n", l1);
%!   fprintf(fd, "l2=%g;\n", l2);
%!   fprintf(fd, "l3=%g;\n", l3);
%!   fprintf(fd, "w=%g;\n", w);
%!   fprintf(fd, "h=%g;\n", h);
%!   fprintf(fd, "dx = %g;\n", dx);
%!   fputs(fd, "Point(1) = {0,0,0};\n");
%!   fputs(fd, "Point(2) = {0,w,0};\n");
%!   fputs(fd, "Point(3) = {0,w,h};\n");
%!   fputs(fd, "Point(4) = {0,0,h};\n");
%!   fputs(fd, "Line(1) = {1,2};\n");
%!   fputs(fd, "Line(2) = {2,3};\n");
%!   fputs(fd, "Line(3) = {3,4};\n");
%!   fputs(fd, "Line(4) = {4,1};\n");
%!   fputs(fd, "Line Loop(5) = {1,2,3,4};\n");
%!   fputs(fd, "Plane Surface(6) = {5};\n");
%!   fputs(fd, "tmp1[] = Extrude {l1,0,0} {\n");
%!   fputs(fd, "  Surface{6}; Layers{Ceil(l1/dx)}; Recombine;\n");
%!   fputs(fd, "};\n");
%!   fputs(fd, "Recombine Surface{6, tmp1[0]};\n");
%!   fputs(fd, "tmp2[] = Extrude {l2,0,0} {\n");
%!   fputs(fd, "  Surface{tmp1[0]}; Layers{Ceil(l2/dx)}; Recombine;\n");
%!   fputs(fd, "};\n");
%!   fputs(fd, "Recombine Surface{tmp1[0], tmp2[0]};\n");
%!   fputs(fd, "tmp3[] = Extrude {l3,0,0} {\n");
%!   fputs(fd, "  Surface{tmp2[0]}; Layers{Ceil(l3/dx)}; Recombine;\n");
%!   fputs(fd, "};\n");
%!   fputs(fd, "Recombine Surface{tmp2[0], tmp3[0]};\n");
%!   fputs(fd, "Physical Surface(\"input\",1) = {6};\n");
%!   fputs(fd, "Physical Surface(\"output\",2) = {tmp3[0]};\n");
%!   fputs(fd, "Physical Surface(\"fsi-bnd1\",3) = {tmp1[0]};\n");
%!   fputs(fd, "Physical Surface(\"fsi-bnd2\",4) = {tmp2[0]};\n");
%!   fputs(fd, "Physical Surface(\"slider1\",8) = {tmp2[2]};\n");
%!   fputs(fd, "Physical Surface(\"slider2\",9) = {tmp2[5]};\n");
%!   fputs(fd, "Physical Volume(\"volume1\",5) = {tmp1[1]};\n");
%!   fputs(fd, "Physical Volume(\"volume2\",6) = {tmp2[1]};\n");
%!   fputs(fd, "Physical Volume(\"volume3\",7) = {tmp3[1]};\n");
%!   fputs(fd, "Mesh.OptimizeThreshold=0.99;\n");
%!   fputs(fd, "ReorientMesh Volume{tmp2[1]};\n");
%!   fputs(fd, "Mesh.SecondOrderIncomplete=1;\n");
%!   fputs(fd, "Mesh.HighOrderOptimize=2;\n");
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
%!   opt_mesh.elem_type = {"iso20", "quad8"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_mesh));
%!   [~] = unlink([filename, ".msh"]);
%!   grp_idx_input = find([mesh.groups.quad8.id] == 1);
%!   grp_idx_output = find([mesh.groups.quad8.id] == 2);
%!   grp_idx_fsi1 = find([mesh.groups.quad8.id] == 3);
%!   grp_idx_fsi2 = find([mesh.groups.quad8.id] == 4);
%!   grp_idx_volume1 = find([mesh.groups.iso20.id] == 5);
%!   grp_idx_volume2 = find([mesh.groups.iso20.id] == 6);
%!   grp_idx_volume3 = find([mesh.groups.iso20.id] == 7);
%!   grp_idx_slider1 = find([mesh.groups.quad8.id] == 8);
%!   grp_idx_slider2 = find([mesh.groups.quad8.id] == 9);
%!   load_case_dof.locked_dof = false(rows(mesh.nodes), 7);
%!   load_case_dof.locked_dof(mesh.groups.quad8(grp_idx_slider1).nodes, 3) = true;
%!   load_case_dof.locked_dof(mesh.groups.quad8(grp_idx_slider2).nodes, 2) = true;
%!   load_case_dof.domain = FEM_DO_FLUID_STRUCT;
%!   mesh.materials.iso20 = zeros(rows(mesh.elements.iso20), 1, "int32");
%!   mesh.materials.iso20(mesh.groups.iso20(grp_idx_volume1).elements) = 1;
%!   mesh.materials.iso20(mesh.groups.iso20(grp_idx_volume2).elements) = 2;
%!   mesh.materials.iso20(mesh.groups.iso20(grp_idx_volume3).elements) = 3;
%!   mesh.elements.acoustic_boundary.quad8 = mesh.elements.quad8([[mesh.groups.quad8([grp_idx_input, grp_idx_output])].elements], :);
%!   mesh.materials.acoustic_boundary.quad8 = zeros(rows(mesh.elements.acoustic_boundary.quad8), 1, "int32");
%!   mesh.materials.acoustic_boundary.quad8(mesh.groups.quad8(grp_idx_input).elements) = 1;
%!   mesh.materials.acoustic_boundary.quad8(mesh.groups.quad8(grp_idx_output).elements) = 3;
%!   mesh.elements.fluid_struct_interface.quad8 = mesh.elements.quad8([[mesh.groups.quad8([grp_idx_fsi1, grp_idx_fsi2])].elements], :);
%!   load_case = struct("particle_velocity", cell(1, 2));
%!   node_idx_input = mesh.groups.iso20(grp_idx_volume1).nodes;
%!   node_idx_output = mesh.groups.iso20(grp_idx_volume3).nodes;
%!   node_idx_constr = [node_idx_input, node_idx_output];
%!   p_constr = [repmat(pinput, 1, numel(node_idx_input)), ...
%!               repmat(poutput, 1, numel(node_idx_output))];
%!   mesh.elements.acoustic_constr = struct("C", mat2cell(ones(1, numel(node_idx_constr)), 1, ones(1, numel(node_idx_constr))), ...
%!                                          "nodes", mat2cell(node_idx_constr, 1, ones(1, numel(node_idx_constr))), ...
%!                                          "scale", mat2cell(repmat(1/omega, 1, numel(node_idx_constr)), 1, ones(1, numel(node_idx_constr))));
%!   load_case = struct("acoustic_constr", cell(1, 2));
%!   load_case(1).acoustic_constr = struct("p", mat2cell(real(p_constr), 1, ones(1, numel(p_constr))));
%!   load_case(2).acoustic_constr = struct("p", mat2cell(imag(p_constr), 1, ones(1, numel(p_constr))));
%!   mesh.material_data = struct("E", {[], E2, []}, ...
%!                               "rho", {rho1, rho2, rho3}, ...
%!                               "nu", {[], nu2, []}, ...
%!                               "c", {c1, [], c3}, ...
%!                               "eta", {eta1, [], eta3}, ...
%!                               "zeta", {zeta1, [], zeta3}, ...
%!                               "alpha", {[], alpha2, []}, ...
%!                               "beta", {[], beta2, []});
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!   [mat_ass.Kfs, ...
%!    mat_ass.Mfs, ...
%!    mat_ass.Dfs_re, ...
%!    mat_ass.Dfs_im, ...
%!    mat_ass.Rfs, ...
%!    mat_ass.coll_Kfs, ...
%!    mat_ass.coll_Mfs, ...
%!    mat_ass.mat_info, ...
%!    mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                        dof_map, ...
%!                                        [FEM_MAT_STIFFNESS_FLUID_STRUCT_RE, ...
%!                                         FEM_MAT_MASS_FLUID_STRUCT_RE, ...
%!                                         FEM_MAT_DAMPING_FLUID_STRUCT_RE, ...
%!                                         FEM_MAT_DAMPING_FLUID_STRUCT_IM, ...
%!                                         FEM_VEC_LOAD_FLUID_STRUCT, ...
%!                                         FEM_VEC_COLL_STIFF_FLUID_STRUCT, ...
%!                                         FEM_VEC_COLL_MASS_FLUID_STRUCT], ...
%!                                        load_case);
%!   opt_sol.number_of_threads = mbdyn_solver_num_threads_default();
%!   opt_sol.solver = "pastix";
%!   opt_sol.refine_max_iter = int32(250);
%!   opt_sol.verbose = int32(0);
%!   opt_sol.pre_scaling = true;
%!   Keff = -omega^2 * mat_ass.Mfs + 1j * omega * complex(mat_ass.Dfs_re, mat_ass.Dfs_im) + mat_ass.Kfs;
%!   Reff = complex(mat_ass.Rfs(:, 1), mat_ass.Rfs(:, 2));
%!   Z = fem_sol_factor(Keff, opt_sol) \ Reff;
%!   Psi = linspace(0, 2 * pi, 37);
%!   sol.t = Psi / omega;
%!   x = mesh.nodes(:, 1);
%!   idxPhi = dof_map.ndof(:, 7);
%!   idxPhi1 = find(idxPhi > 0);
%!   idxPhi = idxPhi(idxPhi1);
%!   sol.p = zeros(rows(mesh.nodes), numel(sol.t));
%!   sol.def = zeros(rows(mesh.nodes), 6, numel(sol.t));
%!   sol.p(idxPhi1, :) = real(-1j * omega * Z(idxPhi, :) .* exp(1j * omega * sol.t));
%!   for j=1:6
%!     idxU = dof_map.ndof(:, j);
%!     idxU1 = find(idxU > 0);
%!     idxU = idxU(idxU1);
%!     sol.def(idxU1, j, :) = real(Z(idxU, :) * exp(1j * omega * sol.t));
%!   endfor
%!   node_idx = mesh.groups.iso20(grp_idx_volume2).nodes;
%!   tol = 1e-5;
%!   for i=1:numel(node_idx)
%!     assert_simple(reshape(sol.def(node_idx(i), 1, :), 1, numel(sol.t)), real(Uref * exp(1j * omega * sol.t)), tol * abs(Uref));
%!   endfor
%!   [~, idx] = sort(mesh.nodes(:, 1));
%!   for i=1:numel(Psi)
%!     figure("visible", "off");
%!     hold on;
%!     plot(mesh.nodes(idx, 1) * unit_second, sol.p(idx, i) * unit_pascal, "-;p;r");
%!     ylim([min(min(sol.p)), max(max(sol.p))] * unit_pascal);
%!     xlabel("x [m]");
%!     ylabel("p [Pa]");
%!     grid on;
%!     grid minor on;
%!     title(sprintf("pressure distribution Psi=%.1fdeg", Psi(i) * 180 / pi));
%!   endfor
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       [~] = unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
