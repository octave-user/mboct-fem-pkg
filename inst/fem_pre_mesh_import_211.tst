## fem_pre_mesh_import.m:211
%!test
%! try
%! ### TEST 211
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
%!   do_plot = false;
%!   if (do_plot)
%!     close all;
%!   endif
%!   unit_meters = 1e-3;
%!   unit_second = 1e3;
%!   unit_kilograms = 1e3;
%!   unit_newton = unit_kilograms * unit_meters / unit_second^2;
%!   unit_pascal = unit_newton / unit_meters^2;
%!   unit_watt = unit_newton * unit_meters / unit_second;
%!   r = 20e-3 / unit_meters;
%!   R = 30e-3 / unit_meters;
%!   H = 20e-3 / unit_meters;
%!   h = 5e-3 / unit_meters;
%!   rho1 = 1.25 / (unit_kilograms / unit_meters^3);
%!   c1 = 340 / (unit_meters / unit_second);
%!   eta1 = 0 / (unit_pascal * unit_second);
%!   zeta1 = 0 / (unit_pascal * unit_second);
%!   E2 = 1000e6 / unit_pascal;
%!   rho2 = 1200 / (unit_kilograms / unit_meters^3);
%!   nu2 = 0.4;
%!   K2 = E2 / (3 * (1 - 2 * nu2));
%!   c2 = sqrt(K2 / rho2);
%!   alpha2 = 0 / (unit_second^-1);
%!   beta2 = 0 / unit_second;
%!   vz = 1e-3 / (unit_meters / unit_second);
%!   omega = 2 * pi * linspace(1, 10000, 20) / unit_second^-1;
%!   fmax = max(omega) / (2 * pi);
%!   lambdamin = c1 / fmax;
%!   dx = lambdamin / 6;
%!   fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!   fprintf(fd, "r=%g;\n", r);
%!   fprintf(fd, "R=%g;\n", R);
%!   fprintf(fd, "H=%g;\n", H);
%!   fprintf(fd, "h=%g;\n", h);
%!   fprintf(fd, "dx=%g;\n", dx);
%!   fputs(fd, "Point(1) = {0,0,0,dx};\n");
%!   fputs(fd, "Point(2) = {r,0,0,dx};\n");
%!   fputs(fd, "Point(3) = {R,0,0,dx};\n");
%!   fputs(fd, "Line(1) = {1,2};\n");
%!   fputs(fd, "Line(2) = {2,3};\n");
%!   fputs(fd, "A1[] = Extrude{0,0,H} {\n");
%!   fputs(fd, "  Line{1}; Layers{Ceil(H / dx)}; Recombine;\n");
%!   fputs(fd, "};\n");
%!   fputs(fd, "A2[] = Extrude{0,0,H} {\n");
%!   fputs(fd, "  Line{2}; Layers{Ceil(H / dx)}; Recombine;\n");
%!   fputs(fd, "};\n");
%!   fputs(fd, "A3[] = Extrude{0,0,-h} {\n");
%!   fputs(fd, "  Line{1}; Layers{Ceil(h / dx)}; Recombine;\n");
%!   fputs(fd, "};\n");
%!   fputs(fd, "A4[] = Extrude{0,0,-h} {\n");
%!   fputs(fd, "  Line{2}; Layers{Ceil(h / dx)}; Recombine;\n");
%!   fputs(fd, "};\n");
%!   fputs(fd, "A5[] = Extrude{0,0,-H} {\n");
%!   fputs(fd, "  Line{A3[0]}; Layers{Ceil(H / dx)}; Recombine;\n");
%!   fputs(fd, "};\n");
%!   fputs(fd, "A6[] = Extrude{0,0,-H} {\n");
%!   fputs(fd, "  Line{A4[0]}; Layers{Ceil(H / dx)}; Recombine;\n");
%!   fputs(fd, "};\n");
%!   fputs(fd, "V1[] = Extrude{{0,0,1},{0,0,0},-Pi/2} {\n");
%!   fputs(fd, "  Surface{A1[1]}; Layers{Ceil(Pi / 2 * r / dx)}; Recombine;\n");
%!   fputs(fd, "};\n");
%!   fputs(fd, "V2[] = Extrude{{0,0,1},{0,0,0},-Pi/2} {\n");
%!   fputs(fd, "  Surface{A2[1]}; Layers{Ceil(Pi / 2 * r / dx)}; Recombine;\n");
%!   fputs(fd, "};\n");
%!   fputs(fd, "V3[] = Extrude{{0,0,1},{0,0,0},-Pi/2} {\n");
%!   fputs(fd, "  Surface{A3[1]}; Layers{Ceil(Pi / 2 * r / dx)}; Recombine;\n");
%!   fputs(fd, "};\n");
%!   fputs(fd, "V4[] = Extrude{{0,0,1},{0,0,0},-Pi/2} {\n");
%!   fputs(fd, "  Surface{A4[1]}; Layers{Ceil(Pi / 2 * r / dx)}; Recombine;\n");
%!   fputs(fd, "};\n");
%!   fputs(fd, "V5[] = Extrude{{0,0,1},{0,0,0},-Pi/2} {\n");
%!   fputs(fd, "  Surface{A5[1]}; Layers{Ceil(Pi / 2 * r / dx)}; Recombine;\n");
%!   fputs(fd, "};\n");
%!   fputs(fd, "V6[] = Extrude{{0,0,1},{0,0,0},-Pi/2} {\n");
%!   fputs(fd, "  Surface{A6[1]}; Layers{Ceil(Pi / 2 * r / dx)}; Recombine;\n");
%!   fputs(fd, "};\n");
%!   fputs(fd, "Coherence;\n");
%!   fputs(fd, "Physical Volume(\"v1\",1)={V3[1]};\n");
%!   fputs(fd, "Physical Volume(\"v2\",2)={V1[1],V2[1],V4[1],V5[1],V6[1]};\n");
%!   fputs(fd, "Physical Surface(\"s1\", 1) = {13};\n");
%!   fputs(fd, "Physical Surface(\"s2\", 2) = {14};\n");
%!   fputs(fd, "Physical Surface(\"s3\", 3) = {2};\n");
%!   fputs(fd, "Physical Surface(\"s4\", 4) = {12};\n");
%!   fputs(fd, "Physical Surface(\"s5\", 5) = {11};\n");
%!   fputs(fd, "Physical Surface(\"s6\", 6) = {3, 8};\n");
%!   fputs(fd, "Physical Surface(\"s7\", 7) = {20, 24};\n");
%!   fputs(fd, "ReverseMesh Surface{2};\n");
%!   fputs(fd, "Mesh.SecondOrderIncomplete=1;\n");
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
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!   [~] = unlink([filename, ".msh"]);
%!   grp_idx_v1_iso20 = find([mesh.groups.iso20.id] == 1);
%!   grp_idx_v1_penta15 = find([mesh.groups.penta15.id] == 1);
%!   grp_idx_v2_iso20 = find([mesh.groups.iso20.id] == 2);
%!   grp_idx_v2_penta15 = find([mesh.groups.penta15.id] == 2);
%!   grp_idx_s1_quad8 = find([mesh.groups.quad8.id] == 1);
%!   grp_idx_s2_quad8 = find([mesh.groups.quad8.id] == 2);
%!   grp_idx_s3_quad8 = find([mesh.groups.quad8.id] == 3);
%!   grp_idx_s3_tria6 = find([mesh.groups.tria6.id] == 3);
%!   grp_idx_s4_quad8 = find([mesh.groups.quad8.id] == 4);
%!   grp_idx_s4_tria6 = find([mesh.groups.tria6.id] == 4);
%!   grp_idx_s5_quad8 = find([mesh.groups.quad8.id] == 5);
%!   grp_idx_s6_quad8 = find([mesh.groups.quad8.id] == 6);
%!   grp_idx_s6_tria6 = find([mesh.groups.tria6.id] == 6);
%!   grp_idx_s7_quad8 = find([mesh.groups.quad8.id] == 7);
%!   grp_idx_s7_tria6 = find([mesh.groups.tria6.id] == 7);
%!   mesh.materials.iso20 = zeros(rows(mesh.elements.iso20), 1, "int32");
%!   mesh.materials.penta15 = zeros(rows(mesh.elements.penta15), 1, "int32");
%!   mesh.materials.iso20(mesh.groups.iso20(grp_idx_v1_iso20).elements) = 2;
%!   mesh.materials.penta15(mesh.groups.penta15(grp_idx_v1_penta15).elements) = 2;
%!   mesh.materials.iso20(mesh.groups.iso20(grp_idx_v2_iso20).elements) = 1;
%!   mesh.materials.penta15(mesh.groups.penta15(grp_idx_v2_penta15).elements) = 1;
%!   mesh.material_data(1).c = c1;
%!   mesh.material_data(1).rho = rho1;
%!   mesh.material_data(1).eta = eta1;
%!   mesh.material_data(1).zeta = zeta1;
%!   mesh.material_data(2).E = E2;
%!   mesh.material_data(2).nu = nu2;
%!   mesh.material_data(2).rho = rho2;
%!   mesh.material_data(2).alpha = alpha2;
%!   mesh.material_data(2).beta = beta2;
%!   elid_fsi_quad8 = [[mesh.groups.quad8([grp_idx_s3_quad8, grp_idx_s4_quad8, grp_idx_s5_quad8])].elements];
%!   elid_fsi_tria6 = [[mesh.groups.tria6([grp_idx_s3_tria6, grp_idx_s4_tria6])].elements];
%!   noid_sliding_s1 = [mesh.groups.quad8(grp_idx_s1_quad8).nodes];
%!   noid_sliding_s2 = [mesh.groups.quad8(grp_idx_s2_quad8).nodes];
%!   elid_impe_quad8 = mesh.groups.quad8(grp_idx_s6_quad8).elements;
%!   elid_impe_tria6 = mesh.groups.tria6(grp_idx_s6_tria6).elements;
%!   elid_vel_quad8 = mesh.groups.quad8(grp_idx_s7_quad8).elements;
%!   elid_vel_tria6 = mesh.groups.tria6(grp_idx_s7_tria6).elements;
%!   noid_outlet = [mesh.groups.quad8(grp_idx_s6_quad8).nodes, mesh.groups.tria6(grp_idx_s6_tria6).nodes];
%!   noid_inlet = [mesh.groups.quad8(grp_idx_s7_quad8).nodes, mesh.groups.tria6(grp_idx_s7_tria6).nodes];
%!   mesh.elements.fluid_struct_interface.quad8 = mesh.elements.quad8(elid_fsi_quad8, :);
%!   mesh.elements.fluid_struct_interface.tria6 = mesh.elements.tria6(elid_fsi_tria6, :);
%!   load_case_dof.locked_dof = false(rows(mesh.nodes), 7);
%!   load_case_dof.locked_dof(noid_sliding_s1, [2,3]) = true;
%!   load_case_dof.locked_dof(noid_sliding_s2, [1,3]) = true;
%!   load_case_dof.domain = FEM_DO_FLUID_STRUCT;
%!   mesh.elements.acoustic_impedance.quad8.nodes = mesh.elements.quad8(elid_impe_quad8, :);
%!   mesh.elements.acoustic_impedance.quad8.z = repmat(rho1 * c1, size(mesh.elements.acoustic_impedance.quad8.nodes));
%!   mesh.materials.acoustic_impedance.quad8 = ones(rows(mesh.elements.acoustic_impedance.quad8.nodes), 1, "int32");
%!   mesh.elements.acoustic_impedance.tria6.nodes = mesh.elements.tria6(elid_impe_tria6, :);
%!   mesh.elements.acoustic_impedance.tria6.z = repmat(rho1 * c1, size(mesh.elements.acoustic_impedance.tria6.nodes));
%!   mesh.materials.acoustic_impedance.tria6 = ones(rows(mesh.elements.acoustic_impedance.tria6.nodes), 1, "int32");
%!   mesh.elements.particle_velocity.quad8.nodes = mesh.elements.quad8(elid_vel_quad8, :);
%!   mesh.materials.particle_velocity.quad8 = ones(rows(mesh.elements.particle_velocity.quad8.nodes), 1, "int32");
%!   load_case(1).particle_velocity.quad8.vn = repmat(-real(vz), size(mesh.elements.particle_velocity.quad8.nodes));
%!   load_case(2).particle_velocity.quad8.vn = repmat(-imag(vz), size(mesh.elements.particle_velocity.quad8.nodes));
%!   mesh.elements.particle_velocity.tria6.nodes = mesh.elements.tria6(elid_vel_tria6, :);
%!   mesh.materials.particle_velocity.tria6 = ones(rows(mesh.elements.particle_velocity.tria6.nodes), 1, "int32");
%!   load_case(1).particle_velocity.tria6.vn = repmat(-real(vz), size(mesh.elements.particle_velocity.tria6.nodes));
%!   load_case(2).particle_velocity.tria6.vn = repmat(-imag(vz), size(mesh.elements.particle_velocity.tria6.nodes));
%!   mesh.elements.acoustic_boundary.quad8 = mesh.elements.quad8([[mesh.groups.quad8([grp_idx_s6_quad8, grp_idx_s7_quad8])].elements], :);
%!   mesh.materials.acoustic_boundary.quad8 = ones(rows(mesh.elements.acoustic_boundary.quad8), 1, "int32");
%!   mesh.elements.acoustic_boundary.tria6 = mesh.elements.tria6([[mesh.groups.tria6([grp_idx_s6_tria6, grp_idx_s7_tria6])].elements], :);
%!   mesh.materials.acoustic_boundary.tria6 = ones(rows(mesh.elements.acoustic_boundary.tria6), 1, "int32");
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!   [mat_ass.Kfs, ...
%!    mat_ass.Mfs, ...
%!    mat_ass.Dfs_re, ...
%!    mat_ass.Dfs_im, ...
%!    mat_ass.Rfs, ...
%!    mat_ass.n, ...
%!    mat_ass.mat_info, ...
%!    mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                        dof_map, ...
%!                                        [FEM_MAT_STIFFNESS_FLUID_STRUCT_RE, ...
%!                                         FEM_MAT_MASS_FLUID_STRUCT_RE, ...
%!                                         FEM_MAT_DAMPING_FLUID_STRUCT_RE, ...
%!                                         FEM_MAT_DAMPING_FLUID_STRUCT_IM, ...
%!                                         FEM_VEC_LOAD_FLUID_STRUCT, ...
%!                                         FEM_VEC_SURFACE_NORMAL_VECTOR], ...
%!                                        load_case);
%!   mat_ass.Dfs = complex(mat_ass.Dfs_re, mat_ass.Dfs_im);
%!   opt_sol.number_of_threads = mbdyn_solver_num_threads_default();
%!   opt_sol.solver = "pastix";
%!   opt_sol.verbose = int32(0);
%!   opt_sol.refine_max_iter = int32(250);
%!   idxPhi = dof_map.ndof(:, 7);
%!   idxPhi1 = find(idxPhi > 0);
%!   idxPhi = idxPhi(idxPhi1);
%!   Reff = complex(mat_ass.Rfs(:, 1), mat_ass.Rfs(:, 2));
%!   P = zeros(1, numel(omega));
%!   z0 = [linspace(-h - H, -h, 20), linspace(0, H, 20)];
%!   x0 = y0 = zeros(size(z0));
%!   for i=1:numel(omega)
%!     Keff = -omega(i)^2 * mat_ass.Mfs + 1j * omega(i) * mat_ass.Dfs + mat_ass.Kfs;
%!     Z = fem_sol_factor(Keff, opt_sol) \ Reff;
%!     if (mod(i, 10) == 1)
%!       fprintf(stderr, "%d/%d %.1fHz %e\n", i, numel(omega), omega(i) / (2 * pi) * unit_second^-1, norm(Keff * Z - Reff) / norm(Keff * Z + Reff));
%!     endif
%!     sol.PhiP = sol.Phi = zeros(rows(mesh.nodes), columns(Z));
%!     sol.Phi(idxPhi1, :) = Z(idxPhi, :);
%!     sol.PhiP(idxPhi1, :) = 1j * omega(i) * Z(idxPhi, :);
%!     sol.p = -sol.PhiP;
%!     sol.def = zeros(rows(mesh.nodes), 6, columns(Z));
%!     for j=1:6
%!       idxU = dof_map.ndof(:, j);
%!       idxU1 = find(idxU > 0);
%!       idxU = idxU(idxU1);
%!       sol.def(idxU1, j, :) = Z(idxU, :);
%!     endfor
%!     [sol.particle_velocity, ...
%!      sol.acoustic_intensity] = fem_ass_matrix(mesh, ...
%!                                               dof_map, ...
%!                                               [FEM_VEC_PARTICLE_VELOCITY_C, ...
%!                                                FEM_SCA_ACOUSTIC_INTENSITY_C], ...
%!                                                load_case, ...
%!                                                sol);
%!     P(i) = sum(sol.acoustic_intensity.P.quad8) + sum(sol.acoustic_intensity.P.tria6);
%!     if (do_plot)
%!       p0 = griddata3(mesh.nodes(:, 1), mesh.nodes(:, 2), mesh.nodes(:, 3), real(sol.p), x0, y0, z0);
%!       figure("visible", "off");
%!       plot(z0 * unit_meters, p0 * unit_pascal, "-;p(z);r");
%!       grid on;
%!       grid minor on;
%!       xlabel("z [m]");
%!       ylabel("p [Pa]");
%!       title(sprintf("pressure distribution at the center %.0fHz", omega(i)/(2 * pi) * unit_second^-1));
%!     endif
%!   endfor
%!   if (do_plot)
%!   figure("visible", "off");
%!   plot(omega / (2 * pi) * unit_second^-1, P * unit_watt, "-;P(f);r");
%!   xlabel("f [Hz]");
%!   ylabel("P [W]");
%!   grid on;
%!   grid minor on;
%!   title("radiated sound power versus frequency");
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
