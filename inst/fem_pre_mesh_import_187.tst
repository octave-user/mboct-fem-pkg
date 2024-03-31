## fem_pre_mesh_import.m:187
%!test
%! ### TEST 187
%! ### 1D wave propagation with two fluid domains and a solid domain between them
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
%!   unit_second = 1e3;
%!   unit_kilograms = 1e3;
%!   unit_newton = unit_kilograms * unit_meters / unit_second^2;
%!   unit_pascal = unit_newton / unit_meters^2;
%!   l1 = 100e-3 / unit_meters;
%!   l2 = 10e-3 / unit_meters;
%!   l3 = 100e-3 / unit_meters;
%!   c1 = 1400 / (unit_meters / unit_second);
%!   rho1 = 1000 / (unit_kilograms / unit_meters^3);
%!   eta1 = 0 / (unit_pascal * unit_second);
%!   zeta1 = 0 / (unit_pascal * unit_second);
%!   E2 = 80000e6 / unit_pascal;
%!   rho2 = 1700 / (unit_kilograms / unit_meters^3);
%!   nu2 = 0.3;
%!   K2 = E2 / (3 * (1 - 2 * nu2));
%!   c2 = sqrt(K2 / rho2);
%!   alpha2 = 0 / (unit_second^-1);
%!   beta2 = 0 / unit_second;
%!   c3 = 220 / (unit_meters / unit_second);
%!   rho3 = 2.1 / (unit_kilograms / unit_meters^3);
%!   eta3 = 0 / (unit_pascal * unit_second);
%!   zeta3 = 0 / (unit_pascal * unit_second);
%!   k1 = 2 * pi / l1;
%!   omega = k1 * c1;
%!   k3 = omega / c3;
%!   f = omega / (2 * pi);
%!   lambda1 = c1 / f;
%!   lambda2 = c2 / f;
%!   lambda3 = c3 / f;
%!   dx1 = lambda1 / 150;
%!   dx2 = lambda2 / 1000;
%!   dx3 = lambda3 / 150;
%!   dx = mean([dx1, dx3]);
%!   w = dx;
%!   h = dx;
%!   vx0 = 1e-6 * (1 + 2j) / (unit_meters / unit_second);
%!   pI = (c1*e^((2*1i*l1*omega)/c1)*rho1*(c3*rho3+c1*rho1+1i*l2*omega*rho2)*vx0)/(c3*e^((2*1i*l1*omega)/c1)*rho3-c3*rho3+c1*e^((2*1i*l1*omega)/c1)*rho1+c1*rho1+1i*l2*omega*e^((2*1i*l1*omega)/c1)*rho2-1i*l2*omega*rho2);
%!   pR = (c1*rho1*(c3*rho3-c1*rho1+1i*l2*omega*rho2)*vx0)/(c3*e^((2*1i*l1*omega)/c1)*rho3-c3*rho3+c1*e^((2*1i*l1*omega)/c1)*rho1+c1*rho1+1i*l2*omega*e^((2*1i*l1*omega)/c1)*rho2-1i*l2*omega*rho2);
%!   pT = (2*c1*c3*e^((1i*l2*omega)/c3+(1i*l1*omega)/c3+(1i*l1*omega)/c1)*rho1*rho3*vx0)/(c3*e^((2*1i*l1*omega)/c1)*rho3-c3*rho3+c1*e^((2*1i*l1*omega)/c1)*rho1+c1*rho1+1i*l2*omega*e^((2*1i*l1*omega)/c1)*rho2-1i*l2*omega*rho2);
%!   U =-(2*1i*c1*e^((1i*l1*omega)/c1)*rho1*vx0)/(omega*(c3*e^((2*1i*l1*omega)/c1)*rho3-c3*rho3+c1*e^((2*1i*l1*omega)/c1)*rho1+c1*rho1+1i*l2*omega*e^((2*1i*l1*omega)/c1)*rho2-1i*l2*omega*rho2));
%!   p1 = @(x, t) pR*exp(1i*((omega*x)/c1+omega*t))+pI*exp(1i*(omega*t-(omega*x)/c1));
%!   p3 = @(x, t) pT*exp(1i*(omega*t-(omega*x)/c3));
%!   p = @(x, t) (p1(x, t) .* (x <= l1)) + (p3(x, t) .* (x >= l1 + l2));
%!   vx1 = @(x, t) ((pI*exp(1i*(omega*t-(omega*x)/c1)))/c1-(pR*exp(1i*((omega*x)/c1+omega*t)))/c1)/rho1;
%!   vx3 = @(x, t) (pT*exp(1i*(omega*t-(omega*x)/c3)))/(c3*rho3);
%!   v = @(x, t) (vx1(x, t) .* (x <= l1)) + (vx3(x, t) .* (x >= l1 + l2)) + (1j * omega * U * exp(1i * omega * t)) .* ((x > l1) & (x < l1 + l2));
%!   fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!   fprintf(fd, "l1=%g;\n", l1);
%!   fprintf(fd, "l2=%g;\n", l2);
%!   fprintf(fd, "l3=%g;\n", l3);
%!   fprintf(fd, "w=%g;\n", w);
%!   fprintf(fd, "h=%g;\n", h);
%!   fprintf(fd, "dx1 = %g;\n", dx1);
%!   fprintf(fd, "dx2 = %g;\n", dx2);
%!   fprintf(fd, "dx3 = %g;\n", dx3);
%!   fputs(fd, "Point(1) = {0,0,0,w};\n");
%!   fputs(fd, "Point(2) = {0,w,0,w};\n");
%!   fputs(fd, "Point(3) = {0,w,h,w};\n");
%!   fputs(fd, "Point(4) = {0,0,h,w};\n");
%!   fputs(fd, "Line(1) = {1,2};\n");
%!   fputs(fd, "Line(2) = {2,3};\n");
%!   fputs(fd, "Line(3) = {3,4};\n");
%!   fputs(fd, "Line(4) = {4,1};\n");
%!   fputs(fd, "Line Loop(5) = {1,2,3,4};\n");
%!   fputs(fd, "Plane Surface(6) = {5};\n");
%!   fputs(fd, "tmp1[] = Extrude {l1,0,0} {\n");
%!   fputs(fd, "  Surface{6};\n");
%!   fputs(fd, "};\n");
%!   fputs(fd, "tmp2[] = Extrude {l2,0,0} {\n");
%!   fputs(fd, "  Surface{tmp1[0]};\n");
%!   fputs(fd, "};\n");
%!   fputs(fd, "tmp3[] = Extrude {l3,0,0} {\n");
%!   fputs(fd, "  Surface{tmp2[0]};\n");
%!   fputs(fd, "};\n");
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
%!   fputs(fd, "ReverseMesh Surface{6};\n");
%!   fputs(fd, "ReverseMesh Surface{tmp1[0]};\n");
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
%!   opt_mesh.elem_type = {"tet10h", "tria6h"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_mesh));
%!   [~] = unlink([filename, ".msh"]);
%!   grp_idx_input = find([mesh.groups.tria6h.id] == 1);
%!   grp_idx_output = find([mesh.groups.tria6h.id] == 2);
%!   grp_idx_fsi1 = find([mesh.groups.tria6h.id] == 3);
%!   grp_idx_fsi2 = find([mesh.groups.tria6h.id] == 4);
%!   grp_idx_volume1 = find([mesh.groups.tet10h.id] == 5);
%!   grp_idx_volume2 = find([mesh.groups.tet10h.id] == 6);
%!   grp_idx_volume3 = find([mesh.groups.tet10h.id] == 7);
%!   grp_idx_slider1 = find([mesh.groups.tria6h.id] == 8);
%!   grp_idx_slider2 = find([mesh.groups.tria6h.id] == 9);
%!   load_case_dof.locked_dof = false(rows(mesh.nodes), 7);
%!   load_case_dof.locked_dof(mesh.groups.tria6h(grp_idx_slider1).nodes, 3) = true;
%!   load_case_dof.locked_dof(mesh.groups.tria6h(grp_idx_slider2).nodes, 2) = true;
%!   load_case_dof.domain = FEM_DO_FLUID_STRUCT;
%!   mesh.materials.tet10h = zeros(rows(mesh.elements.tet10h), 1, "int32");
%!   mesh.materials.tet10h(mesh.groups.tet10h(grp_idx_volume1).elements) = 1;
%!   mesh.materials.tet10h(mesh.groups.tet10h(grp_idx_volume2).elements) = 2;
%!   mesh.materials.tet10h(mesh.groups.tet10h(grp_idx_volume3).elements) = 3;
%!   mesh.elements.acoustic_boundary.tria6h = mesh.elements.tria6h([[mesh.groups.tria6h([grp_idx_input, grp_idx_output])].elements], :);
%!   mesh.materials.acoustic_boundary.tria6h = zeros(rows(mesh.elements.acoustic_boundary.tria6h), 1, "int32");
%!   mesh.materials.acoustic_boundary.tria6h(1:numel(mesh.groups.tria6h(grp_idx_input).elements)) = 1;
%!   mesh.materials.acoustic_boundary.tria6h(numel(mesh.groups.tria6h(grp_idx_output).elements)+1:end) = 3;
%!   mesh.elements.particle_velocity.tria6h.nodes = mesh.elements.tria6h([[mesh.groups.tria6h([grp_idx_input, grp_idx_output])].elements], :);
%!   mesh.materials.particle_velocity.tria6h = zeros(rows(mesh.elements.particle_velocity.tria6h.nodes), 1, "int32");
%!   mesh.materials.particle_velocity.tria6h(1:numel(mesh.groups.tria6h(grp_idx_input).elements)) = int32(1);
%!   mesh.materials.particle_velocity.tria6h(numel(mesh.groups.tria6h(grp_idx_input).elements) + 1:end) = int32(3);
%!   mesh.elements.acoustic_impedance.tria6h.nodes = mesh.elements.tria6h(mesh.groups.tria6h(grp_idx_output).elements, :);
%!   mesh.elements.acoustic_impedance.tria6h.z = repmat(rho3 * c3, size(mesh.elements.acoustic_impedance.tria6h.nodes));
%!   mesh.materials.acoustic_impedance.tria6h = repmat(int32(3), rows(mesh.elements.acoustic_impedance.tria6h.nodes), 1);
%!   mesh.elements.fluid_struct_interface.tria6h = mesh.elements.tria6h([[mesh.groups.tria6h([grp_idx_fsi1, grp_idx_fsi2])].elements], :);
%!   load_case = struct("particle_velocity", cell(1, 2));
%!   for i=1:numel(load_case)
%!     load_case(i).particle_velocity.tria6h.vn = zeros(size(mesh.elements.particle_velocity.tria6h.nodes));
%!   endfor
%!   load_case(1).particle_velocity.tria6h.vn(1:numel(mesh.groups.tria6h(grp_idx_input).elements), :) = -real(vx0);
%!   load_case(2).particle_velocity.tria6h.vn(1:numel(mesh.groups.tria6h(grp_idx_input).elements), :) = -imag(vx0);
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
%!   opt_sol.number_of_threads = mbdyn_solver_num_threads_default();
%!   opt_sol.solver = "pastix";
%!   opt_sol.refine_max_iter = int32(50);
%!   opt_sol.verbose = int32(0);
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
%!   sol.def = sol.vel = zeros(rows(mesh.nodes), 6, numel(sol.t));
%!   sol.p(idxPhi1, :) = real(-1j * omega * Z(idxPhi, :) .* exp(1j * omega * sol.t));
%!   sol.Phi = sol.PhiP = zeros(rows(mesh.nodes), numel(sol.t));
%!   sol.Phi(idxPhi1, :) = real(Z(idxPhi, :) .* exp(1j * omega * sol.t));
%!   sol.PhiP(idxPhi1, :) = real(1j * omega * Z(idxPhi, :) .* exp(1j * omega * sol.t));
%!   [sol.particle_velocity, ...
%!    sol.acoustic_intensity] = fem_ass_matrix(mesh, ...
%!                                             dof_map, ...
%!                                             [FEM_VEC_PARTICLE_VELOCITY, ...
%!                                              FEM_SCA_ACOUSTIC_INTENSITY], ...
%!                                             load_case, ...
%!                                             sol);
%!   for j=1:6
%!     idxU = dof_map.ndof(:, j);
%!     idxU1 = find(idxU > 0);
%!     idxU = idxU(idxU1);
%!     sol.def(idxU1, j, :) = real(Z(idxU, :) * exp(1j * omega * sol.t));
%!     sol.vel(idxU1, j, :) = real(1j * omega * Z(idxU, :) * exp(1j * omega * sol.t));
%!   endfor
%!   vx = zeros(rows(mesh.nodes), numel(sol.t));
%!   for i=1:columns(mesh.elements.tet10h)
%!     vx(mesh.elements.tet10h(:, i), :) = reshape(sol.particle_velocity.v.tet10h(:, i, 1, :), rows(mesh.elements.tet10h), numel(sol.t));
%!   endfor
%!   vx(mesh.groups.tet10h(grp_idx_volume2).nodes, :) = sol.vel(mesh.groups.tet10h(grp_idx_volume2).nodes, 1, :);
%!   node_idx = mesh.groups.tet10h(grp_idx_volume2).nodes;
%!   tol = 1e-5;
%!   [~, idx] = sort(mesh.nodes(:, 1));
%!   pref = real(p(mesh.nodes(idx, 1), sol.t));
%!   vxref = real(v(mesh.nodes(idx, 1), sol.t));
%!   Uref = real(U * exp(1j * omega * sol.t));
%!   for i=1:numel(Psi)
%!     figure("visible", "off");
%!     hold on;
%!     plot(mesh.nodes(idx, 1) * unit_meters, sol.p(idx, i) * unit_pascal, "-;p;r");
%!     plot(mesh.nodes(idx, 1) * unit_meters, pref(:, i) * unit_pascal, "-;pref;k");
%!     ylim([min(min(sol.p)), max(max(sol.p))] * unit_pascal);
%!     xlabel("x [m]");
%!     ylabel("p [Pa]");
%!     grid on;
%!     grid minor on;
%!     title(sprintf("pressure distribution Psi=%.1fdeg", Psi(i) * 180 / pi));
%!   endfor
%!   for i=1:numel(Psi)
%!     figure("visible", "off");
%!     hold on;
%!     plot(mesh.nodes(idx, 1) * unit_meters, vx(idx, i) * unit_meters / unit_second, "-;vx;r");
%!     plot(mesh.nodes(idx, 1) * unit_meters, vxref(:, i) * unit_meters / unit_second, "-;vxref;k");
%!     ylim([min(min(vxref)), max(max(vxref))] * unit_meters / unit_second);
%!     xlabel("x [m]");
%!     ylabel("vx [m/s]");
%!     grid on;
%!     grid minor on;
%!     title(sprintf("velocity distribution Psi=%.1fdeg", Psi(i) * 180 / pi));
%!   endfor
%!   figure("visible", "off");
%!   hold on;
%!   plot(sol.t * unit_second, reshape(mean(sol.def(node_idx, 1, :), 1), 1, numel(sol.t)) * unit_meters, "-;U;r");
%!   plot(sol.t * unit_second, Uref * unit_meters, "-;Uref;k");
%!   xlabel("t [s]");
%!   ylabel("Ux [m]");
%!   grid on;
%!   grid minor on;
%!   title("displacement of solid domain");
%!   tol = 1e-2;
%!   assert_simple(vx(idx, :), vxref, tol * max(max(abs(vxref))));
%!   assert_simple(sol.p(idx, :), pref, tol * max(max(abs(pref))));
%!   assert_simple(reshape(mean(sol.def(node_idx, 1, :), 1), 1, numel(sol.t)), Uref, tol * max(abs(Uref)));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       [~] = unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
