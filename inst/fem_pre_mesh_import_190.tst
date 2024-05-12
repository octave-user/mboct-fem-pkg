## fem_pre_mesh_import.m:190
%!test
%! try
%! ### TEST 190
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
%!   r0 = 200e-3 / unit_meters;
%!   r1 = 300e-3 / unit_meters;
%!   r2 = 310e-3 / unit_meters;
%!   r3 = 400e-3 / unit_meters;
%!   c1 = 1400 / (unit_meters / unit_second);
%!   rho1 = 1000 / (unit_kilograms / unit_meters^3);
%!   eta1 = 0 / (unit_pascal * unit_second);
%!   zeta1 = 0 / (unit_pascal * unit_second);
%!   E2 = 70000e6 / unit_pascal;
%!   rho2 = 1700 / (unit_kilograms / unit_meters^3);
%!   nu2 = 0.3;
%!   K2 = E2 / (3 * (1 - 2 * nu2));
%!   c2 = sqrt(K2 / rho2);
%!   alpha2 = 0 / (unit_second^-1);
%!   beta2 = 0 / unit_second;
%!   c3 = 800 / (unit_meters / unit_second);
%!   rho3 = 500 / (unit_kilograms / unit_meters^3);
%!   eta3 = 0 / (unit_pascal * unit_second);
%!   zeta3 = 0 / (unit_pascal * unit_second);
%!   k1 = 2 * pi / (r1 - r0);
%!   omega = k1 * c1;
%!   k3 = omega / c3;
%!   f = omega / (2 * pi);
%!   lambda1 = c1 / f;
%!   lambda2 = c2 / f;
%!   lambda3 = c3 / f;
%!   dx1 = lambda1 / 20;
%!   dx2 = lambda2 / 20;
%!   dx3 = lambda3 / 20;
%!   dx = min([dx1, dx2, dx3]);
%!   alpha = 2.5 * pi / 180;
%!   v0 = 1e-3 / (unit_meters / unit_second);
%!   ## See Jont Allen
%!   ## THE ACOUSTIC WAVE EQUATION AND SIMPLE SOLUTIONS
%!   ## Chapter 5 spherical waves, equation 5.11.9
%!   z3 = rho3 * c3 * k3 * r3 .* exp(1j * acot(k3 * r3)) ./ sqrt(1 + (k3 * r3).^2);
%!   M = [((k1.*r0-1j).*exp(-1j.*k1.*r0))./(omega.*r0.^2.*rho1), -((k1.*r0+1j).*exp(1j.*k1.*r0))./(omega.*r0.^2.*rho1), 0, 0, -v0;
%!        ((k1.*r1-1j).*exp(-1j.*k1.*r1))./(omega.*r1.^2.*rho1), -((k1.*r1+1j).*exp(1j.*k1.*r1))./(omega.*r1.^2.*rho1), 0, -1j.*omega, 0;
%!        0, 0, ((k3.*r2-1j).*exp(-1j.*k3.*r2))./(omega.*r2.^2.*rho3), -1j.*omega, 0;
%!        -exp(-1j.*k1.*r1)./r1, -exp(1j.*k1.*r1)./r1, exp(-1j.*k3.*r2)./r2, -((r2-r1).*(((nu2-1).*omega.^2.*r2.^2+(2.*nu2-2).*omega.^2.*r1.*r2+(nu2-1).*omega.^2.*r1.^2).*rho2+8.*E2))./((nu2-1).*(r2+r1).^2), 0];
%!   x = -M(:, 1:4) \ M(:, 5);
%!   A1 = x(1);
%!   B1 = x(2);
%!   A3 = x(3);
%!   U2 = x(4);
%!   p1 = @(r, t) (B1*exp(1i*(omega*t+k1*r)))./r+(A1*exp(1i*(omega*t-k1*r)))./r;
%!   p3 = @(r, t) (A3*exp(1i*(omega*t-k3*r)))./r;
%!   v1 = @(r, t) -((B1*k1*r.*exp(2*1i*k1*r)+1i*B1*exp(2*1i*k1*r)-A1*k1*r+1i*A1).*exp(1i*omega*t-1i*k1*r))./(omega*r.^2*rho1);
%!   v3 = @(r, t) (A3*(k3*r-1i).*exp(1i*omega*t-1i*k3*r))./(omega*r.^2*rho3);
%!   v = @(r, t) (r < r1) .* v1(r, t) + (r > r2) .* v3(r, t) + (r >= r1 & r <= r2) .* (1j * omega * U2) * exp(1j * omega * t);
%!   p = @(r, t) (r <= r1) .* p1(r, t) + (r >= r2) .* p3(r, t);
%!   fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!   fprintf(fd, "r0=%g;\n", r0);
%!   fprintf(fd, "r1=%g;\n", r1);
%!   fprintf(fd, "r2=%g;\n", r2);
%!   fprintf(fd, "r3=%g;\n", r3);
%!   fprintf(fd, "alpha=%g;\n", alpha);
%!   fprintf(fd, "dx = %g;\n", dx);
%!   fputs(fd, "Point(1) = {0,0,0,dx};\n");
%!   fputs(fd, "Point(2) = {r0*Cos(-alpha/2),r0*Sin(-alpha/2),0,dx};\n");
%!   fputs(fd, "Point(3) = {r0*Cos(alpha/2),r0*Sin(alpha/2),0,dx};\n");
%!   fputs(fd, "Point(4) = {r1*Cos(-alpha/2),r1*Sin(-alpha/2),0,dx};\n");
%!   fputs(fd, "Point(5) = {r1*Cos(alpha/2),r1*Sin(alpha/2),0,dx};\n");
%!   fputs(fd, "Point(6) = {r2*Cos(-alpha/2),r2*Sin(-alpha/2),0,dx};\n");
%!   fputs(fd, "Point(7) = {r2*Cos(alpha/2),r2*Sin(alpha/2),0,dx};\n");
%!   fputs(fd, "Point(8) = {r3*Cos(-alpha/2),r3*Sin(-alpha/2),0,dx};\n");
%!   fputs(fd, "Point(9) = {r3*Cos(alpha/2),r3*Sin(alpha/2),0,dx};\n");
%!   fputs(fd, "Line(1) = {2,4};\n");
%!   fputs(fd, "Circle(2) = {4,1,5};\n");
%!   fputs(fd, "Line(3) = {5,3};\n");
%!   fputs(fd, "Circle(4) = {3,1,2};\n");
%!   fputs(fd, "Line(5) = {4,6};\n");
%!   fputs(fd, "Circle(6) = {6,1,7};\n");
%!   fputs(fd, "Line(7) = {7,5};\n");
%!   fputs(fd, "Circle(8) = {5,1,4};\n");
%!   fputs(fd, "Line(9) = {6,8};\n");
%!   fputs(fd, "Circle(10) = {8,1,9};\n");
%!   fputs(fd, "Line(11) = {9,7};\n");
%!   fputs(fd, "Circle(12) = {7,1,6};\n");
%!   fputs(fd, "Curve Loop(13) = {1,2,3,4};\n");
%!   fputs(fd, "Curve Loop(14) = {5,6,7,8};\n");
%!   fputs(fd, "Curve Loop(15) = {9,10,11,12};\n");
%!   fputs(fd, "Plane Surface(16) = {13};\n");
%!   fputs(fd, "Plane Surface(17) = {14};\n");
%!   fputs(fd, "Plane Surface(18) = {15};\n");
%!   fputs(fd, "tmp1[] = Extrude{{0,1,0},{0,0,0},-alpha/2} {\n");
%!   fputs(fd, "  Surface{16,17,18}; Layers{Ceil(alpha/2 * r2 / dx)}; Recombine;\n");
%!   fputs(fd, "};\n");
%!   fputs(fd, "tmp2[] = Extrude{{0,1,0},{0,0,0},alpha/2} {\n");
%!   fputs(fd, "  Surface{16,17,18}; Layers{Ceil(alpha/2 * r2 / dx)}; Recombine;\n");
%!   fputs(fd, "};\n");
%!   fputs(fd, "Coherence;\n");
%!   fputs(fd, "ReverseMesh Surface{35,20};\n");
%!   fputs(fd, "Physical Volume(\"v1\",1) = {tmp1[1],tmp2[1]};\n");
%!   fputs(fd, "Physical Volume(\"v2\",2) = {tmp1[7],tmp2[7]};\n");
%!   fputs(fd, "Physical Volume(\"v3\",3) = {tmp1[13],tmp2[13]};\n");
%!   fputs(fd, "Physical Surface(\"s1\", 4) = {37, 22};\n");
%!   fputs(fd, "Physical Surface(\"s2\", 5) = {35, 20};\n");
%!   fputs(fd, "Physical Surface(\"s3\", 6) = {50, 40};\n");
%!   fputs(fd, "Physical Surface(\"s4\", 7) = {54, 45};\n");
%!   fputs(fd, "Physical Surface(\"s5\", 8) = {49, 39};\n");
%!   fputs(fd, "Physical Surface(\"s6\", 9) = {41, 51};\n");
%!   fputs(fd, "Physical Surface(\"s7\", 10) = {52};\n");
%!   fputs(fd, "Physical Surface(\"s8\", 11) = {43};\n");
%!   fputs(fd, "Mesh.OptimizeThreshold=0.9;\n");
%!   fputs(fd, "Mesh.HighOrderOptimize=2;\n");
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
%!   opts.elem_type = {"tria6h", "quad8", "penta15"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opts));
%!   [~] = unlink([filename, ".msh"]);
%!   grp_idx_v1 = find([mesh.groups.penta15.id] == 1);
%!   grp_idx_v2 = find([mesh.groups.penta15.id] == 2);
%!   grp_idx_v3 = find([mesh.groups.penta15.id] == 3);
%!   grp_idx_s1 = find([mesh.groups.quad8.id] == 4);
%!   grp_idx_s2 = find([mesh.groups.quad8.id] == 5);
%!   grp_idx_s3 = find([mesh.groups.quad8.id] == 6);
%!   grp_idx_s4 = find([mesh.groups.quad8.id] == 7);
%!   grp_idx_s5 = find([mesh.groups.quad8.id] == 8);
%!   grp_idx_s6 = find([mesh.groups.quad8.id] == 9);
%!   grp_idx_s7 = find([mesh.groups.tria6h.id] == 10);
%!   grp_idx_s8 = find([mesh.groups.tria6h.id] == 11);
%!   mesh.materials.penta15 = zeros(rows(mesh.elements.penta15), 1, "int32");
%!   mesh.materials.penta15(mesh.groups.penta15(grp_idx_v1).elements) = 1;
%!   mesh.materials.penta15(mesh.groups.penta15(grp_idx_v2).elements) = 2;
%!   mesh.materials.penta15(mesh.groups.penta15(grp_idx_v3).elements) = 3;
%!   empty_cell = cell(1, 3);
%!   mesh.material_data = struct("E", empty_cell, ...
%!                               "nu", empty_cell, ...
%!                               "rho", empty_cell, ...
%!                               "c", empty_cell);
%!   mesh.material_data(1).c = c1;
%!   mesh.material_data(1).rho = rho1;
%!   mesh.material_data(2).E = E2;
%!   mesh.material_data(2).nu = nu2;
%!   mesh.material_data(2).rho = rho2;
%!   mesh.material_data(3).c = c3;
%!   mesh.material_data(3).rho = rho3;
%!   mesh.elements.fluid_struct_interface.quad8 = mesh.elements.quad8([[mesh.groups.quad8([grp_idx_s2, grp_idx_s3])].elements], :);
%!   mesh.elements.particle_velocity.quad8.nodes = mesh.elements.quad8(mesh.groups.quad8(grp_idx_s1).elements, :);
%!   mesh.materials.particle_velocity.quad8 = repmat(int32(1), rows(mesh.elements.particle_velocity.quad8.nodes), 1);
%!   load_case = struct("particle_velocity", cell(1, 2));
%!   load_case(1).particle_velocity.quad8.vn = repmat(-real(v0), size(mesh.elements.particle_velocity.quad8.nodes));
%!   load_case(2).particle_velocity.quad8.vn = repmat(-imag(v0), size(mesh.elements.particle_velocity.quad8.nodes));
%!   mesh.elements.acoustic_impedance.quad8.nodes = mesh.elements.quad8(mesh.groups.quad8(grp_idx_s4).elements, :);
%!   mesh.elements.acoustic_impedance.quad8.z = repmat(z3, size(mesh.elements.acoustic_impedance.quad8.nodes));
%!   mesh.materials.acoustic_impedance.quad8 = repmat(int32(3), rows(mesh.elements.acoustic_impedance.quad8.nodes), 1);
%!   elemid = [[mesh.groups.quad8([grp_idx_s5,grp_idx_s6])].elements];
%!   nodeid = mesh.elements.quad8(elemid, :);
%!   elemid_2 = [[mesh.groups.tria6h([grp_idx_s7,grp_idx_s8])].elements];
%!   nodeid_2 = mesh.elements.tria6h(elemid_2, :);
%!   C = zeros(numel(nodeid) + numel(nodeid_2), 6);
%!   nodes = zeros(numel(nodeid) + numel(nodeid_2), 1, "int32");
%!   idxn1 = int32([2, 3, 4, 1, 2, 3, 4, 1]);
%!   idxn2 = int32([4, 1, 2, 3, 7, 8, 5, 6]);
%!   idxn1_2 = int32([4,5,6,2,3,1]);
%!   idxn2_2 = int32([6,4,5,6,1,5]);
%!   k = int32(0);
%!   jointidx = zeros(rows(mesh.nodes), 4, "int32");
%!   tolortho = 1e-3;
%!   for i=1:rows(nodeid)
%!     X = mesh.nodes(nodeid(i, :), 1:3).';
%!     for j=1:columns(nodeid)
%!       n1 = X(:, idxn1(j)) - X(:, j);
%!       n2 = X(:, idxn2(j)) - X(:, j);
%!       n3 = cross(n1, n2);
%!       n3 /= norm(n3);
%!       fduplicate = false;
%!       for l=1:columns(jointidx)
%!         jointidxprev = jointidx(nodeid(i, j), l);
%!         if (~jointidxprev)
%!           continue;
%!         endif
%!         n4 = C(jointidxprev, 1:3).';
%!         n5 = cross(n3, n4);
%!         if (norm(n5) > tolortho)
%!           n3 = cross(n4, n5);
%!           n3 /= norm(n3);
%!         else
%!           fduplicate = true;
%!         endif
%!       endfor
%!       if (fduplicate)
%!         continue;
%!       endif
%!       C(++k, 1:3) = n3;
%!       nodes(k) = nodeid(i, j);
%!       for l=1:columns(jointidx)
%!         if (~jointidx(nodeid(i, j), l))
%!           jointidx(nodeid(i, j), l) = k;
%!           break;
%!         endif
%!       endfor
%!     endfor
%!   endfor
%!   for i=1:rows(nodeid_2)
%!     X = mesh.nodes(nodeid_2(i, :), 1:3).';
%!     for j=1:columns(nodeid_2)
%!       n1 = X(:, idxn1_2(j)) - X(:, j);
%!       n2 = X(:, idxn2_2(j)) - X(:, j);
%!       n3 = cross(n1, n2);
%!       n3 /= norm(n3);
%!       fduplicate = false;
%!       for l=1:columns(jointidx)
%!         jointidxprev = jointidx(nodeid_2(i, j), l);
%!         if (~jointidxprev)
%!           continue;
%!         endif
%!         n4 = C(jointidxprev, 1:3).';
%!         n5 = cross(n3, n4);
%!         if (norm(n5) > tolortho)
%!           n3 = cross(n4, n5);
%!           n3 /= norm(n3);
%!         else
%!           fduplicate = true;
%!         endif
%!       endfor
%!       if (fduplicate)
%!         continue;
%!       endif
%!       C(++k, 1:3) = n3;
%!       nodes(k) = nodeid_2(i, j);
%!       for l=1:columns(jointidx)
%!         if (~jointidx(nodeid_2(i, j), l))
%!           jointidx(nodeid_2(i, j), l) = k;
%!           break;
%!         endif
%!       endfor
%!     endfor
%!   endfor
%!   C = C(1:k, :);
%!   nodes = nodes(1:k);
%!   mesh.elements.joints = struct("C", mat2cell(C, ones(rows(C), 1), 6), "nodes", mat2cell(nodes, ones(numel(nodes), 1)));
%!   load_case_dof.locked_dof = false(rows(mesh.nodes), 7);
%!   load_case_dof.domain = FEM_DO_FLUID_STRUCT;
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
%!   idxPhi = dof_map.ndof(:, 7);
%!   idxPhi1 = find(idxPhi > 0);
%!   idxPhi = idxPhi(idxPhi1);
%!   opt_sol.number_of_threads = mbdyn_solver_num_threads_default();
%!   opt_sol.solver = "pastix";
%!   opt_sol.refine_max_iter = int32(50);
%!   opt_sol.verbose = int32(0);
%!   Keff = -omega^2 * mat_ass.Mfs + 1j * omega * complex(mat_ass.Dfs_re, mat_ass.Dfs_im) + mat_ass.Kfs;
%!   Reff = complex(mat_ass.Rfs(:, 1), mat_ass.Rfs(:, 2));
%!   Z = fem_sol_factor(Keff, opt_sol) \ Reff;
%!   Psi = linspace(0, 2 * pi, 37);
%!   sol.t = Psi / omega;
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
%!   vr = zeros(rows(mesh.nodes), numel(sol.t));
%!   n = mesh.nodes(:, 1:3);
%!   n = diag(1 ./ norm(n, "rows")) * n;
%!   for i=1:columns(mesh.elements.penta15)
%!     for k=1:numel(sol.t)
%!       vr(mesh.elements.penta15(:, i), k) = 0;
%!       for j=1:3
%!         vr(mesh.elements.penta15(:, i), k) += n(mesh.elements.penta15(:, i), j) .* sol.particle_velocity.v.penta15(:, i, j, k);
%!       endfor
%!     endfor
%!   endfor
%!   for j=1:numel(sol.t)
%!     vr(mesh.groups.penta15(grp_idx_v2).nodes, j) = 0;
%!     for i=1:3
%!       vr(mesh.groups.penta15(grp_idx_v2).nodes, j) += n(mesh.groups.penta15(grp_idx_v2).nodes, i) .* sol.vel(mesh.groups.penta15(grp_idx_v2).nodes, i, j);
%!     endfor
%!   endfor
%!   node_idx = mesh.groups.penta15(grp_idx_v2).nodes;
%!   tol = 1e-5;
%!   [r, idx] = sort(norm(mesh.nodes(:, 1:3), "rows"));
%!   vref = real(v(r, sol.t));
%!   pref = real(p(r, sol.t));
%!   for i=1:numel(Psi)
%!     figure("visible", "off");
%!     hold on;
%!     plot(r * unit_meters, vr(idx, i) * unit_meters / unit_second, "-;vr;r");
%!     plot(r * unit_meters, vref(:, i) * unit_meters / unit_second, "-;vref;k");
%!     xlabel("r [m]");
%!     ylabel("v [m/s]");
%!     grid on;
%!     grid minor on;
%!     title(sprintf("velocity distribution Psi=%.1fdeg", Psi(i) * 180 / pi));
%!   endfor
%!   for i=1:numel(Psi)
%!     figure("visible", "off");
%!     hold on;
%!     plot(r * unit_meters, sol.p(idx, i) * unit_pascal, "-;p;r");
%!     plot(r * unit_meters, pref(:, i) * unit_pascal, "-;pref;k");
%!     xlabel("r [m]");
%!     ylabel("p [Pa]");
%!     grid on;
%!     grid minor on;
%!     title(sprintf("pressure distribution Psi=%.1fdeg", Psi(i) * 180 / pi));
%!   endfor
%!   tol = 2e-2;
%!   idx2 = find(r < r1 | r > r2);
%!   assert_simple(max(max(abs(vr(idx, :) - vref))) < tol * max(max(abs(vref))));
%!   assert_simple(max(max(abs(sol.p(idx(idx2), :) - pref(idx2,:)))) < tol * max(max(abs(pref))));
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
