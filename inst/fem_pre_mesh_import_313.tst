## fem_pre_mesh_import.m:313
%!test
%! ### TEST 313
%! ### rigid body motion of a solid domain between two fluid domains
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
%!     unit_meters = 1e-3;
%!     unit_second = 1;
%!     unit_kilograms = 1e1;
%!     unit_newton = unit_kilograms * unit_meters / unit_second^2;
%!     unit_pascal = unit_newton / unit_meters^2;
%!     g = 9.81 * unit_meters / unit_second^2;
%!     l1 = 100e-3 / unit_meters;
%!     l0 = 10e-3 / unit_meters;
%!     l2 = 100e-3 / unit_meters;
%!     c = 1400 / (unit_meters / unit_second);
%!     rhof = 1000 / (unit_kilograms / unit_meters^3);
%!     eta = 0 / (unit_pascal * unit_second);
%!     zeta = 0 / (unit_pascal * unit_second);
%!     E = 70000e6 / unit_pascal;
%!     rhos = 2700 / (unit_kilograms / unit_meters^3);
%!     nu = 0.3;
%!     alpha = 0 / (unit_second^-1);
%!     beta = 0 / unit_second;
%!     dx = 2e-3 / unit_meters;
%!     f = 10000 / (unit_second^-1);
%!     omega = 2 * pi * f;
%!     k = omega / c;
%!     lambda = c / f;
%!     dx = min(dx, lambda / 6);
%!     w = dx;
%!     h = dx;
%!     A1 = (1j*g*l0*e^(1j*k*l1)*(e^(1j*k*l2)-1)*(e^(1j*k*l2)+1)*rhof*rhos)/(k*l0*e^(2*1j*k*l2+2*1j*k*l1)*rhos-k*l0*e^(2*1j*k*l2)*rhos-k*l0*e^(2*1j*k*l1)*rhos+k*l0*rhos-2*1j*e^(2*1j*k*l2+2*1j*k*l1)*rhof+2*1j*rhof);
%!     B1 = (1j*g*l0*e^(1j*k*l1)*(e^(1j*k*l2)-1)*(e^(1j*k*l2)+1)*rhof*rhos)/(k*l0*e^(2*1j*k*l2+2*1j*k*l1)*rhos-k*l0*e^(2*1j*k*l2)*rhos-k*l0*e^(2*1j*k*l1)*rhos+k*l0*rhos-2*1j*e^(2*1j*k*l2+2*1j*k*l1)*rhof+2*1j*rhof);
%!     A2 = -(1j*g*l0*(e^(1j*k*l1)-1)*(e^(1j*k*l1)+1)*e^(2*1j*k*l2+1j*k*l1+1j*k*l0)*rhof*rhos)/(k*l0*e^(2*1j*k*l2+2*1j*k*l1)*rhos-k*l0*e^(2*1j*k*l2)*rhos-k*l0*e^(2*1j*k*l1)*rhos+k*l0*rhos-2*1j*e^(2*1j*k*l2+2*1j*k*l1)*rhof+2*1j*rhof);
%!     B2 = -(1j*g*l0*(e^(1j*k*l1)-1)*(e^(1j*k*l1)+1)*e^((-1j*k*l1)-1j*k*l0)*rhof*rhos)/(k*l0*e^(2*1j*k*l2+2*1j*k*l1)*rhos-k*l0*e^(2*1j*k*l2)*rhos-k*l0*e^(2*1j*k*l1)*rhos+k*l0*rhos-2*1j*e^(2*1j*k*l2+2*1j*k*l1)*rhof+2*1j*rhof);
%!     U = -(g*k*l0*(e^(1j*k*l1)-1)*(e^(1j*k*l1)+1)*(e^(1j*k*l2)-1)*(e^(1j*k*l2)+1)*rhos)/(omega^2*(k*l0*e^(2*1j*k*l2+2*1j*k*l1)*rhos-k*l0*e^(2*1j*k*l2)*rhos-k*l0*e^(2*1j*k*l1)*rhos+k*l0*rhos-2*1j*e^(2*1j*k*l2+2*1j*k*l1)*rhof+2*1j*rhof));
%!     fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "l0=%g;\n", l0);
%!     fprintf(fd, "l1=%g;\n", l1);
%!     fprintf(fd, "l2=%g;\n", l2);
%!     fprintf(fd, "w=%g;\n", w);
%!     fprintf(fd, "h=%g;\n", h);
%!     fprintf(fd, "dx = %g;\n", dx);
%!     fputs(fd, "Point(1) = {0, -0.5 * w, -0.5 * h};\n");
%!     fputs(fd, "Point(2) = {0,  0.5 * w, -0.5 * h};\n");
%!     fputs(fd, "Point(3) = {0,  0.5 * w,  0.5 * h};\n");
%!     fputs(fd, "Point(4) = {0, -0.5 * w,  0.5 * h};\n");
%!     fputs(fd, "Line(1) = {1,2};\n");
%!     fputs(fd, "Line(2) = {2,3};\n");
%!     fputs(fd, "Line(3) = {3,4};\n");
%!     fputs(fd, "Line(4) = {4,1};\n");
%!     fputs(fd, "Line Loop(5) = {1,2,3,4};\n");
%!     fputs(fd, "Plane Surface(6) = {5};\n");
%!     fputs(fd, "tmp1[] = Extrude {l1,0,0} {\n");
%!     fputs(fd, "  Surface{6}; Layers{Ceil(l1/dx)}; Recombine;\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "Recombine Surface{6, tmp1[0]};\n");
%!     fputs(fd, "tmp2[] = Extrude {l0,0,0} {\n");
%!     fputs(fd, "  Surface{tmp1[0]}; Layers{Ceil(l0/dx)}; Recombine;\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "Recombine Surface{tmp1[0], tmp2[0]};\n");
%!     fputs(fd, "tmp3[] = Extrude {l2,0,0} {\n");
%!     fputs(fd, "  Surface{tmp2[0]}; Layers{Ceil(l2/dx)}; Recombine;\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "Recombine Surface{tmp2[0], tmp3[0]};\n");
%!     fputs(fd, "Physical Surface(\"input\",1) = {6};\n");
%!     fputs(fd, "Physical Surface(\"output\",2) = {tmp3[0]};\n");
%!     fputs(fd, "Physical Surface(\"fsi-bnd1\",3) = {tmp1[0]};\n");
%!     fputs(fd, "Physical Surface(\"fsi-bnd2\",4) = {tmp2[0]};\n");
%!     fputs(fd, "Physical Surface(\"slider1\",8) = {tmp2[2]};\n");
%!     fputs(fd, "Physical Surface(\"slider2\",9) = {tmp2[5]};\n");
%!     fputs(fd, "Physical Volume(\"volume1\",5) = {tmp1[1]};\n");
%!     fputs(fd, "Physical Volume(\"volume2\",6) = {tmp2[1]};\n");
%!     fputs(fd, "Physical Volume(\"volume3\",7) = {tmp3[1]};\n");
%!     fputs(fd, "Mesh.OptimizeThreshold=0.99;\n");
%!     fputs(fd, "ReorientMesh Volume{tmp2[1]};\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete=1;\n");
%!     fputs(fd, "Mesh.HighOrderOptimize=2;\n");
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
%!   mesh.elements.fluid_struct_interface.quad8 = mesh.elements.quad8([[mesh.groups.quad8([grp_idx_fsi1, grp_idx_fsi2])].elements], :);
%!   mesh.material_data = struct("E", {[], E, []}, ...
%!                               "rho", {rhof, rhos, rhof}, ...
%!                               "nu", {[], nu, []}, ...
%!                               "c", {c, [], c}, ...
%!                               "eta", {eta, [], eta}, ...
%!                               "zeta", {zeta, [], zeta}, ...
%!                               "alpha", {[], alpha, []}, ...
%!                               "beta", {[], beta, []});
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!   load_case.g = [g; 0; 0];
%!   [mat_ass.Mfs_re, ...
%!    mat_ass.Mfs_im, ...
%!    mat_ass.Dfs_re, ...
%!    mat_ass.Dfs_im, ...
%!    mat_ass.Kfs_re, ...
%!    mat_ass.Kfs_im, ...
%!    mat_ass.Rfs, ...
%!    mat_ass.mat_info, ...
%!    mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                        dof_map, ...
%!                                        [FEM_MAT_MASS_FLUID_STRUCT_RE, ...
%!                                         FEM_MAT_MASS_FLUID_STRUCT_IM, ...
%!                                         FEM_MAT_DAMPING_FLUID_STRUCT_RE, ...
%!                                         FEM_MAT_DAMPING_FLUID_STRUCT_IM, ...
%!                                         FEM_MAT_STIFFNESS_FLUID_STRUCT_RE, ...
%!                                         FEM_MAT_STIFFNESS_FLUID_STRUCT_IM, ...
%!                                         FEM_VEC_LOAD_FLUID_STRUCT], ...
%!                                        load_case);
%!   opt_sol.number_of_threads = mbdyn_solver_num_threads_default();
%!   opt_sol.refine_max_iter = int32(250);
%!   opt_sol.solver = "pastix";
%!   opt_sol.pre_scaling = true;
%!   opt_sol.verbose = int32(0);
%!   Keff = -omega^2 * complex(mat_ass.Mfs_re, mat_ass.Mfs_im) + 1j * omega * complex(mat_ass.Dfs_re, mat_ass.Dfs_im) + complex(mat_ass.Kfs_re, mat_ass.Kfs_im);
%!   Reff = complex(mat_ass.Rfs);
%!   Phi = fem_sol_factor(Keff, opt_sol) \ Reff;
%!   sol.t = linspace(0, 2 * pi, 18) / omega;
%!   idxp = dof_map.ndof(:, 7);
%!   idxp1 = find(idxp > 0);
%!   idxp = idxp(idxp1);
%!   sol.p = zeros(rows(dof_map.ndof), numel(sol.t));
%!   sol.def = zeros(rows(dof_map.ndof), columns(dof_map.ndof), numel(sol.t));
%!   sol.p(idxp1, :) = real(-1j * omega * Phi(idxp, :) .* exp(1j * omega * sol.t));
%!   for j=1:6
%!     idxU = dof_map.ndof(:, j);
%!     idxU1 = find(idxU > 0);
%!     idxU = idxU(idxU1);
%!     sol.def(idxU1, j, :) = real(Phi(idxU, :) .* exp(1j * omega * sol.t));
%!   endfor
%!   p1 = @(x, t) A1 * exp(1j * (omega * t - k * x)) + B1 * exp(1j * (omega * t + k * x));
%!   p2 = @(x, t) A2 * exp(1j * (omega * t - k * x)) + B2 * exp(1j * (omega * t + k * x));
%!   x1 = linspace(0, l1, 36);
%!   x2 = linspace(l1 + l0, l1 + l0 + l2, 36);
%!   tol = 2e-3;
%!   for i=1:numel(sol.t)
%!     p1ref = real(p1(x1, sol.t(i)));
%!     p2ref = real(p2(x2, sol.t(i)));
%!     p1n = griddata3(mesh.nodes(:, 1), mesh.nodes(:, 2), mesh.nodes(:, 3), sol.p(:, i), x1, zeros(size(x1)), zeros(size(x1)));
%!     p2n = griddata3(mesh.nodes(:, 1), mesh.nodes(:, 2), mesh.nodes(:, 3), sol.p(:, i), x2, zeros(size(x2)), zeros(size(x2)));
%!     Un = mean(sol.def(mesh.groups.iso20(grp_idx_volume2).nodes, 1, i));
%!     if (do_plot)
%!       figure("visible", "off");
%!       hold on;
%!       plot(x1 * unit_meters, p1ref * unit_pascal, "-;p1;r");
%!       plot(x2 * unit_meters, p2ref * unit_pascal, "-;p2;g");
%!       plot(x1 * unit_meters, p1n * unit_pascal, "-;p1n;m");
%!       plot(x2 * unit_meters, p2n * unit_pascal, "-;p2n;c");
%!       ylim([min(min(sol.p)), max(max(sol.p))] * unit_pascal);
%!       xlabel("x [m]");
%!       ylabel("p [Pa]");
%!       grid minor on;
%!       title(sprintf("pressure t=%.2fs", sol.t(i) * unit_second));
%!     endif
%!     assert_simple(p1n.', p1ref, tol * max(abs(p1ref)));
%!     assert_simple(p2n.', p2ref, tol * max(abs(p2ref)));
%!     assert_simple(Un, real(U * exp(1j * omega * sol.t(i))), tol * abs(U));
%!   endfor
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       [~] = unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
