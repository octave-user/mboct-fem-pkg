## fem_pre_mesh_import.m:249
%!test
%! ## TEST 249
%! ####################################################
%! ## Jont Allen
%! ## THE ACOUSTIC WAVE EQUATION AND SIMPLE SOLUTIONS
%! ## Chapter 5.13
%! ####################################################
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
%!     r1 = 25e-3 / unit_meters;
%!     r2 = 2000e-3 / unit_meters;
%!     dx = 5e-3 / unit_meters;
%!     h = dx;
%!     phi = atan(dx / mean([r1, r2]));
%!     c = 340 / (unit_meters / unit_second);
%!     rho = 1.225 / (unit_kilograms / unit_meters^3);
%!     lambda = 40 * dx;
%!     f = c / lambda;
%!     omega = 2 * pi * f;
%!     k = omega / c;
%!     deltaPML = 3 * lambda;
%!     alphaPML = 1;
%!     A = 1;
%!     H1 = @(alpha, x) besselj(alpha, x) + 1j * bessely(alpha, x); ## equation 5.13.8
%!     H2 = @(alpha, x) besselj(alpha, x) - 1j * bessely(alpha, x);
%!     pref = @(r, t) A * H2(0, k * r) .* exp(1j * omega * t); ## equation 5.13.9
%!     vref = @(r, t) -1j * (A / (rho * c)) * H2(1, k * r) .* exp(1j * omega * t); ## equation 5.13.11
%!     vn = vref(r1, 0);
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "r1 = %g;\n", r1);
%!     fprintf(fd, "r2 = %g;\n", r2);
%!     fprintf(fd, "phi = %g;\n", phi);
%!     fprintf(fd, "h = %g;\n", h);
%!     fprintf(fd, "dx = %g;\n", dx);
%!     fprintf(fd, "deltaPML = %g;\n", deltaPML);
%!     fputs(fd, "Point(1) = {r1, 0, 0};\n");
%!     fputs(fd, "Point(2) = {r2, 0, 0};\n");
%!     fputs(fd, "Point(3) = {r2 + deltaPML, 0, 0};\n");
%!     fputs(fd, "Line(1) = {1, 2};\n");
%!     fputs(fd, "Line(2) = {2, 3};\n");
%!     fputs(fd, "tmp1 = Extrude{0, 0, h}{Line{1,2}; };\n");
%!     fputs(fd, "tmp2 = Extrude{{0,0,1},{0, 0, 0},phi}{Surface{tmp1[1],tmp1[5]};};\n");
%!     fputs(fd, "Physical Volume(\"v1\", 1) = {1};\n");
%!     fputs(fd, "Physical Volume(\"v2\", 2) = {2};\n");
%!     fputs(fd, "Physical Surface(\"s1\", 1) = {3};\n");
%!     fputs(fd, "Physical Surface(\"s2\", 2) = {8};\n");
%!     fputs(fd, "MeshSize{PointsOf{Volume{1,2};}} = dx;\n");
%!     fputs(fd, "ReorientMesh Volume{1};\n");
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
%!   opt_mesh.elem_type = {"tria6h", "tet10h"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_mesh));
%!   [~] = unlink([filename, ".msh"]);
%!   grp_idx_v1 = find([mesh.groups.tet10h.id] == 1);
%!   grp_idx_v2 = find([mesh.groups.tet10h.id] == 2);
%!   elem_idx_v1 = mesh.groups.tet10h(grp_idx_v1).elements;
%!   elem_idx_v2 = mesh.groups.tet10h(grp_idx_v2).elements;
%!   grp_idx_s1 = find([mesh.groups.tria6h.id] == 1);
%!   grp_idx_s2 = find([mesh.groups.tria6h.id] == 2);
%!   elem_idx_s1 = mesh.groups.tria6h(grp_idx_s1).elements;
%!   node_idx_s2 = mesh.groups.tria6h(grp_idx_s2).nodes;
%!   mesh.material_data.c = c;
%!   mesh.material_data.rho = rho;
%!   mesh.materials.tet10h = ones(rows(mesh.elements.tet10h), 1, "int32");
%!   mesh.elements.particle_velocity.tria6h.nodes = mesh.elements.tria6h(elem_idx_s1, :);
%!   mesh.materials.particle_velocity.tria6h = ones(rows(mesh.elements.particle_velocity.tria6h.nodes), 1, "int32");
%!   mesh.elements.perfectly_matched_layers.tet10h.f = zeros(3, columns(mesh.elements.tet10h), rows(mesh.elements.tet10h));
%!   mesh.elements.perfectly_matched_layers.tet10h.e1 = zeros(3, columns(mesh.elements.tet10h), rows(mesh.elements.tet10h));
%!   mesh.elements.perfectly_matched_layers.tet10h.e2 = zeros(3, columns(mesh.elements.tet10h), rows(mesh.elements.tet10h));
%!   mesh.elements.perfectly_matched_layers.tet10h.e1(1, :, :) = mesh.nodes(:, 1)(mesh.elements.tet10h.');
%!   mesh.elements.perfectly_matched_layers.tet10h.e1(2, :, :) = mesh.nodes(:, 2)(mesh.elements.tet10h.');
%!   mesh.elements.perfectly_matched_layers.tet10h.e2(3, :, :) = 1;
%!   for i=1:3
%!     mesh.elements.perfectly_matched_layers.tet10h.f(i, :, :) = 1;
%!   endfor
%!   x = mesh.nodes(:, 1)(mesh.elements.tet10h(elem_idx_v2, :).');
%!   y = mesh.nodes(:, 2)(mesh.elements.tet10h(elem_idx_v2, :).');
%!   dr = sqrt(x.^2 + y.^2) - r2;
%!   #####################################################################
%!   ## A parameter-free perfectly matched layer formulation
%!   ## for the finite-element-based solution of the Helmholtz
%!   ## equation
%!   ## Radu Cimpeanu a, Anton Martinsson b , Matthias Heil c
%!   ## a Department
%!   ## of Mathematics, Imperial College London, SW7 2AZ, London, United
%!   ## Kingdom
%!   ## (https://core.ac.uk/download/pdf/77019151.pdf)
%!   #####################################################################
%!   sigmax = alphaPML ./ (deltaPML - dr);
%!   mesh.elements.perfectly_matched_layers.tet10h.f(1, :, elem_idx_v2) = 1 ./ (1 - 1j * sigmax / k);
%!   load_case = struct("particle_velocity", cell(1, 2));
%!   load_case(1).particle_velocity.tria6h.vn = repmat(real(-vn), size(mesh.elements.particle_velocity.tria6h.nodes));
%!   load_case(2).particle_velocity.tria6h.vn = repmat(imag(-vn), size(mesh.elements.particle_velocity.tria6h.nodes));
%!   load_case_dof.locked_dof = false(rows(mesh.nodes), 1);
%!   load_case_dof.locked_dof(node_idx_s2, :) = true;
%!   load_case_dof.domain = FEM_DO_ACOUSTICS;
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!   [mat_ass.Ka_re, ...
%!    mat_ass.Ka_im, ...
%!    mat_ass.Da_re, ...
%!    mat_ass.Da_im, ...
%!    mat_ass.Ma_re, ...
%!    mat_ass.Ma_im, ...
%!    mat_ass.Ra, ...
%!    mat_ass.mat_info, ...
%!    mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                        dof_map, ...
%!                                        [FEM_MAT_STIFFNESS_ACOUSTICS_RE, ...
%!                                         FEM_MAT_STIFFNESS_ACOUSTICS_IM, ...
%!                                         FEM_MAT_DAMPING_ACOUSTICS_RE, ...
%!                                         FEM_MAT_DAMPING_ACOUSTICS_IM, ...
%!                                         FEM_MAT_MASS_ACOUSTICS_RE, ...
%!                                         FEM_MAT_MASS_ACOUSTICS_IM, ...
%!                                         FEM_VEC_LOAD_ACOUSTICS], ...
%!                                        load_case);
%!   opt_sol.number_of_threads = mbdyn_solver_num_threads_default();
%!   opt_sol.solver = "pastix";
%!   opt_sol.refine_max_iter = int32(50);
%!   Keff = -omega^2 * complex(mat_ass.Ma_re, mat_ass.Ma_im) + 1j * omega * complex(mat_ass.Da_re, mat_ass.Da_im) + complex(mat_ass.Ka_re, mat_ass.Ka_im);
%!   Reff = complex(mat_ass.Ra(:, 1), mat_ass.Ra(:, 2));
%!   Phi = fem_sol_factor(Keff, opt_sol) \ Reff;
%!   Psi = linspace(0, 2 * pi, 37);
%!   idx = dof_map.ndof(:, 1);
%!   iact = find(idx > 0);
%!   solC.Phi = zeros(rows(mesh.nodes), columns(Reff));
%!   solC.Phi(iact, :) = Phi(idx(iact), :);
%!   solC.PhiP = 1j * omega * solC.Phi;
%!   solR.t = Psi / omega;
%!   solR.p = real(-solC.PhiP * exp(1j * omega * solR.t));
%!   x = mesh.nodes(:, 1);
%!   y = mesh.nodes(:, 2);
%!   r = sqrt(x.^2 + y.^2);
%!   [r, idx] = sort(r);
%!   idx2 = find(r <= r2);
%!   r = r(idx2);
%!   idx = idx(idx2);
%!   for i=1:numel(solR.t)
%!     figure("visible", "off");
%!     hold on;
%!     plot(r * unit_meters, solR.p(idx, i) * unit_pascal, "-;p(r);r");
%!     plot(r * unit_meters, real(pref(r, solR.t(i))) * unit_pascal, "-;pref;k");
%!     ylim([min(min(solR.p)), max(max(solR.p))] * unit_pascal);
%!     xlabel("r [m]");
%!     ylabel("p [Pa]");
%!     grid on;
%!     grid minor on;
%!     title(sprintf("pressure t=%.2fs", solR.t(i) * unit_second));
%!   endfor
%!   preft = real(pref(r, solR.t));
%!   tol = 5e-4;
%!   assert_simple(solR.p(idx, :), preft, tol * max(max(abs(preft))));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       [~] = unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
