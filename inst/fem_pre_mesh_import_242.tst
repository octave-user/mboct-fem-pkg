## fem_pre_mesh_import.m:242
%!test
%! try
%! ## TEST 242
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
%!     a = 25e-3 / unit_meters;
%!     l = 2000e-3 / unit_meters;
%!     w = 2000e-3 / unit_meters;
%!     dx = 25e-3 / unit_meters;
%!     h = dx;
%!     c = 340 / (unit_meters / unit_second);
%!     rho = 1.225 / (unit_kilograms / unit_meters^3);
%!     lambda = 20 * dx;
%!     f = c / lambda;
%!     omega = 2 * pi * f;
%!     k = omega / c;
%!     deltaPML = 1e-5 / abs(k);
%!     alphaPML = 1;
%!     A = 1;
%!     H1 = @(alpha, x) besselj(alpha, x) + 1j * bessely(alpha, x);
%!     H2 = @(alpha, x) besselj(alpha, x) - 1j * bessely(alpha, x);
%!     pref = @(r, t) A * H2(0, k * r) .* exp(1j * omega * t);
%!     vref = @(r, t) -1j * (A / (rho * c)) * H2(1, k * r) .* exp(1j * omega * t);
%!     vn = vref(a, 0);
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "a = %g;\n", a);
%!     fprintf(fd, "l = %g;\n", l);
%!     fprintf(fd, "w = %g;\n", w);
%!     fprintf(fd, "h = %g;\n", h);
%!     fprintf(fd, "dx = %g;\n", dx);
%!     fprintf(fd, "deltaPML = %g;\n", deltaPML);
%!     fputs(fd, "Point(1) = {0, 0, 0};\n");
%!     fputs(fd, "Point(2) = {a, 0, 0};\n");
%!     fputs(fd, "Point(3) = {0, a, 0};\n");
%!     fputs(fd, "Point(4) = {-a, 0, 0};\n");
%!     fputs(fd, "Point(5) = {0, -a, 0};\n");
%!     fputs(fd, "Point(6) = {0.5 * l, 0.5 * w, 0};\n");
%!     fputs(fd, "Point(7) = {-0.5 * l, 0.5 * w, 0};\n");
%!     fputs(fd, "Point(8) = {-0.5 * l, -0.5 * w, 0};\n");
%!     fputs(fd, "Point(9) = {0.5 * l, -0.5 * w, 0};\n");
%!     fputs(fd, "Point(10) = {0.5 * l + deltaPML, 0.5 * w + deltaPML, 0};\n");
%!     fputs(fd, "Point(11) = {-0.5 * l - deltaPML, 0.5 * w + deltaPML, 0};\n");
%!     fputs(fd, "Point(12) = {-0.5 * l - deltaPML, -0.5 * w - deltaPML, 0};\n");
%!     fputs(fd, "Point(13) = {0.5 * l + deltaPML, -0.5 * w - deltaPML, 0};\n");
%!     fputs(fd, "Point(14) = {0.5 * l, 0.5 * w + deltaPML, 0};\n");
%!     fputs(fd, "Point(15) = {-0.5 * l, 0.5 * w + deltaPML, 0};\n");
%!     fputs(fd, "Point(16) = {-0.5 * l - deltaPML, 0.5 * w, 0};\n");
%!     fputs(fd, "Point(17) = {-0.5 * l - deltaPML, -0.5 * w, 0};\n");
%!     fputs(fd, "Point(18) = {-0.5 * l, -0.5 * w - deltaPML, 0};\n");
%!     fputs(fd, "Point(19) = {0.5 * l, -0.5 * w - deltaPML, 0};\n");
%!     fputs(fd, "Point(20) = {0.5 * l + deltaPML, -0.5 * w, 0};\n");
%!     fputs(fd, "Point(21) = {0.5 * l + deltaPML, 0.5 * w, 0};\n");
%!     fputs(fd, "Circle(1) = {2, 1, 3};\n");
%!     fputs(fd, "Circle(2) = {3, 1, 4};\n");
%!     fputs(fd, "Circle(3) = {4, 1, 5};\n");
%!     fputs(fd, "Circle(4) = {5, 1, 2};\n");
%!     fputs(fd, "Line(5) = {6, 7};\n");
%!     fputs(fd, "Line(6) = {7, 8};\n");
%!     fputs(fd, "Line(7) = {8, 9};\n");
%!     fputs(fd, "Line(8) = {9, 6};\n");
%!     fputs(fd, "Line(9) = {14, 15};\n");
%!     fputs(fd, "Line(10) = {16, 17};\n");
%!     fputs(fd, "Line(11) = {18, 19};\n");
%!     fputs(fd, "Line(12) = {20, 21};\n");
%!     fputs(fd, "Line(13) = {10, 14};\n");
%!     fputs(fd, "Line(14) = {15, 11};\n");
%!     fputs(fd, "Line(15) = {11, 16};\n");
%!     fputs(fd, "Line(16) = {17, 12};\n");
%!     fputs(fd, "Line(17) = {12, 18};\n");
%!     fputs(fd, "Line(18) = {19, 13};\n");
%!     fputs(fd, "Line(19) = {13, 20};\n");
%!     fputs(fd, "Line(20) = {21, 10};\n");
%!     fputs(fd, "Line(21) = {6, 21};\n");
%!     fputs(fd, "Line(22) = {6, 14};\n");
%!     fputs(fd, "Line(23) = {7, 15};\n");
%!     fputs(fd, "Line(24) = {7, 16};\n");
%!     fputs(fd, "Line(25) = {8, 17};\n");
%!     fputs(fd, "Line(26) = {8, 18};\n");
%!     fputs(fd, "Line(27) = {9, 19};\n");
%!     fputs(fd, "Line(28) = {9, 20};\n");
%!     fputs(fd, "Curve Loop(29) = {1, 2, 3, 4};\n");
%!     fputs(fd, "Curve Loop(30) = {5, 6, 7, 8};\n");
%!     fputs(fd, "Curve Loop(31) = {-5, 22, 9, -23};\n");
%!     fputs(fd, "Curve Loop(32) = {24, 10, -25, -6};\n");
%!     fputs(fd, "Curve Loop(33) = {-7, 26, 11, -27};\n");
%!     fputs(fd, "Curve Loop(34) = {12, -21, -8, 28};\n");
%!     fputs(fd, "Curve Loop(35) = {21, 20, 13, -22};\n");
%!     fputs(fd, "Curve Loop(36) = {23, 14, -15, -24};\n");
%!     fputs(fd, "Curve Loop(37) = {25, 16, 17, -26};\n");
%!     fputs(fd, "Curve Loop(38) = {27, 18, -19, -28};\n");
%!     fputs(fd, "Plane Surface(39) = {30, 29};\n");
%!     fputs(fd, "Plane Surface(40) = {31};\n");
%!     fputs(fd, "Plane Surface(41) = {32};\n");
%!     fputs(fd, "Plane Surface(42) = {33};\n");
%!     fputs(fd, "Plane Surface(43) = {34};\n");
%!     fputs(fd, "Plane Surface(44) = {35};\n");
%!     fputs(fd, "Plane Surface(45) = {36};\n");
%!     fputs(fd, "Plane Surface(46) = {37};\n");
%!     fputs(fd, "Plane Surface(47) = {38};\n");
%!     fputs(fd, "Coherence;\n");
%!     fputs(fd, "tmp1 = Extrude{0, 0, h}{Surface{39, 40, 41, 42, 43, 44, 45, 46, 47};};\n");
%!     fputs(fd, "Physical Volume(\"v1\", 1) = {1};\n");
%!     fputs(fd, "Physical Volume(\"v2\", 2) = {3, 5};\n");
%!     fputs(fd, "Physical Volume(\"v3\", 3) = {2, 4};\n");
%!     fputs(fd, "Physical Volume(\"v4\", 4) = {6, 7, 8, 9};\n");
%!     fputs(fd, "Physical Surface(\"s1\", 1) = {53, 54, 55, 52};\n");
%!     gmsh_version = fem_get_software_version("gmsh");
%!     if (gmsh_version(1) == 4 && gmsh_version(2) <= 11)
%!       fputs(fd, "Physical Surface(\"s2\", 2) = {74, 58, 76, 77, 62, 79, 80, 66, 82, 83, 69, 73};\n");
%!     else
%!       fputs(fd, "Physical Surface(\"s2\", 2) = {67, 82, 83, 69, 73, 74, 59, 76, 77, 62, 79, 80};\n");
%!     endif
%!     fputs(fd, "MeshSize{PointsOf{Volume{1,2,3,4,5,6,7,8,9};}} = dx;\n");
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
%!   grp_idx_v3 = find([mesh.groups.tet10h.id] == 3);
%!   grp_idx_v4 = find([mesh.groups.tet10h.id] == 4);
%!   elem_idx_v2 = mesh.groups.tet10h(grp_idx_v2).elements;
%!   elem_idx_v3 = mesh.groups.tet10h(grp_idx_v3).elements;
%!   elem_idx_v4 = mesh.groups.tet10h(grp_idx_v4).elements;
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
%!   mesh.elements.perfectly_matched_layers.tet10h.e1(1, :, :) = 1;
%!   mesh.elements.perfectly_matched_layers.tet10h.e2(2, :, :) = 1;
%!   for i=1:3
%!     mesh.elements.perfectly_matched_layers.tet10h.f(i, :, :) = 1;
%!   endfor
%!   x = mesh.nodes(:, 1);
%!   y = mesh.nodes(:, 2);
%!   dx2 = abs(x(mesh.elements.tet10h(elem_idx_v2, :).')) - 0.5 * l;
%!   dy3 = abs(y(mesh.elements.tet10h(elem_idx_v3, :).')) - 0.5 * w;
%!   dx4 = abs(x(mesh.elements.tet10h(elem_idx_v4, :).')) - 0.5 * l;
%!   dy4 = abs(y(mesh.elements.tet10h(elem_idx_v4, :).')) - 0.5 * w;
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
%!   sigmax2 = alphaPML ./ (deltaPML - dx2);
%!   sigmay3 = alphaPML ./ (deltaPML - dy3);
%!   sigmax4 = alphaPML ./ (deltaPML - dx4);
%!   sigmay4 = alphaPML ./ (deltaPML - dy4);
%!   mesh.elements.perfectly_matched_layers.tet10h.f(1, :, elem_idx_v2) = 1 ./ (1 - 1j * sigmax2 / k);
%!   mesh.elements.perfectly_matched_layers.tet10h.f(2, :, elem_idx_v3) = 1 ./ (1 - 1j * sigmay3 / k);
%!   mesh.elements.perfectly_matched_layers.tet10h.f(1, :, elem_idx_v4) = 1 ./ (1 - 1j * sigmax4 / k);
%!   mesh.elements.perfectly_matched_layers.tet10h.f(2, :, elem_idx_v4) = 1 ./ (1 - 1j * sigmay4 / k);
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
%!   opt_sol.refine_max_iter = int32(50);
%!   opt_sol.symmetric = false; ## FIXME: Workaround suspected issue in PaStiX 6.4.0
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
%!   r = sqrt(x.^2 + y.^2);
%!   [r, idx] = sort(r);
%!   idx2 = find(r <= min([l, w]) / 2);
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
%!   tol = 3e-3;
%!   assert_simple(solR.p(idx, :), preft, tol * max(max(abs(preft))));
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
