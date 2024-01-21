## fem_pre_mesh_unstruct_create.m:10
%!test
%! ## TEST10
%! close all;
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1;
%! SI_unit_kilogram = 1e3;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! geo.D = 148e-3 / SI_unit_meter;
%! geo.d = 110e-3 / SI_unit_meter;
%! geo.h = 10e-3 / SI_unit_meter;
%! geo.x1 = 0 / SI_unit_meter;
%! geo.y1 = 0 / SI_unit_meter;
%! geo.z1 = 100e-3 / SI_unit_meter;
%! geo.x2 = 0 / SI_unit_meter;
%! geo.y2 = 0 / SI_unit_meter;
%! geo.z2 = geo.z1;
%! geo.lsuc = 100e-3 / SI_unit_meter;
%! geo.linlet = 100e-3 / SI_unit_meter;
%! geo.dsuc = 6e-3 / SI_unit_meter;
%! geo.dinlet = 12e-3 / SI_unit_meter;
%! geo.olev = 20e-3 / SI_unit_meter;
%! param.fmin = 0;
%! param.fmax = 2000 / (SI_unit_second^-1);
%! param.maxdef = 10e-3 / SI_unit_meter;
%! param.num_freq_modal = 50000;
%! param.omegaref = 2 * pi * 730 / (SI_unit_second^-1);
%! param.rho1 = 1.33 / (SI_unit_kilogram / SI_unit_meter^3);
%! param.c1 = 225 / (SI_unit_meter / SI_unit_second);
%! param.eta1 = 8.3492e-06 / (SI_unit_pascal * SI_unit_second);
%! param.rho2 = 850 / (SI_unit_kilogram / SI_unit_meter^3);
%! param.c2 = 1680 / (SI_unit_meter / SI_unit_second);
%! param.eta2 = 2.76670261412397E-3 / (SI_unit_pascal * SI_unit_second);
%! param.rho3 = 1.33 / (SI_unit_kilogram / SI_unit_meter^3);
%! param.c3 = 225 / (SI_unit_meter / SI_unit_second);
%! param.s = 0.5 * geo.dinlet * sqrt(param.rho1 * param.omegaref / param.eta1);
%! param.k = param.omegaref * 0.5 * geo.dinlet / param.c1;
%! ## SOUND ATTENUATION IN TUBES DUE TO VISCO-THERMAL EFFECTS
%! ## E. RODARTE, G. SINGH, N.R. MILLER AND P. HRNJAK
%! ## Mechanical and Industrial Engineering Department, Air Conditioning and Refrigeration
%! ## Center, University of Illinois at Urbana-Champaign, 1206 W. Green St.,
%! ## Urbana IL 61801, U.S.A.
%! ## (Received 3 December 1998, and in ,nal form 2 August 1999)
%! param.Gamma3 = 0.0081 + 1.00795j; ## s = 100; k = 0.01; Cp / Cv = 1.1;
%! ## Maxima: subst([GammaRe=real(Gamma3),GammaIm=imag(Gamma3)],rhs(solve(GammaRe+%i*GammaIm=%i/sqrt(1 + %i*omega*tau),[tau])[1]));
%! param.tau3 = -(real(param.Gamma3)^2+2*1i*imag(param.Gamma3)*real(param.Gamma3)-imag(param.Gamma3)^2+1)/((1i*real(param.Gamma3)^2-2*imag(param.Gamma3)*real(param.Gamma3)-1i*imag(param.Gamma3)^2)*param.omegaref);
%! param.eta3 = 3 / 4 * param.tau3 * param.rho3 * param.c3^2;
%! param.tau3_ = 1 / (param.rho3 * param.c3^2) * (4/3 * param.eta3);
%! param.k3 = 1 / sqrt(1 + 1j * param.omegaref * param.tau3_) * param.omegaref / param.c3;
%! assert_simple(1j * param.k3, param.Gamma3 * param.omegaref / param.c3, eps^0.9 * abs(param.Gamma3 * param.omegaref / param.c3));
%! param.m1 = 2.2 / SI_unit_kilogram;
%! param.J1 = 1e-6 * [2.4427674e+03 -3.2393301e+01 -5.3623318e+02
%!                    -3.2393301e+01  3.7506484e+03  7.9745924e+01
%!                    -5.3623318e+02  7.9745924e+01  4.2943071e+03] / (SI_unit_kilogram * SI_unit_meter^2);
%! param.lcg1 = zeros(3, 1);
%! param.m2 = 2.2978223 / SI_unit_kilogram;
%! param.J2 = 1e-6 * [8.6104517e+03, 5.5814351e+01, -3.0103453e+02;
%!                    5.5814351e+01, 1.1905289e+04,  2.0425595e+02;
%!                    -3.0103453e+02, 2.0425595e+02,  1.2109596e+04] / (SI_unit_kilogram * SI_unit_meter^2);
%! param.lcg2 = 1e-3 * [12.940131;
%!                      -6.0192164e-01;
%!                      5.9683120e+01] / SI_unit_meter;
%! param.offset1 = (6.63e-3 + 2.4e-3) / SI_unit_meter;
%! param.N = int32(120);
%! param.use_impedance = false;
%! param.use_pressure_bc = true;
%! param.use_PML = false;
%! helspr1.d = 1.3e-3 / SI_unit_meter;
%! helspr1.D = 12.12e-3 / SI_unit_meter + helspr1.d;
%! helspr1.L = 27.7e-3 / SI_unit_meter;
%! helspr1.n = 5.7;
%! helspr1.ni = 2.7;
%! helspr1.ng = 0.75;
%! helspr1.m = 60;
%! helspr1.material.E = 206000e6 / SI_unit_pascal;
%! helspr1.material.G = 81500e6 / SI_unit_pascal;
%! helspr1.material.nu = helspr1.material.E / (2 * helspr1.material.G) - 1;
%! helspr1.material.rho = 7850 / (SI_unit_kilogram / SI_unit_meter^3);
%! [helspr1.material.alpha, helspr1.material.beta] = fem_pre_mat_rayleigh_damping([1e-2; 0.03e-2], [20, 1000] / (SI_unit_second^-1));
%! helspr1.X = [[45.20e-3,  45.20e-3, -62.5e-3;
%!               43.13e-3, -43.13e-3,   0.0e-3] / SI_unit_meter;
%!              [-helspr1.L, -helspr1.L, -helspr1.L]];
%! section1.A = helspr1.d^2 * pi / 4;
%! section1.Ay = 9 / 10 * section1.A;
%! section1.Az = section1.Ay;
%! section1.Iy = helspr1.d^4 * pi / 64;
%! section1.Iz = section1.Iy;
%! section1.It = section1.Iy + section1.Iz;
%! rubspr1.d = 11.1e-3 / SI_unit_meter;
%! rubspr1.D = 28.6e-3 / SI_unit_meter;
%! rubspr1.L = 10.5e-3 / SI_unit_meter;
%! rubspr1.material.E = 2.6e6 / SI_unit_pascal;
%! rubspr1.material.nu = 0.499;
%! rubspr1.material.rho = 910 / (SI_unit_kilogram / SI_unit_meter^3);
%! [rubspr1.material.alpha, rubspr1.material.beta] = fem_pre_mat_rayleigh_damping([10e-2; 30e-2], [30, 100] / (SI_unit_second^-1));
%! rubspr1.X = [(170e-3 * [0.5, 0.5, -0.5, -0.5] + 8e-3) / SI_unit_meter;
%!              70e-3 * [-0.5, 0.5, -0.5, 0.5] / SI_unit_meter;
%!              -[1, 1, 1, 1] * (helspr1.L + rubspr1.L + param.offset1)];
%! rubspr1.e2 = [1, 0, 0];
%! section2.A = (rubspr1.D^2 - rubspr1.d^2) * pi / 4;
%! section2.Ay = 0.8 * section2.A;
%! section2.Az = section2.Ay;
%! section2.Iy = (rubspr1.D^4 - rubspr1.d^4) * pi / 64;
%! section2.Iz = section2.Iy;
%! section2.It = section2.Iy + section2.Iz;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   unwind_protect
%!     geo_file = [filename, ".geo"];
%!     fd = fopen(geo_file, "w");
%!     if (fd == -1)
%!       error("failed to open file %s", geo_file);
%!     endif
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fputs(fd, "vout = newv;\n");
%!     fputs(fd, "Sphere(vout) = {x2, y2, z2, 0.5 * D};\n");
%!     fputs(fd, "vsuc = newv;\n");
%!     fputs(fd, "Cylinder(vsuc) = {x2, y2, z2, 0.5 * D + lsuc, 0, 0, 0.5 * dsuc};\n");
%!     fputs(fd, "vinlet = newv;\n");
%!     fputs(fd, "Cylinder(vinlet) = {x2 + 0.5 * D + lsuc, y2, z2, linlet, 0, 0, 0.5 * dinlet};\n");
%!     fputs(fd, "vouts = newv;\n");
%!     fputs(fd, "BooleanUnion(vouts) = { Volume{vout}; Delete; }{ Volume{vsuc}; Delete; };\n")
%!     fputs(fd, "vin = newv;\n");
%!     fputs(fd, "Sphere(vin) = {x1, y1, z1, 0.5 * d};\n");
%!     fputs(fd, "v = newv;\n");
%!     fputs(fd, "BooleanDifference(v) = {Volume{vouts}; Delete; }{ Volume{vin}; Delete; };\n");
%!     fputs(fd, "voil = newv;\n");
%!     fputs(fd, "Box(voil) = {x2 - 0.5 * D, y2 - 0.5 * D, z2 - 0.5 * D, D, D, olev};\n");
%!     fputs(fd, "vfrag = BooleanFragments{Volume{v}; Delete; }{ Volume{voil}; Delete;};\n");
%!     fputs(fd, "vfrag2 = BooleanFragments{Volume{vfrag};Delete;}{ Volume{vinlet}; Delete; };\n");
%!     fputs(fd, "Physical Volume(\"gas\", 1) = {8};\n");
%!     fputs(fd, "Physical Volume(\"oil\", 2) = {9};\n");
%!     fputs(fd, "Physical Volume(\"inlet\", 3) = {7};\n");
%!     fputs(fd, "Physical Surface(\"fsi1\", 4) = {12,14};\n");
%!     fputs(fd, "Physical Surface(\"fsi2\", 5) = {8,13,9,24,22};\n");
%!     fputs(fd, "Physical Surface(\"inletsec\", 6) = {23};\n");
%!     fputs(fd, "Physical Surface(\"succon\", 7) = {9};\n");
%!     fputs(fd, "Physical Surface(\"suchose\", 8) = {22};\n");
%!     fputs(fd, "Mesh.HighOrderOptimize = 2;\n");
%!     fputs(fd, "Mesh.OptimizeThreshold = 0.99;\n");
%!     fputs(fd, "MeshSize{PointsOf{Physical Volume{1}; Physical Volume{2}; Physical Volume{3}; } } = h;\n");
%!     fputs(fd, "MeshSize{PointsOf{Physical Surface{7}; } } = dsuc * Pi / 7.;\n");
%!     fputs(fd, "MeshSize{PointsOf{Physical Surface{8}; } } = dinlet * Pi / 7.;\n");
%!     fputs(fd, "ReorientMesh Volume{7};\n");
%!     fputs(fd, "ReorientMesh Volume{8};\n");
%!     fputs(fd, "ReorientMesh Volume{9};\n");
%!     fputs(fd, "ReverseMesh Surface{8,13,9,24,22,12,14};\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   opt.mesh.order = 2;
%!   opt.mesh.elem_type = {"tria6h", "tet10h"};
%!   opt.interactive = false;
%!   mesh = fem_pre_mesh_unstruct_create(geo_file, geo, opt);
%!   mesh = fem_pre_mesh_reorder(mesh);
%!   node_idx_rb1 = int32(rows(mesh.nodes) + 1);
%!   node_idx_rb2 = int32(rows(mesh.nodes) + 2);
%!   grp_idx_v_gas = int32(find([mesh.groups.tet10h.id] == 1));
%!   grp_idx_v_oil = int32(find([mesh.groups.tet10h.id] == 2));
%!   grp_idx_v_inlet = int32(find([mesh.groups.tet10h.id] == 3));
%!   grp_idx_s_fsi1 = int32(find([mesh.groups.tria6h.id] == 4));
%!   grp_idx_s_fsi2 = int32(find([mesh.groups.tria6h.id] == 5));
%!   grp_idx_s_inletsec = int32(find([mesh.groups.tria6h.id] == 6));
%!   mesh.nodes(node_idx_rb1, :) = [0, 0, rubspr1.L + param.offset1 + helspr1.L, zeros(1, 3)];
%!   mesh.nodes(node_idx_rb2, :) = [0, 0, rubspr1.L, zeros(1, 3)];
%!   node_idx_helspr1 = int32(rows(mesh.nodes) + 1);
%!   helspr1.Phi = linspace(0, 2 * pi * helspr1.n, ceil(helspr1.m * helspr1.n)).' + 2 * pi * helspr1.ni;
%!   node_idx_rubspr1 = node_idx_helspr1 + numel(helspr1.Phi) * columns(helspr1.X);
%!   helspr1.x = 0.5 * helspr1.D * cos(helspr1.Phi);
%!   helspr1.y = 0.5 * helspr1.D * sin(helspr1.Phi);
%!   helspr1.z = (helspr1.L - helspr1.d * (2 * (helspr1.ni - helspr1.ng) + 1)) * linspace(0, 1, numel(helspr1.Phi))(:) + helspr1.d * (helspr1.ni - helspr1.ng + 0.5);
%!   helspr1.e2 = [0, 0, 1];
%!   helspr1.Phi0 = [pi/2, -pi/2, pi] - 2 * pi * helspr1.ni;
%!   idx_beam = int32(0);
%!   for i=1:columns(helspr1.X)
%!     R = euler123_to_rotation_matrix([0; 0; helspr1.Phi0(i)]);
%!     idx_node = node_idx_helspr1 - 1 + (1:numel(helspr1.Phi)) + (i - 1) * numel(helspr1.Phi);
%!     elnodes = int32([(1:numel(helspr1.Phi) - 1).', (2:numel(helspr1.Phi)).']) + node_idx_helspr1 - 1 + (i - 1) * numel(helspr1.Phi);
%!     mesh.nodes(idx_node, 1:6) = [[helspr1.x, helspr1.y, helspr1.z] * R.' + helspr1.X(:, i).', zeros(numel(helspr1.Phi), 3)];
%!     mesh.elements.beam2(idx_beam + (1:numel(helspr1.Phi) - 1)) = struct("nodes", mat2cell(elnodes, ones(numel(helspr1.Phi) - 1, 1, "int32"), 2), ...
%!                                                                         "section", mat2cell(repmat(section1, numel(helspr1.Phi) - 1, 1), ones(numel(helspr1.Phi) - 1, 1, "int32")), ...
%!                                                                         "e2", mat2cell(repmat(helspr1.e2, numel(helspr1.Phi) - 1, 1), ones(numel(helspr1.Phi) - 1, 1, "int32"), 3));
%!     idx_beam += numel(helspr1.Phi) - 1;
%!   endfor
%!   for i=1:columns(rubspr1.X)
%!     elnodes = int32([1, 2]) + node_idx_rubspr1 - 1 + (i - 1) * 2;
%!     mesh.nodes(node_idx_rubspr1 - 1 + (1:2) + 2 * (i - 1), 1:6) = [[zeros(1, 3); zeros(1, 2), rubspr1.L] + rubspr1.X(:, i).', zeros(2, 3)];
%!     mesh.elements.beam2(idx_beam + 1) = struct("nodes", mat2cell(elnodes, ones(1, 1, "int32"), 2), ...
%!                                                "section", mat2cell(repmat(section2, 1, 1), ones(1, 1, "int32")), ...
%!                                                "e2", mat2cell(rubspr1.e2, ones(1, 1, "int32"), 3));
%!     ++idx_beam;
%!   endfor
%!   for i=1:columns(helspr1.X)
%!     idx_node = node_idx_helspr1 - 1 + (1:numel(helspr1.Phi)) + (i - 1) * numel(helspr1.Phi);
%!     elnodes = int32([idx_node(1), node_idx_rb2;
%!                      idx_node(end), node_idx_rb1]);
%!     mesh.elements.beam2(idx_beam + (1:rows(elnodes))) = struct("nodes", mat2cell(elnodes, ones(rows(elnodes), 1, "int32"), 2), ...
%!                                                                "section", mat2cell(repmat(section1, rows(elnodes), 1), ones(rows(elnodes), 1, "int32")), ...
%!                                                                "e2", mat2cell(repmat(helspr1.e2, rows(elnodes), 1), ones(rows(elnodes), 1, "int32"), 3));
%!     idx_beam += rows(elnodes);
%!   endfor
%!   for i=1:columns(rubspr1.X)
%!     idx_node = int32([1, 2]) + node_idx_rubspr1 - 1 + (i - 1) * 2;
%!     elnodes = int32([idx_node(2), node_idx_rb2]);
%!     mesh.elements.beam2(idx_beam + (1:rows(elnodes))) = struct("nodes", mat2cell(elnodes, ones(rows(elnodes), 1, "int32"), 2), ...
%!                                                                "section", mat2cell(repmat(section2, rows(elnodes), 1), ones(rows(elnodes), 1, "int32")), ...
%!                                                                "e2", mat2cell(repmat(rubspr1.e2, rows(elnodes), 1), ones(rows(elnodes), 1, "int32"), 3));
%!     idx_beam += rows(elnodes);
%!   endfor
%!   empty_cell = cell(1, 2);
%!   mesh.elements.bodies = struct("m", empty_cell, "J", empty_cell, "lcg", empty_cell, "nodes", empty_cell);
%!   mesh.elements.bodies(1).m = param.m1;
%!   mesh.elements.bodies(1).J = param.J1;
%!   mesh.elements.bodies(1).lcg = param.lcg1;
%!   mesh.elements.bodies(1).nodes = node_idx_rb1;
%!   mesh.elements.bodies(2).m = param.m2;
%!   mesh.elements.bodies(2).J = param.J2;
%!   mesh.elements.bodies(2).lcg = param.lcg2;
%!   mesh.elements.bodies(2).nodes = node_idx_rb2;
%!   empty_cell = cell(1, 3);
%!   mesh.material_data = struct("E", empty_cell, "nu", empty_cell, "rho", empty_cell, "c", empty_cell, "eta", empty_cell);
%!   mesh.material_data(1).c = param.c1;
%!   mesh.material_data(1).rho = param.rho1;
%!   mesh.material_data(1).eta = param.eta1;
%!   mesh.material_data(2).c = param.c2;
%!   mesh.material_data(2).rho = param.rho2;
%!   mesh.material_data(2).eta = param.eta2;
%!   mesh.material_data(3).E = helspr1.material.E;
%!   mesh.material_data(3).nu = helspr1.material.nu;
%!   mesh.material_data(3).beta = helspr1.material.beta;
%!   mesh.material_data(3).alpha = helspr1.material.alpha;
%!   mesh.material_data(3).rho = helspr1.material.rho;
%!   mesh.material_data(4).E = rubspr1.material.E;
%!   mesh.material_data(4).nu = rubspr1.material.nu;
%!   mesh.material_data(4).rho = rubspr1.material.rho;
%!   mesh.material_data(4).beta = rubspr1.material.beta;
%!   mesh.material_data(4).alpha = rubspr1.material.alpha;
%!   mesh.material_data(5).E = helspr1.material.E;
%!   mesh.material_data(5).nu = helspr1.material.nu;
%!   mesh.material_data(5).beta = helspr1.material.beta;
%!   mesh.material_data(5).alpha = helspr1.material.alpha;
%!   mesh.material_data(5).rho = 0;
%!   mesh.material_data(6).c = param.c3;
%!   mesh.material_data(6).rho = param.rho3;
%!   mesh.material_data(6).eta = abs(param.eta3); ## FIXME: complex viscosity not supported yet
%!   load_case_dof.locked_dof = false(rows(mesh.nodes), 7);
%!   if (param.use_pressure_bc)
%!     load_case_dof.locked_dof(mesh.groups.tria6h(grp_idx_s_inletsec).nodes, 7) = true;
%!   endif
%!   load_case_dof.locked_dof(mesh.groups.tria6h(grp_idx_s_fsi1).nodes, 4:6) = true;
%!   load_case_dof.locked_dof(mesh.groups.tria6h(grp_idx_s_fsi2).nodes, 4:6) = true;
%!   for i=1:columns(rubspr1.X)
%!     load_case_dof.locked_dof(node_idx_rubspr1 + 2 * (i - 1), 1:6) = true;
%!   endfor
%!   load_case_dof.domain = FEM_DO_FLUID_STRUCT;
%!   mesh.materials.tet10h = zeros(rows(mesh.elements.tet10h), 1, "int32");
%!   mesh.materials.tet10h(mesh.groups.tet10h(grp_idx_v_gas).elements) = 1;
%!   mesh.materials.tet10h(mesh.groups.tet10h(grp_idx_v_oil).elements) = 2;
%!   mesh.materials.tet10h(mesh.groups.tet10h(grp_idx_v_inlet).elements) = 6;
%!   mesh.materials.beam2 = [repmat(int32(3), columns(helspr1.X) * (numel(helspr1.Phi) - 1), 1);
%!                           repmat(int32(4), columns(rubspr1.X), 1);
%!                           repmat(int32(5), columns(helspr1.X) * 2, 1);
%!                           repmat(int32(5), columns(rubspr1.X), 1)];
%!   mesh.elements.fluid_struct_interface.tria6h = mesh.elements.tria6h([[mesh.groups.tria6h([grp_idx_s_fsi1, grp_idx_s_fsi2])].elements], :);
%!   if (param.use_PML)
%!     mesh.elements.perfectly_matched_layers.tet10h.f = zeros(3, columns(mesh.elements.tet10h), rows(mesh.elements.tet10h));
%!     mesh.elements.perfectly_matched_layers.tet10h.e1 = zeros(3, columns(mesh.elements.tet10h), rows(mesh.elements.tet10h));
%!     mesh.elements.perfectly_matched_layers.tet10h.e2 = zeros(3, columns(mesh.elements.tet10h), rows(mesh.elements.tet10h));
%!     for i=1:3
%!       mesh.elements.perfectly_matched_layers.tet10h.f(i, :, :) = 1;
%!     endfor
%!     mesh.elements.perfectly_matched_layers.tet10h.e1(1, :, :) = 1;
%!     mesh.elements.perfectly_matched_layers.tet10h.e2(2, :, :) = 1;
%!     elem_idx_PML = mesh.groups.tet10h(grp_idx_v_inlet).elements;
%!     node_idx_PML = mesh.elements.tet10h(elem_idx_PML, :);
%!     xPML = mesh.nodes(:, 1)(node_idx_PML) - 0.5 * geo.D - geo.lsuc;
%!     deltaPML = geo.linlet;
%!     sigmax = 1 ./ (deltaPML - xPML);
%!     k = param.omegaref / param.c1;
%!     mesh.elements.perfectly_matched_layers.tet10h.f(1, :, elem_idx_PML) = 1 ./ (1 - 1j * sigmax.' / k);
%!   endif
%!   if (param.use_impedance)
%!     mesh.elements.acoustic_impedance.tria6h.nodes = mesh.elements.tria6h(mesh.groups.tria6h(grp_idx_s_inletsec).elements, :);
%!     mesh.elements.acoustic_impedance.tria6h.z = repmat(param.rho * param.c, size(mesh.elements.acoustic_impedance.tria6h.nodes));
%!     mesh.materials.acoustic_impedance.tria6h = ones(rows(mesh.elements.acoustic_impedance.tria6h.nodes), 1, "int32");
%!   endif
%!   empty_cell = cell(1, numel(mesh.groups.tria6h(grp_idx_s_fsi1).nodes) + numel(mesh.groups.tria6h(grp_idx_s_fsi2).nodes) + 2 * columns(helspr1.X) + columns(rubspr1.X));
%!   mesh.elements.joints = struct("nodes", empty_cell, "C", empty_cell);
%!   idx_joint = int32(0);
%!   for i=1:numel(mesh.groups.tria6h(grp_idx_s_fsi1).nodes)
%!     node_idx_slave = mesh.groups.tria6h(grp_idx_s_fsi1).nodes(i);
%!     lslave = mesh.nodes(node_idx_slave, 1:3) - mesh.nodes(node_idx_rb1, 1:3);
%!     ++idx_joint;
%!     mesh.elements.joints(idx_joint).nodes = int32([node_idx_rb1, node_idx_slave]);
%!     mesh.elements.joints(idx_joint).C = [eye(3), -skew(lslave), -eye(3), zeros(3, 3)];
%!   endfor
%!   for i=1:numel(mesh.groups.tria6h(grp_idx_s_fsi2).nodes)
%!     node_idx_slave = mesh.groups.tria6h(grp_idx_s_fsi2).nodes(i);
%!     lslave = mesh.nodes(node_idx_slave, 1:3) - mesh.nodes(node_idx_rb2, 1:3);
%!     ++idx_joint;
%!     mesh.elements.joints(idx_joint).nodes = int32([node_idx_rb2, node_idx_slave]);
%!     mesh.elements.joints(idx_joint).C = [eye(3), -skew(lslave), -eye(3), zeros(3, 3)];
%!   endfor
%!   for i=1:columns(helspr1.X)
%!     node_idx_slave = node_idx_helspr1 - 1 + numel(helspr1.Phi) * i;
%!     lslave = mesh.nodes(node_idx_slave, 1:3) - mesh.nodes(node_idx_rb1, 1:3);
%!     ++idx_joint;
%!     mesh.elements.joints(idx_joint).nodes = int32([node_idx_rb1, node_idx_slave]);
%!     mesh.elements.joints(idx_joint).C = [     eye(3), -skew(lslave),     -eye(3), zeros(3, 3);
%!                                               zeros(3, 3),        eye(3), zeros(3, 3),     -eye(3)];
%!   endfor
%!   for i=1:columns(helspr1.X)
%!     node_idx_slave = node_idx_helspr1 + numel(helspr1.Phi) * (i - 1);
%!     lslave = mesh.nodes(node_idx_slave, 1:3) - mesh.nodes(node_idx_rb2, 1:3);
%!     ++idx_joint;
%!     mesh.elements.joints(idx_joint).nodes = int32([node_idx_rb2, node_idx_slave]);
%!     mesh.elements.joints(idx_joint).C = [     eye(3), -skew(lslave),     -eye(3), zeros(3, 3);
%!                                               zeros(3, 3),        eye(3), zeros(3, 3),     -eye(3)];
%!   endfor
%!   for i=1:columns(rubspr1.X)
%!     node_idx_slave = node_idx_rubspr1 + 2 * i - 1;
%!     lslave = mesh.nodes(node_idx_slave, 1:3) - mesh.nodes(node_idx_rb2, 1:3);
%!     ++idx_joint;
%!     mesh.elements.joints(idx_joint).nodes = int32([node_idx_rb2, node_idx_slave]);
%!     mesh.elements.joints(idx_joint).C = [     eye(3), -skew(lslave),     -eye(3), zeros(3, 3);
%!                                               zeros(3, 3),        eye(3), zeros(3, 3),     -eye(3)];
%!   endfor
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!   empty_cell = cell(1, 6);
%!   load_case = struct("loads", empty_cell, "loaded_nodes", empty_cell);
%!   idx_load = int32(0);
%!   for i=1:6
%!     ++idx_load;
%!     load_case(idx_load).loaded_nodes = node_idx_rb1;
%!     load_case(idx_load).loads = zeros(1, 6);
%!     load_case(idx_load).loads(i) = 1;
%!   endfor
%!   [mat_ass.Mfs_re, ...
%!    mat_ass.Mfs_im, ...
%!    mat_ass.Kfs_re, ...
%!    mat_ass.Kfs_im, ...
%!    mat_ass.Dfs_re, ...
%!    mat_ass.Dfs_im, ...
%!    mat_ass.Rfs, ...
%!    mat_ass.mat_info, ...
%!    mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                        dof_map, ...
%!                                        [FEM_MAT_MASS_FLUID_STRUCT_RE, ...
%!                                         FEM_MAT_MASS_FLUID_STRUCT_IM, ...
%!                                         FEM_MAT_STIFFNESS_FLUID_STRUCT_RE, ...
%!                                         FEM_MAT_STIFFNESS_FLUID_STRUCT_IM, ...
%!                                         FEM_MAT_DAMPING_FLUID_STRUCT_RE, ...
%!                                         FEM_MAT_DAMPING_FLUID_STRUCT_IM, ...
%!                                         FEM_VEC_LOAD_FLUID_STRUCT], ...
%!                                        load_case);
%!   if (~nnz(mat_ass.Mfs_im))
%!     mat_ass = rmfield(mat_ass, "Mfs_im");
%!   endif
%!   if (~nnz(mat_ass.Dfs_im))
%!     mat_ass = rmfield(mat_ass, "Dfs_im");
%!   endif
%!   if (~nnz(mat_ass.Kfs_im))
%!     mat_ass = rmfield(mat_ass, "Kfs_im");
%!   endif
%!   opt_linsol.solver = "umfpack";
%!   opt_linsol.pre_scaling = true;
%!   opt_linsol.number_of_threads = mbdyn_solver_num_threads_default();
%!   opt_linsol.refine_max_iter = int32(3);
%!   opt_linsol.epsilon_refinement = 1e-10;
%!   opt_linsol.verbose = int32(1);
%!   opt_linsol.weighted_matching = true;
%!   opt_linsol.scaling = true;
%!   opt_eigs.maxit = int32(100);
%!   opt_eigs.disp = int32(0);
%!   opt_eigs.p = 5 * param.N;
%!   [sol_eig, Phi] = fem_sol_modal_fsi(mesh, dof_map, mat_ass, param.N, opt_linsol, opt_eigs);
%!   [Phi, h] = fem_sol_modes_scale(mat_ass.Mfs_re, mat_ass.Kfs_re, sol_eig.lambda, Phi, mat_ass.Rfs);
%!   omega = linspace(2 * pi * param.fmin, 2 * pi * param.fmax, param.num_freq_modal);
%!   U2 = fem_sol_harmonic_modal(h, sol_eig.lambda, Phi(dof_map.ndof(node_idx_rb2, 1:6), :), omega);
%!   for i=1:3
%!     figure("visible", "off");
%!     subplot(2, 1, 1);
%!     hold on;
%!     plot(omega / (2*pi) * (SI_unit_second^-1), 20 * log10(abs(-omega.^2 .* U2(i, :, i)) * (SI_unit_meter / SI_unit_second^2)), "-;modal;r");
%!     xlabel("f [Hz]");
%!     ylabel("a [dB/(1 m/s^2/N)]");
%!     ylim([-80, 10]);
%!     yticks(-100:5:10);
%!     grid on;
%!     grid minor on;
%!     title("frequency response magnitude");
%!     subplot(2, 1, 2);
%!     hold on;
%!     plot(omega / (2*pi) * (SI_unit_second^-1), 180 / pi * arg(-omega.^2 .* U2(i, :, i)), "-;modal;r");
%!     xlabel("f [Hz]");
%!     ylabel("arg(a) [deg]");
%!     ylim([-180,180]);
%!     yticks(-180:30:180);
%!     grid on;
%!     grid minor on;
%!     title("frequency response phase");
%!   endfor
%!   idx_freq = find(sol_eig.f >= 0);
%!   sol_eig.f = sol_eig.f(idx_freq);
%!   sol_eig.D = sol_eig.D(idx_freq);
%!   sol_eig.lambda = sol_eig.lambda(idx_freq);
%!   sol_eig.def = sol_eig.def(:, :, idx_freq);
%!   sol_eig.p = sol_eig.p(:, idx_freq);
%!   for i=1:size(sol_eig.def, 3)
%!     sre = max(max(abs(real(sol_eig.def(:, 1:3, i)))));
%!     sim = max(max(abs(imag(sol_eig.def(:, 1:3, i)))));
%!     if (sim > sre)
%!       sol_eig.def(:, :, i) = imag(sol_eig.def(:, :, i)) * param.maxdef / sim;
%!       sol_eig.p(:, i) = imag(sol_eig.p(:, i)) * param.maxdef / sim;
%!     else
%!       sol_eig.def(:, :, i) = real(sol_eig.def(:, :, i)) * param.maxdef / sre;
%!       sol_eig.p(:, i) = real(sol_eig.p(:, i)) * param.maxdef / sre;
%!     endif
%!   endfor
%!   opt_post.elem_types={"tria6h","beam2"};
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       [~] = unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
