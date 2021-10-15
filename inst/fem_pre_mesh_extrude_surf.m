## Copyright (C) 2021(-2021) Reinhard <octave-user@a1.net>
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} [@var{nodes}, @var{elem}] = fem_pre_mesh_extrude_surf(@var{mesh}, @var{elem_type}, @var{grp_id}, @var{h})
## Extrude a surface mesh along the surface normal in order to generate a volume mesh
##
## @var{mesh} @dots{} existing mesh with surface elements of type @var{elem_type}
##
## @var{elem_type} @dots{} element type of surface elements to be extruded
##
## @var{grp_id} @dots{} group id of surface elements to be extruded
##
## @var{h} @dots{} each element of the vector @var{h} defines the thickness of one extruded layer normal to the surface
##
## @end deftypefn

function [nodes, elem] = fem_pre_mesh_extrude_surf(mesh, elem_type, grp_id, h)
  if (nargin ~= 4 || nargout ~= 2)
    print_usage();
  endif
  
  layers = numel(h);
  
  elem_grp = getfield(mesh.groups, elem_type);

  grp_idx = find([elem_grp.id] == grp_id);

  elem_idx = [[elem_grp(grp_idx)].elements];

  elem_nodes = getfield(mesh.elements, elem_type)(elem_idx, :);

  load_case_dof_n.locked_dof = false(rows(mesh.nodes), 1);
  load_case_dof_n.domain = FEM_DO_ACOUSTICS;

  mesh_n.nodes = mesh.nodes;
  mesh_n.elements.particle_velocity = struct(elem_type, struct("nodes", elem_nodes));
  mesh_n.material_data.c = 0;
  mesh_n.material_data.rho = 0;
  mesh_n.materials.particle_velocity = struct(elem_type, ones(rows(elem_nodes), 1, "int32"));

  dof_map_n = fem_ass_dof_map(mesh_n, load_case_dof_n);

  load_case_n = struct();

  mat_ass.n = fem_ass_matrix(mesh_n, ...
                             dof_map_n, ...
                             [FEM_VEC_SURFACE_NORMAL_VECTOR], ...
                             load_case_n);

  switch (elem_type)
    case {"tria6", "tria6h"}
      corner_nodes = int32([1, 2, 3]);
      num_nodes_elem = int32(15);
      bottom_node_idx = int32([1, 2, 3, 7, 8, 9]);
      top_node_idx = int32([4, 5, 6, 10, 11, 12]);
      interm_node_idx = int32([13, 14, 15]);
    otherwise
      error("elem_type \"%s\" not supported", elem_type);
  endswitch

  top_nodes = unique(elem_nodes(:));
  interm_nodes = unique(elem_nodes(:, corner_nodes)(:));

  elem_n = getfield(mat_ass.n, elem_type);

  nodes = [mesh.nodes;
           zeros(layers * (numel(top_nodes) + numel(interm_nodes)), 6)];

  elem = zeros(layers * rows(elem_nodes), num_nodes_elem, "int32");

  elem(1:rows(elem_nodes), bottom_node_idx) = elem_nodes;

  for k=1:layers
    for i=1:rows(elem_nodes)
      for j=1:columns(elem_nodes)
        node_idx_ij = find(elem_nodes(i, j) == top_nodes) + rows(mesh.nodes) + (k - 1) * (numel(top_nodes) + numel(interm_nodes));
        nodes(node_idx_ij, 1:3) = mesh.nodes(elem_nodes(i, j), 1:3) + sum(h(1:k)) * reshape(elem_n(i, j, :), 1, 3);
        elem(i + (k - 1) * rows(elem_nodes), top_node_idx(j)) = node_idx_ij;
      endfor
      for j=1:numel(corner_nodes)
        node_idx_ij = find(elem_nodes(i, corner_nodes(j)) == interm_nodes) + rows(mesh.nodes) + numel(top_nodes)  + (k - 1) * (numel(top_nodes) + numel(interm_nodes));
        nodes(node_idx_ij, 1:3) = mesh.nodes(elem_nodes(i, corner_nodes(j)), 1:3) + (sum(h(1:k - 1)) + 0.5 * h(k)) * reshape(elem_n(i, corner_nodes(j), :), 1, 3);
        elem(i + (k - 1) * rows(elem_nodes), interm_node_idx(j)) = node_idx_ij;
      endfor
    endfor
  endfor

  for k=2:layers
    elem((1:rows(elem_nodes)) + (k - 1) * rows(elem_nodes), bottom_node_idx) = elem((1:rows(elem_nodes)) + (k - 2) * rows(elem_nodes), top_node_idx);
  endfor
endfunction

%!test
%! ## TEST 1
%! mesh.nodes = [0.0, 0.0, -1.0, 0.0, 0.0, 0.0;
%!               1.0, 0.0, -1.0, 0.0, 0.0, 0.0;
%!               0.0, 1.0, -1.0, 0.0, 0.0, 0.0;
%!               0.5, 0.0, -1.0, 0.0, 0.0, 0.0;
%!               0.5, 0.5, -1.0, 0.0, 0.0, 0.0;
%!               0.0, 0.5, -1.0, 0.0, 0.0, 0.0];
%! e1 = rand(3, 1);
%! e2 = rand(3, 1);
%! e3 = cross(e1, e2);
%! e2 = cross(e3, e1);
%! R = [e1, e2, e3];
%! R /= diag(norm(R, "cols"));
%! mesh.nodes(:, 1:3) *= R.';
%! mesh.elements.tria6h = int32(1:6);
%! mesh.groups.tria6h.id = int32(1);
%! mesh.groups.tria6h.name = "surface";
%! mesh.groups.tria6h.elements = int32(1);
%! mesh.group.tria6h.nodes = unique(mesh.elements.tria6h);
%! h = 2;
%! N = 10;
%! [mesh.nodes, mesh.elements.penta15] = fem_pre_mesh_extrude_surf(mesh, "tria6h", int32(1), repmat(h / N, 1, N));
%! load_case_dof.locked_dof = false(rows(mesh.nodes), 6);
%! mesh.materials.penta15 = ones(rows(mesh.elements.penta15), 1, "int32");
%! mesh.material_data.E = 210000e6;
%! mesh.material_data.nu = 0.3;
%! mesh.material_data.rho = 7850;
%!
%! dof_map = fem_ass_dof_map(mesh, load_case_dof);
%! [mat_ass.K, ...
%!  mat_ass.M, ...
%!  mat_ass.mtot, ...
%!  mat_ass.mat_info, ...
%!  mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_MAT_STIFFNESS, ...
%!                                       FEM_MAT_MASS, ...
%!                                       FEM_SCA_TOT_MASS]);
%! sol_eig = fem_sol_modal(mesh, dof_map, mat_ass, 20);
%! mref = mesh.material_data.rho;
%! rtol = 10 * eps;
%! assert(mat_ass.mtot, mref, rtol * mref);
%! assert(all(sol_eig.f(1:6) < 1e-6 * max(sol_eig.f)));
%! assert(all(sol_eig.f(7:end) > 1e-1 * max(sol_eig.f)));

%!test
%! ## TEST 2
%! mesh.nodes = [-1, -1, -1, 0, 0, 0;
%!                1, -1, -1, 0, 0, 0;
%!                1,  1, -1, 0, 0, 0;
%!               -1,  1, -1, 0, 0, 0;
%!                0, -1, -1, 0, 0, 0;
%!                1,  0, -1, 0, 0, 0;
%!                0,  1, -1, 0, 0, 0;
%!               -1,  0, -1, 0, 0, 0;
%!                0,  0, -1, 0, 0, 0];
%! e1 = rand(3, 1);
%! e2 = rand(3, 1);
%! e3 = cross(e1, e2);
%! e2 = cross(e3, e1);
%! R = [e1, e2, e3];
%! R /= diag(norm(R, "cols"));
%! mesh.nodes(:, 1:3) *= R.';
%! mesh.elements.tria6h = int32([1,2,3,5,6,9;
%!                               1,3,4,9,7,8]);
%! mesh.groups.tria6h.id = int32(1);
%! mesh.groups.tria6h.name = "surface";
%! mesh.groups.tria6h.elements = int32([1, 2]);
%! mesh.group.tria6h.nodes = unique(mesh.elements.tria6h);
%! h = 2;
%! N = 5;
%! [mesh.nodes, mesh.elements.penta15] = fem_pre_mesh_extrude_surf(mesh, "tria6h", int32(1), repmat(h / N, 1, N));
%! load_case_dof.locked_dof = false(rows(mesh.nodes), 6);
%! mesh.materials.penta15 = ones(rows(mesh.elements.penta15), 1, "int32");
%! mesh.material_data.E = 210000e6;
%! mesh.material_data.nu = 0.3;
%! mesh.material_data.rho = 7850;
%!
%! dof_map = fem_ass_dof_map(mesh, load_case_dof);
%! [mat_ass.K, ...
%!  mat_ass.M, ...
%!  mat_ass.mtot, ...
%!  mat_ass.mat_info, ...
%!  mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_MAT_STIFFNESS, ...
%!                                       FEM_MAT_MASS, ...
%!                                       FEM_SCA_TOT_MASS]);
%! sol_eig = fem_sol_modal(mesh, dof_map, mat_ass, 20);
%! mref = 8 * mesh.material_data.rho;
%! rtol = 10 * eps;
%! assert(mat_ass.mtot, mref, rtol * mref);
%! assert(all(sol_eig.f(1:6) < 1e-6 * max(sol_eig.f)));
%! assert(all(sol_eig.f(7:end) > 1e-1 * max(sol_eig.f)));

%!test
%! ### TEST 3
%! ####################################################
%! ## Jont Allen
%! ## THE ACOUSTIC WAVE EQUATION AND SIMPLE SOLUTIONS
%! ## Chapter 5
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
%!     lambda = 10e-3 / unit_meters;
%!     deltaPML = 2 * lambda;
%!     r0 = 10e-3 / unit_meters;
%!     r1 = r0 + 3 * lambda;
%!     r2 = r1 + deltaPML;
%!     dx = lambda / 30;
%!     rho = 1000 / (unit_kilograms / unit_meters^3);
%!     c = 1440 / (unit_meters / unit_second);
%!     layersPML = ceil(deltaPML / dx);
%!     f = c / lambda;
%!     omega = 2 * pi * f;
%!     k = omega / c;
%!     f_enable_plot = true;
%!     Aref = (1 + 0j) / (unit_pascal * unit_meters);
%!     pref = @(r, t) Aref * exp(1j * (omega * t - k * r)) ./ r; ## according equation 5.11.6
%!     zref = @(r) rho * c * k * r .* exp(1j * acot(k * r)) ./ sqrt(1 + (k * r).^2); ## according equation 5.11.9
%!     vnref = @(r, t) Aref * exp(1j * (omega * t - k * r)) ./ (r * zref(r)); ## according equation 5.11.17
%!     alpha = atan2(2 * dx, r0);
%!     beta = atan2(2 * dx, r0);
%!     fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "r0=%g;\n", r0);
%!     fprintf(fd, "r1=%g;\n", r1);
%!     fprintf(fd, "alpha=%g;\n", alpha);
%!     fprintf(fd, "beta=%g;\n", beta);
%!     fprintf(fd, "dx = %g;\n", dx);
%!     fputs(fd, "Point(1) = {r0*Cos(alpha/2),r0*Sin(alpha/2),0,dx};\n");
%!     fputs(fd, "Point(2) = {r1*Cos(alpha/2),r1*Sin(alpha/2),0,dx};\n");
%!     fputs(fd, "Line(1) = {1,2};\n");
%!     fputs(fd, "tmp1[] = Extrude{{-Sin(alpha/2),Cos(alpha/2),0},{0,0,0},-beta/2} {\n");
%!     fputs(fd, "  Line{1};\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "tmp4[] = Extrude{{0,0,1},{0,0,0},-alpha/2} {\n");
%!     fputs(fd, "  Surface{tmp1[1]};\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "Physical Volume(\"v1\",1) = {tmp4[1]};\n");
%!     fputs(fd, "Physical Surface(\"s1\", 2) = {2};\n");
%!     fputs(fd, "Physical Surface(\"s2\", 3) = {3};\n");
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
%!   unlink([filename, ".geo"]);
%!   opt_mesh.elem_type = {"tet10h", "tria6h"};
%!   mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_mesh);
%!   unlink([filename, ".msh"]);
%!   [mesh.nodes, mesh.elements.penta15] = fem_pre_mesh_extrude_surf(mesh, "tria6h", 3, repmat(deltaPML / layersPML, 1, layersPML));
%!   r = sqrt(mesh.nodes(:, 1).^2 + mesh.nodes(:, 2).^2 + mesh.nodes(:, 3).^2);
%!   mesh.groups.tria6h(end + 1).id = 4;
%!   mesh.groups.tria6h(end).nodes = find(r > r2 - 0.1 * deltaPML / layersPML);
%!   grp_idx_s1 = find([mesh.groups.tria6h.id] == 2);
%!   grp_idx_s2 = find([mesh.groups.tria6h.id] == 3);
%!   grp_idx_s3 = find([mesh.groups.tria6h.id] == 4);
%!   mesh.material_data.rho = rho;
%!   mesh.material_data.c = c;
%!   mesh.materials.tet10h = ones(rows(mesh.elements.tet10h), 1, "int32");
%!   mesh.materials.penta15 = ones(rows(mesh.elements.penta15), 1, "int32");
%!   mesh.elements.particle_velocity.tria6h.nodes = mesh.elements.tria6h(mesh.groups.tria6h(grp_idx_s1).elements, :);
%!   mesh.materials.particle_velocity.tria6h = ones(rows(mesh.elements.particle_velocity.tria6h.nodes), 1, "int32");
%!   mesh.elements.acoustic_boundary.tria6h = mesh.elements.tria6h(mesh.groups.tria6h(grp_idx_s2).elements, :);
%!   mesh.materials.acoustic_boundary.tria6h = ones(rows(mesh.elements.acoustic_boundary.tria6h), 1, "int32");
%!   load_case = struct("particle_velocity", cell(1, 2));
%!   load_case(1).particle_velocity.tria6h.vn = repmat(-real(vnref(r0, 0)), size(mesh.elements.particle_velocity.tria6h.nodes));
%!   load_case(2).particle_velocity.tria6h.vn = repmat(-imag(vnref(r0, 0)), size(mesh.elements.particle_velocity.tria6h.nodes));
%!   load_case_dof.domain = FEM_DO_ACOUSTICS;
%!   load_case_dof.locked_dof = false(rows(mesh.nodes), 1);
%!   load_case_dof.locked_dof(mesh.groups.tria6h(grp_idx_s3).nodes, :) = true;
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!   mesh.elements.perfectly_matched_layers.penta15.f = zeros(3, columns(mesh.elements.penta15), rows(mesh.elements.penta15));
%!   mesh.elements.perfectly_matched_layers.penta15.e1 = zeros(3, columns(mesh.elements.penta15), rows(mesh.elements.penta15));
%!   mesh.elements.perfectly_matched_layers.penta15.e2 = zeros(3, columns(mesh.elements.penta15), rows(mesh.elements.penta15));
%!   for j=1:3
%!     mesh.elements.perfectly_matched_layers.penta15.f(j, :, :) = 1;
%!   endfor
%!   xi = mesh.nodes(:, 1)(mesh.elements.penta15.');
%!   yi = mesh.nodes(:, 2)(mesh.elements.penta15.');
%!   zi = mesh.nodes(:, 3)(mesh.elements.penta15.');
%!   ri = sqrt(xi.^2 + yi.^2 + zi.^2) - r1;
%!   mesh.elements.perfectly_matched_layers.penta15.e1(1, :, :) = xi;
%!   mesh.elements.perfectly_matched_layers.penta15.e1(2, :, :) = yi;
%!   mesh.elements.perfectly_matched_layers.penta15.e1(3, :, :) = zi;
%!   mesh.elements.perfectly_matched_layers.penta15.e2(2, :, :) = 1;
%!   sigmax = 1 ./ (deltaPML - ri);
%!   mesh.elements.perfectly_matched_layers.penta15.f(1, :, :) = 1 ./ (1 - 1j * sigmax / k);
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
%!   opt_sol.number_of_threads = int32(4);
%!   opt_sol.solver = "pastix";
%!   opt_sol.refine_max_iter = int32(50);
%!   Keff = -omega^2 * complex(mat_ass.Ma_re, mat_ass.Ma_im) + 1j * omega * complex(mat_ass.Da_re, mat_ass.Da_im) + complex(mat_ass.Ka_re, mat_ass.Ka_im);
%!   Reff = complex(mat_ass.Ra(:, 1), mat_ass.Ra(:, 2));
%!   Phi = fem_sol_factor(Keff, opt_sol) \ Reff;
%!   Psi = linspace(0, 2 * pi, 37);
%!   idx = dof_map.ndof(:, 1);
%!   iact = find(idx > 0);
%!   solC.Phi = solC.PhiP = zeros(rows(mesh.nodes), columns(Phi));
%!   solC.Phi(iact, :) = Phi(idx(iact), :);
%!   solC.PhiP(iact, :) = 1j * omega * Phi(idx(iact), :);
%!   [solC.particle_velocity, ...
%!    solC.acoustic_intensity] = fem_ass_matrix(mesh, ...
%!                                              dof_map, ...
%!                                              [FEM_VEC_PARTICLE_VELOCITY_C, ...
%!                                               FEM_SCA_ACOUSTIC_INTENSITY_C], ...
%!                                              load_case, ...
%!                                              solC);

%!   p = -solC.PhiP(mesh.elements.acoustic_boundary.tria6h);
%!   vx = solC.particle_velocity.vn.tria6h;
%!   R = abs(mean(mean((p - zref(r1) * vx) ./ (p + zref(r1) * vx))));
%!   T = 1 / (1 - R^2);
%!   P = sum(solC.acoustic_intensity.P.tria6h);
%!   P0 = 1e-12 / unit_watt;
%!   fprintf(stderr, "LW=%.2fdB TL=%.2fdB\n", 10 * log10(P/P0), 10 * log10(T));
%!   solR.t = Psi / omega;
%!   solR.p = real(-solC.PhiP * exp(1j * Psi));
%!   rg = sqrt(mesh.nodes(:, 1).^2 + mesh.nodes(:, 2).^2 + mesh.nodes(:, 3).^2);
%!   [rg, idx] = sort(rg);
%!   if (f_enable_plot)
%!     for j=1:columns(solR.p)
%!       figure("visible", "off");
%!       hold on;
%!       plot(rg * unit_meters, solR.p(idx, j) * unit_pascal, "-;p(r);1");
%!       plot(rg * unit_meters, real(pref(rg, solR.t(j))) * unit_pascal, "-;pref(r);0");
%!       xlabel("x [m]");
%!       ylabel("p [Pa]");
%!       grid on;
%!       grid minor on;
%!       title(sprintf("pressure t=%.2fs", solR.t(j)));
%!       ylim([-1, 1] * max(max(abs(solR.p))) * unit_pascal);
%!     endfor
%!   endif
%!   idx2 = find(rg < r1);
%!   rg = rg(idx2);
%!   idx = idx(idx2);
%!   tol = 1e-2;
%!   for j=1:numel(solR.t)
%!     prefj = real(pref(rg, solR.t(j)));
%!     assert(solR.p(idx, j), prefj, tol * max(abs(prefj)));
%!   endfor
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ### TEST 4
%! ####################################################
%! ## Jont Allen
%! ## THE ACOUSTIC WAVE EQUATION AND SIMPLE SOLUTIONS
%! ## Chapter 5
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
%!     lambda = 20e-3 / unit_meters;
%!     deltaPML = 1e-4 * lambda;
%!     r0 = 10e-3 / unit_meters;
%!     r1 = r0 + lambda / 2;
%!     r2 = r1 + 50 * lambda;
%!     r3 = r2 + deltaPML;
%!     dx = lambda / 10;
%!     rho = 1000 / (unit_kilograms / unit_meters^3);
%!     c = 1440 / (unit_meters / unit_second);
%!     layersPML = max([1, ceil(deltaPML / dx)]);
%!     f = c / lambda;
%!     omega = 2 * pi * f;
%!     k = omega / c;
%!     f_enable_plot = true;
%!     Aref = (1 + 0j) / (unit_pascal * unit_meters);
%!     pref = @(r, t) Aref * exp(1j * (omega * t - k * r)) ./ r; ## according equation 5.11.6
%!     zref = @(r) rho * c * k * r .* exp(1j * acot(k * r)) ./ sqrt(1 + (k * r).^2); ## according equation 5.11.9
%!     vnref = @(r, t) Aref * exp(1j * (omega * t - k * r)) ./ (r * zref(r)); ## according equation 5.11.17
%!     alpha = atan2(2 * dx, r0);
%!     beta = atan2(2 * dx, r0);
%!     fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "r0=%g;\n", r0);
%!     fprintf(fd, "r1=%g;\n", r1);
%!     fprintf(fd, "alpha=%g;\n", alpha);
%!     fprintf(fd, "beta=%g;\n", beta);
%!     fprintf(fd, "dx = %g;\n", dx);
%!     fputs(fd, "Point(1) = {r0*Cos(alpha/2),r0*Sin(alpha/2),0,dx};\n");
%!     fputs(fd, "Point(2) = {r1*Cos(alpha/2),r1*Sin(alpha/2),0,dx};\n");
%!     fputs(fd, "Line(1) = {1,2};\n");
%!     fputs(fd, "tmp1[] = Extrude{{-Sin(alpha/2),Cos(alpha/2),0},{0,0,0},-beta/2} {\n");
%!     fputs(fd, "  Line{1};\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "tmp4[] = Extrude{{0,0,1},{0,0,0},-alpha/2} {\n");
%!     fputs(fd, "  Surface{tmp1[1]};\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "Physical Volume(\"v1\",1) = {tmp4[1]};\n");
%!     fputs(fd, "Physical Surface(\"s1\", 2) = {2};\n");
%!     fputs(fd, "Physical Surface(\"s2\", 3) = {3};\n");
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
%!   unlink([filename, ".geo"]);
%!   opt_mesh.elem_type = {"tet10h", "tria6h"};
%!   mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_mesh);
%!   unlink([filename, ".msh"]);
%!   N = ceil((r2 - r1) / dx);
%!   [mesh.nodes, mesh.elements.penta15] = fem_pre_mesh_extrude_surf(mesh, "tria6h", 3, [repmat((r2 - r1) / N, 1, N), r3 - r2]);
%!   r = sqrt(mesh.nodes(:, 1).^2 + mesh.nodes(:, 2).^2 + mesh.nodes(:, 3).^2);
%!   mesh.groups.tria6h(end + 1).id = 4;
%!   mesh.groups.tria6h(end).nodes = find(r >= r3 - deltaPML / 2);
%!   grp_idx_s1 = find([mesh.groups.tria6h.id] == 2);
%!   grp_idx_s2 = find([mesh.groups.tria6h.id] == 3);
%!   grp_idx_s3 = find([mesh.groups.tria6h.id] == 4);
%!   node_idx_constr = mesh.groups.tria6h(grp_idx_s1).nodes;
%!   mesh.material_data.rho = rho;
%!   mesh.material_data.c = c;
%!   mesh.materials.tet10h = ones(rows(mesh.elements.tet10h), 1, "int32");
%!   mesh.materials.penta15 = ones(rows(mesh.elements.penta15), 1, "int32");
%!   mesh.elements.acoustic_boundary.tria6h = mesh.elements.tria6h(mesh.groups.tria6h(grp_idx_s2).elements, :);
%!   mesh.materials.acoustic_boundary.tria6h = ones(rows(mesh.elements.acoustic_boundary.tria6h), 1, "int32");
%!   mesh.elements.acoustic_constr = struct("C", mat2cell(ones(1, numel(node_idx_constr)), 1, ones(1, numel(node_idx_constr))), ...
%!                                          "nodes", mat2cell(node_idx_constr, 1, ones(1, numel(node_idx_constr))), ...
%!                                          "scale", mat2cell(repmat(1/omega, 1, numel(node_idx_constr)), 1, ones(1, numel(node_idx_constr))));
%!   load_case = struct("acoustic_constr", cell(1, 2));
%!   p_constr = repmat(pref(r0, 0), size(node_idx_constr));
%!   load_case(1).acoustic_constr = struct("p", mat2cell(real(p_constr), 1, ones(1, numel(p_constr))));
%!   load_case(2).acoustic_constr = struct("p", mat2cell(imag(p_constr), 1, ones(1, numel(p_constr))));
%!   load_case_dof.domain = FEM_DO_ACOUSTICS;
%!   load_case_dof.locked_dof = false(rows(mesh.nodes), 1);
%!   load_case_dof.locked_dof(mesh.groups.tria6h(grp_idx_s3).nodes, :) = true;
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!   mesh.elements.perfectly_matched_layers.penta15.f = zeros(3, columns(mesh.elements.penta15), rows(mesh.elements.penta15));
%!   mesh.elements.perfectly_matched_layers.penta15.e1 = zeros(3, columns(mesh.elements.penta15), rows(mesh.elements.penta15));
%!   mesh.elements.perfectly_matched_layers.penta15.e2 = zeros(3, columns(mesh.elements.penta15), rows(mesh.elements.penta15));
%!   for j=1:3
%!     mesh.elements.perfectly_matched_layers.penta15.f(j, :, :) = 1;
%!   endfor
%!   xi = mesh.nodes(:, 1)(mesh.elements.penta15.');
%!   yi = mesh.nodes(:, 2)(mesh.elements.penta15.');
%!   zi = mesh.nodes(:, 3)(mesh.elements.penta15.');
%!   vi = sqrt(xi.^2 + yi.^2 + zi.^2) - r2;
%!   [ridx, cidx] = find(vi >= 0);
%!   mesh.elements.perfectly_matched_layers.penta15.e1(1, :, :) = xi;
%!   mesh.elements.perfectly_matched_layers.penta15.e1(2, :, :) = yi;
%!   mesh.elements.perfectly_matched_layers.penta15.e1(3, :, :) = zi;
%!   mesh.elements.perfectly_matched_layers.penta15.e2(2, :, :) = 1;
%!   mesh.elements.perfectly_matched_layers.penta15.f(1, ridx, cidx) = 1j * k * (deltaPML - vi(ridx, cidx));
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
%!   opt_sol.number_of_threads = int32(4);
%!   opt_sol.solver = "pastix";
%!   opt_sol.refine_max_iter = int32(50);
%!   Keff = -omega^2 * complex(mat_ass.Ma_re, mat_ass.Ma_im) + 1j * omega * complex(mat_ass.Da_re, mat_ass.Da_im) + complex(mat_ass.Ka_re, mat_ass.Ka_im);
%!   Reff = complex(mat_ass.Ra(:, 1), mat_ass.Ra(:, 2));
%!   Phi = fem_sol_factor(Keff, opt_sol) \ Reff;
%!   Psi = linspace(0, 2 * pi, 37);
%!   idx = dof_map.ndof(:, 1);
%!   iact = find(idx > 0);
%!   solC.Phi = solC.PhiP = zeros(rows(mesh.nodes), columns(Phi));
%!   solC.Phi(iact, :) = Phi(idx(iact), :);
%!   solC.PhiP(iact, :) = 1j * omega * Phi(idx(iact), :);
%!   [solC.particle_velocity, ...
%!    solC.acoustic_intensity] = fem_ass_matrix(mesh, ...
%!                                              dof_map, ...
%!                                              [FEM_VEC_PARTICLE_VELOCITY_C, ...
%!                                               FEM_SCA_ACOUSTIC_INTENSITY_C], ...
%!                                              load_case, ...
%!                                              solC);

%!   p = -solC.PhiP(mesh.elements.acoustic_boundary.tria6h);
%!   vx = solC.particle_velocity.vn.tria6h;
%!   R = abs(mean(mean((p - zref(r1) * vx) ./ (p + zref(r1) * vx))));
%!   T = 1 / (1 - R^2);
%!   P = sum(solC.acoustic_intensity.P.tria6h);
%!   P0 = 1e-12 / unit_watt;
%!   fprintf(stderr, "LW=%.2fdB TL=%.2fdB\n", 10 * log10(P/P0), 10 * log10(T));
%!   solR.t = Psi / omega;
%!   solR.p = real(-solC.PhiP * exp(1j * Psi));
%!   rg = sqrt(mesh.nodes(:, 1).^2 + mesh.nodes(:, 2).^2 + mesh.nodes(:, 3).^2);
%!   [rg, idx] = sort(rg);
%!   if (f_enable_plot)
%!     for j=1:columns(solR.p)
%!       figure("visible", "off");
%!       hold on;
%!       plot(rg * unit_meters, solR.p(idx, j) * unit_pascal, "-;p(r);1");
%!       plot(rg * unit_meters, real(pref(rg, solR.t(j))) * unit_pascal, "-;pref(r);0");
%!       xlabel("x [m]");
%!       ylabel("p [Pa]");
%!       grid on;
%!       grid minor on;
%!       title(sprintf("pressure t=%.2fs", solR.t(j)));
%!       ylim([-1, 1] * max(max(abs(solR.p))) * unit_pascal);
%!     endfor
%!   endif
%!   idx2 = find(rg < r1);
%!   rg = rg(idx2);
%!   idx = idx(idx2);
%!   tol = 2e-2;
%!   pref = real(pref(rg, solR.t));
%!   assert(solR.p(idx, :), pref, tol * max(max(abs(pref))));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST 5
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
%!     d1 = 800e-3 / unit_meters;
%!     d2 = 120e-3 / unit_meters;
%!     d3 = 60e-3 / unit_meters;
%!     t = 3.5e-3 / unit_meters;
%!     r = 1500e-3 / unit_meters;
%!     h1 = 10 * t;
%!     h2 = 200e-3 / unit_meters;
%!     E1 = 210000e6 / unit_pascal;
%!     rho1 = 7800 / (unit_kilograms / unit_meters^3);
%!     nu1 = 0.3;
%!     alpha1 = 0.1826 / (unit_second^-1);
%!     beta1 = 5.0125e-6 / unit_second;
%!     c2 = 340 / (unit_meters / unit_second);
%!     rho2 = 1.225 / (unit_kilograms / unit_meters^3);
%!     Fz = (1 + 0j) / unit_newton;
%!     f = (100) / (unit_second^-1);
%!     solver = "precond";
%!     f_enable_PML = true;
%!     f_enable_plot = false;
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "d1 = %.16g;\n", d1);
%!     fprintf(fd, "d2 = %.16g;\n", d2);
%!     fprintf(fd, "d3 = %.16g;\n", d3);
%!     fprintf(fd, "t = %.16g;\n", t);
%!     fprintf(fd, "r = %.16g;\n", r);
%!     fprintf(fd, "h1 = %.16g;\n", h1);
%!     fprintf(fd, "h2 = %.16g;\n", h2);
%!     fputs(fd, "v1 = newv;\n");
%!     fputs(fd, "Cylinder(v1) = {0, 0, -0.5 * t, 0, 0, t, 0.5 * d1};\n");
%!     fputs(fd, "v2 = newv;\n");
%!     fputs(fd, "Cylinder(v2) = {0, 0, -0.5 * t, 0, 0, t, 0.5 * d2};\n");
%!     fputs(fd, "v3 = newv;\n");
%!     fputs(fd, "Cylinder(v3) = {0, 0, -0.5 * t, 0, 0, t, 0.5 * d3};\n");
%!     fputs(fd, "v4 = newv;\n");
%!     fputs(fd, "Sphere(v4) = {0, 0, 0, r};\n");
%!     fputs(fd, "v5 = BooleanFragments{Volume{v4};Delete;}{Volume{v1,v2,v3};Delete;};\n");
%!     fputs(fd, "Physical Volume(\"solid\", 1) = {9, 10};\n");
%!     fputs(fd, "Physical Volume(\"fluid\", 2) = {7, 8};\n");
%!     fputs(fd, "Physical Surface(\"fluid-struct\", 1) = {7, 14, 12, 13, 15, 11};\n");
%!     fputs(fd, "Physical Surface(\"clamp\", 2) = {14, 15};\n");
%!     fputs(fd, "Physical Surface(\"fluid-boundary\", 3) = {10};\n");
%!     fputs(fd, "Physical Point(\"load\", 1) = {9};\n");
%!     fputs(fd, "ReverseMesh Surface{7, 14, 12, 13, 15, 11};\n");
%!     fputs(fd, "Mesh.HighOrderOptimize = 1;\n");
%!     fputs(fd, "MeshSize{PointsOf{Volume{7, 8}; } } = h2;\n");
%!     fputs(fd, "MeshSize{PointsOf{Volume{9,10}; } } = h1;\n");
%!     fputs(fd, "Mesh.HighOrderIterMax = 100;\n");
%!     fputs(fd, "Mesh.HighOrderPassMax = 100;\n");
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
%!   unlink([filename, ".geo"]);
%!   mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%!   unlink([filename, ".msh"]);
%!   grp_idx_v_solid = find([mesh.groups.tet10.id] == 1);
%!   grp_idx_v_fluid = find([mesh.groups.tet10.id] == 2);
%!   grp_idx_s_fsi = find([mesh.groups.tria6.id] == 1);
%!   grp_idx_s_clamp = find([mesh.groups.tria6.id] == 2);
%!   grp_idx_s_boundary = find([mesh.groups.tria6.id] == 3);
%!   num_nodes_non_PML = rows(mesh.nodes);
%!   node_idx_constr_PML = num_nodes_non_PML + 1:rows(mesh.nodes);
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
%!   mesh.materials.tet10(mesh.groups.tet10(grp_idx_v_solid).elements) = 1;
%!   mesh.materials.tet10(mesh.groups.tet10(grp_idx_v_fluid).elements) = 2;
%!   mesh.elements.fluid_struct_interface.tria6 = mesh.elements.tria6(mesh.groups.tria6(grp_idx_s_fsi).elements, :);
%!   mesh.elements.acoustic_boundary.tria6 = mesh.elements.tria6(mesh.groups.tria6(grp_idx_s_boundary).elements, :);
%!   mesh.materials.acoustic_boundary.tria6 = repmat(int32(2), rows(mesh.elements.acoustic_boundary.tria6), 1);
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
%!   P = T = zeros(1, numel(f));
%!   ITER = repmat(intmax(), 1, 2);
%!   MAXITER = 20;
%!   TOL = eps^0.5;
%!   RESTART = 10;
%!   MAXITERS = int32(100);
%!   TOLS = eps^0.5;
%!   perm = [];
%!   cpu_factor = 0;
%!   cpu_solve = 0;
%!   FLAG = -1;
%!   alpha = 0.75;
%!   unwind_protect
%!     for i=1:numel(f)
%!       omega = 2 * pi * f(i);
%!       k = omega / c2;
%!       if (f_enable_PML)
%!         lambda = c2 / f(i);
%!         N = 0;
%!         NPML = 1;
%!         deltaPML = 1 / k;
%!         h = lambda / 10;
%!         rPML = r + h * N;
%!         [mesh.nodes, mesh.elements.penta15] = fem_pre_mesh_extrude_surf(setfield(mesh, "nodes", mesh.nodes(1:num_nodes_non_PML, :)), "tria6", 3, [repmat(h, 1, N), repmat(deltaPML / NPML, 1, NPML)]);
%!         mesh.materials.penta15 = repmat(int32(2), rows(mesh.elements.penta15), 1);
%!         mesh.elements.perfectly_matched_layers.penta15.f = zeros(3, columns(mesh.elements.penta15), rows(mesh.elements.penta15));
%!         mesh.elements.perfectly_matched_layers.penta15.e1 = zeros(3, columns(mesh.elements.penta15), rows(mesh.elements.penta15));
%!         mesh.elements.perfectly_matched_layers.penta15.e2 = zeros(3, columns(mesh.elements.penta15), rows(mesh.elements.penta15));
%!         xi = mesh.nodes(:, 1)(mesh.elements.penta15.');
%!         yi = mesh.nodes(:, 2)(mesh.elements.penta15.');
%!         zi = mesh.nodes(:, 3)(mesh.elements.penta15.');
%!         vi = sqrt(xi.^2 + yi.^2 + zi.^2) - rPML;
%!         mesh.elements.perfectly_matched_layers.penta15.e1(1, :, :) = xi;
%!         mesh.elements.perfectly_matched_layers.penta15.e1(2, :, :) = yi;
%!         mesh.elements.perfectly_matched_layers.penta15.e1(3, :, :) = zi;
%!         mesh.elements.perfectly_matched_layers.penta15.e2(1, :, :) = 1;
%!         for j=1:3
%!           mesh.elements.perfectly_matched_layers.penta15.f(j, :, :) = 1;
%!         endfor
%!         idxPML = find(mean(vi, 1) >= 0);
%!         sigmax = 1 ./ (deltaPML - vi(:, idxPML)) - 1 / deltaPML;
%!         mesh.elements.perfectly_matched_layers.penta15.f(1, :, idxPML) = 1 ./ (1 - 1j * sigmax / k);
%!       endif
%!       load_case_dof.locked_dof = false(rows(mesh.nodes), 7);
%!       load_case_dof.locked_dof(mesh.groups.tria6(grp_idx_s_clamp).nodes, 1:3) = true;
%!       if (f_enable_PML)
%!         load_case_dof.locked_dof(node_idx_constr_PML, 7) = true;
%!       endif
%!       load_case_dof.domain = FEM_DO_FLUID_STRUCT;
%!       dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!       dof_map.parallel.threads_ass = int32(4);
%!       Phi = zeros(dof_map.totdof, 1);
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
%!       opt_sol.number_of_threads = int32(6);
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
%!       z = rho2 * c2 * k * r .* exp(1j * acot(k * r)) ./ sqrt(1 + (k * r).^2); ## according Jont Allen equation 5.11.9
%!       p = sol.p(mesh.elements.acoustic_boundary.tria6);
%!       vx = sol.particle_velocity.vn.tria6;
%!       R = abs(mean(mean((p - z * vx) ./ (p + z * vx))));
%!       T(i) = 1 ./ (1 - R.^2);
%!       P(i) = sum(sol.acoustic_intensity.P.tria6);
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
%!         Zc = linspace(0, r + deltaPML , 300);
%!         Xc = Yc = zeros(size(Zc));
%!         deltaC = 0.4 / unit_meters;
%!         delta0 = 0.1 / unit_meters;
%!         idxc = find((mesh.nodes(:, 3) >= 0) & (sqrt(mesh.nodes(:, 1).^2 + mesh.nodes(:, 2).^2) < abs(mesh.nodes(:, 3)) * (deltaC / (r + deltaPML)) + delta0));
%!         pc = zeros(numel(Zc), columns(solt.p));
%!         for j=1:columns(pc)
%!           pc(:, j) = griddata3(mesh.nodes(idxc, 1), mesh.nodes(idxc, 2), mesh.nodes(idxc, 3), solt.p(idxc, j), Xc, Yc, Zc);
%!         endfor
%!         for j=1:columns(pc)
%!           figure("visible","off");hold on
%!           plot(Zc * unit_meters, pc(:, j) * unit_pascal, "-;p(X=0,Y=0,Z);1");
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
%!   plot(f * unit_second^-1, 10 * log10(P / P0), "-;P(f);1");
%!   xlabel("f [Hz]");
%!   ylabel("Lw [dB]");
%!   grid on;
%!   grid minor on;
%!   title("sound power level");
%!   figure("visible", "off");
%!   hold on;
%!   plot(f * unit_second^-1, 10 * log10(T), "-;TL(f);1");
%!   xlabel("f [Hz]");
%!   ylabel("TL [dB]");
%!   title("transmission loss");
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
