## Copyright (C) 2018(-2023) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{load_case}, @var{bearing_surf}, @var{idx_group}] = fem_ehd_pre_comp_mat_load_case(@var{mesh}, @var{bearing_surf})
## Build pressure load cases required for fem_ehd_pre_comp_mat_unstruct.
##
## @var{mesh} @dots{} Return value from fem_pre_mesh_import or fem_pre_mesh_unstruct_create.
##
## @var{bearing_surf} @dots{} Struct array describing all bearing surfaces of this @var{mesh}.
##
## @var{bearing_surf}.group_idx @dots{} Index of a group of triangular elements in @var{mesh}.groups.tria6 used for applying pressure loads.
##
## @var{bearing_surf}.relative_tolerance @dots{} Optional field. If present, it will be checked if the nodes present in the mesh are covering the full bearing width.
##
## @var{bearing_surf}.absolute_tolerance @dots{} Optional field. See also @var{bearing_surf}.relative_tolerance.
##
## @var{bearing_surf}.X0 @dots{} Optional centre of the cylindrical bearing surface. See also @var{bearing_surf}.relative_tolerance.
##
## @var{bearing_surf}.R @dots{} Optional orientation of the cylindrical bearing surface. R(:, 3) will be the axis of the cylinder. See also @var{bearing_surf}.relative_tolerance.
##
## @var{bearing_surf}.w @dots{} Optional axial width of the cylindrical bearing surface. See also @var{bearing_surf}.relative_tolerance.
##
## @var{bearing_surf}.options.mesh_size @dots{} Mesh size used for the hydrodynamic bearing element.
##
## @var{bearing_surf}.options.reference_pressure @dots{} Define the unit pressure applied to bearing surfaces.
##
## @var{bearing_surf}.options.bearing_model @dots{} This value may be "EHD/FE" or "EHD/FD". In case of "EHD/FE" the number of output nodes will valid for a quadratic mesh.
##
## @seealso{fem_ehd_pre_comp_mat_unstruct}
## @end deftypefn

function [load_case, bearing_surf, idx_group] = fem_ehd_pre_comp_mat_load_case(mesh, bearing_surf, options)
  if (nargin < 2 || nargin > 3 || nargout < 2 || nargout > 3)
    print_usage();
  endif

  if (nargin < 3)
    options = struct();
  endif

  if (~isfield(options, "elem_type"))
    options.elem_type = "tria6";
  endif

  [bearing_surf, idx] = fem_ehd_pre_comp_mat_grid(mesh, bearing_surf, options);

  load_case = fem_pre_load_case_create_empty(idx);

  idx_group = zeros(numel(bearing_surf) + 1, 1, "int32");
  idx = int32(0);

  for i=1:numel(bearing_surf)
    idx_group(i) = idx + 1;
    elements = getfield(mesh.elements, options.elem_type)(getfield(mesh.groups, options.elem_type)(bearing_surf(i).group_idx).elements, :);

    N = [numel(bearing_surf(i).grid_x), numel(bearing_surf(i).grid_z)];

    idx_dst_per = [       0,  N(1), N(1) + 1] + 1;
    idx_src_per = [N(1) - 1,     1,        2] + 1;

    Xn = bearing_surf(i).R.' * (mesh.nodes(:, 1:3).' - bearing_surf(i).X0);
    zi = Xn(3, :)(elements);
    xi = mod(atan2(Xn(2, :)(elements), Xn(1, :)(elements)), 2 * pi) * bearing_surf(i).r;

    xi_S = [bearing_surf(i).grid_x(end - 1) - 2 * pi * bearing_surf(i).r, ...
            bearing_surf(i).grid_x, ...
            bearing_surf(i).grid_x(2) + 2 * pi * bearing_surf(i).r];

    zi_S = [2 * bearing_surf(i).grid_z(1) - bearing_surf(i).grid_z(2), ...
            bearing_surf(i).grid_z, ...
            2 * bearing_surf(i).grid_z(end) - bearing_surf(i).grid_z(end - 1)];

    p_S = zeros(numel(zi_S), numel(xi_S));

    for k=1:N(1) - 1
      for j=1:N(2)
        p_S(:, :) = 0;
        p_S(j + 1, k + 1) = bearing_surf(i).options.reference_pressure;
        p_S(j + 1, idx_dst_per) = p_S(j + 1, idx_src_per);
        p_U = interp2(xi_S, zi_S, p_S, xi, zi, "linear");

        idx_press = any(p_U, 2);
        load_case(++idx).pressure = struct(options.elem_type, struct("elements", elements(idx_press, :), "p", p_U(idx_press, :)));
      endfor
    endfor
  endfor

  idx_group(end) = ++idx;

  for i=1:numel(bearing_surf)
    bearing_surf(i).idx_load_case = idx_group(i):idx_group(i + 1) - 1;
  endfor
endfunction

function [bearing_surf, idx] = fem_ehd_pre_comp_mat_grid(mesh, bearing_surf, options)
  for i=1:numel(bearing_surf)
    nodes = getfield(mesh.groups, options.elem_type)(bearing_surf(i).group_idx).nodes;
    elements = getfield(mesh.elements, options.elem_type)(getfield(mesh.groups, options.elem_type)(bearing_surf(i).group_idx).elements, :);

    if (isfield(bearing_surf(i), "relative_tolerance") && isfield(bearing_surf(i), "absolute_tolerance"))
      z = bearing_surf(i).R(:, 3).' * (mesh.nodes(nodes, 1:3).' - bearing_surf(i).X0);

      dz = abs([0.5 * bearing_surf(i).w - max(z);
                0.5 * bearing_surf(i).w + min(z)]);

      tol = bearing_surf(i).relative_tolerance * bearing_surf(i).w + bearing_surf(i).absolute_tolerance;

      if (any(dz > tol))
        error("z-coordinate of nodes at bearing %d does not cover the complete bearing width", i);
      endif
    endif

    clear z dz tol elements nodes;
  endfor

  idx = int32(0);

  for i=1:numel(bearing_surf)
    if (isfield(bearing_surf(i).options, "mesh_size"))
      dx = bearing_surf(i).options.mesh_size;
      Nx = round(2 * pi * bearing_surf(i).r / dx);
      Nz = round(bearing_surf(i).w / dx);
    else
      Nx = bearing_surf(i).options.number_of_nodes_x;
      Nz = bearing_surf(i).options.number_of_nodes_z;
    endif

    N = [max(4, Nx), max(3, Nz)];

    if (isfield(bearing_surf(i).options, "bearing_model") && ischar(bearing_surf(i).options.bearing_model))
      switch (bearing_surf(i).options.bearing_model)
        case "EHD/FE"
          for j=1:numel(N)
            if (mod(N(j), 2) == 0)
              ++N(j);
            endif
          endfor
      endswitch
    endif

    bearing_surf(i).grid_x = linspace(0, 2 * pi * bearing_surf(i).r, N(1));
    bearing_surf(i).grid_z = linspace(-0.5 * bearing_surf(i).w, 0.5 * bearing_surf(i).w, N(2));
    idx += (numel(bearing_surf(i).grid_x) - 1) * numel(bearing_surf(i).grid_z);
    clear dx N Nx Nz;
  endfor
endfunction

%!test
%! close all;
%! fd = -1;
%! filename = "";
%! unwind_protect
%! filename = tempname();
%! if (ispc())
%!   filename(filename == "\\") = "/";
%! endif
%! unwind_protect
%! [fd, msg] = fopen([filename, ".geo"], "w");
%! if (fd == -1)
%!   error("failed to open file \"%s.geo\"", filename);
%! endif
%! ri = 8e-3;
%! ro = 10e-3;
%! h = 12e-3;
%! c = 2e-3;
%! b = h - 2 * c;
%! scale_def = 5e-3;
%! mesh_size = 5e-3;
%! f_post_pro = false;
%! fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%! fprintf(fd, "ri = %g;\n", ri);
%! fprintf(fd, "ro = %g;\n", ro);
%! fprintf(fd, "h = %g;\n", h);
%! fprintf(fd, "c = %g;\n", c);
%! fputs(fd, "Point(1) = {ri,0.0,0.0};\n");
%! fputs(fd, "Point(2) = {ro,0.0,0.0};\n");
%! fputs(fd, "Point(3) = {ro,0.0,c};\n");
%! fputs(fd, "Point(4) = {ro,0.0,h - c};\n");
%! fputs(fd, "Point(5) = {ro,0.0,h};\n");
%! fputs(fd, "Point(6) = {ri,0.0,h};\n");
%! fputs(fd, "Point(7) = {ri,0.0,h - c};\n");
%! fputs(fd, "Point(8) = {ri,0.0,c};\n");
%! fputs(fd, "Line(1) = {1,2};\n");
%! fputs(fd, "Line(2) = {2,3};\n");
%! fputs(fd, "Line(3) = {3,4};\n");
%! fputs(fd, "Line(4) = {4,5};\n");
%! fputs(fd, "Line(5) = {5,6};\n");
%! fputs(fd, "Line(6) = {6,7};\n");
%! fputs(fd, "Line(7) = {7,8};\n");
%! fputs(fd, "Line(8) = {8,1};\n");
%! fputs(fd, "Line Loop(5) = {1,2,3,4,5,6,7,8};\n");
%! fputs(fd, "Plane Surface(6) = {5};\n");
%! fputs(fd, "tmp[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} { Surface{6}; };\n");
%! fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%! fputs(fd, "Physical Surface(\"clamp\",1) = {tmp[2]};\n");
%! fputs(fd, "Physical Surface(\"load1\",2) = {tmp[4]};\n");
%! fputs(fd, "Physical Surface(\"load2\",3) = {tmp[8]};\n");
%! unwind_protect_cleanup
%!  if (fd ~= -1)
%!    fclose(fd);
%!  endif
%! end_unwind_protect
%! fprintf(stderr, "meshing ...\n");

%! pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-clmin", sprintf("%g", 0.75 * mesh_size), "-clmax", sprintf("%g", 1.25 *mesh_size), [filename, ".geo"]});
%! status = spawn_wait(pid);
%! if (status ~= 0)
%!  error("gmsh failed with status %d", status);
%! endif

%! fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%! mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%! fprintf(stderr, "%d nodes\n", rows(mesh.nodes));
%! grp_id_clamp = find([[mesh.groups.tria6].id] == 1);
%! grp_id_p1 = find([[mesh.groups.tria6].id] == 3);
%! grp_id_p2 = find([[mesh.groups.tria6].id] == 2);
%! bearing_surf(1).group_idx = grp_id_p1;
%! bearing_surf(1).options.reference_pressure = 1e6;
%! bearing_surf(1).options.mesh_size = 20e-3;
%! bearing_surf(1).r = ri;
%! bearing_surf(1).w = b;
%! bearing_surf(1).X0 = [0; 0; b/2 + c];
%! bearing_surf(1).R = eye(3);
%! bearing_surf(1).relative_tolerance = 0;
%! bearing_surf(1).absolute_tolerance = sqrt(eps) * ri;
%! bearing_surf(2).group_idx = grp_id_p2;
%! bearing_surf(2).options.reference_pressure = 1e6;
%! bearing_surf(2).options.mesh_size = 20e-3;
%! bearing_surf(2).r = ro;
%! bearing_surf(2).w = b;
%! bearing_surf(2).X0 = [0; 0; b/2 + c];
%! bearing_surf(2).R = eye(3);
%! bearing_surf(2).relative_tolerance = 0;
%! bearing_surf(2).absolute_tolerance = sqrt(eps) * ri;
%! [load_case, bearing_surf, idx_group] = fem_ehd_pre_comp_mat_load_case(mesh, bearing_surf);
%! load_case(1).locked_dof = false(rows(mesh.nodes), 6);
%! load_case(1).locked_dof(mesh.groups.tria6(grp_id_clamp).nodes, :) = true;
%! mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%! E = 210000e6;
%! nu = 0.3;
%! mesh.material_data.rho = 7850;
%! mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%! dof_map = fem_ass_dof_map(mesh, load_case(1));
%! opt_solver.refine_max_iter = int32(0);
%! fprintf(stderr, "assembling matrices ...\n");
%! [mat_ass.K, ...
%!  mat_ass.R, ...
%!  mat_ass.Rlumped] = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_MAT_STIFFNESS, ...
%!                                     FEM_VEC_LOAD_CONSISTENT, ...
%!                                     FEM_VEC_LOAD_LUMPED], ...
%!                                    load_case);
%! fprintf(stderr, "solving for static deflection (consistent load) ...\n");
%! sol_stat = fem_sol_static(mesh, dof_map, mat_ass, opt_solver);
%! fprintf(stderr, "solving for static deflection (lumped load) ...\n");
%! sol_stat_lumped = fem_sol_static(mesh, dof_map, setfield(mat_ass, "R", mat_ass.Rlumped), opt_solver);
%! if (f_post_pro)
%! post_pro_file = [filename, "_post_pro.geo"];
%! for j=1:numel(idx_group) - 1
%!  max_def = 0;
%!  for i=idx_group(j):idx_group(j + 1) - 1
%!    max_def = max(max_def, max(norm(sol_stat.def(:, 1:3, i), "rows")));
%!  endfor
%! for i=0:idx_group(j + 1) - idx_group(j) - 1
%!   deformation_file = sprintf("%s_%d_%03d_post_pro.msh", filename, j, i + 1);
%!   unwind_protect
%!   fem_post_sol_step_export(deformation_file, sol_stat, idx_group(j) + i, i + 1, i + 1, 1);
%!   fd = -1;
%!   unwind_protect
%!     [fd, msg] = fopen(post_pro_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\": %s", post_pro_file, msg);
%!     endif
%!     fprintf(fd, "Merge \"%s\";\n", [filename, ".msh"]);
%!     fputs(fd, "Mesh.SurfaceEdges = 0;\n");
%!     fputs(fd, "Mesh.SurfaceFaces = 0;\n");
%!     fputs(fd, "Mesh.SurfaceNumbers = 0;\n");
%!     fputs(fd, "Mesh.VolumeEdges = 0;\n");
%!     fputs(fd, "Mesh.VolumeFaces = 0;\n");
%!     fprintf(fd, "Merge \"%s\";\n", deformation_file);
%!     fputs(fd, "View[0].Type = 1;\n");
%!     fputs(fd, "View[0].VectorType = 5;\n");
%!     fputs(fd, "View[0].Visible = 1;\n");
%!     fprintf(fd, "View[0].DisplacementFactor = %g;\n", scale_def / max_def);
%!     fputs(fd, "View[0].ShowTime = 1;\n");
%!     fputs(fd, "View[0].ShowElement = 1;\n");
%!     fputs(fd, "View[0].IntervalsType = 3;\n");
%!     fputs(fd, "View[0].NbIso = 20;\n");
%!     fputs(fd, "View[0].DrawSkinOnly = 1;\n");
%!     fputs(fd, "General.Trackball = 0;\n");
%!     fputs(fd, "General.RotationX = -60;\n");
%!     fputs(fd, "General.RotationY = 10;\n");
%!     fputs(fd, "General.RotationZ = 40;\n");
%!     fputs(fd, "General.Axes = 3;\n");
%!     fputs(fd, "General.Orthographic = 1;\n");
%!     fputs(fd, "General.RotationCenterGravity = 1;\n");
%!     fprintf(fd, "View[0].TimeStep = %d;\n", i);
%!     fprintf(fd, "Print \"%s_%d_%03d.jpg\";\n", filename, j, i + 1);
%!     fputs(fd, "Exit;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   pid = spawn("gmsh", {post_pro_file});
%!   status = spawn_wait(pid);
%!   unwind_protect_cleanup
%!     [~] = unlink(post_pro_file);
%!     [~] = unlink(deformation_file);
%!   end_unwind_protect
%!   if (0 ~= status)
%!     error("gmsh failed with status %d", status);
%!   endif
%! endfor
%! jpg_filenames = dir(sprintf("%s_%d_*.jpg", filename, j));
%! for i=1:numel(jpg_filenames)
%!   [img, map, alpha] = imread(fullfile(jpg_filenames(i).folder, jpg_filenames(i).name));
%!   figure("visible", "off");
%!   imshow(img, map);
%!   title(sprintf("load case %d", i));
%! endfor
%! endfor
%! figure_list();
%! endif
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       status = unlink(fullfile(fn(i).folder, fn(i).name));
%!       if (status ~= 0)
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect
