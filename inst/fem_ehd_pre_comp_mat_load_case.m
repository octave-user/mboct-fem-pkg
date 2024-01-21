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

