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
## @deftypefn{Function File} [@var{mesh_split}] = fem_pre_mesh_elem_split(@var{mesh})
## @deftypefnx{} [@dots{}] = fem_pre_mesh_elem_split(@var{mesh}, @var{options})
## Convert eight node hexahedra to six node pentahedra by splitting,
## if the maximum corner angle is above certain limit.
##
## @var{mesh} @dots{} Original mesh with distorted elements.
##
## @var{options}.iso8.max_angle @dots{} Limiting value for corner angles which causes splitting.
##
## @var{mesh_split} @dots{} New mesh with reduced corner angles.
## @seealso{fem_pre_mesh_elem_angle}
## @end deftypefn

function mesh_split = fem_pre_mesh_elem_split(mesh, options)
  if (nargin < 1 || nargin > 2 || nargout > 1)
    print_usage();
  endif

  if (nargin < 2)
    options = struct();
  endif

  if (~isfield(options, "iso8"))
    options.iso8 = struct();
  endif

  if (~isfield(options.iso8, "max_angle"))
    options.iso8.max_angle = 175 * pi / 180;
  endif

  mesh_split = mesh;
  angle_data = fem_pre_mesh_elem_angle(mesh);

  cond_split = any(mean(angle_data.iso8, 3) > options.iso8.max_angle, 2);
  idx_elem_keep = find(~cond_split);
  idx_elem_split = find(cond_split);

  inum_elem = numel(idx_elem_keep);
  mesh_split.elements.iso8 = zeros(numel(idx_elem_keep) + 2 * numel(idx_elem_split), columns(mesh.elements.iso8), "int32");
  mesh_split.materials.iso8 = zeros(numel(idx_elem_keep) + 2 * numel(idx_elem_split), 1, "int32");

  if (numel(idx_elem_keep))
    mesh_split.elements.iso8(1:inum_elem, :) = mesh.elements.iso8(idx_elem_keep, :);
    mesh_split.materials.iso8(1:inum_elem, :) = mesh.materials.iso8(idx_elem_keep);
  endif

  for i=1:numel(idx_elem_split)
    Phi = mean(angle_data.iso8(idx_elem_split(i), :, :), 3);

    idx_max_angle = find(Phi == max(Phi));

    if (numel(idx_max_angle))
      switch (idx_max_angle(1))
        case {1, 3}
          idx_nodes = [1, 2, 3, 3, 5, 6, 7, 7;
                       1, 3, 4, 4, 5, 7, 8, 8];
        case {2, 4}
          idx_nodes = [2, 3, 4, 4, 6, 7, 8, 8;
                       1, 2, 4, 4, 5, 6, 8, 8];
        case {5, 7}
          idx_nodes = [5, 6, 2, 2, 8, 7, 3, 3;
                       5, 2, 1, 1, 8, 3, 4, 4];
        case {6, 8}
          idx_nodes = [6, 2, 1, 1, 7, 3, 4, 4;
                       5, 6, 1, 1, 8, 7, 4, 4];
        case {9, 11}
          idx_nodes = [5, 1, 4, 4, 6, 2, 3, 3;
                       5, 4, 8, 8, 6, 3, 7, 7];
        case {10, 12}
          idx_nodes = [1, 4, 8, 8, 2, 3, 7, 7;
                       1, 8, 5, 5, 2, 7, 6, 6];
        otherwise
          error("invalid node index %d detected", idx_max_angle(1));
      endswitch

      for j=1:rows(idx_nodes)
        mesh_split.elements.iso8(++inum_elem, :) = mesh.elements.iso8(idx_elem_split(i), idx_nodes(j, :));
        mesh_split.materials.iso8(inum_elem) = mesh.materials.iso8(idx_elem_split(i));
      endfor
      continue
    endif
  endfor
endfunction

