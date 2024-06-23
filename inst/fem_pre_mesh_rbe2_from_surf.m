## Copyright (C) 2018(-2024) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{joints}] = fem_pre_mesh_rbe2_from_surf(@var{mesh}, @var{group_id}, @var{master_node_idx})
## @deftypefnx {} [@dots{}] = fem_pre_mesh_rbe2_from_surf(@dots{}, @var{elem_type})
##
## Builds joint elements from specified groups of tria6 elements
##
## @var{mesh} @dots{} Finite element mesh data structure
##
## @var{group_id} @dots{} array of group id-numbers used to identify surface elements with slave nodes
##
## @var{master_node_idx} @dots{} array of master node indices for rbe3 elements
##
## @var{elem_type} @dots{} the element type addressed by @var{group_id} (e.g. "tria6" or "iso4")
##
## @end deftypefn

function joints = fem_pre_mesh_rbe2_from_surf(mesh, group_id, master_node_idx, elem_type)
  if (nargin < 3 || nargin > 4 || nargout > 1)
    print_usage();
  endif

  if (nargin < 4)
    elem_type = "tria6";
  endif

  if (~iscell(elem_type))
    elem_type = {elem_type};
  endif

  rbe2 = fem_pre_mesh_rbe3_from_surf(mesh, group_id, master_node_idx, elem_type);

  X = mesh.nodes(rbe2.nodes, 1:3).';
  Xm = X(:, 1);
  Xs = X(:, 2:end);
  ls = Xs - Xm;

  empty_cell = cell(1, columns(Xs));

  joints = struct("nodes", empty_cell, "C", empty_cell);

  for i=1:numel(joints)
    joints(i).nodes = rbe2.nodes([1, i + 1]);
    Phi = [eye(3),      -skew(ls(:, i))];

    joints(i).C = [Phi, -eye(3), zeros(3, 3)];
  endfor
endfunction
