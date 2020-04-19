## Copyright (C) 2018(-2020) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} @var{dof_in_use} = fem_cms_dof_active(@var{mesh})
## Return a boolean matrix @var{dof_in_use} which indicates active degrees of freedom used by elements.
##
## @var{mesh} @dots{} Finite element mesh data structure
##
## @end deftypefn

function dof_in_use = fem_cms_dof_active(mesh)
  if (nargin ~= 1 || nargout > 1)
    print_usage();
  endif
  
  dof_in_use = false(rows(mesh.nodes), columns(mesh.nodes));

  persistent elem_types = {"iso8", "tet10"};

  for i=1:numel(elem_types)
    if (isfield(mesh.elements, elem_types{i}))
      elem_nodes = getfield(mesh.elements, elem_types{i});
      for j=1:columns(elem_nodes)
        dof_in_use(elem_nodes(:, j), 1:3) = true;
      endfor
    endif
  endfor

  if (isfield(mesh.elements, "rbe3"))
    for i=1:numel(mesh.elements.rbe3)
      dof_in_use(mesh.elements.rbe3(i).nodes(1), 1:6) = true;
      dof_in_use(mesh.elements.rbe3(i).nodes(2:end), 1:3) = true;
    endfor
  endif
endfunction
