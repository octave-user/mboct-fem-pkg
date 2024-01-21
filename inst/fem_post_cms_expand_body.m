## Copyright (C) 2019(-2021) Reinhard <octave-user@a1.net>
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
## @deftypefn{Function File} @var{def} = fem_post_cms_expand_body(@var{mesh}, @var{dof_map}, @var{mat_ass}, @var{q})
## Convert a modal solution @var{Q} to a nodal solution @var{def}.
##
## @var{mesh} @dots{} Finite element mesh data structure
##
## @var{dof_map} @dots{} Degree of freedom mapping.
##
## @var{mat_ass}.Tred @dots{} Selected mode shapes of the reduced order model
##
## @var{q} @dots{} Matrix of modal displacements
##
## @var{def} @dots{} Three dimensional array of nodal displacements
##
## @seealso{fem_post_cms_expand}
##
## @end deftypefn

function def = fem_post_cms_expand_body(mesh, dof_map, mat_ass, q)
  if (nargin ~= 4 || nargout > 1)
    print_usage();
  endif
  
  U = zeros(dof_map.totdof, columns(q));
  U(dof_map.idx_node, :) = mat_ass.Tred * q;
  def = fem_post_def_nodal(mesh, dof_map, U);
endfunction

