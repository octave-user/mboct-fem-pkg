## Copyright (C) 2018(-2021) Reinhard <octave-user@a1.net>
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

  persistent elem_types = {"iso8", "iso20", "tet10", "penta15", "tet10h"};

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

%!demo
%! close all;
%! SI_unit_m = 1e-3;
%! SI_unit_kg = 1e3;
%! SI_unit_s = 1e-1;
%! SI_unit_N = SI_unit_kg * SI_unit_m / SI_unit_s^2;
%! SI_unit_Pa = SI_unit_N / SI_unit_m^2;
%! a = 150e-3 / SI_unit_m;
%! b = 20e-3 / SI_unit_m;
%! c = 45e-3 / SI_unit_m;
%! d = 10e-3 / SI_unit_m;
%! e = 10e-3 / SI_unit_m;
%! X = [ 0.5 * a,  0.5 * b,  0.5 * c;  ##  1
%!             0,  0.5 * b,  0.5 * c;  ##  2
%!             0, -0.5 * b,  0.5 * c;  ##  3
%!       0.5 * a, -0.5 * b,  0.5 * c;  ##  4
%!       0.5 * a,  0.5 * b, -0.5 * c;  ##  5
%!             0,  0.5 * b, -0.5 * c;  ##  6
%!             0, -0.5 * b, -0.5 * c;  ##  7
%!       0.5 * a, -0.5 * b, -0.5 * c,  ##  8
%!             a,  0.5 * b,  0.5 * c;  ##  9
%!             a, -0.5 * b,  0.5 * c;  ## 10
%!             a,  0.5 * b, -0.5 * c;  ## 11
%!             a, -0.5 * b, -0.5 * c,  ## 12
%!         a + d,        0,        0;  ## 13
%!            -e,        0,        0;  ## 14
%!         3 * a,    2 * b,  1.5 * c]; ## 15
%! mesh.nodes = [X, zeros(rows(X), 3)];
%! mesh.elements.iso8 = int32([1:8; 9, 1, 4, 10, 11, 5, 8, 12]);
%! mesh.materials.iso8 = int32([1; 1]);
%! mesh.elements.rbe3(1).nodes = int32([13, 9, 10, 11, 12]);
%! mesh.elements.rbe3(1).weight = ones(1, 4);
%! mesh.elements.rbe3(2).nodes = int32([14, 2, 3, 6, 7]);
%! mesh.elements.rbe3(2).weight = ones(1, 4);
%! dof_status = fem_cms_dof_active(mesh);
%! assert(all(all(dof_status(1:14, 1:3))));
%! assert(all(all(dof_status(13:14, 1:6))));
%! assert(~any(dof_status(15, :)));
