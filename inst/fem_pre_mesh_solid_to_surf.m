## Copyright (C) 2011(-2024) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{mesh}] = fem_pre_mesh_solid_to_surf(@var{mesh})
## Create a surface mesh from solid elements
##
## @var{mesh} @dots{} on input the mesh must containt solid elements, on output the mesh will contain surface elements at all external boundaries
##
## @end deftypefn

function [mesh] = fem_pre_mesh_solid_to_surf(mesh)
  quad8 = zeros(rows(mesh.elements.penta15) * 3, 8);

  idx_node = int32([1, 2, 5, 4, 7, 14, 10, 13;
                    2, 3, 6, 5, 14, 8, 15, 11;
                    3, 1, 4, 6, 9, 13, 12, 15]);

  for i=1:rows(mesh.elements.penta15)
    for j=1:rows(idx_node)
      quad8(rows(idx_node) * (i - 1) + j, :) = mesh.elements.penta15(i, idx_node(j, :));
    endfor
  endfor

  [quad8sort1] = sort(quad8, 2);
  [quad8sort2, idx] = sortrows(quad8sort1, 1:columns(quad8));

  duplicate = false(rows(quad8), 1);

  for i=2:rows(quad8sort2)
    if (all(quad8sort2(i, :) == quad8sort2(i - 1, :)))
      duplicate(idx([i, i - 1])) = true;
    endif
  endfor

  mesh.elements.quad8 = quad8(~duplicate, :);


  tria6 = zeros(rows(mesh.elements.penta15) * 3, 6);

  idx_node = int32([4, 5, 6, 10, 11, 12;
                    1, 3, 2,  9,  8,  7]);

  for i=1:rows(mesh.elements.penta15)
    for j=1:rows(idx_node)
      tria6(rows(idx_node) * (i - 1) + j, :) = mesh.elements.penta15(i, idx_node(j, :));
    endfor
  endfor

  [tria6sort1] = sort(tria6, 2);
  [tria6sort2, idx] = sortrows(tria6sort1, 1:columns(tria6));

  duplicate = false(rows(tria6), 1);

  for i=2:rows(tria6sort2)
    if (all(tria6sort2(i, :) == tria6sort2(i - 1, :)))
      duplicate(idx([i, i - 1])) = true;
    endif
  endfor

  mesh.elements.tria6 = tria6(~duplicate, :);
endfunction
