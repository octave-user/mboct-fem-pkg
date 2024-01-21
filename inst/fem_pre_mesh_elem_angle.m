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
## @deftypefn{Function File} @var{angle_data} = fem_pre_mesh_elem_angle(@var{mesh})
## Compute corner angles for eight node hexahedrons
##
## @var{mesh} @dots{} Finite element mesh data structure containing iso8 elements
##
## @var{angle_data} @dots{} scalar struct contains computed corner angles
## @seealso{fem_pre_mesh_elem_split}
## @end deftypefn

function angle_data = fem_pre_mesh_elem_angle(mesh)
  if (nargin ~= 1 || nargout > 1)
    print_usage();
  endif

  angle_data = struct();
  
  if (isfield(mesh.elements, "iso8"))
    ## angle:a        b         #  edge index
    persistent edges = int32([4, 1, 2, 8, 5, 6;  #  1
                              1, 2, 3, 5, 6, 7;  #  2
                              2, 3, 4, 6, 7, 8;  #  3
                              3, 4, 1, 7, 8, 5;  #  4
                              1, 5, 6, 4, 8, 7;  #  5
                              5, 6, 2, 8, 7, 3;  #  6
                              6, 2, 1, 7, 3, 4;  #  7
                              2, 1, 5, 3, 4, 8;  #  8
                              8, 5, 1, 7, 6, 2;  #  9
                              5, 1, 4, 6, 2, 3;  # 10
                              1, 4, 8, 2, 3, 7;  # 11
                              4, 8, 5, 3, 7, 6]); # 12
    
    angle_data.iso8 = zeros(rows(mesh.elements.iso8), rows(edges), 2);

    X = mesh.nodes(:, 1:3).';
    
    for j=1:rows(edges)      
      for k=1:2
        X1 = X(:, mesh.elements.iso8(:, edges(j, (k - 1) * 3 + 2)));
        X2 = X(:, mesh.elements.iso8(:, edges(j, (k - 1) * 3 + 3)));
        X3 = X(:, mesh.elements.iso8(:, edges(j, (k - 1) * 3 + 1)));

        e1 = X2 - X1;
        e2 = X3 - X1;
        e3 = cross(e1, e2);
        e2 = cross(e3, e1);

        e1 /= diag(norm(e1, "cols"));
        e2 /= diag(norm(e2, "cols"));

        e21_R1 = zeros(1, rows(mesh.elements.iso8));
        e22_R1 = zeros(1, rows(mesh.elements.iso8));

        dX = X3 - X1;
        
        for i=1:rows(e1)
          e21_R1 += e1(i, :) .* dX(i, :);
          e22_R1 += e2(i, :) .* dX(i, :);
        endfor
        
        Phi = atan2(e22_R1, e21_R1);
        angle_data.iso8(:, j, k) = Phi;
      endfor
    endfor
  endif
endfunction

