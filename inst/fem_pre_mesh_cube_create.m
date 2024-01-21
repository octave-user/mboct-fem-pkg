## Copyright (C) 2011(-2021) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{mesh}, @var{load_case}] = fem_pre_mesh_cube_create(@var{geometry}, @var{mesh_size}, @var{material}, @var{f})
## Create a very simple cube mesh of a cantilever beam clamped at x=0 and loaded at x=l
##
## @var{geometry}.l @dots{} length of the cube
##
## @var{geometry}.w @dots{} width of the cube
##
## @var{geometry}.h @dots{} height of the cube
##
## @var{mesh_size}.num_elem_l @dots{} number of elements in l-direction
##
## @var{mesh_size}.num_elem_w @dots{} number of elements in w-direction
##
## @var{mesh_size}.num_elem_h @dots{} number of elements in h-direction
##
## @var{material} @dots{} material data for the complete mesh
##
## @var{f} @dots{} load vector to be applied at the free end of the cantilever beam
##
## @var{mesh} @dots{} Finite element mesh data structure
##
## @var{load_case} @dots{} load case data structure
## @seealso{fem_post_sol_plot, fem_tests}
## @end deftypefn

function [mesh, load_case] = fem_pre_mesh_cube_create(geometry, mesh_size, material, f)
  if (nargin ~= 4 || nargout > 2)
    print_usage();
  endif
  
  N_l = mesh_size.num_elem_l + 1;
  N_w = mesh_size.num_elem_w + 1;
  N_h = mesh_size.num_elem_h + 1;
  num_nodes = N_l * N_w * N_h;
  num_elements = (N_l - 1) * (N_w - 1) * (N_h - 1);
  mesh.nodes = zeros(num_nodes, 6);
  load_case.locked_dof = false(num_nodes, 6);
  load_case.loads = zeros(N_w * N_h, 3);
  load_case.loaded_nodes = zeros(N_w * N_h, 1, "int32");
  mesh.elements.iso8 = zeros(num_elements, 8, "int32");
  node_idx = 0;
  load_idx = 0;

  for x_node_idx=1:N_l
    x = (x_node_idx - 1) / (N_l - 1) * geometry.l;
    for y_node_idx = 1:N_w
      y = (y_node_idx - 1) / (N_w - 1) * geometry.w;
      for z_node_idx = 1:N_h
        z = (z_node_idx - 1) / (N_h - 1) * geometry.h;
        ++node_idx;
        mesh.nodes(node_idx, 1:3) = [x, y, z];
        load_case.locked_dof(node_idx,1:3) = (x == 0);
        if (x == geometry.l)
          ++load_idx;
          load_case.loads(load_idx,1:3) = f / (N_w * N_h);
          load_case.loaded_nodes(load_idx) = node_idx;
        endif
      endfor
    endfor
  endfor

  element_idx = 0;

  for x_elem_idx=1:N_l-1
    for y_elem_idx = 1:N_w-1
      for z_elem_idx = 1:N_h-1
        ++element_idx;
        mesh.elements.iso8(element_idx,1) = (x_elem_idx    ) * N_w * N_h + (y_elem_idx    ) * N_h + z_elem_idx + 1;
        mesh.elements.iso8(element_idx,2) = (x_elem_idx - 1) * N_w * N_h + (y_elem_idx    ) * N_h + z_elem_idx + 1;
        mesh.elements.iso8(element_idx,3) = (x_elem_idx - 1) * N_w * N_h + (y_elem_idx - 1) * N_h + z_elem_idx + 1;
        mesh.elements.iso8(element_idx,4) = (x_elem_idx    ) * N_w * N_h + (y_elem_idx - 1) * N_h + z_elem_idx + 1;

        mesh.elements.iso8(element_idx,5) = (x_elem_idx    ) * N_w * N_h + (y_elem_idx    ) * N_h + z_elem_idx;
        mesh.elements.iso8(element_idx,6) = (x_elem_idx - 1) * N_w * N_h + (y_elem_idx    ) * N_h + z_elem_idx;
        mesh.elements.iso8(element_idx,7) = (x_elem_idx - 1) * N_w * N_h + (y_elem_idx - 1) * N_h + z_elem_idx;
        mesh.elements.iso8(element_idx,8) = (x_elem_idx    ) * N_w * N_h + (y_elem_idx - 1) * N_h + z_elem_idx;
      endfor
    endfor
  endfor

  mesh.material_data = material;
  mesh.materials.iso8 = ones(rows(mesh.elements.iso8), 1, "int32");
endfunction

