## Copyright (C) 2011(-2020) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} @var{def} = fem_post_def_nodal(@var{mesh}, @var{dof_map}, @var{U})
## Map a solution vector or load vector @var{U} back to nodal displacements @var{def}
##
## @var{mesh} @dots{} Finite element mesh data structure
##
## @var{dof_map} @dots{} Mapping for degrees of freedom
##
## @var{U} @dots{} Solution vector or load vector
##
## @end deftypefn

function def = fem_post_def_nodal(mesh, dof_map, U)
  if (nargin ~= 3 || nargout > 1)
    print_usage();
  endif

  def = zeros(rows(dof_map.ndof), columns(dof_map.ndof), columns(U));
  
  for i=1:columns(dof_map.ndof)
    dof_idx = dof_map.ndof(:, i);
    idx_dof_idx_ne_0 = find(dof_idx > 0);
    def(idx_dof_idx_ne_0, i, :) = U(dof_idx(idx_dof_idx_ne_0), :);
  endfor
endfunction

%!demo
%! close all;
%! a = 70e-3;
%! b = 20e-3;
%! c = 10e-3;
%! d = 140e-3;
%! scale_def = 100e-3;
%! tol = eps^0.7;
%! X = [ 0.5 * a,  0.5 * b,  0.5 * c;  #  1
%!             0,  0.5 * b,  0.5 * c;  #  2
%!             0, -0.5 * b,  0.5 * c;  #  3
%!       0.5 * a, -0.5 * b,  0.5 * c;  #  4
%!       0.5 * a,  0.5 * b, -0.5 * c;  #  5
%!             0,  0.5 * b, -0.5 * c;  #  6
%!             0, -0.5 * b, -0.5 * c;  #  7
%!       0.5 * a, -0.5 * b, -0.5 * c,  #  8
%!             a,  0.5 * b,  0.5 * c;  #  9
%!             a, -0.5 * b,  0.5 * c;  # 10
%!             a,  0.5 * b, -0.5 * c;  # 11
%!             a, -0.5 * b, -0.5 * c,  # 12
%!             d,        0,        0]; # 13

%! mesh.nodes = [X, zeros(rows(X), 3)];
%! mesh.elements.iso8 = int32([1:8;
%!                             9, 1, 4, 10, 11, 5, 8, 12]);
%! mesh.materials.iso8 = int32([1; 1]);
%! mesh.elements.rbe3.nodes = int32([13, 9, 10, 11, 12]);
%! mesh.elements.rbe3.weight = ones(1, 4);
%! n1 = [0.5; 0.5; 0];
%! n1 /= norm(n1);
%! n2 = [1; 1; 1];
%! n3 = cross(n1, n2);
%! n2 = cross(n1, n3);
%! n2 /= norm(n2);
%! n3 /= norm(n3);
%! mesh.elements.joints(1).nodes = int32([13]);
%! mesh.elements.joints(1).C = [[1, 0, 0; 
%!                               0, 0, 1],    zeros(2, 3);
%!                               zeros(2, 3), [n2.'; n3.']];
%! load_case.joints(1).U = [0; 0; 0; 0];
%! for i=[2, 3, 6, 7]
%!   mesh.elements.joints(end + 1).nodes = int32(i);
%!   mesh.elements.joints(end).C = [eye(3), zeros(3, 3)];
%!   load_case.joints(end + 1).U = [0; 0; 0];
%! endfor
%! E = 210000e6;
%! nu = 0.3;
%! mesh.material_data.rho = 7850;
%! mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%! load_case.locked_dof = false(rows(mesh.nodes), 6);

%! # load_case.locked_dof([2, 3, 6, 7], 1:6) = true;
%! load_case.loaded_nodes = int32([1; 4; 5; 8]);
%! load_case.loads = repmat([0, 0, 10,  0,   0, 0], 4, 1);
%! dof_map = fem_ass_dof_map(mesh, load_case);  
%! [mat_ass.K, ...
%!  mat_ass.R, ...
%!  mat_ass.mat_info] = fem_ass_matrix(mesh, ...
%!                                     dof_map, ...
%!                                     [FEM_MAT_STIFFNESS, ...
%!                                      FEM_VEC_LOAD_CONSISTENT], ...
%!                                     load_case);
%! [Tred, Kred, ~, Rred] = fem_cms_constr_elim(mesh, dof_map, mat_ass);
%! [sol, U] = fem_sol_static(mesh, dof_map, mat_ass);
%! Ured = zeros(columns(mat_ass.K), columns(mat_ass.R));
%! Ured(dof_map.idx_node, :) = Tred * (Kred \ Rred);
%! sol_red.def = fem_post_def_nodal(mesh, dof_map, Ured);
%! assert(sol_red.def, sol.def, sqrt(eps) * max(max(max(abs(sol.def)))));
%! figure("visible", "off");
%! fem_post_sol_plot(mesh, sol_red, scale_def / max(norm(sol_red.def(:, 1:3), "cols")));
%! xlabel("x [m]");
%! ylabel("y [m]");
%! zlabel("z [m]");
%! title("deformed mesh");
%! figure_list();
%! view(30, 5);
