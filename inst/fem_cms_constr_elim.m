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
## @deftypefn {Function File} [@var{Tred}, @var{Kred}, @var{Mred}, @var{Rred}] = fem_cms_constr_elim(@var{mesh}, @var{dof_map}, @var{mat_ass})
##
## Eliminate algebraic constraint equations from the stiffness matrix and return matrices for a reduced set of equations
## @seealso{fem_cms_create}
## @end deftypefn

function [Tred, Kred, Mred, Rred] = fem_cms_constr_elim(mesh, dof_map, mat_ass)
  if (nargin ~= 3 || nargout > 4)
    print_usage();
  endif
  
  C = mat_ass.K(dof_map.idx_lambda, dof_map.idx_node) / mat_ass.mat_info.beta;

  cond = sum(abs(C), 1);
  idx_constr = find(cond);
  idx_unconstr = find(~cond);
  
  Tred1 = fem_cms_constr_mat(C(:, idx_constr));

  [ridx1, cidx1, data1] = find(Tred1);
  
  ridx = [idx_constr(ridx1.'), idx_unconstr];
  cidx = [cidx1.', columns(Tred1) + (1:numel(idx_unconstr))];
  data = [data1; ones(numel(idx_unconstr), 1)];
  
  Tred = sparse(ridx, cidx, data, rows(Tred1) + numel(idx_unconstr), columns(Tred1) + numel(idx_unconstr));
  
  if (nargout >= 2)
    Kred = fem_transf_mat(mat_ass.K, dof_map.idx_node(idx_constr), dof_map.idx_node(idx_unconstr), Tred1);
  endif

  if (nargout >= 3)
    if (isfield(mat_ass, "M"))
      Mred = fem_transf_mat(mat_ass.M, dof_map.idx_node(idx_constr), dof_map.idx_node(idx_unconstr), Tred1);
    else
      Mred = [];
    endif
  endif

  if (nargout >= 4)
    Rred = Tred.' * mat_ass.R(dof_map.idx_node, :);
  endif
endfunction

function Ared = fem_transf_mat(A, idx1, idx2, T1)
  A11 = A(idx1, idx1);
  A12 = A(idx1, idx2);
  A22 = A(idx2, idx2);

  Ared11 = fem_cms_matrix_trans(T1, A11, "Full");
  Ared12 = T1.' * A12;
  
  Ared = [Ared11, Ared12;
          Ared12.', A22];
endfunction

%!test ##demo
%! ## Build a simple mesh made of a single hexahedron
%! X = [ 1,  1,  1;
%!      -1,  1,  1;
%!      -1, -1,  1;
%!       1, -1,  1;
%!       1,  1, -1;
%!      -1,  1, -1;
%!      -1, -1, -1;
%!       1, -1, -1];
%! mesh.nodes = [X, zeros(rows(X), 3)];
%! mesh.elements.iso8 = int32(1:8);
%! mesh.materials.iso8 = int32(1);
%! mesh.material_data.E = 210000e6;
%! mesh.material_data.nu = 0.3;
%! mesh.material_data.rho = 7850;
%! mesh.material_data.C = fem_pre_mat_isotropic(mesh.material_data.E, mesh.material_data.nu);
%! load_case.locked_dof = false(rows(mesh.nodes), 6);
%! ## Clamp all bottom nodes to the ground
%! load_case.locked_dof(5:8, :) = true;
%! ## Apply several loads at the top nodes
%! load_case.loaded_nodes = int32(1:4).';
%! load_case.loads = [1e10, 0, 0;
%!                    0, 2e10, 0;
%!                    0, 0, -3e10;
%!                    0, 1e10, 0];
%! ## Impose a constraint to the top nodes in a way,
%! ## that the displacement for all four nodes is identical in x, y and z-direction
%! mesh.elements.joints.nodes = int32(1:4);
%! mesh.elements.joints.C = [  eye(3, 6),  -eye(3, 6), zeros(3, 6), zeros(3, 6);
%!                           zeros(3, 6),   eye(3, 6),  -eye(3, 6), zeros(3, 6);
%!                           zeros(3, 6), zeros(3, 6),   eye(3, 6),  -eye(3, 6)];
%! ## Build the degree of freedom table   
%! [dof_map] = fem_ass_dof_map(mesh, load_case);
%! ## Assemble global matrices
%! [mat_ass.K, ...
%!  mat_ass.M, ...
%!  mat_ass.R, ...
%!  mat_ass.mat_info, ...
%!  mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_MAT_STIFFNESS, ...
%!                                       FEM_MAT_MASS, ...
%!                                       FEM_VEC_LOAD_CONSISTENT], ...
%!                                      load_case);
%! ## Eliminate redundant equations
%! [Tred, Kred, Mred, Rred] = fem_cms_constr_elim(mesh, dof_map, mat_ass);
%! ## Solve the reduced set of equations
%! Ured = Tred * fem_sol_linsolve(Kred, Rred);
%! sol_red.def = zeros(size(mesh.nodes));
%! ## Compute nodal displacements
%! for i=1:columns(dof_map.ndof)
%!   idx = find(dof_map.ndof(:, i) > 0);
%!   sol_red.def(idx, i) = Ured(dof_map.ndof(idx, i));
%! endfor
%! tol = eps^0.9;
%! ## Solve the full set of equations
%! [sol_ref] = fem_sol_static(mesh, dof_map, mat_ass);
%! ## Check if we got the same result with reduced set of equations and full set of equations
%! assert(sol_red.def, sol_ref.def, tol * norm(sol_ref.def));
%! ## Check if constraint we imposed for the top nodes is satisfied
%! assert(sol_red.def(1:4, :), repmat(sol_red.def(1, :), 4, 1), tol * norm(sol_red.def(1, :)));

