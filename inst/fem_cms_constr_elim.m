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
## @deftypefn {Function File} [@var{Tred}, @var{Kred}, @var{Mred}, @var{Rred}] = fem_cms_constr_elim(@var{mesh}, @var{dof_map}, @var{mat_ass})
##
## Eliminate algebraic constraint equations from the stiffness matrix and return matrices for a reduced set of equations
## @seealso{fem_cms_create}
## @end deftypefn

function [Tred, Kred, Mred, Rred] = fem_cms_constr_elim(mesh, dof_map, mat_ass)
  if (nargin ~= 3 || nargout > 4)
    print_usage();
  endif

  diagK = abs(diag(mat_ass.K));
  max_diagK = full(max(diagK));
  min_diagK = full(min(diagK));

  ## According to Code_Aster r3.03.08
  beta = 0.5 * (min_diagK + max_diagK);

  C = mat_ass.K(dof_map.idx_lambda, dof_map.idx_node) / beta;

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

