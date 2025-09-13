## Copyright (C) 2019(-2023) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{lambda_elem}] = fem_fem_post_lambda_elem(@var{dof_map}, @var{mat_ass}, @var{U})
## Extract Lagrange multipliers @var{lambda_elem} from a Finite Element solution vector @var{U}.
## @end deftypefn

function lambda_elem = fem_post_lambda_elem(dof_map, mat_ass, U)
  lambda_elem = struct();

  if (isfield(dof_map, "edof") && isfield(mat_ass, "mat_info"))
    mat_type_beta = int32(-1);

    switch (dof_map.domain)
      case FEM_DO_STRUCTURAL
        mat_type_beta = [FEM_MAT_STIFFNESS, FEM_MAT_STIFFNESS_SYM, FEM_MAT_STIFFNESS_SYM_L];
      case FEM_DO_THERMAL
        mat_type_beta = FEM_MAT_THERMAL_COND;
      case FEM_DO_ACOUSTIC
        mat_type_beta = MAT_DAMPING_ACOUSTICS_RE;
      case FEM_DO_FLUID_STRUCT
        mat_type_beta = FEM_MAT_STIFFNESS_FLUID_STRUCT_RE;
      otherwise
        error("unknown domain %d", dof_map.domain);
    endswitch

    for i=1:numel(mat_type_beta)
      idx_beta = find(mat_ass.mat_info.mat_type == mat_type_beta(i));

      if (~isempty(idx_beta))
        break;
      endif
    endfor

    if (numel(idx_beta) ~= 1)
      error("failed to detect matrix type based on mat_ass.mat_info");
    endif

    beta = mat_ass.mat_info.beta(idx_beta);

    elem_lambda = {"rbe2", "rbe3", "joints"};

    for i=1:numel(elem_lambda)
      if (isfield(dof_map.edof, elem_lambda{i}))
        idx_lambda = getfield(dof_map.edof, elem_lambda{i});
        lambda = zeros(rows(idx_lambda), columns(idx_lambda), columns(U));
        for j=1:columns(U)
          for k=1:columns(idx_lambda)
            lambda(idx_lambda(:, k) > 0, k, j) = beta * U(idx_lambda(idx_lambda(:, k) > 0, k), j);
          endfor
        endfor
        lambda_elem = setfield(lambda_elem, elem_lambda{i}, lambda);
      endif
    endfor
  endif
endfunction
