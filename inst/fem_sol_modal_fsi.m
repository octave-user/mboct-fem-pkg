## Copyright (C) 2023(-2023) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{sol}, @var{Phi}] = fem_sol_modal_fsi(@var{mesh}, @var{dof_map}, @var{mat_ass}, @var{N}, @var{varargin})
## Compute damped natural frequencies and mode shapes.
##
## @var{mesh} @dots{} Finite element mesh data structure
##
## @var{dof_map} @dots{} Degree of freedom mapping
##
## @var{mat_ass}.Kfs_re @dots{} Nodal stiffness matrix real part
##
## @var{mat_ass}.Kfs_im @dots{} Nodal stiffness matrix imaginary part
##
## @var{mat_ass}.Dfs_re @dots{} Nodal damping matrix real part
##
## @var{mat_ass}.Dfs_im @dots{} Nodal damping matrix imaginary part
##
## @var{mat_ass}.Mfs_re @dots{} Nodal mass matrix real part
##
## @var{mat_ass}.Mfs_im @dots{} Nodal mass matrix imaginary part
##
## @var{N} @dots{} Number of modes to compute
##
## @var{sol}.def @dots{} nodal solution for mode shapes
##
## @var{sol}.f @dots{} undamped natural frequencies
##
## @var{Phi} @dots{} mode shapes
##
## @seealso{fem_sol_eigsd}
## @end deftypefn

function [sol, Phi] = fem_sol_modal_fsi(mesh, dof_map, mat_ass, N, varargin)
  if (nargin < 4 || nargout > 3)
    print_usage();
  endif

  switch (dof_map.domain)
    case FEM_DO_FLUID_STRUCT
    otherwise
      error("invalid value for dof_map.domain");
  endswitch

  if (isfield(mat_ass, "Mfs_im"))
    M = complex(mat_ass.Mfs_re, mat_ass.Mfs_im);
  else
    M = mat_ass.Mfs_re;
  endif

  if (isfield(mat_ass, "Dfs_im"))
    D = complex(mat_ass.Dfs_re, mat_ass.Dfs_im);
  else
    D = mat_ass.Dfs_re;
  endif
  
  if (isfield(mat_ass, "Kfs_im"))
    K = complex(mat_ass.Kfs_re, mat_ass.Kfs_im);
  else
    K = mat_ass.Kfs_re;
  endif
  
  [Phi, lambda] = fem_sol_eigsd(K, D, M, N, varargin{:});

  [sol.lambda, sol.f, sol.D, Phi, sol.def] = fem_sol_eigsd_post_proc(mesh, dof_map, lambda, Phi);

  sol.p = zeros(rows(mesh.nodes), columns(Phi));
  idx_act_press_dof = dof_map.ndof(:, 7);
  idx_act_press_node = find(idx_act_press_dof > 0);
  
  for i=1:columns(Phi)
    sol.p(idx_act_press_node, i) = Phi(idx_act_press_dof(idx_act_press_node), i);
  endfor
endfunction
