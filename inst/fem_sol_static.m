## Copyright (C) 2011(-2019) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{sol}] = fem_sol_static(@var{mesh}, @var{dof_map}, @var{mat_ass})
## @deftypefnx {} [@dots{}, @var{U}] = fem_sol_static(@dots{}, @var{options})
##
## Solve a linear static problem and return a structure with nodal displacements
##
## @var{mesh} @dots{} finite element mesh data structure
##
## @var{dof_map} @dots{} degree of freedom mapping
##
## @var{mat_ass} @dots{} data structure containing the global stiffness matrix @var{K} and load vectors @var{R}
##
## @var{options} @dots{} options passed to the linear solver fem_sol_linsolve
##
## @end deftypefn

function [sol, U] = fem_sol_static(mesh, dof_map, mat_ass, options)
  if (nargin < 3 || nargin > 4 || nargout > 2)
    print_usage();
  endif

  if (nargin < 4)
    options = struct();
  endif

  U = fem_sol_linsolve(mat_ass.K, mat_ass.R, options);
  sol.def = fem_post_def_nodal(mesh, dof_map, U);
endfunction

%!demo
%! close all;
%! material.E = 210000e6;
%! material.nu = 0.3;
%! material.rho = 7850;
%! Fz = 15000;
%! h = 10e-3 / 2;
%! geometry.l = 1000e-3;
%! geometry.w = 10e-3;
%! geometry.h = 50e-3;
%! mesh_size.num_elem_l = ceil(geometry.l / h);
%! mesh_size.num_elem_w = ceil(geometry.w / h);
%! mesh_size.num_elem_h = ceil(geometry.h / h);
%! f = [0; 0; Fz];
%! [mesh, load_case] = fem_pre_mesh_cube_create(geometry, mesh_size, material, f);
%! [dof_map] = fem_ass_dof_map(mesh, load_case);
%! [mat_ass.K, ...
%!  mat_ass.R] = fem_ass_matrix(mesh, ...
%!                              dof_map, ...
%!                              [FEM_MAT_STIFFNESS, ...
%!                               FEM_VEC_LOAD_CONSISTENT], ...
%!                              load_case);
%! sol = fem_sol_static(mesh, dof_map, mat_ass);
