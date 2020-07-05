## Copyright (C) 2019(-2020) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} @var{load_case} = fem_pre_load_case_create_empty(@var{num_load_cases})
## Allocate a struct array of size <@var{num_load_cases}> of empty load cases
## @end deftypefn

function load_case = fem_pre_load_case_create_empty(num_load_cases)
  if (nargin ~= 1 || nargout > 1)
    print_usage();
  endif
  
  empty_data = cell(1, num_load_cases);
  
  load_case = struct("pressure", empty_data, ...
                     "locked_dof", empty_data, ...
                     "joints", empty_data, ...
                     "loads", empty_data, ...
                     "loaded_nodes", empty_data);
endfunction

%!test ##demo
%! close all;
%! a = 70e-3;
%! b = 20e-3;
%! c = 10e-3;
%! d = 140e-3;
%! scale_def = 20e-3;
%! Fx = 100;
%! Fy = 150;
%! Fz = -230;
%! X = [ a,  0.5 * b,  0.5 * c;
%!       0,  0.5 * b,  0.5 * c;
%!       0, -0.5 * b,  0.5 * c;
%!       a, -0.5 * b,  0.5 * c;
%!       a,  0.5 * b, -0.5 * c;
%!       0,  0.5 * b, -0.5 * c;
%!       0, -0.5 * b, -0.5 * c;
%!       a, -0.5 * b, -0.5 * c];
%! mesh.nodes = [X, zeros(rows(X), 3)];
%! mesh.elements.iso8 = int32(1:8);
%! mesh.materials.iso8 = int32(1);
%! mesh.material_data.E = 210000e6;
%! mesh.material_data.nu = 0.3;
%! mesh.material_data.rho = 7850;
%! mesh.material_data.C = fem_pre_mat_isotropic(mesh.material_data.E, mesh.material_data.nu);
%! load_case = fem_pre_load_case_create_empty(3);
%! load_case(1).locked_dof = false(rows(mesh.nodes), 6);
%! load_case(1).locked_dof([2; 3; 6; 7], 1:3) = true;
%! load_case(1).loaded_nodes = int32([1; 4; 5; 8]);
%! load_case(1).loads = Fx * [0.25, 0, 0, 0, 0, 0;
%!                            0.25, 0, 0, 0, 0, 0;
%!                            0.25, 0, 0, 0, 0, 0;
%!                            0.25, 0, 0, 0, 0, 0];
%! load_case(2).loaded_nodes = int32([1; 4; 5; 8]);
%! load_case(2).loads = Fy * [0, 0.25, 0, 0, 0, 0;
%!                            0, 0.25, 0, 0, 0, 0;
%!                            0, 0.25, 0, 0, 0, 0;
%!                            0, 0.25, 0, 0, 0, 0];
%! load_case(3).loaded_nodes = int32([1; 4; 5; 8]);
%! load_case(3).loads = Fz * [0, 0, 0.25, 0, 0, 0;
%!                            0, 0, 0.25, 0, 0, 0;
%!                            0, 0, 0.25, 0, 0, 0;
%!                            0, 0, 0.25, 0, 0, 0];
%! dof_map = fem_ass_dof_map(mesh, load_case(1));  
%! [mat_ass.K, ...
%!  mat_ass.R] = fem_ass_matrix(mesh, ...
%!                              dof_map, ...
%!                              [FEM_MAT_STIFFNESS, ...
%!                               FEM_VEC_LOAD_CONSISTENT], ...
%!                              load_case);
%! sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%! for i=1:size(sol_stat.def, 3)
%!   figure("visible", "off");
%!   fem_post_sol_plot(mesh, sol_stat, scale_def / max(max(abs(sol_stat.def(:, :, i)))), i);
%!   xlabel("x [m]");
%!   ylabel("y [m]");
%!   zlabel("z [m]");
%!   grid on;
%!   grid minor on;
%!   title(sprintf("load case %d", i));
%!   view(30, 30);
%! endfor
%! figure_list();
