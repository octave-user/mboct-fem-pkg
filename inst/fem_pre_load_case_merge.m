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
## @deftypefn {Function File} @var{load_case} = fem_pre_load_case_merge(@var{varargin})
## Combine several struct arrays of load cases into a single struct array
## @end deftypefn

function load_case = fem_pre_load_case_merge(varargin)
  load_cases = varargin;

  fn = {};
  
  for i=1:numel(load_cases)
    fn = {fn{:}, fieldnames(load_cases{i}){:}};
  endfor
  
  fn = unique(fn);
  
  for i=1:numel(load_cases)
    for j=1:numel(fn)
      if (~isfield(load_cases{i}, fn{j}))
        eval(sprintf("load_cases{i}(1).%s = [];", fn{j}));
      endif
    endfor
  endfor

  load_case = [load_cases{:}];
endfunction

%!test
%! load_case1.locked_dof = false(10, 6);
%! load_case2.pressure.tria6 = zeros(5, 6, "int32");
%! load_case3(1).joints.U = zeros(3, 1);
%! load_case3(2).joints.U = ones(3, 1);
%! load_case = fem_pre_load_case_merge(load_case1, load_case2, load_case3);
%! assert(numel(load_case), 4);
%! assert(isfield(load_case, "locked_dof"));
%! assert(isfield(load_case, "pressure"));
%! assert(isfield(load_case, "joints"));
%! assert(all(all(load_case(1).locked_dof == load_case1.locked_dof)));
%! assert(all(all(load_case(2).pressure.tria6 == load_case2.pressure.tria6)));
%! assert(all(all(load_case(3).joints.U == load_case3(1).joints.U)));
%! assert(all(all(load_case(4).joints.U == load_case3(2).joints.U)));

%!test
%! load_case = fem_pre_load_case_merge();

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
%! load_case1.locked_dof = false(rows(mesh.nodes), 6);
%! load_case1.locked_dof([2; 3; 6; 7], 1:3) = true;
%! load_case1.loaded_nodes = int32([1; 4; 5; 8]);
%! load_case1.loads = Fx * [0.25, 0, 0, 0, 0, 0;
%!                          0.25, 0, 0, 0, 0, 0;
%!                          0.25, 0, 0, 0, 0, 0;
%!                          0.25, 0, 0, 0, 0, 0];
%! load_case2.loaded_nodes = int32([1; 4; 5; 8]);
%! load_case2.loads = Fy * [0, 0.25, 0, 0, 0, 0;
%!                          0, 0.25, 0, 0, 0, 0;
%!                          0, 0.25, 0, 0, 0, 0;
%!                          0, 0.25, 0, 0, 0, 0];
%! load_case3.loaded_nodes = int32([1; 4; 5; 8]);
%! load_case3.loads = Fz * [0, 0, 0.25, 0, 0, 0;
%!                          0, 0, 0.25, 0, 0, 0;
%!                          0, 0, 0.25, 0, 0, 0;
%!                          0, 0, 0.25, 0, 0, 0];
%! load_case = fem_pre_load_case_merge(load_case1, load_case2, load_case3);
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
