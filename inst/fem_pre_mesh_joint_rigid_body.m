## Copyright (C) 2019(-2026) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{joint}] = fem_pre_mesh_joint_rigid_body(@var{mesh}, @var{node_id}, @var{R}, @var{dof_idx})
## Build a joint which imposes a rigid body constraint between node_id(1) and node_id(1)
## @end deftypefn

function joint = fem_pre_mesh_joint_rigid_body(mesh, node_id, R, dof_idx)
  if (nargin < 2 || nargout > 1)
    print_usage();
  endif
  
  if (nargin < 3)
    R = eye(3);
  endif

  if (nargin < 4)
    dof_idx = 1:6;
  endif
  
  joint.nodes = node_id;
  l21 = mesh.nodes(node_id(2), 1:3) - mesh.nodes(node_id(1), 1:3);
  joint.C = [        R.', -R.' * skew(l21),        -R.', zeros(3, 3);
             zeros(3, 3),              R.', zeros(3, 3),        -R.'];
  joint.C = joint.C(dof_idx, :);
endfunction

%!test
%! Fx = 1;
%! Fy = 2;
%! Fz = 3;
%! Mx = 4;
%! My = 5;
%! Mz = 6;
%! lx = 20e-3;
%! ly = 10e-3;
%! lz = 5e-3;
%! Ux = 1e-3;
%! Uy = 2.5e-3;
%! Uz = -7e-3;
%! Phix = 10 * pi / 180;
%! Phiy = 5 * pi / 180;
%! Phiz = -7 * pi / 180;
%! R = [0, 1, 0;
%!      0, 0, 1;
%!      1, 0, 0];
%! mesh.nodes = [ 0,  0,  0, 0, 0, 0;
%!               lx, ly, lz, 0, 0, 0];
%! mesh.material_data = struct()([]);
%! mesh.materials = struct();
%! mesh.elements.joints(1) = fem_pre_mesh_joint_rigid_body(mesh, int32([1, 2]), R);
%! mesh.elements.joints(2).nodes = int32(1);
%! mesh.elements.joints(2).C = eye(6);
%! load_case_dof.locked_dof = false(size(mesh.nodes));
%! load_case_dof.dof_in_use = true(2, 6);
%! load_case.joints(1).U = zeros(6, 1);
%! load_case.joints(2).U = [Ux; Uy; Uz; Phix; Phiy; Phiz];
%! load_case.loads = [Fx, Fy, Fz, Mx, My, Mz];
%! load_case.loaded_nodes = int32(2);
%! dof_map = fem_ass_dof_map(mesh, load_case_dof);
%! [mat_ass.K, ...
%!  mat_ass.R, ...
%!  mat_ass.mat_info, ...
%!  mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_MAT_STIFFNESS, ...
%!                                       FEM_VEC_LOAD_CONSISTENT], ...
%!                                      load_case);
%! sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%! tol = eps;
%! assert_simple(sol_stat.def(1, 1:6), [Ux, Uy, Uz, Phix, Phiy, Phiz], tol * max(abs([Phix, Phiy, Phiz])));
%! assert_simple(sol_stat.def(2, 1:6), [[Ux, Uy, Uz] - (skew([lx;ly;lz])*[Phix;Phiy;Phiz]).', Phix, Phiy, Phiz], tol * max(abs([Phix, Phiy, Phiz])));
%! assert_simple(sol_stat.lambda.joints(1, :), [(-R.' * [Fx; Fy; Fz]).', (-R.' * [Mx; My; Mz]).'], tol  * max(abs([Fx, Fy, Fz, Mx, My, Mz])));
%! assert_simple(sol_stat.lambda.joints(2, :), [Fx,Fy,Fz,([Mx;My;Mz]+skew([lx;ly;lz])*[Fx;Fy;Fz]).'], tol * max(abs([Fx, Fy, Fz, Mx, My, Mz])));
