## Copyright (C) 2025(-2025) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{X}, @var{XP}, @var{XPP}] = fem_post_cms_sol_nodal(@var{sol_dyn}, @var{cms_data}, @var{node_idx})
## Evaluate the time history of position, velocity and acceleration of an arbitrary node in a reduced order component mode synthesis model.
##
## @var{sol_dyn} @dots{} Modal solution imported from MBDyn
##
## @var{cms_data} @dots{} Finite element data including mesh, options, and mode shapes used to create modal element data for MBDyn
##
## @seealso{fem_cms_export, fem_post_cms_sol_import, fem_post_cms_expand, fem_post_cms_sol_merge, fem_post_mesh_merge, fem_post_mesh_export}
## @end deftypefn

function [X, XP, XPP] = fem_post_cms_sol_nodal(sol_dyn, cms_data, node_idx)
  if (nargin ~= 3 || nargout > 3)
    print_usage();
  endif

  idx_body = int32(-1);

  for i=1:numel(sol_dyn.bodies)
    switch (sol_dyn.bodies(i).element.name)
      case cms_data.cms_opt.element.name
        idx_body = i;
        break;
    endswitch
  endfor

  if (idx_body == -1)
    error("body %s not found in argument sol_dyn", cms_data.cms_opt.element.name);
  endif

  dof_idx = cms_data.dof_map.ndof(node_idx, 1:3);
  idx_active_dof = find(dof_idx > 0);
  dof_idx = dof_idx(idx_active_dof);
  idx_Tred = zeros(numel(dof_idx), 1, "int32");

  for i=1:numel(dof_idx)
    idx_Tred(i) = find(dof_idx(i) == cms_data.dof_map.idx_node);
  endfor

  Tred = zeros(3, columns(cms_data.mat_ass.Tred));

  Tred(idx_active_dof, :) = cms_data.mat_ass.Tred(idx_Tred, :);

  X0 = sol_dyn.bodies(idx_body).X;
  XP0 = sol_dyn.bodies(idx_body).XP;
  XPP0 = sol_dyn.bodies(idx_body).XPP;
  R0 = sol_dyn.bodies(idx_body).R;
  W0 = sol_dyn.bodies(idx_body).W;
  WP0 = sol_dyn.bodies(idx_body).WP;
  f = cms_data.mesh.nodes(node_idx, 1:3).' - cms_data.mesh.nodes(cms_data.cms_opt.nodes.modal.number, 1:3).';
  q = sol_dyn.bodies(idx_body).q;
  qdot = sol_dyn.bodies(idx_body).qdot;
  qddot = sol_dyn.bodies(idx_body).qddot;

  X = XP = XPP = nan(3, columns(X0));

  for i=1:columns(X)
    X(:, i) = X0(:, i) + R0(:, :, i) * (f + Tred * q(:, i));
    XP(:, i) = XP0(:, i) + cross(W0(:, i), R0(:, :, i) * (f + Tred * q(:, i))) + R0(:, :, i) * Tred * qdot(:, i);
    if (~isempty(XPP0))
      XPP(:, i) = XPP0(:, i) + (skew(WP0(:, i)) + skew(W0(:, i)) * skew(W0(:, i))) * R0(:, :, i) * (f + Tred * q(:, i)) ...
                  + cross(W0(:, i), R0(:, :, i) * Tred * qdot(:, i)) + R0(:, :, i) * Tred * qddot(:, i);
    endif
  endfor
endfunction

%!test
%! pkg load mbdyn_util_oct;
%! a = [1; 2; 3];
%! v0 = [0.1; 0.5; 0.6];
%! s0 = [0.05; 0.06; 0.07];
%! Phi0 = 0;
%! omega0 = 0;
%! alpha = 0;
%! cms_data.cms_opt.element.name = "rotor1";
%! cms_data.cms_opt.nodes.modal.number = int32(2);
%! cms_data.mat_ass.Tred = eye(6);
%! cms_data.dof_map.ndof = int32([1, 2, 3, 4, 5, 6]);
%! cms_data.dof_map.idx_node = int32([1, 2, 3, 4, 5, 6]);
%! cms_data.dof_map.totdof = int32(6);
%! cms_data.mesh.nodes = [0.5, 0.6, 0.7, 0, 0, 0;
%!                        0.5, 0.6, 0.7, 0, 0, 0];
%! sol_dyn.t = [0, 0.01];
%! sol_dyn.bodies(1).element.name = "rotor1";
%! sol_dyn.bodies(1).q = zeros(6, 2);
%! sol_dyn.bodies(1).qdot = zeros(6, 2);
%! sol_dyn.bodies(1).qddot = zeros(6, 2);
%! sol_dyn.bodies(1).XPP = repmat(a, 1, 2);
%! sol_dyn.bodies(1).XP = v0 + a * sol_dyn.t;
%! sol_dyn.bodies(1).X = s0 + v0 * sol_dyn.t + 0.5 * a * sol_dyn.t.^2;
%! sol_dyn.bodies(1).WP = [zeros(2, 2);
%!                         repmat(alpha, 1, 2)];
%! sol_dyn.bodies(1).W = [zeros(2, 2);
%!                        omega0 + alpha * sol_dyn.t];
%! sol_dyn.bodies(1).R = euler123_to_rotation_matrix([zeros(2, 2);
%!                                                    Phi0 + omega0 * sol_dyn.t + 0.5 * alpha * sol_dyn.t.^2]);
%! node_idx = int32(1);
%! [X, XP, XPP] = fem_post_cms_sol_nodal(sol_dyn, cms_data, node_idx);
%! assert_simple(X, s0 + v0 * sol_dyn.t + 0.5 * a * sol_dyn.t.^2, eps * norm(s0));
%! assert_simple(XP, v0 + a * sol_dyn.t, eps * norm(v0));
%! assert_simple(XPP, a * ones(1, 2), eps * norm(a));

%!test
%! pkg load mbdyn_util_oct;
%! a = [0; 0; 0];
%! v0 = [0; 0; 0];
%! s0 = [0; 0; 0];
%! Phi0 = 0;
%! omega0 = pi / 0.01;
%! alpha = 0;
%! cms_data.cms_opt.element.name = "rotor1";
%! cms_data.cms_opt.nodes.modal.number = int32(2);
%! cms_data.mat_ass.Tred = eye(6);
%! cms_data.dof_map.ndof = int32([1, 2, 3, 4, 5, 6]);
%! cms_data.dof_map.idx_node = int32([1, 2, 3, 4, 5, 6]);
%! cms_data.dof_map.totdof = int32(6);
%! cms_data.mesh.nodes = [1.5, 2.6, 3.7, 0, 0, 0;
%!                        0.5, 0.6, 0.7, 0, 0, 0];
%! sol_dyn.t = [0, 0.01];
%! sol_dyn.bodies(1).element.name = "rotor1";
%! sol_dyn.bodies(1).q = zeros(6, 2);
%! sol_dyn.bodies(1).qdot = zeros(6, 2);
%! sol_dyn.bodies(1).qddot = zeros(6, 2);
%! sol_dyn.bodies(1).XPP = repmat(a, 1, 2);
%! sol_dyn.bodies(1).XP = v0 + a * sol_dyn.t;
%! sol_dyn.bodies(1).X = s0 + v0 * sol_dyn.t + 0.5 * a * sol_dyn.t.^2;
%! sol_dyn.bodies(1).WP = [zeros(2, 2);
%!                         repmat(alpha, 1, 2)];
%! sol_dyn.bodies(1).W = [zeros(2, 2);
%!                        omega0 + alpha * sol_dyn.t];
%! sol_dyn.bodies(1).R = euler123_to_rotation_matrix([zeros(2, 2);
%!                                                    Phi0 + omega0 * sol_dyn.t + 0.5 * alpha * sol_dyn.t.^2]);
%! node_idx = int32(1);
%! [X, XP, XPP] = fem_post_cms_sol_nodal(sol_dyn, cms_data, node_idx);
%! l = (cms_data.mesh.nodes(1, 1:3) - cms_data.mesh.nodes(2, 1:3)).';
%! Xref = Vref = Aref = zeros(3, 2);
%! for i=1:numel(sol_dyn.t)
%!   Xref(:, i) = sol_dyn.bodies(1).R(:, :, i) * l;
%!   Vref(:, i) = skew([0; 0; omega0]) * sol_dyn.bodies(1).R(:, :, i) * l;
%!   Aref(:, i) = skew([0; 0; omega0]) * skew([0; 0; omega0]) * sol_dyn.bodies(1).R(:, :, i) * l;
%! endfor
%! assert(X, Xref, eps * norm(Xref));
%! assert(XP, Vref, eps * norm(Vref));
%! assert(XPP, Aref, eps * norm(Aref));
