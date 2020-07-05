## Copyright (C) 2018(-2020) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{Kred}, @var{Tred}] = fem_cms_red_guyan(@var{K}, @var{keep_dof})
## Build a reduced order model by static condensation.
##
## @var{K} @dots{} Nodal stiffness matrix
##
## @var{keep_dof} @dots{} Index vector for master degrees of freedom
##
## @var{Kred} @dots{} Reduced stiffness matrix
##
## @var{Tred} @dots{} Static mode shapes
##
## @end deftypefn

function [Kred, Tred] = fem_cms_red_guyan(K, keep_dof)
  if (nargin ~= 2 || nargout > 2)
    print_usage();
  endif
  
  [K, dof_list, keep_dof, condense_dof] = fem_reorder_equations_of_motion(K, keep_dof);

  KNN = K(condense_dof, condense_dof);
  KNH = K(condense_dof, keep_dof);

  KNNfact = fem_sol_factor(KNN);

  PHI = KNNfact \ KNH;
  
  T = [eye(length(keep_dof));
       -PHI];

  Kred = fem_cms_matrix_trans(T, K, "Upper");
  
  Tred = eye(length(K))(:, dof_list) * T;
endfunction

function [Kr, dof_list, keep_dof, condense_dof] = fem_reorder_equations_of_motion(K, keep_dof)
  N = length(K);
  
  if (any(keep_dof < 1 | keep_dof > N))
    error("keep_dof not in range [%d:%d]", 1, N);
  endif
  
  condense_dof = true(1, N);
  
  condense_dof(keep_dof) = false;
  
  condense_dof = find(condense_dof);
  
  N_keep = length(keep_dof);
  N_condense = length(condense_dof);
  
  dof_list = [keep_dof, condense_dof];
  
  Kr = K(dof_list, dof_list);
  
  keep_dof = 1:N_keep;
  condense_dof = N_keep + (1:N_condense);
endfunction

%!test
%! E = 210000e6;
%! A = 100e-6;
%! l = 50e-3;
%! s = A * E / l;
%! K = [ s,    -s,     0,      0;
%!      -s, 2 * s,    -s,      0;
%!       0,    -s, 2 * s,     -s;
%!       0,     0,    -s,  2 * s];
%!
%! [Kred, Tred] = fem_cms_red_guyan(K, [1, columns(K)]);
%!
%! Kred2 = [ s / 3,    -s / 3;
%!          -s / 3,     s + s / 3];
%! assert(Kred, Kred2, eps * norm(Kred2));

%!test ##demo
%! E = 210000e6;
%! A = 100e-6;
%! l = 50e-3;
%! s = A * E / l;
%! K = [ s,    -s,     0,     0;
%!      -s, 2 * s,    -s,     0;
%!       0,    -s, 2 * s,    -s;
%!       0,     0,    -s,     s];
%!
%! [Kred, Tred] = fem_cms_red_guyan(K, [1, columns(K)]);
%!
%! s2 = A * E / (3 * l);
%! Kred2 = [ s2,    -s2;
%!          -s2,     s2];
%! assert(Kred, Kred2, eps * norm(Kred2));
