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
## @deftypefn {Function File} [@var{lambda}, @var{f}, @var{D}, @var{Phi}, @var{def}] = fem_sol_eigsd_post_proc(@var{mesh}, @var{dof_map}, @var{lambda}, @var{Phi})
## Post-processing support for fem_sol_modal_damped and fem_sol_modal_fsi
##
## @var{mesh} @dots{} Finite element mesh
##
## @var{dof_map} @dots{} Degree of freedom mapping
##
## @var{lambda} @dots{} Eigenvalues
##
## @var{Phi} @dots{} Mode shapes
##
## @var{f} @dots{} Natural frequencies
##
## @var{D} @dots{} Damping ratio
##
## @var{def} @dots{} Nodal deformation
##
## @end deftypefn

function [lambda, f, D, Phi, def] = fem_sol_eigsd_post_proc(mesh, dof_map, lambda, Phi)
  [~, idx] = sortrows([imag(lambda)(:), real(lambda)(:)], [1, 2]);

  lambda = lambda(idx);

  if (nargout >= 2)
    f = imag(lambda) / (2 * pi);
  endif

  if (nargout >= 3)
    D = -real(lambda) ./ abs(lambda);
  endif

  if (nargin >= 4 && nargout >= 4)
    Phi = Phi(:, idx);
  endif

  if (nargout >= 5 && nargin >= 4)
    def = fem_post_def_nodal(mesh, dof_map, Phi * diag(1 ./ norm(Phi, "cols")));
  endif
endfunction
