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
## @deftypefn {Function File} [@var{alpha}, @var{beta}] = fem_pre_mat_rayleigh_damping(@var{D}, @var{f})
## Compute Rayleigh damping coefficients
##
## @var{D} @dots{} modal damping ratio
##
## @var{f} @dots{} frequency
##
## @var{alpha} @dots{} mass damping coefficient
##
## @var{beta} @dots{} stiffness damping coefficient
##
## @end deftypefn

function [alpha, beta] = fem_pre_mat_rayleigh_damping(D, f)
  if (nargin ~= 2 && nargout ~= 2)
    print_usage();
  endif

  if (numel(D) ~= 2 || numel(f) ~= 2)
    print_usage();
  endif
  
  A = [1 / (4 * pi * f(1)), pi * f(1);
       1 / (4 * pi * f(2)), pi * f(2)];
  
  x = A \ D(:);
  
  alpha = x(1);
  beta = x(2);

  if (alpha < 0)
    warning("negative value for alpha detected");
  endif

  if (beta < 0)
    warning("negative value for beta detected");
  endif
endfunction

%!test
%! D = [1e-2;  0.1e-2];
%! f = [10; 1000];
%! omega = 2 * pi * f;
%! [alpha, beta] = fem_pre_mat_rayleigh_damping(D, f);
%! Dref = 1/2* (alpha ./ omega + beta * omega);
%! assert_simple(D, Dref, eps * norm(Dref));
