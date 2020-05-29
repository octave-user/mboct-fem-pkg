## Copyright (C) 2020(-2020) Reinhard <octave-user@a1.net>
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
## @deftypefn {} @var{Asym} = fem_mat_sym(@var{A})
## Create a symmetric matrix which stores only the lower
## or upper triangular part of the matrix
## in order to reduce memory consumption.
## Only a subset of operations is available for this matrix type. 
## @end deftypefn

function Asym = fem_mat_sym(A)
  if (~(nargin == 1 && ismatrix(A) && issquare(A)))
    print_usage();
  endif
  
  Asym.A = A;
  Asym = class(Asym, "fem_mat_sym");
endfunction

%!test
%! A = rand(10, 10);
%! A += A.';
%! [r, c, d] = find(A);
%! idx = find(r >= c);
%! Asym = fem_mat_sym(sparse(r(idx), c(idx), d(idx), rows(A), columns(A)));
%! b = rand(columns(Asym), 10);
%! x1 = A * b;
%! x2 = Asym * b;
%! assert(x1, x2, eps^0.8 * norm(b));
%! x1 = Asym(2:4, 2:4) * b(2:4, :);
%! x2 = A(2:4, 2:4) * b(2:4, :);
%! assert(x1, x2, eps * norm(b(2:4, :)));
