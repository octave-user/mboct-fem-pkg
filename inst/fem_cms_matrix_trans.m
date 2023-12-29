## Copyright (C) 2011(-2021) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} @var{TAT} = fem_cms_matrix_trans(@var{T}, @var{A}, @var{mat_type})
##
## Effectively compute @var{TAT} = @var{T}.' * @var{A} * @var{T} but more efficient in time and memory for large symmetric sparse matrices @var{A}
##
## @end deftypefn

function TAT = fem_cms_matrix_trans(T, A, mat_type)
  if (nargin ~= 3 || nargout > 1)
    print_usage();
  endif
  
  switch (mat_type)
    case "Upper"
      TAT = zeros(columns(T), columns(T));
      for j=1:columns(T)
        AT = A * T(:, j);
        for i=1:j
          TAT(i, j) = T(:, i).' * AT;
          if (i ~= j)
            TAT(j, i) = TAT(i, j);
          endif
        endfor
      endfor
    case "Lower"
      TAT = zeros(columns(T), columns(T));
      for j=1:columns(T)
        AT = A * T(:, j);
        for i=j:columns(T)
          TAT(i, j) = T(:, i).' * AT;
          if (i ~= j)
            TAT(j, i) = TAT(i, j);
          endif
        endfor
      endfor
    otherwise
      ## unsymmetrical case
      TAT = T.' * A * T;
  endswitch
endfunction

%!test
%! N = 10;
%! for t={"Upper", "Lower"}
%! for i=1:100
%! A = rand(N, N);
%! A += A.';
%! T = eye(N, N)(:, 1:floor(0.75*N));
%! TAT = fem_cms_matrix_trans(T, A, t);
%! assert_simple(TAT, T.' * A * T, eps);
%! assert_simple(issymmetric(TAT));
%! endfor
%! endfor

%!test
%! for t={"Upper", "Lower"}
%! N = 100;
%! for i=1:100
%! A = sprand(N, N, 0.1);
%! A += A.';
%! T = eye(N, N)(:, 1:floor(0.9*N));
%! TAT = fem_cms_matrix_trans(T, A, t);
%! assert_simple(TAT, T.' * A * T, eps);
%! assert_simple(issymmetric(TAT));
%! endfor
%! endfor

%!test
%! N = 1000;
%! for j=0:1
%! for i=1:100
%! A = sprand(N, N, 0.01) + abs(diag(rand(N, 1)));
%! A *= A.';
%! [r, c, d] = find(A);
%! if j
%! idx = find(r <= c);
%! mat_type = "Upper";
%! else
%! idx = find(r >= c);
%! mat_type = "Lower";
%! endif
%! r = r(idx);
%! c = c(idx);
%! d = d(idx);
%! Asym = sparse(r, c, d, N, N);
%! switch matrix_type(Asym)
%! case "Upper"
%! assert_simple(j, 1);
%! case "Lower"
%! assert_simple(j, 0);
%! otherwise
%! assert_simple(false);
%! endswitch
%! T = eye(N, N)(:, 1:floor(0.1 * N));
%! TAT = fem_cms_matrix_trans(T, Asym, mat_type);
%! assert_simple(TAT, T.' * A * T, eps);
%! assert_simple(issymmetric(TAT));
%! endfor
%! endfor

%!test
%! for i=1:100
%! A = [rand(10, 10), rand(10, 5);
%!      rand(5, 10), zeros(5, 5)];
%! A += A.';
%! T = rand(rows(A), floor(0.75 * columns(A)));
%! TAT = fem_cms_matrix_trans(T, A, "Lower");
%! assert_simple(TAT, T.' * A * T, eps^0.7);
%! endfor

%!test
%! A = rand(10, 10);
%! A *= A.';
%! T = rand(10, 3);
%! TAT1 = fem_cms_matrix_trans(T, A, "Upper");
%! TAT2 = fem_cms_matrix_trans(T, A, "Lower");
%! assert_simple(isdefinite(A));
%! assert_simple(isdefinite(TAT1));
%! assert_simple(TAT2, TAT1, eps * norm(TAT1));
