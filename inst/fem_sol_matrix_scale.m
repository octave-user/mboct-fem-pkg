## Copyright (C) 2021(-2021) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{AS}, @var{D1}, @var{D2}] = fem_sol_matrix_scale(@var{A}, @var{tol}, @var{maxiter})
## Compute row scaling factors @var{D1} and column scaling factors @var{D2} to improve the condition number of matrix @var{A}
##
## @var{A} @dots{} input matrix
##
## @var{tol} @dots{} tolerance for the maximum difference between each row/column norm from unity
##
## @var{maxiter} @dots{} maximum number of scaling iterations
##
## @var{AS} @dots{} scaled output matrix
##
## @var{D1} @dots{} row scaling factors
##
## @var{D2} @dots{} column scaling factors
##
## @end deftypefn

function [AS, D1, D2] = fem_sol_matrix_scale(A, tol, maxiter)
  if (nargin < 3)
    print_usage();
  endif
  
  D1 = ones(rows(A), 1);
  D2 = ones(1, columns(A));
  
  AS = A;
    
  iter = int32(0);

  max_rn = max_cn = realmax;
  
  do
    DR = sqrt(norm(AS, inf, "rows"));
    DC = sqrt(norm(AS, inf, "cols"));

    D1 ./= max(1, DR);
    D2 ./= max(1, DC);
    
    AS = diag(D1) * A * diag(D2);

    max_rn_prev = max_rn;
    max_cn_prev = max_cn;
    
    max_rn = max(abs(1 - DR));
    max_cn = max(abs(1 - DC));

    if (abs(max_rn - max_rn_prev) < tol && abs(max_cn - max_cn_prev) < tol)
      break;
    endif
    
    if (++iter > maxiter)
      error("maximum number of iterations exceeded when scaling the matrix");
    endif
  until (max_rn <= tol && max_cn <= tol)
endfunction
