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

