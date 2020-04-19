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
## @deftypefn {} @var{Afact} = fem_fact_lu(@var{A}, @var{opts})
## Create a factor object to solve @var{A} * @var{x} = @var{b} by means of Octave's builtin lu function via  @var{x} = @var{Afact} \ @var{b}
## @seealso{lu}
## @end deftypefn

function Afact = fem_fact_lu(A, opts)
  if (~(nargin == 2 && ismatrix(A) && issquare(A) && isstruct(opts)))
    print_usage();
  endif
  
  if (issparse(A))
    [Afact.L, Afact.U, Afact.P, Afact.Q, Afact.R] = lu(A);
  else
    [Afact.L, Afact.U, Afact.P] = lu(A);
    Afact.Q = eye(size(A));
    Afact.R = eye(size(A));
  endif
  
  Afact = class(Afact, "fem_fact_lu");
endfunction
