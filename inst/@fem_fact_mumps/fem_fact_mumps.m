## Copyright (C) 2019(-2021) Reinhard <octave-user@a1.net>
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
## @deftypefn {} @var{Afact} = fem_fact_mumps(@var{A}, @var{opts})
## Create a factor object which uses MUMPS to solve @var{A} * @var{x} = @var{b} by means of @var{x} = @var{Afact} \ @var{b}
## @seealso{mumps}
## @end deftypefn

function Afact = fem_fact_mumps(A, opts)
  if (~(nargin == 2 && ismatrix(A) && issquare(A) && isstruct(opts)))
    print_usage();
  endif
  
  Afact.mumpsobj = mumps(A, opts);

  Afact = class(Afact, "fem_fact_mumps");
endfunction
