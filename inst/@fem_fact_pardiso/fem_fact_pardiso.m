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
## @deftypefn {} @var{fact} = fem_fact_pardiso(@var{A}, @var{opts})
## Create a factor object which uses Intel MKL Pardiso to solve @var{A} * x = b via @var{x} = @var{Afact} \ @var{b}
## @seealso{pastix}
## @end deftypefn

function fact = fem_fact_pardiso(A, opts)
  if (~(nargin == 2 && ismatrix(A) && issquare(A) && isstruct(opts)))
    print_usage();
  endif
    
  fact.parobj = pardiso(A, opts);

  fact = class(fact, "fem_fact_pardiso");
endfunction
