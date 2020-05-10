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
## @deftypefn {Function File} @var{X} = mtimes(@var{Asym}, @var{b})
## Compute @var{X} = @var{Asym} * @var{b} by accessing only the lower or upper triangular part of @var{Asym}.
## @end deftypefn

function X = mtimes(Asym, b)
  if (nargin ~= 2)
    print_usage();
  endif                
  
  X = sp_sym_mtimes(Asym.A, b);
endfunction
