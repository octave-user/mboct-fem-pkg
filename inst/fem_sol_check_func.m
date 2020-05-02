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
## @deftypefn {Function File} @var{flag} = fem_sol_check_func(@var{funcname})
##
## Check if the solver <@var{funcname}> is available.
##
## @seealso{fem_sol_linsolve}
## @end deftypefn

function flag = fem_sol_check_func(func)
  if (nargin ~= 1 || nargout > 1)
    print_usage();
  endif
  
  flag = false;
  
  switch (exist(func))
    case {3, 5}
    otherwise
      return;
  endswitch

  ## If the function was not compiled, but the PKG_ADD file
  ## was generated from comments in the source code by Octave's package manager,
  ## then exist incorrectly reports, that the function exists.

  fname = which(func);

  if (numel(fname))
    flag = true;
  endif
endfunction
