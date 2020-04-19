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
## @deftypefn {Function File} @var{load_case} = fem_pre_load_case_create_empty(@var{num_load_cases})
## Allocate a struct array of size <@var{num_load_cases}> of empty load cases
## @end deftypefn

function load_case = fem_pre_load_case_create_empty(num_load_cases)
  if (nargin ~= 1 || nargout > 1)
    print_usage();
  endif
  
  empty_data = cell(1, num_load_cases);
  
  load_case = struct("pressure", empty_data, ...
                     "locked_dof", empty_data, ...
                     "joints", empty_data, ...
                     "loads", empty_data, ...
                     "loaded_nodes", empty_data);
endfunction
