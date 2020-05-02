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
## @deftypefn {Function File} @var{sol} = fem_post_def_import(@var{filename}, @var{mesh}, @var{num_load_case})
## Import deformations from Nastran pch files
##
## @var{filename} @dots{} Nastran pch filename
##
## @var{mesh} @dots{} Finite element mesh data structure
##
## @var{sol} @dots{} Finite element solution data structure
## @end deftypefn

function sol = fem_post_def_import(filename, mesh, num_load_case)
  if (nargin < 2 || nargin > 3 || nargout > 1)
    print_usage();
  endif

  if (nargin < 3)
    num_load_case = 1;
  endif
  
  fd = -1;
  
  unwind_protect
    [fd] = fopen(filename, "rt");

    if (fd == -1)
      error("failed to open file \"%s\"", filename);
    endif

    sol.def = zeros(rows(mesh.nodes), 6, num_load_case);
    
    subcase_idx = int32(0);
    node_idx = int32(0);
    
    while (true)
      line = fgetl(fd);

      if (~ischar(line))
        break;
      endif

      [subcase_id, count] = sscanf(line, "$SUBCASE ID = %d\n", "C");

      if (count == 1)
        ++subcase_idx;
        node_idx = int32(0);
        continue;
      endif

      [node_id, type, u, v, w, count] = sscanf(line, " %d %s %g %g %g", "C");

      if (count == 5)
        sol.def(++node_idx, :, subcase_idx) = [u, v, w, zeros(1, 3)];
      endif
    endwhile
  unwind_protect_cleanup
    if (fd ~= -1)
      fclose(fd);
    endif
  end_unwind_protect
endfunction
