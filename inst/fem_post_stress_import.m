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
## @deftypefn {Function File} @var{sol} = fem_post_stress_import(@var{filename}, @var{mesh}, @var{num_load_cases})
## Import stresses from a Nastran PCH file.
##
## @var{filename} @dots{} Nastran PCH filename.
##
## @var{mesh} @dots{} Finite element mesh data structure
##
## @var{sol} @dots{} Finite element solution data structure
##
## @end deftypefn

function sol = fem_post_stress_import(filename, mesh, num_load_case)
  if (nargin < 2 || nargin > 3 || nargout > 1)
    print_usage();
  endif

  if (nargin < 3)
    num_load_case = 1;
  endif

  fd = -1;

  unwind_protect
    [fd, msg] = fopen(filename);

    if (fd == -1)
      error("failed to open file \"%s\"", filename);
    endif

    taun = zeros(rows(mesh.nodes), num_load_case);
    subcase_idx = int32(0);
    node_idx = int32(0);
    
    while (true)
      line = fgetl(fd);

      if (~ischar(line))
        break;
      endif

      [subcase_id, count] = sscanf(line, "LASTFALL %d\n", "C");

      if (count == 1)
        ++subcase_idx;
        node_idx = int32(0);
        continue;
      endif

      [node_id, tau, count] = sscanf(line, " %d %g", "C");

      if (count == 2)
        taun(++node_idx, subcase_idx) = tau;
      endif
    endwhile

    elem_types = {"iso4", "iso8", "tet10", "tria6"};

    sol.stress.vmis = struct();
    
    for i=1:numel(elem_types)
      if (isfield(mesh.elements, elem_types{i}))
        elem_nodes = getfield(mesh.elements, elem_types{i});
        vmis = zeros(rows(elem_nodes), columns(elem_nodes), columns(taun));
        for j=1:columns(taun)
          vmis(:, :, j) = taun(:, j)(elem_nodes);
        endfor
        sol.stress.vmis = setfield(sol.stress.vmis, elem_types{i}, vmis);
        clear vmis;
      endif
    endfor
  unwind_protect_cleanup
    if (fd ~= -1)
      fclose(fd);
    endif
  end_unwind_protect
endfunction
