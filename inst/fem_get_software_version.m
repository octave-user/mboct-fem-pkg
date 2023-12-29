## Copyright (C) 2023(-2023) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{version}] = fem_get_software_version(@var{command})
## Parse the software version of an external pre/post-processor (e.g. gmsh) and return the version number as a matrix.
##
## @end deftypefn

function [version] = fem_get_software_version(command)
  [status, version] = shell(sprintf("\"%s\" --version", command));

  if (status)
    error("command \"%s\" failed with status %d", command, status);
  endif

  [ver_major, ver_minor, ver_sub, count] = sscanf(version, "%d.%d.%d.%d.%d", "C");

  if (count < 2)
    error("failed to parse software version");
  endif

  version = [ver_major, ver_minor, ver_sub];
endfunction

%!test
%! version = fem_get_software_version("gmsh");
%! assert(size(version), [1, 3]);
%! assert(version(1) >= 4);
