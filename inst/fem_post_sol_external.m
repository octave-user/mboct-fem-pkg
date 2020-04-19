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
## @deftypefn{Function File} fem_post_sol_external(@var{mesh}, @var{sol}, @var{varargin})
## Do postprocessing of a finite element solution within Gmsh.
##
## @var{mesh} @dots{} Finite element mesh data structure.
##
## @var{sol} @dots{} Finite element solution data structure.
##
## @var{varargin} @dots{} Additional arguments passed to fem_post_sol_export
##
## @seealso{fem_post_sol_export}
## @end deftypefn

function fem_post_sol_external(mesh, sol, varargin)
  if (nargin < 2 || nargin > 3 || nargout > 0)
    print_usage();
  endif
  
  prefix = "";

  unwind_protect
    prefix = tempname();
    
    post_pro_geo = fem_post_sol_export(prefix, mesh, sol, varargin{:});

    pid = spawn("gmsh", {post_pro_geo});

    status = spawn_wait(pid);
    
    fprintf(stderr, "gmsh returned with status %d\n", status);

  unwind_protect_cleanup
    if (numel(prefix))
      files = dir([prefix, "*"]);

      for i=1:numel(files)
        [err, msg] = unlink(fullfile(files(i).folder, files(i).name));

        if (err ~= 0)
          warning("failed to remove file \"%s\": %d (%s)", files(i).name, err, msg);
        endif
      endfor
    endif
  end_unwind_protect
endfunction
