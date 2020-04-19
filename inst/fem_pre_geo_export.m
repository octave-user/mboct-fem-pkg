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
## @deftypefn {Function File} fem_pre_geo_export(@var{param_file}, @var{geo_file}, @var{param_val})
## Create a .geo file which can be loaded into Gmsh (http://gmsh.info).
## All Octave variables in <@var{param_val}> are written to <@var{param_file}> and <@var{geo_file}> is included at the end of <@var{param_file}>.
##
## @var{param_file} @dots{} Output filename which can be processed by Gmsh.
##
## @var{geo_file} @dots{} Gmsh source file which uses variables defined in <@var{param_val}>.
##
## @var{param_val} @dots{} Scalar struct with variables needed by <@var{geo_file}>.
## @seealso{fem_pre_mesh_unstruct_create}
## @end deftypefn

function fem_pre_geo_export(param_file, geo_file, param_val)
  if (nargin ~= 3 || nargout > 0)
    print_usage();
  endif
  
  [fd, msg] = fopen(param_file, "wt");

  if fd == -1
    error("failed to open file \"%s\": %s", param_file, msg);
  endif

  unwind_protect
    fn = fieldnames(param_val);

    for i=1:numel(fn)
      fval = getfield(param_val, fn{i});

      if (isscalar(fval) || ischar(fval))
        array_spec = "";
      elseif (ismatrix(fval) || iscell(fval))
        array_spec = "[]";
      else
        error("unsupported data type for parameter \"%s\"", fn{i});
      endif

      if (iscell(fval) && numel(fval) && ischar(fval{1}))
        prefix_str = "Str("; ## Required for Gmsh
        postfix_str = ")";
      else
        prefix_str = "";
        postfix_str = "";
      endif
      
      fprintf(fd, "%s%s = %s", fn{i}, array_spec, prefix_str);

      if (ischar(fval))
        fprintf(fd, "\"%s\"", fval);
      elseif (isscalar(fval))
        fprintf(fd, fem_format_gmsh(fn{i}, fval), fval);
      elseif (ismatrix(fval) || iscell(fval))
        fprintf(fd, "{");
        
        for j=1:numel(fval)
          if (iscell(fval))
            fvalj = fval{j};
          else
            fvalj = fval(j);
          endif
          
          fprintf(fd, fem_format_gmsh(fn{i}, fvalj), fvalj);
          
          if (j < numel(fval))
            fputs(fd, ", ");
          endif
        endfor
        fputs(fd, "}");
      else
        error("unsupported data type for parameter \"%s\"", fn{i});
      endif

      fprintf(fd, "%s;\n", postfix_str);
    endfor

    fprintf(fd, "Include \"%s\";\n", geo_file);
  unwind_protect_cleanup
    fclose(fd);
    fd = -1;
  end_unwind_protect
endfunction

function format = fem_format_gmsh(fname, fval)
  if (ischar(fval))
    format = "\"%s\"";
  elseif (isinteger(fval))
    format = "%d";
  elseif (isreal(fval))
    format = "%.16g";
  else
    error("unsupported data type for parameter \"%s\": %s", fname, class(fval));
  endif
endfunction
