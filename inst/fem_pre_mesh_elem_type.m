## Copyright (C) 2019(-2024) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{eltype}] = fem_pre_mesh_elem_type()
## Return a list of supported elements
## @end deftypefn

function eltype_out = fem_pre_mesh_elem_type()
  if (nargin > 0 || nargout > 1)
    print_usage();
  endif

  persistent eltype = [];

  if (isempty(eltype))
    empty_cell = cell(1, 38);

    eltype = struct("dim", empty_cell, ...
                    "id", empty_cell, ...
                    "name", empty_cell, ...
                    "norder", empty_cell, ...
                    "nordernonp", empty_cell, ...
                    "promote", empty_cell, ...
                    "default_import", mat2cell(false(size(empty_cell)), 1, ones(size(empty_cell))));

    idx = int32(0);

    ++idx;

    eltype(idx).dim = 3;
    eltype(idx).id = 5;
    eltype(idx).name = "iso8";
    eltype(idx).norder = [5,6,7,8,1,2,3,4];
    eltype(idx).default_import = true;

    ++idx;

    eltype(idx).dim = 3;
    eltype(idx).id = 5;
    eltype(idx).name = "iso8upc";
    eltype(idx).norder = [5,6,7,8,1,2,3,4];

    ++idx;

    eltype(idx).dim = 3;
    eltype(idx).id = 5;
    eltype(idx).name = "iso8f";
    eltype(idx).norder = [5,6,7,8,1,2,3,4];

    ++idx;

    eltype(idx).dim = 3;
    eltype(idx).id = 17;
    eltype(idx).name = "iso20";
    eltype(idx).norder = [5,6,7,8,1,2,3,4,17,19,20,18,9,12,14,10,11,13,15,16];
    eltype(idx).default_import = true;

    ++idx;

    eltype(idx).dim = 3;
    eltype(idx).id = 17;
    eltype(idx).name = "iso20upc";
    eltype(idx).norder = [5,6,7,8,1,2,3,4,17,19,20,18,9,12,14,10,11,13,15,16];

    ++idx;

    eltype(idx).dim = 3;
    eltype(idx).id = 17;
    eltype(idx).name = "iso20f";
    eltype(idx).norder = [5,6,7,8,1,2,3,4,17,19,20,18,9,12,14,10,11,13,15,16];

    ++idx;

    eltype(idx).dim = 3;
    eltype(idx).id = 17;
    eltype(idx).name = "iso20r";
    eltype(idx).norder = [1,2,3,4,5,6,7,8,9,12,14,10,17,19,20,18,11,13,15,16];

    ++idx;

    eltype(idx).dim = 3;
    eltype(idx).id = 17;
    eltype(idx).name = "iso20fr";
    eltype(idx).norder = [1,2,3,4,5,6,7,8,9,12,14,10,17,19,20,18,11,13,15,16];

    ++idx;

    eltype(idx).dim = 3;
    eltype(idx).id = 17;
    eltype(idx).name = "iso20upcr";
    eltype(idx).norder = [1,2,3,4,5,6,7,8,9,12,14,10,17,19,20,18,11,13,15,16];

    ++idx;

    eltype(idx).dim = 3;
    eltype(idx).id = 12;
    eltype(idx).name = "iso27";
    eltype(idx).norder = [1,2,3,4,5,6,7,8,9,12,14,10,11,13,15,16,17,19,20,18,21,22,24,25,23,26,27];
    eltype(idx).default_import = true;

    ++idx;

    eltype(idx).dim = 3;
    eltype(idx).id = 12;
    eltype(idx).name = "iso27f";
    eltype(idx).norder = [1,2,3,4,5,6,7,8,9,12,14,10,11,13,15,16,17,19,20,18,21,22,24,25,23,26,27];

    ++idx;

    eltype(idx).dim = 3;
    eltype(idx).id = 6;
    eltype(idx).name = "penta6";
    eltype(idx).norder = [6,4,4,5,3,1,1,2];
    eltype(idx).promote = "iso8";
    eltype(idx).default_import = true;

    ++idx;

    eltype(idx).dim = 3;
    eltype(idx).id = 6;
    eltype(idx).name = "penta6f";
    eltype(idx).norder = [6,4,4,5,3,1,1,2];
    eltype(idx).promote = "iso8f";

    ++idx;

    eltype(idx).dim = 3;
    eltype(idx).id = 6;
    eltype(idx).name = "penta6upc";
    eltype(idx).norder = [6,4,4,5,3,1,1,2];
    eltype(idx).promote = "iso8upc";

    ++idx;

    eltype(idx).dim = 3;
    eltype(idx).id = 18;
    eltype(idx).name = "penta15";
    eltype(idx).norder = [1,2,3,4,5,6,7,10,8,13,15,14,9,11,12];
    eltype(idx).default_import = true;

    ++idx;

    eltype(idx).dim = 3;
    eltype(idx).id = 18;
    eltype(idx).name = "penta15f";
    eltype(idx).norder = [1,2,3,4,5,6,7,10,8,13,15,14,9,11,12];

    ++idx;

    eltype(idx).dim = 3;
    eltype(idx).id = 18;
    eltype(idx).name = "penta15upc";
    eltype(idx).norder = [1,2,3,4,5,6,7,10,8,13,15,14,9,11,12];

    ++idx;

    eltype(idx).dim = 3;
    eltype(idx).id = 7;
    eltype(idx).name = "pyra5";
    eltype(idx).norder = [5,5,5,5,1,2,3,4];
    eltype(idx).nordernonp = [1,2,3,4,5];
    eltype(idx).promote = "iso8";
    eltype(idx).default_import = true;

    ++idx;

    eltype(idx).dim = 3;
    eltype(idx).id = 7;
    eltype(idx).name = "pyra5f";
    eltype(idx).norder = [5,5,5,5,1,2,3,4];
    eltype(idx).nordernonp = [1,2,3,4,5];
    eltype(idx).promote = "iso8f";

    ++idx;

    eltype(idx).dim = 3;
    eltype(idx).id = 4;
    eltype(idx).name = "tet4";
    eltype(idx).norder = [4,4,4,4,1,2,3,3];
    eltype(idx).promote = "iso8";
    eltype(idx).default_import = true;

    ++idx;

    eltype(idx).dim = 3;
    eltype(idx).id = 4;
    eltype(idx).name = "tet4f";
    eltype(idx).norder = [4,4,4,4,1,2,3,3];
    eltype(idx).promote = "iso8f";

    ++idx;

    eltype(idx).dim = 3;
    eltype(idx).id = 11;
    eltype(idx).name = "tet10";
    eltype(idx).norder = [1,2,3,4,5,6,7,8,10,9];
    eltype(idx).default_import = true;

    ++idx;

    eltype(idx).dim = 3;
    eltype(idx).id = 11;
    eltype(idx).name = "tet10h";
    eltype(idx).norder = [1,2,3,4,5,6,7,8,10,9];

    ++idx;

    eltype(idx).dim = 3;
    eltype(idx).id = 11;
    eltype(idx).name = "tet10hf";
    eltype(idx).norder = [1,2,3,4,5,6,7,8,10,9];

    ++idx;

    eltype(idx).dim = 3;
    eltype(idx).id = 11;
    eltype(idx).name = "tet10upc";
    eltype(idx).norder = [1,2,3,4,5,6,7,8,10,9];

    ++idx;

    eltype(idx).dim = 3;
    eltype(idx).id = 29;
    eltype(idx).name = "tet20";
    eltype(idx).norder = [1,5,6,2,7,8,3,9,10,17,16,20,14,19,12,18,15,13,11,4];
    eltype(idx).default_import = true;

    ++idx;

    eltype(idx).dim = 2;
    eltype(idx).id = 9;
    eltype(idx).name = "tria6";
    eltype(idx).norder = [1,2,3,4,5,6];
    eltype(idx).default_import = true;

    ++idx;

    eltype(idx).dim = 2;
    eltype(idx).id = 2;
    eltype(idx).name = "tria3";
    eltype(idx).norder = [1,2,3,3];
    eltype(idx).nordernonp = [1,2,3];
    eltype(idx).promote = "iso4";
    eltype(idx).default_import = true;

    ++idx;

    eltype(idx).dim = 2;
    eltype(idx).id = 3;
    eltype(idx).name = "iso4";
    eltype(idx).norder = [1,2,3,4];
    eltype(idx).default_import = true;

    ++idx;

    eltype(idx).dim = 2;
    eltype(idx).id = 16;
    eltype(idx).name = "quad8";
    eltype(idx).norder = [1,2,3,4,5,6,7,8];
    eltype(idx).default_import = true;

    ++idx;

    eltype(idx).dim = 2;
    eltype(idx).id = 16;
    eltype(idx).name = "quad8r";
    eltype(idx).norder = [3,4,1,2,7,8,5,6];

    ++idx;

    eltype(idx).dim = 2;
    eltype(idx).id = 10;
    eltype(idx).name = "quad9";
    eltype(idx).norder = [1,2,3,4,5,6,7,8,9];
    eltype(idx).default_import = true;

    ++idx;

    eltype(idx).dim = 2;
    eltype(idx).id = 9;
    eltype(idx).name = "tria6h";
    eltype(idx).norder = [1,2,3,4,5,6];

    ++idx;

    eltype(idx).dim = 2;
    eltype(idx).id = 21;
    eltype(idx).name = "tria10";
    eltype(idx).norder = [1,2,3,4,5,6,7,8,9,10];
    eltype(idx).default_import = true;

    ++idx;

    eltype(idx).dim = 1;
    eltype(idx).id = 1;
    eltype(idx).name = "line2";
    eltype(idx).norder = [1,2];
    eltype(idx).default_import = true;

    ++idx;

    eltype(idx).dim = 1;
    eltype(idx).id = 8;
    eltype(idx).name = "line3";
    eltype(idx).norder = [1,2,3];
    eltype(idx).default_import = true;

    ++idx;

    eltype(idx).dim = 1;
    eltype(idx).id = 26;
    eltype(idx).name = "line4";
    eltype(idx).norder = [1,2,3,4];
    eltype(idx).default_import = true;

    ++idx;

    eltype(idx).dim = 0;
    eltype(idx).id = 15;
    eltype(idx).name = "point1";
    eltype(idx).norder = 1;
    eltype(idx).default_import = true;

    for i=1:numel(eltype)
      if (isempty(eltype(i).promote))
        elem_idx_promote_to = int32(-1);
      else
        elem_idx_promote_to = fem_pre_mesh_elem_type_index({eltype.name}, eltype(i).promote);
      endif

      eltype(i).promote = elem_idx_promote_to;
    endfor
  endif

  eltype_out = eltype;
endfunction

%!test
%! eltype = fem_pre_mesh_elem_type();
%! assert_simple(numel(eltype), 38);
%! for i=1:numel(eltype)
%!   switch (eltype(i).dim)
%!   case {0, 1, 2, 3}
%!   otherwise
%!     assert(false);
%!   endswitch
%!   if (eltype(i).default_import)
%!     default_flag = "*";
%!   else
%!     default_flag = " ";
%!   endif
%!   printf("%2d(%2d)%-1s%-10s->", i, eltype(i).id, default_flag, eltype(i).name);
%!   assert(~isempty(eltype(i).promote));
%!   elem_name_promote = "";
%!   if (eltype(i).promote ~= -1)
%!     assert(eltype(i).promote >= 1 && eltype(i).promote <= numel(eltype));
%!     elem_name_promote = eltype(eltype(i).promote).name;
%!   else
%!     assert(isempty(eltype(i).nordernonp));
%!   endif
%!   printf("%-10s", elem_name_promote);
%!   printf(" %d", eltype(i).norder);
%!   printf("\n");
%!   assert(~isempty(eltype(i).id));
%!   assert(~isempty(eltype(i).name));
%!   assert(~isempty(eltype(i).dim));
%!   assert(~isempty(eltype(i).norder));
%! endfor
