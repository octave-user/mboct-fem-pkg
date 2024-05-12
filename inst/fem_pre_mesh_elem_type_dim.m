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
## @deftypefn {Function File} [@var{eltype}] = fem_pre_mesh_elem_type_dim(@var{idim})
## Return a list of supported elements
## @end deftypefn

function eltype = fem_pre_mesh_elem_type_dim(idim)
  eltype = fem_pre_mesh_elem_type();

  idx = false(1, numel(eltype));

  for i=1:numel(idim)
    idx |= [eltype.dim] == idim(i);
  endfor

  eltype = eltype(idx);
endfunction

%!test
%! try
%! for idim=[0,1,2,3]
%!   eltype = fem_pre_mesh_elem_type_dim(idim);
%!   assert_simple(~isempty(eltype));
%!   for j=1:numel(eltype)
%!     assert_simple(eltype(j).dim == idim);
%!   endfor
%! endfor

%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
%!test
%! try
%! eltype = fem_pre_mesh_elem_type_dim([2, 3]);
%! for i=1:numel(eltype)
%!   assert(eltype(i).dim == 2 || eltype(i).dim == 3);
%! endfor

%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
%!test
%! try
%! eltype = fem_pre_mesh_elem_type_dim([1, 2]);
%! for i=1:numel(eltype)
%!   assert(eltype(i).dim == 1 || eltype(i).dim == 2);
%! endfor
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
