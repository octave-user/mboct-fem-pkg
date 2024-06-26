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
## @deftypefn {Function File} [@var{idx}] = fem_pre_mesh_elem_type_index(@var{elem_names}, @var{name})
## Return the value @var{idx} for which strcmp(@var{elem_names}@{@var{idx}@}, @var{name}) returns true.
## @end deftypefn

function idx = fem_pre_mesh_elem_type_index(elem_names, name)
  idx = int32(find(cellfun(@(elem_name) strcmp(elem_name, name), elem_names, "UniformOutput", true)));
endfunction

%!test
%! try
%! eltype = fem_pre_mesh_elem_type();
%! for i=1:numel(eltype)
%!   idx = fem_pre_mesh_elem_type_index({eltype.name}, eltype(i).name);
%!   assert_simple(isscalar(idx));
%!   assert_simple(~isempty(idx));
%!   assert_simple(idx == i);
%! endfor
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
