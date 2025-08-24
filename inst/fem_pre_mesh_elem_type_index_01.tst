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
