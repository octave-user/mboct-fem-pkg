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
%!   assert_simple(eltype(i).dim == 2 || eltype(i).dim == 3);
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
%!   assert_simple(eltype(i).dim == 1 || eltype(i).dim == 2);
%! endfor
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
