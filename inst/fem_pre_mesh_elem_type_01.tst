%!test
%! try
%! eltype = fem_pre_mesh_elem_type();
%! assert_simple(numel(eltype), 39);
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
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
