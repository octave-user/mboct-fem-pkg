## fem_get_software_version.m:01
%!test
%! try
%! version = fem_get_software_version("gmsh");
%! assert(size(version), [1, 3]);
%! assert(version(1) >= 4);
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
