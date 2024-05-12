## fem_cms_red_guyan.m:02
%!test
%! try
%! E = 210000e6;
%! A = 100e-6;
%! l = 50e-3;
%! s = A * E / l;
%! K = [ s,    -s,     0,     0;
%!      -s, 2 * s,    -s,     0;
%!       0,    -s, 2 * s,    -s;
%!       0,     0,    -s,     s];
%!
%! [Kred, Tred] = fem_cms_red_guyan(K, [1, columns(K)]);
%!
%! s2 = A * E / (3 * l);
%! Kred2 = [ s2,    -s2;
%!          -s2,     s2];
%! assert_simple(Kred, Kred2, eps * norm(Kred2));
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
