## fem_pre_mesh_elem_angle.m:01
%!test
%! try
%! pkg load mbdyn_util_oct;
%! a = 70e-3;
%! b = 20e-3;
%! c = 10e-3;
%! tol = eps;
%! X = [ 0.5 * a,  0.5 * b,  0.5 * c;  #  1
%!             0,  0.5 * b,  0.5 * c;  #  2
%!             0, -0.5 * b,  0.5 * c;  #  3
%!       0.5 * a, -0.5 * b,  0.5 * c;  #  4
%!       0.5 * a,  0.5 * b, -0.5 * c;  #  5
%!             0,  0.5 * b, -0.5 * c;  #  6
%!             0, -0.5 * b, -0.5 * c;  #  7
%!       0.5 * a, -0.5 * b, -0.5 * c,  #  8
%!             a,  0.5 * b,  0.5 * c;  #  9
%!             a, -0.5 * b,  0.5 * c;  # 10
%!             a,  0.5 * b, -0.5 * c;  # 11
%!             a, -0.5 * b, -0.5 * c]; # 12
%! Phi1 = [0, 20, 45, 330] * pi / 180;
%! Phi2 = [0, 60, 270, 285] * pi / 180;
%! Phi3 = [0, -15, 30, -195] * pi / 180;
%! for j=1:length(Phi1)
%!   R1 = euler123_to_rotation_matrix([Phi1(j); Phi2(j); Phi3(j)]);
%!   mesh.nodes = [X * R1.', zeros(rows(X), 3)];
%!   mesh.elements.iso8 = int32([1:8;
%!                             9, 1, 4, 10, 11, 5, 8, 12]);
%!   mesh.materials.iso8 = int32([1; 1]);
%!   angle_data(j) = fem_pre_mesh_elem_angle(mesh);
%! endfor
%! for j=1:numel(angle_data)
%!   assert_simple(all(all(all(abs(angle_data(j).iso8 - pi / 2) < 2 * pi * tol))));
%! endfor
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
