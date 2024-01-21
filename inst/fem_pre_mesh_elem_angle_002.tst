## fem_pre_mesh_elem_angle.m:02
%!test
%! a = 70e-3;
%! b = 20e-3;
%! c = 10e-3;
%! tol = eps;
%!
%! alpha = [atan(b / a) + pi / 2, ...
%!          pi / 2, ...
%!          pi / 2, ...
%!          pi / 2 - atan(b / a)];
%!
%! X = [0.5 * b, 0.5 * a, c;
%!            0, 0.5 * a, c;
%!            0,       0, c;
%!            b,       0, c;
%!      0.5 * b, 0.5 * a, 0;
%!            0, 0.5 * a, 0;
%!            0,       0, 0;
%!            b,       0, 0];
%! Phi1 = [0, 20, 45, 330] * pi / 180;
%! Phi2 = [0, 60, 270, 285] * pi / 180;
%! Phi3 = [0, -15, 30, -195] * pi / 180;
%! for j=1:length(Phi1)
%!   R1 = euler123_to_rotation_matrix([Phi1(j); Phi2(j); Phi3(j)]);
%!   mesh.nodes = [X * R1.', zeros(rows(X), 3)];
%!   mesh.elements.iso8 = int32([1:8]);
%!   mesh.materials.iso8 = int32([1]);
%!   angle_data(j) = fem_pre_mesh_elem_angle(mesh);
%! endfor
%! for i=1:numel(angle_data)
%!   for j=1:rows(angle_data(j).iso8)
%!     for k=1:2
%!       assert_simple(all(abs(angle_data(j).iso8(j, 1:4, k) - alpha) < 2 * pi * tol));
%!       assert_simple(all(abs(angle_data(j).iso8(j, 5:end, k) - pi / 2) < 2 * pi * tol));
%!     endfor
%!   endfor
%! endfor
