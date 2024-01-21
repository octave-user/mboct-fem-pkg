## fem_pre_mesh_reorder.m:01
%!test
%! mesh.nodes = [0, 0, 0;
%!               0, 1, 0;
%!               0, 2, 0;
%!               1, 0, 0;
%!               1, 1, 0;
%!               1, 2, 0;
%!               2, 0, 0;
%!               2, 1, 0;
%!               2, 2, 0];
%!
%! mesh.elements.iso4 = int32([1,2,5,4;
%!                             2,3,6,5;
%!                             4,5,8,7;
%!                             5,6,9,8]);
%! mesh.groups.iso4(1).elements = int32([1, 2]);
%! mesh.groups.iso4(2).elements = int32([3, 4]);
%! mesh.groups.iso4(1).id = 1;
%! mesh.groups.iso4(2).id = 2;
%! mesh.groups.iso4(1).name = "group1";
%! mesh.groups.iso4(2).name = "group2";
%!
%! for i=1:numel(mesh.groups.iso4)
%!   mesh.groups.iso4(i).nodes = sort(unique(mesh.elements.iso4(mesh.groups.iso4(i).elements, :)(:)));
%! endfor
%!
%! x = rand(rows(mesh.nodes), 10);
%! s = zeros(rows(mesh.elements.iso4), 1);
%!
%! for i=1:rows(mesh.elements.iso4)
%!   s(i) = sum(x(mesh.elements.iso4(i, :)));
%! endfor
%!
%! g = zeros(numel(mesh.groups.iso4), 1);
%!
%! for i=1:numel(mesh.groups.iso4)
%!   g(i) += sum(x(mesh.groups.iso4(i).nodes));
%! endfor
%!
%! algorithm = {"metis", "amd", "symrcm", "csymamd", "symamd"};
%!
%! for j=1:numel(algorithm)
%!   [mesh2, perm, iperm] = fem_pre_mesh_reorder(mesh, struct("algorithm", algorithm{j}));
%!
%!   x2 = x(perm, :);
%!   s2 = zeros(rows(mesh2.elements.iso4), 1);
%!
%!   for i=1:rows(mesh2.elements.iso4)
%!     s2(i) = sum(x2(mesh2.elements.iso4(i, :)));
%!   endfor
%!
%!   A = zeros(rows(mesh.nodes), rows(mesh.nodes));
%!
%!   for i=1:rows(mesh.elements.iso4)
%!     A(mesh.elements.iso4(i, :), mesh.elements.iso4(i, :)) += x(mesh.elements.iso4(i, :));
%!   endfor
%!
%!   A2 = zeros(rows(mesh2.nodes), rows(mesh2.nodes));
%!
%!   for i=1:rows(mesh2.elements.iso4)
%!     A2(mesh2.elements.iso4(i, :), mesh2.elements.iso4(i, :)) += x2(mesh2.elements.iso4(i, :));
%!   endfor
%!
%!   assert_simple(s2, s);
%!   assert_simple(mesh2.nodes(iperm,:), mesh.nodes)
%!   assert_simple(A2(iperm, iperm), A);
%!   assert_simple(A2, A(perm, perm));
%!   for i=1:numel(mesh.groups.iso4)
%!     assert_simple(mesh2.groups.iso4(i).elements, mesh.groups.iso4(i).elements);
%!     assert_simple(mesh2.groups.iso4(i).name, mesh.groups.iso4(i).name);
%!     assert_simple(perm(mesh2.groups.iso4(i).nodes), mesh.groups.iso4(i).nodes);
%!   endfor
%!   g2 = zeros(numel(mesh2.groups.iso4), 1);
%!   for i=1:numel(mesh2.groups.iso4)
%!    g2(i) += sum(x2(mesh2.groups.iso4(i).nodes));
%!   endfor
%! endfor
