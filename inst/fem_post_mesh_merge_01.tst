%!test
%! mesh_data(1).mesh.nodes = ones(15, 6) * 1000;
%! mesh_data(1).mesh.elements.penta15 = int32([1:15]);
%! mesh_data(1).mesh.materials.penta15 = int32(1);
%! mesh_data(1).mesh.material_data = struct("E", 100000e6, "nu", 0.3, "rho", 7850);
%! mesh_data(1).mesh.groups.penta15.id = int32(1000);
%! mesh_data(1).mesh.groups.penta15.nodes = unique(mesh_data(1).mesh.elements.penta15);
%! mesh_data(1).mesh.groups.penta15.elements = int32(1);
%! mesh_data(1).mesh.groups.penta15.name = "penta15 group 1";
%! mesh_data(2).mesh.nodes = ones(20+15, 6) * 2000;
%! mesh_data(2).mesh.elements.iso20 = int32([1:20]);
%! mesh_data(2).mesh.elements.penta15 = int32([1:15] + 20);
%! mesh_data(2).mesh.materials.penta15 = int32(1);
%! mesh_data(2).mesh.materials.iso20 = int32(1);
%! mesh_data(2).mesh.material_data = struct("E", 200000e6, "nu", 0.3, "rho", 7850);

%! mesh_data(2).mesh.groups.penta15.id = int32(2000);
%! mesh_data(2).mesh.groups.penta15.nodes = unique(mesh_data(2).mesh.elements.penta15);
%! mesh_data(2).mesh.groups.penta15.elements = int32(1);
%! mesh_data(2).mesh.groups.penta15.name = "penta15 group 2";

%! mesh_data(2).mesh.groups.iso20.id = int32(2100);
%! mesh_data(2).mesh.groups.iso20.nodes = unique(mesh_data(2).mesh.elements.iso20);
%! mesh_data(2).mesh.groups.iso20.elements = int32(1);
%! mesh_data(2).mesh.groups.iso20.name = "iso20 group 2";

%! mesh_data(3).mesh.nodes = ones(20 + 15, 6) * 3000;
%! mesh_data(3).mesh.elements.iso20 = int32([1:20]);
%! mesh_data(3).mesh.elements.penta15 = int32([1:15] + 20);
%! mesh_data(3).mesh.materials.penta15 = int32(1);
%! mesh_data(3).mesh.materials.iso20 = int32(1);
%! mesh_data(3).mesh.material_data = struct("E", 300000e6, "nu", 0.3, "rho", 7850);

%! mesh_data(3).mesh.groups.penta15.id = int32(3000);
%! mesh_data(3).mesh.groups.penta15.nodes = unique(mesh_data(3).mesh.elements.penta15);
%! mesh_data(3).mesh.groups.penta15.elements = int32(1);
%! mesh_data(3).mesh.groups.penta15.name = "penta15 group 3";

%! mesh_data(3).mesh.groups.iso20.id = int32(3100);
%! mesh_data(3).mesh.groups.iso20.nodes = unique(mesh_data(3).mesh.elements.iso20);
%! mesh_data(3).mesh.groups.iso20.elements = int32(1);
%! mesh_data(3).mesh.groups.iso20.name = "iso20 group 3";

%! mesh_data(4).mesh.nodes = ones(10, 6) * 4000;
%! mesh_data(4).mesh.material_data = struct()([]);
%! mesh_data(4).mesh.elements.point1 = int32(1:10)(:);
%! for i=1:10
%!   mesh_data(4).mesh.groups.point1(i).id = 4000 + i;
%!   mesh_data(4).mesh.groups.point1(i).nodes = i;
%!   mesh_data(4).mesh.groups.point1(i).elements = i;
%!   mesh_data(4).mesh.groups.point1(i).name = sprintf("point1 %d", i);
%! endfor

%! mesh1 = fem_post_mesh_merge(mesh_data);
%! mesh2 = fem_post_mesh_merge(mesh_data([4, 3, 2, 1]));
%! mesh3 = fem_post_mesh_merge(mesh_data([2, 1, 4, 3]));
%! mesh4 = fem_post_mesh_merge(mesh_data([4, 3, 2, 1]));
%! mesh5 = fem_post_mesh_merge(mesh_data([1, 3, 4, 2]));
%! mesh6 = fem_post_mesh_merge(mesh_data([2, 3, 4, 1]));
%! mesh7 = fem_post_mesh_merge(mesh_data([3, 1, 2, 4]));
