%!test
%! mesh_data(1).mesh.nodes = repmat(1000, 3, 6);
%! mesh_data(1).mesh.elements.iso8 = repmat(int32(1000), 3, 8);
%! mesh_data(1).mesh.materials.iso8 = int32([1;2;2]);
%! mesh_data(1).mesh.elements.line2 = repmat(int32(1100), 2, 2);
%! mesh_data(1).mesh.material_data(1).E = 1000;
%! mesh_data(1).mesh.material_data(1).nu = 0.1;
%! mesh_data(1).mesh.material_data(1).rho = 10000;
%! mesh_data(1).mesh.material_data(2).E = 1100;
%! mesh_data(1).mesh.material_data(2).nu = 0.11;
%! mesh_data(1).mesh.material_data(2).rho = 11000;
%! mesh_data(1).mesh.groups.iso8(1).id = 1000;
%! mesh_data(1).mesh.groups.iso8(1).name = "iso8 1000";
%! mesh_data(1).mesh.groups.iso8(1).nodes = int32(1:3)(:);
%! mesh_data(1).mesh.groups.iso8(1).elements = int32(1:3)(:);
%! mesh_data(1).mesh.groups.iso8(2).id = 1100;
%! mesh_data(1).mesh.groups.iso8(2).name = "iso8 1100";
%! mesh_data(1).mesh.groups.iso8(2).nodes = int32(1:2)(:);
%! mesh_data(1).mesh.groups.iso8(2).elements = int32(1:2)(:);

%! mesh_data(2).mesh.nodes = repmat(2000, 2, 6);
%! mesh_data(2).mesh.elements.iso8 = repmat(int32(2000), 2, 8);
%! mesh_data(2).mesh.materials.iso8 = int32([1;2]);
%! mesh_data(2).mesh.elements.line2 = repmat(int32(2100), 4, 2);
%! mesh_data(2).mesh.material_data(1).E = 2000;
%! mesh_data(2).mesh.material_data(1).nu = 0.2;
%! mesh_data(2).mesh.material_data(1).rho = 20000;
%! mesh_data(2).mesh.material_data(2).E = 2200;
%! mesh_data(2).mesh.material_data(2).nu = 0.22;
%! mesh_data(2).mesh.material_data(2).rho = 22000;

%! mesh_data(2).mesh.groups.iso8(1).id = 2000;
%! mesh_data(2).mesh.groups.iso8(1).name = "iso8 2000";
%! mesh_data(2).mesh.groups.iso8(1).nodes = int32(1:2)(:);
%! mesh_data(2).mesh.groups.iso8(1).elements = int32(1:2)(:);
%! mesh_data(2).mesh.groups.iso8(2).id = 2100;
%! mesh_data(2).mesh.groups.iso8(2).name = "iso8 2100";
%! mesh_data(2).mesh.groups.iso8(2).nodes = int32(1)(:);
%! mesh_data(2).mesh.groups.iso8(2).elements = int32(1)(:);

%! mesh_data(3).mesh.nodes = repmat(3000, 4, 6);
%! mesh_data(3).mesh.elements.iso8 = repmat(int32(3000), 2, 8);
%! mesh_data(3).mesh.materials.iso8 = int32([1;1]);
%! mesh_data(3).mesh.elements.line2 = repmat(int32(3100), 5, 2);
%! mesh_data(3).mesh.material_data.E = 3000;
%! mesh_data(3).mesh.material_data.nu = 0.3;
%! mesh_data(3).mesh.material_data.rho = 30000;
%! mesh_data(3).mesh.groups.iso8(1).id = 3000;
%! mesh_data(3).mesh.groups.iso8(1).name = "iso8 3000";
%! mesh_data(3).mesh.groups.iso8(1).nodes = int32(1:2)(:);
%! mesh_data(3).mesh.groups.iso8(1).elements = int32(1:2)(:);

%! mesh_data(4).mesh.nodes = repmat(4000, 4, 6);
%! mesh_data(4).mesh.elements.iso20 = repmat(int32(4000), 2, 20);
%! mesh_data(4).mesh.materials.iso20 = int32([1;1]);
%! mesh_data(4).mesh.elements.line2 = repmat(int32(4100), 5, 2);
%! mesh_data(4).mesh.material_data.E = 4000;
%! mesh_data(4).mesh.material_data.nu = 0.3;
%! mesh_data(4).mesh.material_data.rho = 40000;
%! mesh_data(4).mesh.groups.iso20(1).id = 4000;
%! mesh_data(4).mesh.groups.iso20(1).name = "iso8 4000";
%! mesh_data(4).mesh.groups.iso20(1).nodes = int32(1:2)(:);
%! mesh_data(4).mesh.groups.iso20(1).elements = int32(1:2)(:);

%! [mesh, offset] = fem_pre_mesh_merge(mesh_data);

%! grp_idx_1000 = find([mesh.groups.iso8.id] == 1000);
%! assert(all(mesh.elements.iso8(mesh.groups.iso8(grp_idx_1000).elements) == 1000));
%! grp_idx_2000 = find([mesh.groups.iso8.id] == 2000);
%! assert(all(mesh.elements.iso8(mesh.groups.iso8(grp_idx_2000).elements) == 2000 + rows(mesh_data(1).mesh.nodes)));
%! grp_idx_3000 = find([mesh.groups.iso8.id] == 3000);
%! assert(all(mesh.elements.iso8(mesh.groups.iso8(grp_idx_3000).elements) == 3000 + rows(mesh_data(1).mesh.nodes) + rows(mesh_data(2).mesh.nodes)));
%! grp_idx_4000 = find([mesh.groups.iso20.id] == 4000);
%! assert(all(mesh.elements.iso20(mesh.groups.iso8(grp_idx_4000).elements) == 4000 + rows(mesh_data(1).mesh.nodes) + rows(mesh_data(2).mesh.nodes) + rows(mesh_data(3).mesh.nodes)));
