## fem_ehd_pre_comp_mat_load_case.m:01
%!test
%! close all;
%! fd = -1;
%! filename = "";
%! unwind_protect
%! filename = tempname();
%! if (ispc())
%!   filename(filename == "\\") = "/";
%! endif
%! unwind_protect
%! [fd, msg] = fopen([filename, ".geo"], "w");
%! if (fd == -1)
%!   error("failed to open file \"%s.geo\"", filename);
%! endif
%! ri = 8e-3;
%! ro = 10e-3;
%! h = 12e-3;
%! c = 2e-3;
%! b = h - 2 * c;
%! scale_def = 5e-3;
%! mesh_size = 5e-3;
%! f_post_pro = false;
%! fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%! fprintf(fd, "ri = %g;\n", ri);
%! fprintf(fd, "ro = %g;\n", ro);
%! fprintf(fd, "h = %g;\n", h);
%! fprintf(fd, "c = %g;\n", c);
%! fputs(fd, "Point(1) = {ri,0.0,0.0};\n");
%! fputs(fd, "Point(2) = {ro,0.0,0.0};\n");
%! fputs(fd, "Point(3) = {ro,0.0,c};\n");
%! fputs(fd, "Point(4) = {ro,0.0,h - c};\n");
%! fputs(fd, "Point(5) = {ro,0.0,h};\n");
%! fputs(fd, "Point(6) = {ri,0.0,h};\n");
%! fputs(fd, "Point(7) = {ri,0.0,h - c};\n");
%! fputs(fd, "Point(8) = {ri,0.0,c};\n");
%! fputs(fd, "Line(1) = {1,2};\n");
%! fputs(fd, "Line(2) = {2,3};\n");
%! fputs(fd, "Line(3) = {3,4};\n");
%! fputs(fd, "Line(4) = {4,5};\n");
%! fputs(fd, "Line(5) = {5,6};\n");
%! fputs(fd, "Line(6) = {6,7};\n");
%! fputs(fd, "Line(7) = {7,8};\n");
%! fputs(fd, "Line(8) = {8,1};\n");
%! fputs(fd, "Line Loop(5) = {1,2,3,4,5,6,7,8};\n");
%! fputs(fd, "Plane Surface(6) = {5};\n");
%! fputs(fd, "tmp[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} { Surface{6}; };\n");
%! fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%! fputs(fd, "Physical Surface(\"clamp\",1) = {tmp[2]};\n");
%! fputs(fd, "Physical Surface(\"load1\",2) = {tmp[4]};\n");
%! fputs(fd, "Physical Surface(\"load2\",3) = {tmp[8]};\n");
%! unwind_protect_cleanup
%!  if (fd ~= -1)
%!    fclose(fd);
%!  endif
%! end_unwind_protect
%! fprintf(stderr, "meshing ...\n");
%! pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-clmin", sprintf("%g", 0.75 * mesh_size), "-clmax", sprintf("%g", 1.25 *mesh_size), [filename, ".geo"]});
%! status = spawn_wait(pid);
%! if (status ~= 0)
%!  error("gmsh failed with status %d", status);
%! endif
%! fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%! mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%! fprintf(stderr, "%d nodes\n", rows(mesh.nodes));
%! grp_id_clamp = find([[mesh.groups.tria6].id] == 1);
%! grp_id_p1 = find([[mesh.groups.tria6].id] == 3);
%! grp_id_p2 = find([[mesh.groups.tria6].id] == 2);
%! bearing_surf(1).group_idx = grp_id_p1;
%! bearing_surf(1).options.reference_pressure = 1e6;
%! bearing_surf(1).options.mesh_size = 20e-3;
%! bearing_surf(1).r = ri;
%! bearing_surf(1).w = b;
%! bearing_surf(1).X0 = [0; 0; b/2 + c];
%! bearing_surf(1).R = eye(3);
%! bearing_surf(1).relative_tolerance = 0;
%! bearing_surf(1).absolute_tolerance = sqrt(eps) * ri;
%! bearing_surf(2).group_idx = grp_id_p2;
%! bearing_surf(2).options.reference_pressure = 1e6;
%! bearing_surf(2).options.mesh_size = 20e-3;
%! bearing_surf(2).r = ro;
%! bearing_surf(2).w = b;
%! bearing_surf(2).X0 = [0; 0; b/2 + c];
%! bearing_surf(2).R = eye(3);
%! bearing_surf(2).relative_tolerance = 0;
%! bearing_surf(2).absolute_tolerance = sqrt(eps) * ri;
%! [load_case, bearing_surf, idx_group] = fem_ehd_pre_comp_mat_load_case(mesh, bearing_surf);
%! load_case(1).locked_dof = false(rows(mesh.nodes), 6);
%! load_case(1).locked_dof(mesh.groups.tria6(grp_id_clamp).nodes, :) = true;
%! mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%! E = 210000e6;
%! nu = 0.3;
%! mesh.material_data.rho = 7850;
%! mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%! dof_map = fem_ass_dof_map(mesh, load_case(1));
%! opt_solver.refine_max_iter = int32(0);
%! fprintf(stderr, "assembling matrices ...\n");
%! [mat_ass.K, ...
%!  mat_ass.R, ...
%!  mat_ass.Rlumped] = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_MAT_STIFFNESS, ...
%!                                     FEM_VEC_LOAD_CONSISTENT, ...
%!                                     FEM_VEC_LOAD_LUMPED], ...
%!                                    load_case);
%! fprintf(stderr, "solving for static deflection (consistent load) ...\n");
%! sol_stat = fem_sol_static(mesh, dof_map, mat_ass, opt_solver);
%! fprintf(stderr, "solving for static deflection (lumped load) ...\n");
%! sol_stat_lumped = fem_sol_static(mesh, dof_map, setfield(mat_ass, "R", mat_ass.Rlumped), opt_solver);
%! if (f_post_pro)
%! post_pro_file = [filename, "_post_pro.geo"];
%! for j=1:numel(idx_group) - 1
%!  max_def = 0;
%!  for i=idx_group(j):idx_group(j + 1) - 1
%!    max_def = max(max_def, max(norm(sol_stat.def(:, 1:3, i), "rows")));
%!  endfor
%! for i=0:idx_group(j + 1) - idx_group(j) - 1
%!   deformation_file = sprintf("%s_%d_%03d_post_pro.msh", filename, j, i + 1);
%!   unwind_protect
%!   fem_post_sol_step_export(deformation_file, sol_stat, idx_group(j) + i, i + 1, i + 1, 1);
%!   fd = -1;
%!   unwind_protect
%!     [fd, msg] = fopen(post_pro_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\": %s", post_pro_file, msg);
%!     endif
%!     fprintf(fd, "Merge \"%s\";\n", [filename, ".msh"]);
%!     fputs(fd, "Mesh.SurfaceEdges = 0;\n");
%!     fputs(fd, "Mesh.SurfaceFaces = 0;\n");
%!     fputs(fd, "Mesh.SurfaceNumbers = 0;\n");
%!     fputs(fd, "Mesh.VolumeEdges = 0;\n");
%!     fputs(fd, "Mesh.VolumeFaces = 0;\n");
%!     fprintf(fd, "Merge \"%s\";\n", deformation_file);
%!     fputs(fd, "View[0].Type = 1;\n");
%!     fputs(fd, "View[0].VectorType = 5;\n");
%!     fputs(fd, "View[0].Visible = 1;\n");
%!     fprintf(fd, "View[0].DisplacementFactor = %g;\n", scale_def / max_def);
%!     fputs(fd, "View[0].ShowTime = 1;\n");
%!     fputs(fd, "View[0].ShowElement = 1;\n");
%!     fputs(fd, "View[0].IntervalsType = 3;\n");
%!     fputs(fd, "View[0].NbIso = 20;\n");
%!     fputs(fd, "View[0].DrawSkinOnly = 1;\n");
%!     fputs(fd, "General.Trackball = 0;\n");
%!     fputs(fd, "General.RotationX = -60;\n");
%!     fputs(fd, "General.RotationY = 10;\n");
%!     fputs(fd, "General.RotationZ = 40;\n");
%!     fputs(fd, "General.Axes = 3;\n");
%!     fputs(fd, "General.Orthographic = 1;\n");
%!     fputs(fd, "General.RotationCenterGravity = 1;\n");
%!     fprintf(fd, "View[0].TimeStep = %d;\n", i);
%!     fprintf(fd, "Print \"%s_%d_%03d.jpg\";\n", filename, j, i + 1);
%!     fputs(fd, "Exit;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   pid = spawn("gmsh", {post_pro_file});
%!   status = spawn_wait(pid);
%!   unwind_protect_cleanup
%!     [~] = unlink(post_pro_file);
%!     [~] = unlink(deformation_file);
%!   end_unwind_protect
%!   if (0 ~= status)
%!     error("gmsh failed with status %d", status);
%!   endif
%! endfor
%! jpg_filenames = dir(sprintf("%s_%d_*.jpg", filename, j));
%! for i=1:numel(jpg_filenames)
%!   [img, map, alpha] = imread(fullfile(jpg_filenames(i).folder, jpg_filenames(i).name));
%!   figure("visible", "off");
%!   imshow(img, map);
%!   title(sprintf("load case %d", i));
%! endfor
%! endfor
%! figure_list();
%! endif
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       status = unlink(fullfile(fn(i).folder, fn(i).name));
%!       if (status ~= 0)
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect
