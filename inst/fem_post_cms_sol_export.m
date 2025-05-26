## Copyright (C) 2011(-2025) Reinhard <octave-user@a1.net>
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} fem_post_cms_sol_export(@var{cms_data}, @var{input_file}, @var{output_file}, @var{opt_scale}, @var{opt_post})
## Expand and export the solution from MBDyn models with modal elements.
## @seealso{fem_post_cms_sol_import}
## @end deftypefn

function fem_post_cms_sol_export(cms_data, input_file, output_file, opt_scale, opt_post)
  if (~isfield(opt_scale, "scale_type"))
    opt_scale.scale_type = "least square";
  endif
  
  if (~isfield(opt_scale, "scale"))
    opt_scale.scale = 1;
  endif
  
  if (~isfield(opt_post, "print_and_exit"))
    opt_post.print_and_exit = false;
  endif
  
  if (~isfield(opt_scale, "output_stress"))
    opt_scale.output_stress = FEM_SCA_STRESS_VMIS;
  endif
  
  if (~isfield(opt_post, "step_start"))
    opt_post.step_start = -1;
  endif
  
  if (~isfield(opt_post, "step_inc"))
    opt_post.step_inc = -1;
  endif
  
  if (~isfield(opt_post, "step_end"))
    opt_post.step_end = intmax;
  endif
  
  opt_mesh.group_id = "unique";
  res.log_dat = mbdyn_post_load_log(input_file);
  res.sol_dyn = fem_post_cms_sol_import(input_file, cms_data);

  num_steps_tot = numel(res.sol_dyn.t);
  opt_post.step_end = min(num_steps_tot, opt_post.step_end);

  if (opt_post.step_start > 0 || opt_post.step_inc > 0)
    if (opt_post.step_start < 0)
      opt_post.step_start = 1;
    endif

    if (opt_post.step_inc < 0)
      opt_post.step_inc = 1;
    endif

    idx_t = opt_post.step_start:opt_post.step_inc:opt_post.step_end;
  elseif (isfield(res.log_dat.vars, "output_steps_per_revolution"))
    num_steps_rev = res.log_dat.vars.output_steps_per_revolution;
    opt_post.step_inc = max([1, floor(num_steps_rev / 36)]);
    opt_post.step_start = max([1, opt_post.step_end - num_steps_rev]);
    idx_t = opt_post.step_start:opt_post.step_inc:opt_post.step_end;
  else
    idx_t = 1:opt_post.step_end;
  endif

  chunk_size = 5;
  [res.mesh_combined, res.dof_map_combined] = fem_post_mesh_merge(cms_data, opt_mesh);

  opt_post.result_filenames = cell(1, numel(idx_t));

  k = int32(0);

  opt_post.stress_output = {};
  opt_post.num_time_steps = numel(idx_t);
  opt_post.script = [output_file, "_struct.geo"];

  if (opt_post.print_and_exit)
    opt_post.print_to_file = [output_file, "_struct_step"];
  endif

  opt_post.mesh_filename = [output_file, "_struct.msh"];
  fem_post_mesh_export(opt_post.mesh_filename, res.mesh_combined, opt_post);

  for i=1:chunk_size:numel(idx_t)
    i0 = i;
    i1 = min([numel(idx_t), i0 + chunk_size - 1]);
    fprintf(stderr, "computing deformation and stress %d:%d/%d ...\n", i0, i1, numel(idx_t));
    
    sol_scaled = fem_post_cms_expand(res.sol_dyn, cms_data, idx_t(i0:i1), opt_scale);
    sol_combined = fem_post_cms_sol_merge(res.mesh_combined, res.dof_map_combined, sol_scaled);
    
    clear sol_scaled;
    
    if (isfield(sol_combined, "stress"))
      opt_post.stress_output = fieldnames(sol_combined.stress);
    endif
    
    for j=1:numel(i0:i1)
      ++k;
      opt_post.result_filenames{k} = sprintf("%s_struct_%03d.msh", output_file, k);      
      fprintf(stderr, "creating file %d/%d: %s\n", k, numel(opt_post.result_filenames), opt_post.result_filenames{k});      
      fem_post_sol_step_export(opt_post.result_filenames{k}, sol_combined, j, k, res.sol_dyn.t(idx_t(i0 + j - 1)), 1);
    endfor

    clear sol_combined;
  endfor

  fem_post_sol_script(opt_post.script, opt_post);

  if (opt_post.print_and_exit)
    pid = spawn("gmsh", {opt_post.script});
    status = spawn_wait(pid);

    if (status ~= 0)
      warning("gmsh failed with status %d", status);
    endif
  endif
endfunction
