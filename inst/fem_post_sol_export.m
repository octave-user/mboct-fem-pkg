## Copyright (C) 2019(-2020) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} @var{post_pro_geo} = fem_post_sol_export(@var{prefix}, @var{mesh}, @var{sol}, @var{options})
## Export a finite element mesh and related solution to a set of files for post-processing with Gmsh.
##
## @var{prefix} @dots{} Prefix of output files
##
## @var{mesh} @dots{} Finite element mesh data structure
##
## @var{sol} @dots{} Finite element solution data structure
##
## @var{options}.elem_types @dots{} Cell array of element names to be exported
##
## @var{post_pro_geo} @dots{} Gmsh script filename which must be executed in order to load the results into Gmsh.
##
## @seealso{fem_post_sol_external, fem_post_mesh_export, fem_post_sol_step_export, fem_post_sol_script}
## @end deftypefn

function post_pro_geo = fem_post_sol_export(prefix, mesh, sol, options)
  if (nargin < 3 || nargin > 4 || nargout > 1 || ~isscalar(sol) || ~isstruct(sol) || ~isscalar(mesh) || ~isstruct(mesh))
    print_usage();
  endif

  if (nargin < 4)
    options = struct();
  endif

  if (~isfield(options, "elem_types"))
    eltypes = fem_pre_mesh_elem_type();
    ## Do not export surface elements by default
    ## because of issues related to the "Skin" plugin in Gmsh
    options.elem_types = {eltypes([eltypes.dim] == 3).name, "beam2", "beam3"};
  endif

  elem_types = fieldnames(mesh.elements);

  for i=1:numel(elem_types)
    switch (elem_types{i})
      case options.elem_types
      otherwise
        if (isfield(mesh.elements, elem_types{i}))
          mesh.elements = rmfield(mesh.elements, elem_types{i});
        endif
        
        if (isfield(mesh, "groups") && isfield(mesh.groups, elem_types{i}))
          mesh.groups = rmfield(mesh.groups, elem_types{i});
        endif
    endswitch
  endfor
  
  options.mesh_filename = [prefix, ".msh"];
  fem_post_mesh_export(options.mesh_filename, mesh);

  if (isfield(sol, "def"))
    num_steps = size(sol.def, 3);
  elseif (isfield(sol, "theta"))
    num_steps = size(sol.theta, 2);
  elseif (isfield(sol, "p"))
    num_steps = size(sol.p, 2);
  else
    error("there are no nodal fields in argument sol");
  endif
  
  if (isfield(sol, "t"))
    t = sol.t;
    options.show_time = 1;
  elseif (isfield(sol, "f"))
    t = sol.f;
    options.show_time = 6;
  else
    t = 1:num_steps;
    options.show_time = 1;
  endif

  if (~isfield(options, "output_step_idx"))
    options.output_step_idx = 1:num_steps;
  endif

  options.num_time_steps = numel(options.output_step_idx);

  options.result_filenames = cell(1, options.num_time_steps);
  
  for j=1:options.num_time_steps
    options.result_filenames{j} = sprintf("%s_%03d.msh", prefix, j);
    fem_post_sol_step_export(options.result_filenames{j}, sol, options.output_step_idx(j), j, t(options.output_step_idx(j)), 1);
  endfor

  post_pro_geo = [prefix, ".geo"];
  
  if (isfield(sol, "stress"))
    options.stress_output = fieldnames(sol.stress);
  else
    options.stress_output = {};
  endif
  
  fem_post_sol_script(post_pro_geo, options);
endfunction

