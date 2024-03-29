## Copyright (C) 2019(-2021) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} @var{load_case} = fem_pre_load_case_merge(@var{varargin})
## Combine several struct arrays of load cases into a single struct array
## @end deftypefn

function load_case = fem_pre_load_case_merge(varargin)
  load_cases = varargin;

  fn = {};
  
  for i=1:numel(load_cases)
    fn = {fn{:}, fieldnames(load_cases{i}){:}};
  endfor
  
  fn = unique(fn);
  
  for i=1:numel(load_cases)
    for j=1:numel(fn)
      if (~isfield(load_cases{i}, fn{j}))
        eval(sprintf("load_cases{i}(1).%s = [];", fn{j}));
      endif
    endfor
  endfor

  load_case = [load_cases{:}];
endfunction

