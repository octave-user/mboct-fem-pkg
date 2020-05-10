## Copyright (C) 2020(-2020) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} @var{Asub} = subsref(@var{Asym}, @var{varargin})
## Return @var{Asub} = @var{Asym}(@var{varargin}).
## @end deftypefn

function Asub = subsref(Asym, s)
  if (nargin ~= 2)
    print_usage();
  endif

  if (~(isstruct(s) && isfield(s, "subs") && iscell(s.subs) && numel(s.subs) == 2))
    error("index operation not valid for class fem_mat_sym");
  endif
  
  i = s.subs{1};
  j = s.subs{2};

  if (numel(i) ~= numel(j) || any(i ~= j))
    error("row and column index must be the same for class fem_mat_sym");
  endif
  
  Asub = fem_mat_sym(Asym.A(i, j));
endfunction