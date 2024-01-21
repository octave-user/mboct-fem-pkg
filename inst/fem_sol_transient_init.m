## Copyright (C) 2021(-2021) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{data}] = fem_sol_transient_init(@var{M}, @var{D}, @var{K}, @var{dt}, @var{options})
## Initialize the transient solver
##
## @var{M} @dots{} Mass matrix
##
## @var{D} @dots{} Damping matrix
##
## @var{K} @dots{} Stiffness matrix
##
## @var{dt} @dots{} Time step
##
## @var{options}.delta @dots{} Coefficient for Newmark method
##
## @var{options}.alpha @dots{} Coefficient for Newmark method
##
## @var{options}.solver @dots{} String name of the linear solver
##
## @var{options}.number_of_threads @dots{} Number of threads for the linear solver
##
## @end deftypefn

function data = fem_sol_transient_init(M, D, K, dt, options)
  if (nargin ~= 5 || nargout ~= 1)
    print_usage();
  endif
  
  if (~isfield(options, "delta"))
    options.delta = 1/2;
  endif

  if (~isfield(options, "alpha"))
    options.alpha = (options.delta + 1/2)^2 / 4;
  endif
  
  data.a0 = 1 / ( options.alpha * dt^2 );
  data.a1 = options.delta / ( options.alpha * dt );
  data.a2 = 1 / ( options.alpha * dt );
  data.a3 = 1 / ( 2 * options.alpha ) - 1;
  data.a4 = options.delta / options.alpha - 1;
  data.a5 = dt / 2 * ( options.delta / options.alpha - 2 );
  data.a6 = dt * ( 1 - options.delta );
  data.a7 = options.delta * dt;

  data.M = M;
  data.D = D;
  data.Keff = K + data.a0 * M + data.a1 * D;
  data.Keff_fact = fem_sol_factor(data.Keff, options);
endfunction

