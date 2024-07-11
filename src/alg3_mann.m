function [V, converged, i, stats] = alg3_mann(Ybus, Sbus, V0, ref, pv, pq, mpopt)
%NEWTONPF  Solves power flow using Mann's iteration method (power/polar)
%   [V, CONVERGED, I] = NEWTONPF(YBUS, SBUS, V0, REF, PV, PQ, MPOPT)
%
%   Solves for bus voltages using a full Newton-Raphson method, using nodal
%   power balance equations and polar coordinate representation of
%   voltages, given the following inputs:
%       YBUS  - full system admittance matrix (for all buses)
%       SBUS  - handle to function that returns the complex bus power
%               injection vector (for all buses), given the bus voltage
%               magnitude vector (for all buses)
%       V0    - initial vector of complex bus voltages
%       REF   - bus index of reference bus (voltage ang reference & gen slack)
%       PV    - vector of bus indices for PV buses
%       PQ    - vector of bus indices for PQ buses
%       MPOPT - (optional) MATPOWER option struct, used to set the
%               termination tolerance, maximum number of iterations, and
%               output options (see MPOPTION for details).
%
%   The bus voltage vector contains the set point for generator
%   (including ref bus) buses, and the reference angle of the swing
%   bus, as well as an initial guess for remaining magnitudes and
%   angles.
%
%   Returns the final complex voltages, a flag which indicates whether it
%   converged or not, and the number of iterations performed.
%
%   See also RUNPF, NEWTONPF_S_CART, NEWTONPF_I_POLAR, NEWTONPF_I_CART.

%   MATPOWER
%   Copyright (c) 1996-2019, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%% default arguments
if nargin < 7
    mpopt = mpoption;
end
lastwarn('')
%% options
tol         = mpopt.pf.tol;
max_it      = mpopt.pf.nr.max_it;
lin_solver  = mpopt.pf.nr.lin_solver;
%---- implemented options -----------------
rtol = mpopt.fmincon.tol_x;
atol = mpopt.fmincon.tol_f;
uround=eps; fac1=0.2; fac2=mpopt.mips.max_it; f_savety=0.9;  hmin = 1e5*eps/2; facmax = fac2;
hmax=1e7; 
step_control=mpopt.mips.step_control;
stats = struct('ni',0,'timecost',0,'h',zeros(max_it,1),'nfeval',0, ...
    'dev',zeros(max_it+1,1),'nJeval',0,'nHzeval',0,'ndecompose',0);
stats.m = 'mann';
%% initialize
converged = 0;
i = 0;
V = V0;
Va = angle(V);
Vm = abs(V);
h = mpopt.opf.violation;
%% set up indexing for updating V
npv = length(pv);
npq = length(pq);
j1 = 1;         j2 = npv;           %% j1:j2 - V angle of pv buses
j3 = j2 + 1;    j4 = j2 + npq;      %% j3:j4 - V angle of pq buses
j5 = j4 + 1;    j6 = j4 + npq;      %% j5:j6 - V mag of pq buses
%% set up identity
nv = 2*npq+npv; % number of total unknown variables
%% evaluate G(y)
mis = V .* conj(Ybus * V) - Sbus(Vm);
G = [   real(mis([pv; pq]));
    imag(mis(pq))   ];

%% check tolerance
normF = norm(G, inf);
if mpopt.verbose > 1
    fprintf('\n it    max P & Q mismatch (p.u.)');
    fprintf('\n----  ---------------------------');
    fprintf('\n%3d        %10.3e', i, normF);
end
if normF < tol
    converged = 1;
    if mpopt.verbose > 1
        fprintf('\nConverged!\n');
    end
end

%% attempt to pick fastest linear solver, if not specified
if isempty(lin_solver)
    nx = length(G);
    if nx <= 10 || have_feature('octave')
        lin_solver = '\';       %% default \ operator
    else    %% MATLAB and nx > 10 or Octave and nx > 2000
        lin_solver = 'LU3';     %% LU decomp with 3 output args, AMD ordering
    end
end

%% do Newton iterations
stats.dev(1)=normF;
mu = 0.5;
tic
while (~converged && i < max_it)
    %% update iteration counter, that is, RODAS steps
    i = i + 1;

    %% evaluate Jacobian
    J = makeJ(Ybus, V, Vm, Sbus, pv, pq);
    stats.nJeval=stats.nJeval+1;
    dx0 = mplinsolve(J, -G, lin_solver);
    stats.ndecompose = stats.ndecompose + 1;
    dx = mu * dx0;

    mu = min([norm(dx0, inf)^(-1) 1]);
    %% update voltage (y)
    if npv
        Va(pv) = Va(pv) + dx(j1:j2);
    end
    if npq
        Va(pq) = Va(pq) + dx(j3:j4);
        Vm(pq) = Vm(pq) + dx(j5:j6);
    end
    V = Vm .* exp(1j * Va);
    Vm = abs(V);            %% update Vm and Va again in case
    Va = angle(V);          %% we wrapped around with a negative Vm
    %% evalute G(y)
    mis = V .* conj(Ybus * V) - Sbus(Vm);
    G = [   real(mis([pv; pq]));
        imag(mis(pq))   ];

    %% check for convergence
    normF = norm(G, inf);
    if mpopt.verbose > 1
        fprintf('\n%3d        %10.3e', i, normF);
    end
    if normF < tol
        converged = 1;
        if mpopt.verbose
            fprintf('\n Mann converged in %d iterations.\n', i);
        end
    end
    %% record V and deviation
    stats.dev(i+1)=normF;
end
stats.timecost=toc;
stats.ni=i;
if mpopt.verbose
    if ~converged
        fprintf('\n Mann power flow (power balance, polar) did not converge in %d iterations.\n', i);
    end
end
