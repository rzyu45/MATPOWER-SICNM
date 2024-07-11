function [V, converged, i, stats] = alg4_romberg(Ybus, Sbus, V0, ref, pv, pq, mpopt)
%NEWTONPF  Solves power flow using romberg's iteration method (power/polar)
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
stats = struct('ni',0,'timecost',0,'h',zeros(max_it,1),'nfeval',0, ...
    'dev',zeros(max_it+1,1),'nJeval',0,'nHzeval',0,'ndecompose',0);
stats.m = 'romberg';
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
% declare F, J, phi
F = @(y) F_romberg(y, pv, pq, Vm, Va, Ybus, Sbus, j1, j2, j3, j4, j5, j6);
J = @(y) J_romberg(y, Ybus, Vm, Va, Sbus, pv, pq, j1, j2, j3, j4, j5, j6);
phi = @(u, v) mplinsolve(J(u), -F(v), lin_solver);
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
x = makey(Va, Vm, pv, pq);
r = 2;
psi = 8;
N = 4;
tic
while (~converged && i < max_it)
    %% update iteration counter, that is, RODAS steps
    i = i + 1;
    
    %% evaluate
    dx = phi(x, x);
    % update H
    H = max([0.1, min([norm(dx, inf)^(-1), 1])]);
    % update b
    b = x + H * dx;
    % romberg loop
    c = zeros(size(x, 1), N);
    h = H;
    for j = 1:N

        if j == 1
            c(:, j) = x;
        else
            c(:, j) = c(:, j-1);
        end

        h = h/2^(j-1);
        for l = 0:j-1
            if l == 0
                c(:, j) = c(:, j) + h*phi(x, x);
            else
                c(:, j) = c(:, j) + h*phi(x, c(:, l));
            end
        end

    end
    
    xnew = (r^psi*c(:, end) - b)/(r^psi - 1);

    %% update voltage (y)
    if npv
        Va(pv) = xnew(j1:j2);
    end
    if npq
        Va(pq) = xnew(j3:j4);
        Vm(pq) = xnew(j5:j6);
    end
    V = Vm .* exp(1j * Va);
    Vm = abs(V);            %% update Vm and Va again in case
    Va = angle(V);          %% we wrapped around with a negative Vm
    %% evalute G(y)
    mis = V .* conj(Ybus * V) - Sbus(Vm);
    G = [   real(mis([pv; pq]));
        imag(mis(pq))   ];
    x=xnew;
    N = max([N-1, 1]);
    %% check for convergence
    normF = norm(G, inf);
    if mpopt.verbose > 1
        fprintf('\n%3d        %10.3e', i, normF);
    end
    if normF < tol
        converged = 1;
        if mpopt.verbose
            fprintf('\n Romberg converged in %d iterations.\n', i);
        end
    end
    if any(isnan(x)) || any(x>1000)
        break
    end
    %% record V and deviation
    stats.dev(i+1)=normF;
end
stats.timecost=toc;
stats.ni=i;
if mpopt.verbose
    if ~converged
        fprintf('\n Romberg power flow (power balance, polar) did not converge in %d iterations.\n', i);
    end
end
