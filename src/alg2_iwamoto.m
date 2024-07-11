function [V, converged, i, stats] = alg2_iwamoto(Ybus, Sbus, V0, ref, pv, pq, mpopt)
%NEWTONPF_S_CART  Solves power flow using Iwamoto's method (power/cartesian)
% [1] IWAMOTO S, TAMURA Y. A Load Flow Calculation Method for Ill-Conditioned Power Systems[J/OL]. 
% IEEE Transactions on Power Apparatus and Systems, 1981, PAS-100(4): 1736-1743. DOI:10.1109/TPAS.1981.316511.

%   [V, CONVERGED, I] = alg2_iwamoto(YBUS, SBUS, V0, REF, PV, PQ, MPOPT)
%
%   Solves for bus voltages using a full Newton-Raphson method, using nodal
%   power balance equations and cartesian coordinate representation of
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
%   See also RUNPF, NEWTONPF, NEWTONPF_I_POLAR, NEWTONPF_I_CART.

%   MATPOWER
%   Copyright (c) 1996-2019, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Baljinnyam Sereeter, Delft University of Technology
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%% default arguments
if nargin < 7
    mpopt = mpoption;
end

%% options
tol         = mpopt.pf.tol;
max_it      = mpopt.pf.nr.max_it;
lin_solver  = mpopt.pf.nr.lin_solver;
stats = struct('ni',0,'timecost',0,'h',zeros(max_it,1),'nfeval',0, ...
    'dev',zeros(max_it+1,1),'nJeval',0,'nHzeval',0,'ndecompose',0);
stats.m = 'IWAMOTO';
%% initialize
converged = 0;
i = 0;
V = V0;
Vm = abs(V);
Vmpv = Vm(pv);

%% set up indexing for updating V
nb = length(V);
npv = length(pv);
npq = length(pq);
j1 = 1;         j2 = npq;           %% j1:j2 - Vr of pq buses
j3 = j2 + 1;    j4 = j2 + npv;      %% j3:j4 - Vr of pv buses
j5 = j4 + 1;    j6 = j4 + npq;      %% j5:j6 - Vi of pq buses
j7 = j6 + 1;    j8 = j6 + npv;      %% j7:j8 - Vi of pv buses

%% evaluate F(x0)
mis = V .* conj(Ybus * V) - Sbus(Vm);
F = [   real(mis([pq; pv]));
        imag(mis(pq));
        V(pv) .* conj(V(pv)) - Vmpv.^2  ];
stats.nfeval = stats.nfeval + 1;
%% check tolerance
normF = norm(F, inf);
stats.dev(1)=normF;
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
    nx = length(F);
    if nx <= 10 || have_feature('octave')
        lin_solver = '\';       %% default \ operator
    else    %% MATLAB and nx > 10 or Octave and nx > 2000
        lin_solver = 'LU3';     %% LU decomp with 3 output args, AMD ordering
    end
end

%% do Newton iterations
while (~converged && i < max_it)
    %% update iteration counter
    i = i + 1;

    %% evaluate Jacobian
    [dSbus_dVr, dSbus_dVi] = dSbus_dV(Ybus, V, 1);
    dV2_dVr = sparse(1:npv, npq+(1:npv), 2*real(V(pv)), npv, npv+npq);
    dV2_dVi = sparse(1:npv, npq+(1:npv), 2*imag(V(pv)), npv, npv+npq);

    %% handling of derivatives for voltage dependent loads
    %% (not yet implemented) goes here

    j11 = real(dSbus_dVr([pq; pv], [pq; pv]));
    j12 = real(dSbus_dVi([pq; pv], [pq; pv]));
    j21 = imag(dSbus_dVr(pq, [pq; pv]));
    j22 = imag(dSbus_dVi(pq, [pq; pv]));
    j31 = dV2_dVr;
    j32 = dV2_dVi;

    J = [   j11 j12;
            j21 j22;
            j31 j32;    ];

    %% compute update step
    dx = mplinsolve(J, -F, lin_solver);
    stats.ndecompose=stats.ndecompose+1;
    %% compute c0, c1, c2, c3
    c0 = -F; %F(x)
    c1 = -J*dx; %J(x)@dx
    c2 = -Fdx(Ybus,dx, pv, pq, nb, npv, npq, j1, j2, j3, j4,j5,j6,j7,j8);%F(dx)
    stats.nfeval = stats.nfeval + 1;
    %% compute g0, g1, g2, g3
    g0 = (c0.' * c1);
    g1 = c1.' * c1 + 2*c0.'*c2;
    g2 = 3*(c1.'*c2);
    g3 = 2*(c2.'*c2);
    r = roots([g3, g2, g1, g0]);
%     if real(r(3))>0 && real(r(3))<=1
        mu = real(r(3));
%         fprintf('iwamoto muliplier: %f\n', mu)
%     else
%         fprintf('Iwamoto mu roots are [%.2f + 1j*%.2f, %.2f + 1j*%.2f, %.2f + 1j*%.2f]\n', ...
%             real(r(1)), imag(r(1)), ...
%             real(r(2)), imag(r(2)), ...
%             real(r(3)), imag(r(3)))
%         mu = 1;
%     end
    %% compute mu
    dx = dx*mu;
    %% update voltage
    if npv
        V(pv) = V(pv) + dx(j3:j4) + 1j * dx(j7:j8);
    end
    if npq
        V(pq) = V(pq) + dx(j1:j2) + 1j * dx(j5:j6);
    end

    %% evalute F(x)
    mis = V .* conj(Ybus * V) - Sbus(Vm);
    F = [   real(mis([pq; pv]));
            imag(mis(pq));
            V(pv) .* conj(V(pv)) - Vmpv.^2  ];
    stats.nfeval = stats.nfeval + 1;
    %% check for convergence
    normF = norm(F, inf);
    if mpopt.verbose > 1
        fprintf('\n%3d        %10.3e', i, normF);
    end
    if normF < tol
        converged = 1;
        if mpopt.verbose
            fprintf('\nIwamoto power flow (power balance, cartesian) converged in %d iterations.\n', i);
        end
    end
    %% record V and deviation
    stats.dev(i+1)=normF;
end
stats.timecost=toc;
stats.ni=i;
if mpopt.verbose
    if ~converged
        fprintf('\nIwamoto power flow (power balance, cartesian) did not converge in %d iterations.\n', i);
    end
end
end
function Fdx_ = Fdx(Ybus,dx, pv, pq, nb, npv, npq, j1, j2, j3, j4,j5,j6,j7,j8)
% do not consider ZIP load case
    dV = zeros(nb,1);
    if npv
        dV(pv) = dx(j3:j4) + 1j * dx(j7:j8);
    end
    if npq
        dV(pq) = dx(j1:j2) + 1j * dx(j5:j6);
    end
    mis_ = dV .* conj(Ybus * dV);
    Fdx_ = [real(mis_([pq; pv]));
            imag(mis_(pq));
            dV(pv) .* conj(dV(pv))];
end