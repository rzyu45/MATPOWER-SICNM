function [V, converged, i, stats] = sicnm(Ybus, Sbus, V0, ref, pv, pq, mpopt)
%NEWTONPF  Solves power flow using semi-implicit continuous Newton's method (power/polar)
%   [V, CONVERGED, I] = NEWTONPF(YBUS, SBUS, V0, REF, PV, PQ, MPOPT)
%
%   Solves for bus voltages using an LICN method, using nodal
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
Identity=speye(nv);
%% set up zero
Zero = sparse(nv,nv);
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

%% initialize z
% evaluate Jacobian
J0 = makeJ(Ybus, V, Vm, Sbus, pv, pq);
y = makey(Va,Vm,pv,pq);
z = mplinsolve(J0, -G, lin_solver);
M = sparse(1:nv,1:nv,1,2*nv,2*nv);%mass matrix
F = @(y,z) [z;tilde_g(y, z, pv, pq, Vm, Va, Ybus, Sbus, j1, j2, j3, j4, j5, j6)];
stats.dev(1)=normF;
%% initialize compute parameters alpha, b, gamma
if mpopt.mips.sc.red_it==1
    [stage,gamma,b,bd,~,alpha,gammatilde,~,~,~]=coeff_rodas;
    stats.m='rodas4';
    pord = 4;
elseif mpopt.mips.sc.red_it==2
    [stage,gamma,b,bd,alpha,gammatilde]=coeff_rodas3d;
    stats.m='rodas3d';
    pord = 3;
end
%% do Newton iterations
tic
while (~converged && i < max_it)
    %% update iteration counter, that is, RODAS steps
    i = i + 1;

    %% evaluate Jacobian
    J = makeJ(Ybus, V, Vm, Sbus, pv, pq);
    stats.nJeval=stats.nJeval+1;
    %% evaluate Hessian
    Hz = makeHz(Ybus, V, pv, pq, z, j1, j2, j3, j4, j5, j6);
    stats.nHzeval=stats.nHzeval+1;
    A=Hz+J;
    %% do RODAS
    err = 2;
    while err>1
        %% evaluate tilde_E, which contains h
        E2=-h*gamma*(J+h*gamma*A);
        q=amd(E2);
        [L1, U1, p] = lu(E2(q,q), 1.0, 'vector');
        stats.ndecompose=stats.ndecompose+1;
        L=[Identity,Zero;
            -h*gamma*A(q(p),:),L1];
        U=[Identity, -h*gamma*Identity(:,q);
            Zero, U1];
        p1=[1:nv,p+nv];
        q1=[1:nv,q+nv];
        K=zeros(2*nv,stage);
        %% compute update step
        rhs=F(y,z);
        stats.nfeval=stats.nfeval+1;
        K(q1,1)=U \ ( L \ rhs(q1(p1)) );
        [warnMsg, warnId] = lastwarn;
        if ~isempty(warnMsg)
            if warnId=="MATLAB:nearlySingularMatrix"
                errorStruct.message = 'Failed to Converge!';
                errorStruct.identifier = 'MATLAB:nearlySingularMatrix';
                error(errorStruct)
            end
        end
        for j=2:stage
            sum_1=K*alpha(:,j);sum_2=K*gammatilde(:,j);
            y1=y+h*sum_1(1:nv);
            z1=z+h*sum_1(nv+1:end);
            rhs=F(y1,z1)+M*sum_2;
            stats.nfeval=stats.nfeval+1;
            K(q1,j)=U \ ( L \ rhs(q1(p1)) );
            [warnMsg, warnId] = lastwarn;
            if ~isempty(warnMsg)
                if warnId=="MATLAB:nearlySingularMatrix"
                    errorStruct.message = 'Failed to Converge!';
                    errorStruct.identifier = 'MATLAB:nearlySingularMatrix';
                    error(errorStruct)
                end
            end
            K(:,j)=K(:,j)-sum_2;
        end
        sum_1 = K*(h*b); ynew = y + sum_1(1:nv); znew=z+sum_1(nv+1:end);
        sum_2 = K*(h*bd);
        %% ---- error test -----------
        if step_control
            SK = atol + rtol.*abs([ynew; znew]);
            err = max(abs((sum_1-sum_2)./SK));  %-- L-infinity Norm
            if ~isempty(find(~isfinite([ynew; znew])))
                err=1.0e6; disp('Warning Rodas5P: NaN or Inf occurs');
            end
            err = max(err,1.0e-6);
            if err<=1
                stats.h(i)=h;
            end
            %% update h
            fac = f_savety/err^(1/pord); fac=min(facmax,max(fac1,fac));
            hnew=h*fac;
            h=min(hnew,hmax);
            if h<hmin
                errorStruct.message = sprintf("h=%.4e < hmin\n", h);
                errorStruct.identifier = 'Matlab:TooSmallStep';
                error(errorStruct)
            end
        else
            err=0.1;
        end
    end
    %% update voltage (y)
    if npv
        Va(pv) = ynew(j1:j2);
    end
    if npq
        Va(pq) = ynew(j3:j4);
        Vm(pq) = ynew(j5:j6);
    end
    V = Vm .* exp(1j * Va);
    Vm = abs(V);            %% update Vm and Va again in case
    Va = angle(V);          %% we wrapped around with a negative Vm
    %% evalute G(y)
    mis = V .* conj(Ybus * V) - Sbus(Vm);
    G = [   real(mis([pv; pq]));
        imag(mis(pq))   ];
    y=ynew;
    z=znew;
    %% check for convergence
    normF = norm(G, inf);
    if mpopt.verbose > 1
        fprintf('\n%3d        %10.3e', i, normF);
    end
    if normF < tol
        converged = 1;
        if mpopt.verbose
            fprintf('\nSICNM %s power flow (power balance, polar) converged in %d iterations.\n', stats.m, i);
        end
    end
    %% record V and deviation
    stats.dev(i+1)=normF;
end
stats.timecost=toc;
stats.ni=i;
if mpopt.verbose
    if ~converged
        fprintf('\nSICNM %s power flow (power balance, polar) did not converge in %d iterations.\n', stats.m, i);
    end
end
