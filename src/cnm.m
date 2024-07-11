function [V, converged, i, stats] = cnm(Ybus, Sbus, V0, ref, pv, pq, mpopt)
%NEWTONPF  Solves power flow using continuous Newton method (power/polar)
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
mtype=mpopt.mips.sc.red_it;
jmax=mpopt.mips.max_it;
hmin=1e5*eps/2;
stats = struct('ni',0,'timecost',0,'h',zeros(max_it,1),'nj',zeros(max_it,1), ...
    'dev',zeros(max_it+1,1));
switch mtype
    case 1
        stats.m='ICNM-JH';
    case 2
        stats.m='ICNM-J';
    case 3
        stats.m='ICNM-J1';
    case 4
        stats.m='ICNM-J0';
    case 5
        stats.m='CNM-ERK4';
        stats.nfeval=0;
end
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
%% evaluate G(x0)
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
stats.dev(1)=normF;
%% attempt to pick fastest linear solver, if not specified
if isempty(lin_solver)
    nx = length(G);
    if nx <= 10 || have_feature('octave')
        lin_solver = '\';       %% default \ operator
    else    %% MATLAB and nx > 10 or Octave and nx > 2000
        lin_solver = 'LU3';     %% LU decomp with 3 output args, AMD ordering
    end
end
%% initialize
y = makey(Va,Vm,pv,pq);
ynew=y;
nv = 2*npq+npv; % number of total unknown variables
if mtype==4
    if npv
        Va(pv) = ynew(j1:j2);
    end
    if npq
        Va(pq) = ynew(j3:j4);
        Vm(pq) = ynew(j5:j6);
    end
    V=Vm .* exp(1j * Va);
    J0=makeJ(Ybus,V,Vm,Sbus,pv,pq);
    dA0=decomposition(J0);
end
if mtype==5 || mtype==6
    F = @(y) f_cnm_e(y, pv, pq, Vm, Va, Ybus, Sbus, j1, j2, j3, j4, j5, j6);
end
%% do Newton iterations
tic
while (~converged && i < max_it)
    %% update iteration counter
    i = i + 1;

    dy=ones(nv,1);
    switch mtype
        case 1 % ICNM-JH
            for j = 1:jmax
                % phi
                if npv
                    Va(pv) = ynew(j1:j2);
                end
                if npq
                    Va(pq) = ynew(j3:j4);
                    Vm(pq) = ynew(j5:j6);
                end
                V=Vm .* exp(1j * Va);
                J=makeJ(Ybus,V,Vm,Sbus,pv,pq);
                mis = V .* conj(Ybus * V) - Sbus(Vm);
                G = [   real(mis([pv; pq]));
                    imag(mis(pq))   ];
                phi=J*(ynew-y)+h*G;
                %phiy
                Hydy=makeHz(Ybus, V, pv, pq, ynew-y, j1, j2, j3, j4, j5, j6);
                phiy=Hydy+(1+h)*J;
                dA=decomposition(phiy);
                dy=dA\phi;
                [warnMsg, warnId] = lastwarn;
                if ~isempty(warnMsg)
                    if warnId=="MATLAB:nearlySingularMatrix"
                        errorStruct.message = 'Failed to Converge!';
                        errorStruct.identifier = 'MATLAB:nearlySingularMatrix';
                        error(errorStruct)
                    end
                end
                ynew=ynew-dy;
                if norm(dy,inf)<mpopt.pf.tol && j<=jmax
                    stats.nj(i)=j;
                    stats.h(i)=h;
                    break
                end
                if mod(j,5)==0
                    h=h*0.75;
                end
                if h<hmin
                    errorStruct.message = sprintf("h=%.4e < hmin\n", h);
                    errorStruct.identifier = 'Matlab:TooSmallStep';
                    error(errorStruct)
                end
            end
            if j>0 && j<4
                h=h*1.25;
            end
            if j==jmax
                stats.h(i)=h;
                stats.nj(i)=j;
            end
        case 2  %ICNM-J
            for j = 1:jmax
                if npv
                    Va(pv) = ynew(j1:j2);
                end
                if npq
                    Va(pq) = ynew(j3:j4);
                    Vm(pq) = ynew(j5:j6);
                end
                V=Vm .* exp(1j * Va);
                J=makeJ(Ybus,V,Vm,Sbus,pv,pq);
                mis = V .* conj(Ybus * V) - Sbus(Vm);
                G = [   real(mis([pv; pq]));
                    imag(mis(pq))   ];
                dA=decomposition(J);
                dy=(ynew-y+h*(dA\G))/(1+h);
                [warnMsg, warnId] = lastwarn;
                if ~isempty(warnMsg)
                    if warnId=="MATLAB:nearlySingularMatrix"
                        errorStruct.message = 'Failed to Converge!';
                        errorStruct.identifier = 'MATLAB:nearlySingularMatrix';
                        error(errorStruct)
                    end
                end
                ynew=ynew-dy;
                if norm(dy,inf)<mpopt.pf.tol && j<=jmax
                    stats.nj(i)=j;
                    stats.h(i)=h;
                    break
                end
                if mod(j,5)==0
                    h=h*0.75;
                end
                if h<hmin
                    errorStruct.message = sprintf("h=%.4e < hmin\n", h);
                    errorStruct.identifier = 'Matlab:TooSmallStep';
                    error(errorStruct)
                end
            end
            if j>0 && j<4
                h=h*1.25;
            end
            if j==jmax
                stats.h(i)=h;
                stats.nj(i)=j;
            end
        case 3  %ICNM-J1
            if npv
                Va(pv) = ynew(j1:j2);
            end
            if npq
                Va(pq) = ynew(j3:j4);
                Vm(pq) = ynew(j5:j6);
            end
            V=Vm .* exp(1j * Va);
            J=makeJ(Ybus,V,Vm,Sbus,pv,pq);
            mis = V .* conj(Ybus * V) - Sbus(Vm);
            G = [   real(mis([pv; pq]));
                imag(mis(pq))   ];
            dA=decomposition(J);
            dy=(ynew-y+h*(dA\G))/(1+h);
            [warnMsg, warnId] = lastwarn;
            if ~isempty(warnMsg)
                if warnId=="MATLAB:nearlySingularMatrix"
                    errorStruct.message = 'Failed to Converge!';
                    errorStruct.identifier = 'MATLAB:nearlySingularMatrix';
                    error(errorStruct)
                end
            end
            ynew=ynew-dy;
        case 4  %ICNM-J0
            mis = V .* conj(Ybus * V) - Sbus(Vm);
            G = [   real(mis([pv; pq]));
                imag(mis(pq))   ];
            dy=(ynew-y+h*(dA0\G))/(1+h);
            [warnMsg, warnId] = lastwarn;
            if ~isempty(warnMsg)
                if warnId=="MATLAB:nearlySingularMatrix"
                    errorStruct.message = 'Failed to Converge!';
                    errorStruct.identifier = 'MATLAB:nearlySingularMatrix';
                    error(errorStruct)
                end
            end
            ynew=ynew-dy;
        case 5 % CNM-ERK4
                k1=F(y);
                k2=F(y+0.5*h*k1);
                k3=F(y+0.5*h*k2);
                k4=F(y+h*k3);
                stats.nfeval=stats.nfeval+4;
                dy= h*(k1+2*k2+2*k3+k4)/6;
                ynew=y+dy;
                epsilon=max(abs(k2-ynew));
                stats.h(i)=h;
                if epsilon>0.01 
                    h=max(0.985*h,0.75);
                else
                    h=min(1.015*h,0.75);
                end
        case 6
        err = 2;
            while err>1
                %
                K=zeros(nv,stage);
                f=zeros(nv,stage+1);
                % compute update step
                % stage=1;
                K(:,1)=F(y);
                f(:,1)=K(:,1);
                stats.nfeval=stats.nfeval+1;
                for s=2:stage
                    y1=y+h*K*B(:,s-1);
                    K(:,s)=F(y1);
                    f(:,s)=K(:,s);
                    stats.nfeval=stats.nfeval+1;
                    [warnMsg, warnId] = lastwarn;
                    if ~isempty(warnMsg)
                        if warnId=="MATLAB:nearlySingularMatrix"
                            errorStruct.message = 'Failed to Converge!';
                            errorStruct.identifier = 'MATLAB:nearlySingularMatrix';
                            error(errorStruct)
                        end
                    end
                end
                ynew = y + h*K*B(:,end);
                f(:,stage+1)=F(ynew);
                % ---- error test -----------
                SK = atol + rtol.*abs(ynew);
%                 err = sqrt(sum( ((fE)./SK).^2 )/nv);  %-- L_2 Norm
                %      err = sum(abs((sum_1-sum_2)./SK))/neq;  %-- L-1 Norm
                err = max(abs((f*E)./SK));  %-- L-infinity Norm
%                 err = max(abs((fE)./SK));  %-- L-infinity Norm
                if ~isempty(find(~isfinite(ynew)))
                    err=1.0e6; disp('Warning Rodas5P: NaN or Inf occurs');
                end
                err = max(err,1.0e-6);
                if err<=1
                    stats.h(i)=h;
                end
                %% update h
                fac = f_savety/err^(1/pow); fac=min(facmax,max(fac1,fac));
                hnew=h*fac;
                h=min(hnew,hmax);
                if h<hmin
                    errorStruct.message = sprintf("h=%.4e < hmin\n", h);
                    errorStruct.identifier = 'Matlab:TooSmallStep';
                    error(errorStruct)
                end
            end
    end
    %% update voltage
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

    %% evalute G(x)
    mis = V .* conj(Ybus * V) - Sbus(Vm);
    G = [   real(mis([pv; pq]));
        imag(mis(pq))   ];

    y=ynew;
    %% check for convergence
    normF = norm(G, inf);
    if mpopt.verbose > 1
        fprintf('\n%3d        %10.3e', i, normF);
    end
    if normF < tol
        converged = 1;
        if mpopt.verbose
            switch mtype
                case 1
                    m='ICNM-JH';
                case 2
                    m='ICNM-J';
                case 3
                    m='ICNM-J1';
                case 4
                    m='ICNM-J0';
                case 5
                    m='CNM-ERK4';
                case 6
                    m='CNM-DOPRI5(4)';
            end
            fprintf('\n%s power flow (power balance, polar) converged in %d iterations.\n', m, i);
        end
    end
    stats.dev(i+1)=normF;
end
stats.timecost=toc;
stats.ni=i;
if mpopt.verbose
    if ~converged
        switch mtype
            case 1
                m='ICNM-JH';
            case 2
                m='ICNM-J';
            case 3
                m='ICNM-J1';
            case 4
                m='ICNM-J0';
            case 5
                m='CNM-ERK4';
            case 6
                m='CNM-DOPRI5(4)';
        end
        fprintf('\n%s power flow (power balance, polar) did not converge in %d iterations.\n', m, i);
    end
end
