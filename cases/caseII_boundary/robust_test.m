clc
clear
%% options
ntest=1000;
tol=1e-5;
mpoptnr = mpoption('verbose', 0,'out.all',0,'opf.start',1,'pf.nr.max_it', 40);
% mips.sc.red_it for cnm type
% opf.violation is used as setter for h
% fmincon.tol_x is used as setter for rtol
% fmincon.tol_f is used as setter for atol
% mips.max_it for fac2
% mips.max_it for maximum inner loop iteration
mpoptcnm = mpoption('pf.alg','cnm','verbose', 0,'out.all',0,'opf.start',1, ...
    'opf.violation', 0.1,'pf.nr.max_it', 40,...
    'mips.max_it', 30);
% opf.violation is used as setter for h
% fmincon.tol_x is used as setter for rtol
% fmincon.tol_f is used as setter for atol
% mips.max_it for fac2
% mips.step_control for step size control

mpoptrodas4 = mpoption('pf.alg', 'sicnm', 'verbose', 0,'out.all',0, ...
    'pf.nr.max_it', 40, 'opf.violation', 0.1,...
    'fmincon.tol_x', 1e-1,...
    'fmincon.tol_f', 1e-1, ...
    'mips.max_it',6.0,...
    'mips.step_control',1, ...
    'mips.sc.red_it',1);
mpoptrodas3d = mpoption('pf.alg', 'sicnm', 'verbose', 0,'out.all',0, ...
    'pf.nr.max_it', 40, 'opf.violation', 0.1,...
    'fmincon.tol_x', 1e-1,...
    'fmincon.tol_f', 1e-1, ...
    'mips.max_it',6.0,...
    'mips.step_control',1, ...
    'mips.sc.red_it',2);

mpoptIWA = mpoption('pf.alg', 'IWAMOTO', 'verbose', 0,'out.all',0, 'pf.nr.max_it', 40);
mpoptMANN = mpoption('pf.alg', 'MANN', 'verbose', 0,'out.all',0, 'pf.nr.max_it', 40);
mpoptROMBERG = mpoption('pf.alg', 'ROMBERG', 'verbose', 0,'out.all',0, 'pf.nr.max_it', 40);
mpopt2S4 = mpoption('pf.alg', '2S4', 'verbose', 0,'out.all',0, 'pf.nr.max_it', 40);

%%
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

mpc=case9241pegase;name="9241";

% convert to internal indexing
mpc = ext2int(mpc);

% feed 1/2 load with random initials
% [baseMVA, bus, gen, branch] = deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch);
% [ref, pv, pq] = bustypes(mpc.bus, mpc.gen);
% npv=length(pv);npq=length(pq);nref=length(ref);
% p_id=[pv;pq];
% p_size=round((npv+npq)/2);
% p_part=randi([1 npv+npq], p_size,1);
% Va0 = zeros(size(mpc.bus(:, VA), 1), ntest);
% for i = 1:ntest
%     Va0(:, i) = mpc.bus(:, VA);
%     Va0(p_id(p_part), i) = Va0(p_id(p_part), i)+(-0.005+0.01*randn(p_size,1))/pi*180;
% end
load("Va0_1-2_1000.mat")

nmethod = 12;
profile=cell(ntest,nmethod);
status=zeros(ntest,nmethod);
mpoptcnm.pf.tol=tol;
mpoptnr.pf.tol=tol;
mpoptrodas4.pf.tol=tol;
mpoptrodas3d.pf.tol=tol;
mpoptIWA.pf.tol=tol;
mpoptMANN.pf.tol=tol;
mpoptROMBERG.pf.tol=tol;
mpopt2S4.pf.tol=tol;
%%

for i = 398
    mpc.bus(:,VA)=Va0(:, i);
    fprintf('\n-★-★-★-★-★test %d-★-★-★-★-★-\n', i);
    %% newton
    mid = 1;
    fprintf("\n-----Newton-Raphson-----\n")
    try %nr method
        lastwarn('')
        pfnr=runpfa(mpc, mpoptnr);
        pfnr.stats.success=pfnr.success;
        profile{i,mid}=pfnr.stats;
        if pfnr.success
            status(i,mid)=1;
            fprintf("\nConverged within %i iterations\n", pfnr.iterations)
        end
    catch ME
        if ME.identifier=="MATLAB:nearlySingularMatrix"
            errid=1;
        elseif ME.identifier=="Matlab:TooSmallStep"
            errid=0;
        end
        profile{i,mid}=errid;
    end
    %% licnmpf-rodas4
    mid = 11;
    fprintf("\n-----SICNM-RODAS4-----\n")
    try
        lastwarn('')
        pflir=runpfa(mpc, mpoptrodas4);
        pflir.stats.success=pflir.success;
        profile{i,mid}=pflir.stats;
        if pflir.success
            status(i,mid)=1;
            fprintf("\nConverged within %i iterations\n", pflir.stats.ni)
        end
    catch ME
        if ME.identifier=="MATLAB:nearlySingularMatrix"
            errid=1;
        elseif ME.identifier=="Matlab:TooSmallStep"
            errid=0;
        end
        profile{i,mid}=errid;
    end
    %% licnmpf-rodas3d
    mid = 12;
    fprintf("\n-----SICNM-RODAS3D-----\n")
    try
        lastwarn('')
        pfliy=runpfa(mpc, mpoptrodas3d);
        pfliy.stats.success=pfliy.success;
        profile{i,mid}=pfliy.stats;
        if pfliy.success
            status(i,mid)=1;
            fprintf("\nConverged within %i iterations\n", pfliy.stats.ni)
        end
    catch ME
        if ME.identifier=="MATLAB:nearlySingularMatrix"
            errid=1;
        elseif ME.identifier=="Matlab:TooSmallStep"
            errid=0;
        end
        profile{i,mid}=errid;
    end
    %% ICNM-JH
    mid = 7;
    mpoptcnm.mips.sc.red_it=1;
    fprintf("\n-----ICNM-JH-----\n")
    try
        lastwarn('')
        pf_icnmjh=runpfa(mpc, mpoptcnm);
        pf_icnmjh.stats.success=pf_icnmjh.success;
        profile{i,mid}=pf_icnmjh.stats;
        if pf_icnmjh.success
            status(i,mid)=1;
            fprintf("\nConverged within %i iterations\n", pf_icnmjh.stats.ni)
        end
    catch ME
        if ME.identifier=="MATLAB:nearlySingularMatrix"
            errid=1;
        elseif ME.identifier=="Matlab:TooSmallStep"
            errid=0;
        end
        profile{i,mid}=errid;
    end
    %% ICNM-J
    mid = 8;
    mpoptcnm.mips.sc.red_it=2;
    fprintf("\n-----ICNM-J-----\n")
    try
        lastwarn('')
        pf_icnmj=runpfa(mpc, mpoptcnm);
        pf_icnmj.stats.success=pf_icnmj.success;
        profile{i,mid}=pf_icnmj.stats;
        if pf_icnmj.success
            status(i,mid)=1;
            fprintf("\nConverged within %i iterations\n", pf_icnmj.stats.ni)
        end
    catch ME
        if ME.identifier=="MATLAB:nearlySingularMatrix"
            errid=1;
        elseif ME.identifier=="Matlab:TooSmallStep"
            errid=0;
        end
        profile{i,mid}=errid;
    end
    %% ICNM-J1
    mid = 9;
    mpoptcnm.mips.sc.red_it=3;
    fprintf("\n-----ICNM-J1-----\n")
    try
        lastwarn('')
        pf_icnmj1=runpfa(mpc, mpoptcnm);
        pf_icnmj1.stats.success=pf_icnmj1.success;
        profile{i,mid}=pf_icnmj1.stats;
        if pf_icnmj1.success
            status(i,mid)=1;
            fprintf("\nConverged within %i iterations\n", pf_icnmj1.stats.ni)
        end
    catch ME
        if ME.identifier=="MATLAB:nearlySingularMatrix"
            errid=1;
        elseif ME.identifier=="Matlab:TooSmallStep"
            errid=0;
        end
        profile{i,mid}=errid;
    end
    %% ICNM-J0
    mid = 10;
    mpoptcnm.mips.sc.red_it=4;
    fprintf("\n-----ICNM-J0-----\n")
    try
        lastwarn('')
        pf_icnmj1=runpfa(mpc, mpoptcnm);
        pf_icnmj1.stats.success=pf_icnmj1.success;
        profile{i,mid}=pf_icnmj1.stats;
        if pf_icnmj1.success
            status(i,mid)=1;
            fprintf("\nConverged within %i iterations\n", pf_icnmj1.stats.ni)
        end
    catch ME
        if ME.identifier=="MATLAB:nearlySingularMatrix"
            errid=1;
        elseif ME.identifier=="Matlab:TooSmallStep"
            errid=0;
        end
        profile{i,mid}=errid;
    end
    %% CNM-ERK4
    mid = 3;
    mpoptcnm.mips.sc.red_it=5;
    fprintf("\n-----CNM-ERK4-----\n")
    try
        lastwarn('')
        pf_cnme=runpfa(mpc, mpoptcnm);
        pf_cnme.stats.success=pf_cnme.success;
        profile{i,mid}=pf_cnme.stats;
        if pf_cnme.success
            status(i,mid)=1;
            fprintf("\nConverged within %i iterations\n", pf_cnme.stats.ni)
        end
    catch ME
        if ME.identifier=="MATLAB:nearlySingularMatrix"
            errid=1;
        elseif ME.identifier=="Matlab:TooSmallStep"
            errid=0;
        end
        profile{i,mid}=errid;
    end
    %% IWAMOTO
    mid = 2;
    fprintf("\n-----IWAMOTO-----\n")
    try
        lastwarn('')
        pf_iwa=runpfa(mpc, mpoptIWA);
        pf_iwa.stats.success=pf_iwa.success;
        profile{i,mid}=pf_iwa.stats;
        if pf_iwa.success
            status(i,mid)=1;
            fprintf("\nConverged within %i iterations\n", pf_iwa.stats.ni)
        end
    catch ME
        if ME.identifier=="MATLAB:nearlySingularMatrix"
            errid=1;
        elseif ME.identifier=="Matlab:TooSmallStep"
            errid=0;
        end
        profile{i,mid}=errid;
    end
    %% ROMBERG
    mid = 4;
    fprintf("\n-----ROMBERG-----\n")
    try
        lastwarn('')
        pf_romb=runpfa(mpc, mpoptROMBERG);
        pf_romb.stats.success=pf_romb.success;
        profile{i,mid}=pf_romb.stats;
        if pf_romb.success
            status(i,mid)=1;
            fprintf("\nConverged within %i iterations\n", pf_romb.stats.ni)
        end
    catch ME
        if ME.identifier=="MATLAB:nearlySingularMatrix"
            errid=1;
        elseif ME.identifier=="Matlab:TooSmallStep"
            errid=0;
        end
        profile{i,mid}=errid;
    end
    %% 2S4
    mid = 5;
    fprintf("\n-----2S4-----\n")
    try
        lastwarn('')
        pf_2s4=runpfa(mpc, mpopt2S4);
        pf_2s4.stats.success=pf_2s4.success;
        profile{i,mid}=pf_2s4.stats;
        if pf_2s4.success
            status(i,mid)=1;
            fprintf("\nConverged within %i iterations\n", pf_2s4.stats.ni)
        end
    catch ME
        if ME.identifier=="MATLAB:nearlySingularMatrix"
            errid=1;
        elseif ME.identifier=="Matlab:TooSmallStep"
            errid=0;
        end
        profile{i,mid}=errid;
    end
    %% MANN
    mid = 6;
    fprintf("\n-----MANN-----\n")
    try
        lastwarn('')
        pf_mann=runpfa(mpc, mpoptMANN);
        pf_mann.stats.success=pf_mann.success;
        profile{i,mid}=pf_mann.stats;
        if pf_mann.success
            status(i,mid)=1;
            fprintf("\nConverged within %i iterations\n", pf_mann.stats.ni)
        end
    catch ME
        if ME.identifier=="MATLAB:nearlySingularMatrix"
            errid=1;
        elseif ME.identifier=="Matlab:TooSmallStep"
            errid=0;
        end
        profile{i,mid}=errid;
    end
end
