function [Gaa_ap, Gva_vp, Gav_ap, Gvv_vp]=...
    makeGxx_xp(Ybus, V, pv, pq, z, j1, j2, j3, j4, j5, j6)
    %% initialize
    n=length(V);
    Va=angle(V);
    Ibus = Ybus * V;
    diagV = sparse(1:n, 1:n, V, n, n);
    diagI = sparse(1:n, 1:n, Ibus, n, n);
    X=sparse(1:n,1:n,exp(1j*Va),n,n);
    Y=sparse(1:n,1:n,exp(-1j*Va),n,n);
    A=conj(diagI)-conj(Ybus)*conj(diagV);
    %% LAMBDA
    Vap=zeros(n,1);%Va'
    Vmp=zeros(n,1);%Vp'
    Vap(pv)=z(j1:j2);%pv
    Vap(pq)=z(j3:j4);%pq
    Vmp(pq)=z(j5:j6);%pq

    %% derive Gaa_ap, Gav_ap
    diagvap = sparse(1:n, 1:n, Vap, n, n);
    Z=diagvap*conj(Ybus);
    W=conj(Ybus)*diagvap;
    B=Z-W;
    diagAlam=sparse(1:n,1:n,A*Vap,n,n);
    Gaa_ap = diagV*B*conj(diagV)-diagAlam*diagV;
    Gav_ap = 1j*diagV*B*Y+1j*diagAlam*X;
    
    %% derive Gva_vp, Gvv_vp
    diagvmp = sparse(1:n, 1:n, Vmp, n, n);
    Z=diagvmp*conj(Ybus);
    W=conj(Ybus)*diagvmp;
    C=sparse(1:n,1:n,conj(Ybus)*Y*Vmp,n,n);

    Gva_vp = 1j*(X*diagvmp*A-diagV*W*Y+C*diagV);
    Gvv_vp = X*Z*Y+C*X;
end


