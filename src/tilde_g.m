function dtilde_g = tilde_g(y, z, pv, pq, Vm, Va, Ybus, Sbus, j1, j2, j3, j4, j5, j6)
    Va(pv) = y(j1:j2);
    Va(pq) = y(j3:j4);
    Vm(pq) = y(j5:j6);
    V = Vm .* exp(1j * Va);
    mis = V .* conj(Ybus * V) - Sbus(Vm);
    G = [   real(mis([pv; pq]));
            imag(mis(pq))   ];
    J = makeJ(Ybus, V, Vm, Sbus, pv, pq);
    dtilde_g = J*z+G;
end