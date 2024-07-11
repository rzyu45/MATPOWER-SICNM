function J_romberg = J_romberg(y, Ybus, Vm, Va, Sbus, pv, pq, j1, j2, j3, j4, j5, j6)
    Va(pv) = y(j1:j2);
    Va(pq) = y(j3:j4);
    Vm(pq) = y(j5:j6);
    V = Vm .* exp(1j * Va);

    [dSbus_dVa, dSbus_dVm] = dSbus_dV(Ybus, V);
    [dummy, neg_dSd_dVm] = Sbus(Vm); % this is for zip load
    dSbus_dVm = dSbus_dVm - neg_dSd_dVm;

    j11 = real(dSbus_dVa([pv; pq], [pv; pq]));
    j12 = real(dSbus_dVm([pv; pq], pq));
    j21 = imag(dSbus_dVa(pq, [pv; pq]));
    j22 = imag(dSbus_dVm(pq, pq));

    J_romberg = [   j11 j12;
            j21 j22;    ];
end