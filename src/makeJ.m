function J = makeJ(Ybus, V, Vm, Sbus, pv, pq)
    [dSbus_dVa, dSbus_dVm] = dSbus_dV(Ybus, V);
    [dummy, neg_dSd_dVm] = Sbus(Vm); % this is for zip load
    dSbus_dVm = dSbus_dVm - neg_dSd_dVm;

    j11 = real(dSbus_dVa([pv; pq], [pv; pq]));
    j12 = real(dSbus_dVm([pv; pq], pq));
    j21 = imag(dSbus_dVa(pq, [pv; pq]));
    j22 = imag(dSbus_dVm(pq, pq));

    J = [   j11 j12;
            j21 j22;    ];
end