function Hz = makeHz(Ybus, V, pv, pq, z0, j1, j2, j3, j4, j5, j6)
    [Gaa_ap, Gva_vp, Gav_ap, Gvv_vp]=...
            makeGxx_xp(Ybus, V, pv, pq, z0, j1, j2, j3, j4, j5, j6);
    temp1 = Gaa_ap+Gva_vp;
    temp2 = Gav_ap+Gvv_vp;
    % Hz = H(y)\otimes z
    Hz11 = real(temp1([pv; pq], [pv; pq]));
    Hz12 = real(temp2([pv; pq], pq));
    Hz21 = imag(temp1(pq, [pv; pq]));
    Hz22 = imag(temp2(pq, pq));
    Hz = [ Hz11 Hz12;
           Hz21 Hz22];
end