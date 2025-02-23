function eulerCompliant = checkEulerBuckling(directStress, E, b, h, A_str, L_eff)

    I_str = (b * h^3) / 12;
    sigma_euler = (pi^2 * E * I_str) / (A_str * L_eff^2);
    sigma_max_compression = max(abs(directStress)) * 1e6;
    eulerCompliant = sigma_max_compression <= sigma_euler;

end
