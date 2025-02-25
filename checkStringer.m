function [stressCompliant,eulerCompliant] = checkStringer(directStress, tensileYieldStress,E, h, A_str, L, t_s)
    stressCompliant = all(abs(directStress) <= tensileYieldStress);

    a = L;
    b = h;
    c = L-(t_s);
    t = t_s;


    A_web = t_s * h;
    I_web_local = (t_s * h^3) / 12;  
    y_web_centroid = h/2;          


    A_flange = t_s * L;
    I_flange_local = (L * t_s^3) / 12;  
    y_flange_centroid = h + t_s/2;    

    A_total = A_web + A_flange;
    y_NA = (A_web*y_web_centroid + A_flange*y_flange_centroid) / A_total;

    d_web    = abs(y_web_centroid    - y_NA);
    d_flange = abs(y_flange_centroid - y_NA);

    I_web    = I_web_local    + A_web    * d_web^2;
    I_flange = I_flange_local + A_flange * d_flange^2;
    
    I_str = I_web + I_flange;

    sigma_euler = (pi^2 * E * I_str) / (A_str * L^2);
    sigma_max_compression = max(abs(directStress)) * 1e6;
    eulerCompliant = sigma_max_compression <= sigma_euler;

    %increase b and h
    %decrease stringer area
    %decrease L_eff
    %decrease direct stress

end
