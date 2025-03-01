function results = OptFunctionManufacture(t_s_vec, b1_vec, h_vec, M_pos, sf375, T375, c_root, c_tip, y_wing)

    E = 72e9;               % Young's Modulus (Pa)
    rho = 2710;             % Density of Al (kg/m^3) 
    M = M_pos;              % Bending moment (N/m)
    c = linspace(c_root, c_tip, length(y_wing));  % Wing chord (m)
    b2 = 0.0758 .* c;       % Wing box height (m)
    width = 0.6 .* c;       % Wing box width (m)
    
    tau_tresca = 175e6;                          % TRESCA shear yield stress (Pa)
    Limit_sigma_compressive = 400e6;             % Compressive yield strength (Pa)
    
    results = struct();
    
    index = 1;
    for t_s = t_s_vec
        for b1 = b1_vec
            for h = h_vec
                
                n = width ./ b1;              % Number of stringers
                d_to_h = 0.3;                 % Width-to-height ratio
                d = h * d_to_h;               % Stringer width (m)
                A_s = t_s * h + 2 * t_s * d;  % Z-stringer area (m^2)
                
                t_stringers = ((M .* b1 .^ 2) ./ (3.62 .* E .* width .* b2)).^(1/3);
                N = M ./ b2 ./ width;   
                sigma = N ./ t_stringers;  % Critical global buckling stress (Pa)
                
                t_eff = t_stringers + A_s ./ b1;  % Effective skin thickness (m)
                t_eff = round(t_eff, 3);  % Round to nearest mm
                
                % Catchpole correction factor
                thickness_ratio = t_s ./ t_stringers;
                area_ratio = A_s ./ (b1 .* t_stringers);
                [stress_ratio, ~] = catchpole_interpolation(thickness_ratio, area_ratio);
                sigma_crit = stress_ratio .* sigma;
                [efficiency, ~] = Farrar_interpolation(thickness_ratio, area_ratio);
                
                rib_spacing_skin = ((efficiency ./ sigma_crit) .^ 2) .* N .* E;
                
                % Spar design
                S = sf375;
                T = -T375;
                q_torque = T ./ (2 .* width .* b2);
                q_shear = -S ./ (2 .* b2);
                q_fs = abs(q_torque + q_shear);
                q_rs = abs(q_shear - q_torque);
                
                Ks = 13.4;
                t_front = (abs((q_fs .* b2) ./ (Ks * E))).^(1/3);
                t_rear = (abs((q_rs .* b2) ./ (Ks * E))).^(1/3);
                tau_front = q_fs ./ t_front;
                tau_rear = q_rs ./ t_rear;
                
                tau_skin = q_torque ./ t_stringers;
                tau_skin(end) = tau_skin(end - 1);
                
                combination = (tau_skin ./ tau_tresca) .^ 2 + sigma ./ sigma_crit;
                
                % Rib Design
                I_skin_panel = (width .* (t_eff) .^ 3) ./ 12 + width .* t_eff .* (b2 ./ 2) .^ 2;
                Crush_load = 0.5 .* ((M .^ 2) .* rib_spacing_skin .* b2 .* (t_eff .* width)) ./ (E .* I_skin_panel .^ 2);
                t_rib = ((Crush_load .* b2 .^ 2) ./ (3.62 .* E .* width)).^(1/3);
                t_rib = round(t_rib, 3);  % Round to nearest mm
                
                sigma_rib = Crush_load ./ (t_rib .* width);
                
                obj = t_eff + (b2 .* t_rib) ./ rib_spacing_skin;
                
                results(index).t_s = t_s;
                results(index).b1 = b1;
                results(index).h = h;
                results(index).t_eff = t_eff;  % Store rounded thickness
                results(index).t_rib = t_rib;  % Store rounded rib thickness
                results(index).obj = obj;
                results(index).max_tau = max(tau_skin);
                results(index).max_sigma_crit = max(sigma_crit);
                results(index).max_combination = max(combination);
                results(index).min_efficiency = min(efficiency);
                
                index = index + 1;
            end
        end
    end
end
