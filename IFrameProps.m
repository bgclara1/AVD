function [A, Ixx] = IFrameProps(b, h, t)
    % I-beam cross-sectional area
    A = 2*(b*t) + h*(b-t);

    % Moment of Inertia (Ixx)
    I_flange = (b * t^3) / 12;  
    A_flange = b * t;  
    d_flange = (h/2) + (t/2);  
    I_flange_total = 2 * (I_flange + A_flange * d_flange^2);  

    I_web = ((b - t) * h^3) / 12;  

    Ixx = I_flange_total + I_web;
end
