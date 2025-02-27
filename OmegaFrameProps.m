function [A, Ixx] = OmegaFrameProps(h, b, c, t)
    A_web = t * h;            
    A_long_flange = t * c;    
    A_short_flange = t * b;   

    y_web = h / 2;             
    y_long_flange = 0;          
    y_short_flange = h;       

    A = 2 * A_web + A_long_flange + 2 * A_short_flange;
    y_centroid = (2 * A_web * y_web + A_long_flange * y_long_flange + 2 * A_short_flange * y_short_flange) / A;

    I_web = (t * h^3) / 12;        
    I_long_flange = (c * t^3) / 12;   
    I_short_flange = (b * t^3) / 12;  

    % Step 5: Apply the parallel axis theorem for each rectangle
    d_web = y_web - y_centroid;                     % Distance for vertical webs
    d_long_flange = y_long_flange - y_centroid;     % Distance for long horizontal flange
    d_short_flange = y_short_flange - y_centroid;   % Distance for short horizontal flanges

    Ixx_web = 2 * (I_web + A_web * d_web^2);                % Contribution from vertical webs
    Ixx_long_flange = I_long_flange + A_long_flange * d_long_flange^2; % Contribution from long horizontal flange
    Ixx_short_flange = 2 * (I_short_flange + A_short_flange * d_short_flange^2); % Contribution from short horizontal flanges

    Ixx = Ixx_web + Ixx_long_flange + Ixx_short_flange;




    