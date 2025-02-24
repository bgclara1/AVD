%% Z-Stringer Area Calculation
function area = ZStringerArea(t_s, h, L)
 %   area = t_s * (2*L + h);
    area = L*t_s*2+(h-2*t_s)*t_s;
end
