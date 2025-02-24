function [Np, Sp, Mp, Nq, Sq, Mq, Nt, St, Mt, Nf, Sf, Mf] = calculateFrameLoads(Pt, Qt, Tt, R, phi)


    Np = (Pt ./ (2*pi)) .* ((pi - phi).*cos(phi) - 0.5*sin(phi));
    Sp = (Pt ./ (2*pi)) .* (1 + 0.5*cos(phi) - (pi - phi).*sin(phi));
    Mp = (Pt .* R ./ (2*pi)) .* ((pi - phi).*(1 - cos(phi)) - (3/2)*sin(phi));

    Nq = (Qt ./ (2*pi)) .* ((3/2)*cos(phi) + (pi - phi).*sin(phi));
    Sq = (Qt ./ (2*pi)) .* ((pi - phi).*cos(phi) - 0.5*sin(phi));
    Mq = (Qt .* R ./ (2*pi)) .* (1 + 0.5*cos(phi) - (pi - phi).*sin(phi));

    Nt = (Tt ./ (pi .* R)) .* sin(phi);
    St = (Tt ./ (2*pi .* R)) .* (1 + 2*cos(phi));
    Mt = (Tt ./ (2*pi)) .* ((pi - phi) - 2*sin(phi));

    Nf = Np + Nq + Nt;
    Sf = Sp + Sq + St;
    Mf = Mp + Mq + Mt;
    
end
