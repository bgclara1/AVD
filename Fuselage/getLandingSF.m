function SFland = getLandingSF(xDiscr, total)

     SFland = total(1);
    for i = (2:length(total))
        %SF_3_75(i) = SF_3_75(i-1) + aero(i) + total(i);
       SFland(i) = SFland(i-1) + total(i);
    end

end



