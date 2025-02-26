function SF_OEI = getOEI(xDiscr, total) 

    SF_OEI = total(1);
    for i = (2:length(total))
       SF_OEI(i) = SF_OEI(i-1) + total(i);
    end


end