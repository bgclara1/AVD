function SF_1_5 = getSF_1_5(xDiscr, aero15,total)


    SF_1_5 = total(1);
    for i = (2:length(total))
        SF_1_5(i) = SF_1_5(i-1) + total(i);
    end


end