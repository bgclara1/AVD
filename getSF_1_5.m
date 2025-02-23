function SF_1_5 = getSF_1_5(xDiscr, aero15,total)

    total = total*-1;
    SF_1_5 = total(1);
    for i = (2:length(total))
        SF_1_5(i) = SF_1_5(i-1) + aero15(i) + total(i);
    end

    figure;
    plot(xDiscr, SF_1_5)
    title('Fuselage + Aero SF at load factor -1.5')


end