function SF_3_75 = getSF_3_75(total)

    SF_3_75 = total(1);
    for i = (2:length(total))
        %SF_3_75(i) = SF_3_75(i-1) + aero(i) + total(i);
       SF_3_75(i) = SF_3_75(i-1) + total(i);
    end

    % figure;
    % plot(xDiscr, SF_3_75)
    % title('Fuselage + Aero SF at load factor 3.75')


end