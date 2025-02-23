function plotFinalSF(xDiscr, SF_3_75, SF_1_5, SF_OEI, SFland)
    
    figure;
    plot(xDiscr, SF_3_75)
    %hold on;
    %plot(xDiscr,SF_1_5)
    hold on;
    plot(xDiscr, SF_OEI)
    hold on;
    plot(xDiscr, SFland)
    title('SF Plot')
    legend('n=3.75','OEI', 'landing')


end