function plotFinalSF(xDiscr, SF_3_75, SF_1_5, SF_OEI, SFland)
    
    figure;
    stairs(xDiscr, SF_3_75)
    %hold on;
    %plot(xDiscr,SF_1_5)
    hold on;
    stairs(xDiscr, SF_OEI)
    hold on;
    stairs(xDiscr, SFland)
    title('SF Plot')
    legend('n=3.75','OEI', 'landing')


end