function plotFinalSF(xDiscr, SF_3_75, SF_1_5, SF_OEI)
    
    figure;
    plot(xDiscr, SF_3_75)
    hold on;
    plot(xDiscr,SF_1_5)
    %hold on;
   % plot(xDiscr, SF_OEI)
    title('SF Plot')
    legend('n=3.75','n=-1.5')


end