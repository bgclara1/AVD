function plotFinalSF(xDiscr, SF_3_75, SF_1_5, SF_OEI, SFland)
    set(0,'defaultfigurecolor',[1 1 1])
set(groot,'defaultAxesFontSize',15)
set(groot,'defaulttextfontsize',15)
set(groot,'defaultLineMarkerSize',4)
set(groot,'defaultLineLineWidth',1.5)
set(groot,'defaultAxesXGrid','on')
set(groot,'defaultAxesYGrid','on')
set(groot, 'DefaultAxesBox', 'off')
set(groot, 'defaultAxesFontName','Cambria Math')

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