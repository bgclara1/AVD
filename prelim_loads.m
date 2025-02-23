close all;

% ======== set default params for plotting ================

close all;

set(0,'defaultfigurecolor',[1 1 1])
set(groot,'defaultAxesFontSize',15)
set(groot,'defaulttextfontsize',15)
set(groot,'defaultLineMarkerSize',4)
set(groot,'defaultLineLineWidth',3)
set(groot,'defaultAxesXGrid','on')
set(groot,'defaultAxesYGrid','on')
set(groot, 'DefaultAxesBox', 'off')
set(groot, 'defaultAxesFontName','Cambria Math')

% ======== get inertial distrobution ====================

[xDiscr,inertialDistro] = getInertialDistro();

% ======== get air load plot ============================

[xDiscr,aero,aero15] = getAirLoad();

% ======== get air SF plot ==============================

getAirSF(xDiscr,aero)

% ======= get fusealge SF plot ==========================

total = getFuselageSF(xDiscr,inertialDistro);

            % you get a bump at 31 that looks diff to reports but our aero
            % load is more disproportional so when u add it on it makes the
            % thing shoot above zero


SF_3_75 = getSF_3_75(xDiscr,aero,total)    %combines the air and fuselage SF at 3.75 loading

SF_1_5 = getSF_1_5(xDiscr, aero15,total)

SF_OEI = getOEI(xDiscr,total)

plotFinalSF(xDiscr, SF_3_75, SF_1_5, SF_OEI)



%getLandingSF(xDiscr, total)









