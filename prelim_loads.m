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

[xDiscr,aero] = getAirLoad();

% ======== get air SF plot ==============================

getAirSF(xDiscr,aero)

% ======= get fusealge SF plot ==========================

getFuselageSF(xDiscr,inertialDistro)

            % you get a bump at 31 that looks diff to reports but our aero
            % load is more disproportional so when u add it on it makes the
            % thing shoot above zero

            





