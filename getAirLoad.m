
function [xDiscr,aero,aero15] = getAirLoad()

    lTot = 80;

    for i = 1:lTot
        xDiscr(i) = i-1;
    end

    Sref = 115.1508; %HT
    SrefWing = 469.44;
    frontSparWP = 0.2;
    backSparWP = 0.8;
    frontSparTP = 0.2;
    chordHT = 8.25;
    
    wingWidth = 13.6; %chord
    WingLE = 28.4779;
    xPosFSW = WingLE+0.2*wingWidth;
    xPosBSW = WingLE + 0.8*wingWidth;
    xPosFST = 74.86; % Nadir
    
    Vd = 192.4;
    Va = 92.4;
    
    rho = 1.225;
    rhoC = 0.4592; % C cruise
    
    MAC = 8.75;
    Cm = -0.14;
    
    %------------------- 3.75 --------------------------------
    MoW = 3.75*Cm*(0.5*rhoC*Vd^2*SrefWing*MAC);
    LHT = (-330800*9.81*(38.8831-30.86)+MoW)/(71.1-30.86);
    
    A = [ xPosFSW	xPosBSW ;
            1	1];
    B = [xPosFST*LHT*-1; LHT*-1];
    
    
    X = linsolve(A,B);
    forceRF = X(1);
    forceRR = X(2);
    
    
    momentRF = forceRF*xPosFSW;
    momentRR = forceRR*xPosBSW;

    clear aero
    aero = zeros(1,80);
    aero(31) = forceRF;
    aero(39) = forceRR;
    aero(71) = LHT;

    %------------ - 1.5 -------------------------------------

    
    MoW15 = -1.5*Cm*(0.5*rhoC*Vd^2*SrefWing*MAC);
    LHT15 = (-330800*9.81*(38.8831-30.86)+MoW15)/(71.1-30.86);
    
    A = [ xPosFSW	xPosBSW ;
            1	1];
    B = [xPosFST*LHT15; LHT15];
    
    X15 = linsolve(A,B);
    forceRF15 = X15(1);
    forceRR15 = X15(2);
    
    aero15 = zeros(1,80);
    aero15(31) = forceRF15;
    aero15(39) = forceRR15;
    aero15(71) = LHT*-1;
    % 
    % figure;
    % bar(xDiscr,aero);
    % hold on;
    % bar(xDiscr,aero15);


end




