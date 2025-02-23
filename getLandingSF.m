function getLandingSF(xDiscr, total)

        
    lTot = 79;
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


    MaxLandingMassFrac = 0.85;
    MoW15 = 2.18*Cm*(0.5*rhoC*Va^2*SrefWing*MAC);
    LHTland = (MaxLandingMassFrac*-330800*9.81*(38.8831-30.86)+MoW15)/(71.1-30.86);
    
    A = [ 31	39 ;
            1	1];
    B = [71*LHTland*-1; LHTland*-1];
    
    Xland = linsolve(A,B);
    forceRFland = Xland(1);
    forceRRland = Xland(2);
    
    aeroland = zeros(lTot+1);
    aeroland(39) = forceRRland; % back wing spar
    aeroland(31) = forceRFland; %front wing spar
    aeroland(71) = LHTland; % HT ac
    aeroland = aeroland(:,1);
    
    
    
    SFaeroland = [];
    SFaeroland(1) = aeroland(1);
    for i = 2:length(aeroland)
        SFaeroland(i) = SFaeroland(i-1) + aeroland(i);
    end
    
    SFland = [];
    SFland(1) = SFaeroland(1) + total(1);
    for i = 2:length(SFaeroland)
        SFland(i) = SFaeroland(i) + total(i);
    end


    figure;
    plot(xDiscr,SFland)
end