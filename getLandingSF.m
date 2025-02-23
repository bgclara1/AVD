function SFland = getLandingSF(xDiscr, total)

        
    lTot = 79;
    SrefWing = 469.44;
    
    wingWidth = 13.6; %chord
    WingLE = 28.4779;
    
    Va = 92.4;
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
    
    aeroland = zeros(1,79);
    aeroland(39) = forceRRland; % back wing spar
    aeroland(31) = forceRFland; %front wing spar
    aeroland(71) = LHTland; % HT ac
    
    
    
    SFland = [];
    SFland(1) = total(1);
    for i = 2:length(aeroland)
        SFland(i) = SFland(i-1) + aeroland(i) + total(i);
    end

    figure;
    plot(xDiscr,SFland)
end