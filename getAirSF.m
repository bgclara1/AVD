function getAirSF(xDiscr,aero)

    aeroSF = 0;

    for i = (2:length(aero))
        aeroSF(i) = aeroSF(i-1) + aero(i);
    end

    figure
    plot(xDiscr,aeroSF)
    title('Aero Only SF')


end