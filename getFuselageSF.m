function getFuselageSF(xDiscr,inertialDistro)
    
    inertialDistro = inertialDistro/9.81;

    totalLoad = sum(inertialDistro);

    for i = 1:length(xDiscr)
        inertialMoments(i) = inertialDistro(i) * xDiscr(i);
    end

    totalMoment = sum(inertialMoments);

    fuselageCG = totalMoment/totalLoad;

    A = [ 31	39 ;
            1	1];
    B = [totalMoment; totalLoad];
    
    
    X = linsolve(A,B);
    forceRF = X(1);
    forceRR = X(2);

    momentRF = forceRF*31;
    momentRR = forceRR*39;

    inertialDistro = inertialDistro*9.81;
    
    sparReaction = zeros(1,79);
    sparReaction(31) = forceRF*9.81*-1;
    sparReaction(39) = forceRR*9.81*-1;

    for i = 1:length(xDiscr)
        total(i) = inertialDistro(i) + sparReaction(i);
    end
    
    fuselageSF = total(1);
    for i = (2:length(total))
        fuselageSF(i) = fuselageSF(i-1) + total(i);
    end

    figure;
    plot(xDiscr, fuselageSF)


end