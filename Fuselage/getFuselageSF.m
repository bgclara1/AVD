function [fuselageSF,total] = getFuselageSF(xDiscr,inertialDistro)
    
    inertialDistro = inertialDistro; 

    totalLoad = sum(inertialDistro);

    for i = 1:length(xDiscr)
        inertialMoments(i) = inertialDistro(i) * xDiscr(i);
    end

    totalMoment = sum(inertialMoments);

%    fuselageCG = totalMoment/totalLoad;

    A = [ 31	39 ;
            1	1];
    B = [totalMoment; totalLoad];
    
    
    X = linsolve(A,B);
    forceRF = X(1);
    forceRR = X(2);
    
    sparReaction = zeros(1,80);
    sparReaction(31) = forceRF*-1;
    sparReaction(39) = forceRR*-1;

    for i = 1:length(xDiscr)
        total(i) = inertialDistro(i) + sparReaction(i);
    end
    
    fuselageSF = total(1);
    for i = (2:length(total))
        fuselageSF(i) = fuselageSF(i-1) + total(i);
    end

    % figure;
    % plot(xDiscr, fuselageSF)
    % title('Fuselage Only SF')


end