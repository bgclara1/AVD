function frameThicknessMatrix = frameThickness(I_xx,bRange,hRange)


    for i = 1:length(bRange)
        for j = 1:length(hRange)
            frameThicknessMatrix(i,j) = I_xx/(bRange(i)*hRange(j)^2/2+hRange(j)^3/12);
        end
    end

end