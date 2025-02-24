function frameAreaMatrix = frameArea(frameThicknessMatrix,bRange,hRange)


    for i = 1:length(bRange)
        for j = 1:length(hRange)
            frameAreaMatrix(i,j) = (2*bRange(i)+hRange(j))*frameThicknessMatrix(i,j);
        end
    end

end