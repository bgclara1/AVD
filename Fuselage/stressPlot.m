function directStress = stressPlot(skinThickness,StringerArea,y,numStringers,density)

    r = 3.195;   
    M = 1e+07;
    C = 2 * pi * r;

  %  StringerArea = 0.005;
    equivBoomArea = 15 * skinThickness^2;
    SingleBoomArea = StringerArea + equivBoomArea;

    areaxarm = SingleBoomArea .* y;
    Iarray = SingleBoomArea .* y.^2;

    I = sum(Iarray);
    c = max(abs(y));
    
    sigma = (M*c)/(I*10^6);
    
    A = StringerArea*numStringers+C*skinThickness;
    
   % density = 2.78e+03; %change
 %   FuselageWeightPUL = density*A;
    
    sigma = (M * c) ./ (I * 1e6);
    directStress = sigma .* y / c; 
    
    % figure;
    % hold on;
    % grid on;
    % plot3(x, y, zeros(size(x)), 'r-', 'LineWidth', 2);
    % quiver3(x, y, zeros(size(x)), zeros(size(x)), zeros(size(x)), directStress, 'b', 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
    % xlabel('z');
    % ylabel('y');
    % zlabel('Direct stress \sigma_z (MPa)');
    % title('Direct stress distribution around fuselage');
    % view(3);
    % legend('Fuselage cross-section', 'Direct stress at each stringer');
    % axis equal;

end