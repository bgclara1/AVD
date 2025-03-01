function plotOmegaFrame(h, b, c, t)
    figure;
    hold on;
    

    x = [0,0.5*c,0.5*c,0.5*c+b-t,0.5*c+b-t,0.5*c-t,0.5*c-t,-1*(0.5*c-t),-1*(0.5*c-t),-1*(0.5*c+b-t),-1*(0.5*c+b-t),-1*0.5*c,-1*0.5*c,0];
    y = [0,0,h-t,h-t,h,h,t,t,h,h,h-t,h-t,0,0];
    

    x = [-c/2, -b/2, -b/2, b/2, b/2, c/2, c/2, b/2, b/2, -b/2, -b/2, -c/2, -c/2];
    y = [h+t, h+t, t, t, h+t, h+t, h, h, 0, 0, h, h, h+t];


    plot(x, y);

    xlabel('Width (m)');
    ylabel('Height (m)');
    title('Omega Section Geometry');

    grid on;
    hold off;
end









