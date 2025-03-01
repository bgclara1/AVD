function  [A, Ixx] = IFrameProps(currentT, currentH, currentL)
    t = currentT;
    b = currentL;
    d = currentH;
    h = d-2*t;
    A = b*d - h*(b-t);

    Ixx = (b*d^3-h^3*(b-t))/12;
    
end