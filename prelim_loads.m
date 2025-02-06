clmax = 1.5;
MTOW = 729289.163; %lbs
%nmax = clmax/cd0; % it might also just be 2.5 from regulations - max possible load for any maneouvre so from dive to pull up inst
nlim = 2*1+(24000/(MTOW+1000)); %mtow lbs gives 2.0329 but cs25 says no less than 2.5 so we'll take 2.5
nlim = 2.5;
vlim = 0.83*303.1; %speed of sound at 30,000ft
rho0 = 1.227;
veaslim = sqrt(0.4583/rho0)*vlim;
sref = 469.44;
WoverS = MTOW/sref; %W/s
n=[];
veas=[];
x=[];
a = 1;

for i = -1:0.1:2.5
    n(a) = i;
    a = a + 1;
end

a = 1;
for j = 0:veaslim/length(n):veaslim
    veas(a) = j;
    a = a+1;
end 
for k = 1:length(n)
    x(k) = ((rho0/2)*veas(k)^2*clmax/WoverS);
end

figure;
plot(n,x)
xline(veaslim)










