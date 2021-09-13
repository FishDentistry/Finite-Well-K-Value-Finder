function [allowedE] = Asym_FiniteWell_Fxn(Lnm,UL,UR,state)
Lnm = Lnm * 10^-9;
UL = UL*(1.6*10^-19);
UR = UR*(1.6*10^-19);
m = 9.109*10^-31;
hBar = (6.626*10^-34)/(2*pi); %Define and convert input values to correct form
if UL > UR
    temp = UR;
elseif UR > UL %Find limiting side of well
    temp = UL;
end
maxE = temp; 
maxKval = sqrt((2*m*maxE)/(hBar)^2); %Find limiting kvalue
allowedKvals = [];
%disp(maxE/(1.6021*10^-19));
%disp(maxKval);
alphaL = @(k) sqrt(((2*m*UL)/(hBar^2))-k.^2); %Define functions for both alphas
alphaR = @(k) sqrt(((2*m*UR)/(hBar^2))-k.^2);
quantization = @(k) cot(k.*Lnm)+((alphaL(k)./alphaR(k)).*cot(k*Lnm)) - (k./alphaR(k)) + (alphaL(k)./k); %Create quantization function

kHold = 0;
count = 1;
count2 = 0;
min = 0.1;%Set beginning minimum and maximum
max = (pi/Lnm);
%disp(max)
numintervals = round(maxKval/(pi/Lnm)); %Find number of intervals to look over
for j=1:numintervals-1
    x0 = [min+0.1 max-0.1];
    allowedKvals(j) = fzero(quantization,x0); %Add kvalues
    min = max;
    max = max + (pi/Lnm);   
end
%disp(allowedKvals)
x0 = [min+0.1 maxKval-0.1];
allowedKvals(j+1) = fzero(quantization,x0); %Add final kvalue
K = allowedKvals(state); %Set kvalue determined by state parameter
%disp(K)
alphaL = alphaL(K);
alphaR = alphaR(K); %Set alphas based on chosen k
part1 = 1/(2*alphaL);
part2numerator = ((alphaL*sin(K*Lnm)/K)+cos(K*Lnm))^2;
part2 = (part2numerator)/(2*alphaR);
part3 = ((alphaL^2*Lnm)/(2*K^2))-((alphaL^2*sin(2*K*Lnm))/(4*K^3))+((alphaL*(sin(K*Lnm))^2)/(K^2))+(Lnm/2)+((sin(2*K*Lnm))/(4*K));
full = part1 + part2 + part3;
C= sqrt(1/(full)); %Find C from normalization conditions
A = (alphaL*C)/K;
B=C; 
G = ((alphaL*C/K)*sin(K*Lnm)+C*cos(K*Lnm))/exp(-alphaR*Lnm); %Define other coefficients from C
dx = 0.000000000001;
x=[0:dx:Lnm];
xl = [-Lnm*2.5:dx:0];
xr = [Lnm:dx:Lnm*2.5];
psimiddle = A*sin(K*x)+C*cos(K*x); %Create functions and bounds for x axis
psileft = C*exp(alphaL*xl);
psiright = G*exp(-alphaR*xr);


figure(1);
hold on
plot(x,psimiddle);
plot(xl,psileft);
plot(xr,psiright); %Plot spatial 
xline(0,'k--');
xline(Lnm,'k--');
title(['Asym Finite Well: L=',num2str((Lnm/10^-9)),' nm,',' UL/UR=',num2str(UL/(1.6*10^-19)),'/',num2str(UR/(1.6*10^-19)),' eV,',' n=',num2str(state)])
xlabel('Position,x(m)') 
ylabel('Spatial part of wave function, {\psi}(x)')
hold off

figure(2);
hold on
plot(x,psimiddle.^2);
plot(xl,psileft.^2);
plot(xr,psiright.^2); %Plot probability density
xline(0,'k--');
xline(Lnm,'k--');
title(['Asym Finite Well: L=',num2str((Lnm/10^-9)),' nm,',' UL/UR=',num2str(UL/(1.6*10^-19)),'/',num2str(UR/(1.6*10^-19)),' eV,',' n=',num2str(state)])
xlabel('Position,x(m)') 
ylabel('Probability density, |{\psi}(x)|^2')
hold off

allowedEnergies = []; %Calculate allowed energies
for c= 1:length(allowedKvals)
    tempvar = (allowedKvals(c)^2*hBar^2)/(2*m);
    allowedEnergies(c) = tempvar / (1.602*10^-19);
end
disp(["The ",num2str(length(allowedEnergies))," allowed energies ", "are: ", num2str(allowedEnergies)," eV"])
riemann = sum([psimiddle.^2,psiright.^2,psileft.^2]*dx); %Calculate riemann sum
disp(["The Riemann sum is:", num2str(riemann)])
end



