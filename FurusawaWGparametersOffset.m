%Furusawa waveguide, variance vs phase with offset. 

dB_to_pow = @(d) 10.^(d./10);
pow_to_dB = @(p) 10.*log10(p);

%Waveguide parameters taken from "Over-8-dB squeezed light generation by a broadband
%waveguide optical parametric amplifier toward fault-tolerant
%ultra-fast quantum computers"
P_in = 660e-3; %Input power
L = 45e-3;
eta_norm = 8.2./L.^2; %1/W*m^2
eta_tot = 0.88; 
rad = 14; 

%Sellmeier coefficients
a1 = 4.9048;
a2 = 0.11775;
a3 = 0.21802;
a4 = 0.027153;
b1 = 2.2314e-6;
b2 = -2.9671e-8;
b3 = 2.1429e-8;

B = sqrt(P_in).*tanh(L.*sqrt(P_in.*eta_norm)); %Output power square root
nB = sqrt(eta_norm).*L.*B; %sqrt(Output power) times sqrt(eta_norm)*L

%Determines what values to plot between
lim = pi/2;
phi = linspace(-lim,lim,1000); 

T = 31; %Celsius
T_0 = 30; %Celsius
dT = T-T_0;
f = (T-T_0).*(T+T_0+546);
dfdT = 2.*T + 546;

%Note that n assume lambda to be in microns
n = @(lambda)  sqrt(a1 + b3.*f + (a2+b1.*f)./(lambda.^2 - (a3 + b2.*f).^2) - a4.*lambda.^2);
dndT = @(lambda) dfdT.*(2.*b2.*(a2+b1.*f).*(a3+b2.*f)./(lambda.^2-(a3+b2.*f).^2).^2 + b1./(lambda.^2-(a3+b2.*f).^2) +b3)./(2*n(lambda));

alfa = 7.5e-6; %Thermal expansion coefficient
lambda_0 = 1550e-9; %Input wavelength
lambda_SHG = 775e-9; %Wavelength of generated SHG light

delta_k = (T-T_0).*pi.*4./lambda_0.*abs((dndT(lambda_SHG*1e3) - dndT(lambda_0*1e3)) + alfa.*(n(lambda_SHG) - n(lambda_0))); %Phase mismatch

pm = delta_k.*L;

%Squeezing
V_asqz = @(pm) eta_tot.*abs(sin(pm).*sinh(nB)).^2 + abs(cos(pm).*sinh(nB) + cosh(nB)).^2 + (1-eta_tot);
V_sqz = @(pm) eta_tot.*abs(sin(pm).*sinh(nB)).^2 + abs(cos(pm).*sinh(nB) - cosh(nB)).^2 + (1-eta_tot);

V_no_offset = V_sqz(phi).*cos(rad./1000).^2 + V_asqz(phi).*sin(rad./1000).^2;
V_offset = V_sqz(phi+pm).*cos(rad./1000).^2 + V_asqz(phi+pm).*sin(rad./1000).^2;

figure(1) %Creates the first plot

hold on
plot(phi,pow_to_dB(V_no_offset), 'Displayname', '\Delta T = 0')
plot(phi,pow_to_dB(V_offset), 'Displayname', strcat('\Delta T = ',num2str(dT.*1e3),' mK'))
plot(phi,phi.*0,'k','Displayname', 'Shot noise')

lgd = legend;
lgd.FontSize = 10;
lgd.Location = 'best';

ax=gca;
ax.FontSize = 12;
set(ax, 'xlim', [-lim lim]);
set(ax,'XTick',-lim:lim/2:lim) 
set(ax,'XTickLabel',{strcat(num2str(-lim./pi),'\pi'),strcat(num2str(-lim/(2*pi)),'\pi'), '0',strcat(num2str(lim/(2*pi)),'\pi'), strcat(num2str(lim/pi),'\pi')})

xlabel('SH phase, \phi','FontSize',16);
ylabel('Variance (dB)','Fontsize',16);
title('Variance vs phase with temperature offset','Fontsize',14);
grid
%% Power fluctuations

dB_to_pow = @(d) 10.^(d./10);
pow_to_dB = @(p) 10.*log10(p);

%Waveguide parameters taken from "Over-8-dB squeezed light generation by a broadband
%waveguide optical parametric amplifier toward fault-tolerant
%ultra-fast quantum computers"
P_in = 660e-3; %Input power
L = 45e-3;
eta_norm = 8.2./L.^2; %1/W*m^2
eta_tot = 0.88; 
rad = 14; 

%Sellmeier coefficients
a1 = 4.9048;
a2 = 0.11775;
a3 = 0.21802;
a4 = 0.027153;
b1 = 2.2314e-6;
b2 = -2.9671e-8;
b3 = 2.1429e-8;

alfa = 7.5e-6; %Thermal expansion coefficient
lambda_0 = 1550e-9; %Input wavelength
lambda_SHG = 775e-9; %Wavelength of generated SHG light

B = sqrt(P_in).*tanh(L.*sqrt(P_in.*eta_norm)); %Output power square root
nB = sqrt(eta_norm).*L.*B; %sqrt(Output power) times sqrt(eta_norm)*L

delta_P = 0.1; %Input fluctuations in delta_P_in/P_in
P_in_eff = (1+delta_P).*P_in;

tau = 9.533e-4; %Relaxation time (s), Assuming a nonlinear interaction radius of 25 microns
V = 1.296e-11; %Mode volume (m^3)
C = 633; %Specific heat capacity, J/(kg*k)
rho = 4648; %Density, kg/m^3
alfa_LiNbO = 0.01e2; %Absorption coefficient (1/m)

OMEGA = 0; %Sideband frequency 
s = 1i.*OMEGA; 

B_eff = sqrt(P_in_eff).*tanh(L.*sqrt(P_in_eff.*eta_norm)); %Output power square root
nB_eff = sqrt(eta_norm).*L.*B_eff; %sqrt(Output power) times sqrt(eta_norm)*L

delta_b = (B_eff-B)./B; %Re(delta_b/B), assuming same phase
P_abs = alfa_LiNbO.*L.*abs(B_eff).^2;% Absorbed power

dT = 1./(1+s.*tau).*2.*tau.*P_abs./(C.*rho.*V).*delta_b; %Temperature fluctuation

T_0 = 30;%Celsius
f = dT.*(dT + 2.*T_0 +546);
dfdT = 2.*(T_0+dT) + 546;

%Note that n assumes lambda to be in microns
n = @(lambda)  sqrt(a1 + b3.*f + (a2+b1.*f)./(lambda.^2 - (a3 + b2.*f).^2) - a4.*lambda.^2);
dndT = @(lambda) dfdT.*(2.*b2.*(a2+b1.*f).*(a3+b2.*f)./(lambda.^2-(a3+b2.*f).^2).^2 + b1./(lambda.^2-(a3+b2.*f).^2) +b3)./(2*n(lambda));

%Determines what values to plot between
lim = pi/4;
phi = linspace(-lim,lim,1000); 

delta_k = dT.*pi.*4./lambda_0.*abs((dndT(lambda_SHG*1e3) - dndT(lambda_0*1e3)) + alfa.*(n(lambda_SHG) - n(lambda_0))); %Phase mismatch

pm = delta_k.*L;

%Squeezing
V_asqz = @(pm,nB) eta_tot.*abs(sin(pm).*sinh(nB)).^2 + abs(cos(pm).*sinh(nB) + cosh(nB)).^2 + (1-eta_tot);
V_sqz = @(pm,nB) eta_tot.*abs(sin(pm).*sinh(nB)).^2 + abs(cos(pm).*sinh(nB) - cosh(nB)).^2 + (1-eta_tot);

V_no_offset = V_sqz(phi,nB).*cos(rad./1000).^2 + V_asqz(phi,nB).*sin(rad./1000).^2;
V_offset = V_sqz(phi+pm,nB_eff).*cos(rad./1000).^2 + V_asqz(phi+pm,nB_eff).*sin(rad./1000).^2;

figure(2) %Creates the second plot

hold on
plot(phi,pow_to_dB(V_no_offset), 'Displayname', '\Delta P_{in} = 0')
plot(phi,pow_to_dB(V_offset), 'Displayname', strcat('\Delta P_{in}/P_{in} = ',num2str(delta_P)))
plot(phi,phi.*0,'k','Displayname', 'Shot noise')

lgd = legend;
lgd.FontSize = 10;
lgd.Location = 'best';

ax=gca;
ax.FontSize = 12;
set(ax, 'xlim', [-lim lim]);
set(ax,'XTick',-lim:lim/2:lim) 
set(ax,'XTickLabel',{strcat(num2str(-lim./pi),'\pi'),strcat(num2str(-lim/(2*pi)),'\pi'), '0',strcat(num2str(lim/(2*pi)),'\pi'), strcat(num2str(lim/pi),'\pi')})

xlabel('SH phase, \phi','FontSize',16);
ylabel('Variance (dB)','Fontsize',16);
title('Variance vs phase with amplitude offset','Fontsize',14);
grid