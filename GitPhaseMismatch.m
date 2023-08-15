%WARNING: RUN EACH SECTION SEPERATELY FOR CLEARER RESULTS
%PhaseMismatch for GitHub

%% Squeezing vs temperature fluctuations

%Phase mismatch causing phase noise fluctuation
dB_to_pow = @(d) 10.^(d./10);
pow_to_dB = @(p) 10.*log10(p);

P_in = 660e-3; %Input power
eta_norm = 0.4e4; %1/W*m^2
L = 45e-3;

%Sellmeier coefficients
a1 = 4.9048;
a2 = 0.11775;
a3 = 0.21802;
a4 = 0.027153;
b1 = 2.2314e-6;
b2 = -2.9671e-8;
b3 = 2.1429e-8;


T = linspace(28,32,1000); %Celsius
T_0 = 30; %Celsius
f = (T-T_0).*(T+T_0+546);
dfdT = 2.*T + 546; %Derivative of f

%Formula for refractive index
n = @(lambda)  sqrt(a1 + b3.*f + (a2+b1.*f)./(lambda.^2 - (a3 + b2.*f).^2) - a4.*lambda.^2);

%Deriative w.r.t T of the above formula
dndT = @(lambda) dfdT.*(2.*b2.*(a2+b1.*f).*(a3+b2.*f)./(lambda.^2-(a3+b2.*f).^2).^2 + b1./(lambda.^2-(a3+b2.*f).^2) +b3)./(2*n(lambda));

alfa = 7.5e-6; %Thermal expansion coefficient
lambda_0 = 1550e-9; %Input wavelength
lambda_SHG = 775e-9; %Wavelength of generated SHG light

%Delta k, note that the formula for n inputs the wavelength in microns
delta_k = (T-T_0).*pi.*4./lambda_0.*abs((dndT(lambda_SHG*1e3) - dndT(lambda_0*1e3)) + alfa.*(n(lambda_SHG) - n(lambda_0))); %Phase mismatch

pm = delta_k.*L; %Total phase mismatch

B = sqrt(P_in).*tanh(L.*sqrt(P_in.*eta_norm)); %Output power square root
nB = sqrt(eta_norm).*L.*B; %sqrt(Output power) times sqrt(eta_norm)*L

%Analytical solution, several assumptions made
V_an_X = (sin(pm).*sinh(nB)).^2 + (cos(pm).*sinh(nB) + cosh(nB)).^2;
V_an_Y = (sin(pm).*sinh(nB)).^2 + (cos(pm).*sinh(nB) - cosh(nB)).^2;

hold on
plot(T,pow_to_dB(V_an_Y), 'Displayname', 'Squeezing')
plot(T,T.*0,'k','Displayname', 'Shot noise')

lgd = legend;
lgd.FontSize = 10;
lgd.Location = 'best';

xlabel('Temperature ( ^oC)','FontSize',16);
ylabel('Variance (dB)','Fontsize',16);

%% Squeezing vs amplitude fluctuations

dB_to_pow = @(d) 10.^(d./10);
pow_to_dB = @(p) 10.*log10(p);

P_in = 660e-3; %Input power
eta_norm = 0.4e4; %1/W*m^2
L = 45e-3;

tau = 9.5e-7; %Relaxation time (s)
V = 2.4e-14; %Mode volume (m^3)
C = 633; %Specific heat capacity, J/(kg*k)
rho = 4648; %Density, kg/m^3
alfa_KTP = 0.01e2; %Absorption coefficient (1/m)

OMEGA = 0; %Sideband frequency 
s = 1i.*OMEGA; 

delta_P = linspace(-0.1,0.1,1000); %Input fluctuations in delta_P_in/P_in
P_in_eff = (1+delta_P).*P_in;

B = sqrt(P_in).*tanh(L.*sqrt(P_in.*eta_norm)); %Output power square root
nB = sqrt(eta_norm).*L.*B; %sqrt(Output power) times sqrt(eta_norm)*L

B_eff = sqrt(P_in_eff).*tanh(L.*sqrt(P_in_eff.*eta_norm)); %Output power square root
nB_eff = sqrt(eta_norm).*L.*B_eff; %sqrt(Output power) times sqrt(eta_norm)*L

delta_b = (B_eff-B)./B; %Re(delta_b/B), assuming same phase
P_abs = alfa_KTP.*L.*abs(B_eff).^2;% Absorbed power

dT = 1./(1+s.*tau).*2.*tau.*P_abs./(C.*rho.*V).*delta_b; %Temperature fluctuation

%Sellmeier coefficients
a1 = 4.9048;
a2 = 0.11775;
a3 = 0.21802;
a4 = 0.027153;
b1 = 2.2314e-6;
b2 = -2.9671e-8;
b3 = 2.1429e-8;


T_0 = 30; %Celsius
f = dT.*(dT+2*T_0+546);
dfdT = 2.*(dT+T_0) + 546;

n = @(lambda)  sqrt(a1 + b3.*f + (a2+b1.*f)./(lambda.^2 - (a3 + b2.*f).^2) - a4.*lambda.^2);

dndT = @(lambda) dfdT.*(2.*b2.*(a2+b1.*f).*(a3+b2.*f)./(lambda.^2-(a3+b2.*f).^2).^2 + b1./(lambda.^2-(a3+b2.*f).^2) +b3)./(2*n(lambda));

alfa = 7.5e-6; %Thermal expansion coefficient
lambda_0 = 1550e-9; %Input wavelength
lambda_SHG = 775e-9; %Wavelength of generated SHG light

%Delta k, note that the formula for n inputs the wavelength in microns
delta_k = dT.*pi.*4./lambda_0.*abs((dndT(lambda_SHG*1e3) - dndT(lambda_0*1e3)) + alfa.*(n(lambda_SHG) - n(lambda_0))); %Phase mismatch

pm = delta_k.*L;

%Analytical solution, several assumptions made
V_an_X = (sin(pm).*sinh(nB)).^2 + (cos(pm).*sinh(nB) + cosh(nB)).^2;
V_an_Y = (sin(pm).*sinh(nB)).^2 + (cos(pm).*sinh(nB) - cosh(nB)).^2;

hold on
plot(delta_P,pow_to_dB(V_an_Y), 'Displayname', 'Squeezing')
plot(delta_P,delta_P.*0,'k','Displayname', 'Shot noise')

lgd = legend;
lgd.FontSize = 10;
lgd.Location = 'best';

xlabel('\delta P_{in}/P_{in}','FontSize',16);
ylabel('Variance (dB)','Fontsize',16);

%% Squeezing vs frequency fluctuations

dB_to_pow = @(d) 10.^(d./10);
pow_to_dB = @(p) 10.*log10(p);

P_in = 660e-3; %Input power
eta_norm = 0.4e4; %1/W*m^2
L = 45e-3;

B = sqrt(P_in).*tanh(L.*sqrt(P_in.*eta_norm)); %Output power square root
nB = sqrt(eta_norm).*L.*B; %sqrt(Output power) times sqrt(eta_norm)*L

%Sellmeier coefficients
a1 = 4.9048;
a2 = 0.11775;
a3 = 0.21802;
a4 = 0.027153;
b1 = 2.2314e-6;
b2 = -2.9671e-8;
b3 = 2.1429e-8;

T = 30; %Celsius
T_0 = 30; %Celsius
dT = T-T_0;
f = (T-T_0).*(T+T_0+546);

%Formula for n and its derivative. Note that lambda here is assumed to be 
%in microns.
n = @(lambda)  sqrt(a1 + b3.*f + (a2+b1.*f)./(lambda.^2 - (a3 + b2.*f).^2) - a4.*lambda.^2);
dndlambda = @(lambda) -(2.*lambda.*(a2+b1.*f)./(lambda.^2 - (a3 + b2.*f).^2).^2 +2.*a4.*lambda)./(2.*n(lambda));

lambda_0 = 1550e-9;
lambda_SHG = lambda_0/2; 


d_omega = linspace(-1e-5,1e-5,1000); %Frequency fluctuations
d_lambda = lambda_0./(1+d_omega); %Convert to wavelength

%Formula for delta k
delta_k = 2.*pi.*(2.*dndlambda(lambda_0.*1e3).*1e6-dndlambda(lambda_SHG.*1e3).*1e6).*(1-lambda_0./d_lambda); %Frequency fluctuation of input pump

pm = delta_k.*L; 

%Analytical solution, several assumptions made
V_an_X = (sin(pm).*sinh(nB)).^2 + (cos(pm).*sinh(nB) + cosh(nB)).^2;
V_an_Y = (sin(pm).*sinh(nB)).^2 + (cos(pm).*sinh(nB) - cosh(nB)).^2;

hold on
plot(1e6.*d_omega,pow_to_dB(V_an_Y), 'Displayname', 'Squeezing')
plot(1e6.*d_omega,d_omega.*0,'k','Displayname', 'Shot noise')

lgd = legend;
lgd.FontSize = 10;
lgd.Location = 'best';

xlabel('\delta \omega /\omega \cdot 10^{6}','FontSize',16);
ylabel('Variance (dB)','Fontsize',16);