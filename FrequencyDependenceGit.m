%Frequency dependency script for GitHub
%Run each section separately for separate plots! 
%% Frequency dependence squeezer WG

sideband_freq = linspace(1e6,1e9,10000);
carrier_freq = 1.9355*2*pi*1e14; %Assuming an input wavelength of 1550 nm
L = 45e-3;
P_in = 660e-3;

eta_norm = 0.4e4; %1/(W*m^2) Normalized efficiency

eta_norm = eta_norm.*((carrier_freq + 2.*pi.*sideband_freq)./carrier_freq).^2; %Converting to conversion effiency and including dependence on sideband freq.

%Requires to have sqzWG in the same folder! Script exists on GitHub
[s,as] = sqzWG(L,P_in,eta_norm,1,0.88,1,0,14e-3); 

figure(1) %Creates the first plot

semilogx(sideband_freq,s,'Displayname', 'Squeezing')
hold on 
semilogx(sideband_freq,as,'Displayname', 'Anti-squeezing')
semilogx(sideband_freq,sideband_freq.*0,'k','Displayname', 'Shot noise') %Plots shot noise level

lgd = legend;
lgd.FontSize = 10;
lgd.Location = 'best';

xlabel('Sideband Frequency (Hz)','FontSize',16);
ylabel('Variance (dB)','Fontsize',16);

title('Squeezing vs sideband frequency, carrier wavelength = 1550 nm','Fontsize',14)
grid
%% SHG efficiency for a certain wavelength 

sinc = @(x) sin(x)./x; %sinc function

lambda_0 = 1550e-9; %1550 nm input wavelength
dn_0 = -0.033912e6; %Derivative of refractive index at lambda_0 (LiNbO_3)
dn_SHG = -0.14181e6; %Derivative of refractive index at lambda_0(LiNbO_3)

L = 45e-3; %10 mm length

lambda = linspace(1548e-9,1552e-9,1000);
delta_k = 2.*pi.*(2.*dn_0-dn_SHG).*(1-lambda_0./lambda);

eff = (sinc(delta_k.*L./2)).^2;

figure(2) %Creates the second plot

plot(lambda.*1e9,eff)

xlabel('Wavelength (nm)','FontSize',16);
ylabel('SHG efficiency','Fontsize',16);

%Calculating FWHM. 1.39155 is the HWHM of sinc^2(x)
FWHM = (2.*lambda_0.*(1.39155)./(pi.*L.*(2.*dn_0-dn_SHG)))./(1-1.39155./(pi.*L.*(2.*dn_0-dn_SHG))); 

t = annotation('textbox', [0.65, 0.8, 0.1, 0.1], 'String', strcat('FWHM = ', num2str(round(FWHM*1e9,2)), ' nm'));
t.FontSize = 14;
grid
%% SHG efficiency vs sideband frequency 

sinc = @(x) sin(x)./x;

c = 3e8;
lambda_0 = 1550e-9; %1550 nm input wavelength
carrier_freq = c./lambda_0; %Carrier frequency, f
lambda_SHG = lambda_0/2; %Half the wavelength is the SHG output
dn_0 = -0.033912e6; %Derivative of refractive index at lambda_0 (LiNbO_3)
dn_SHG = -0.14181e6; %Derivative of refractive index at lambda_0(LiNbO_3)

OMEGA = 1e9; %Sideband frequency
n = 10;
freq_ls = linspace(-n*OMEGA,n*OMEGA,10000);

L = 45e-3; %10 mm length

delta_k = 2.*pi.*(2.*dn_0-dn_SHG).*freq_ls./carrier_freq;

eff = (sinc(delta_k.*L./2)).^2;

figure(3) %Creates the third plot

plot(freq_ls./OMEGA,eff)

xlabel('Sideband frequencies (\Omega)','FontSize',16);
ylabel('SHG efficiency','Fontsize',16);

ax=gca;
ax.FontSize = 12;
set(ax, 'xlim', [-10 10]);
set(ax,'XTick',-10:2:10) 
x_tick_label_ls = cell(1,11);
for i=1:11
    num = -12 + 2*i; 
    x_tick_label_ls{i} = strcat(num2str(num), '\Omega'); 
end
set(ax,'XTickLabel',x_tick_label_ls)

title(strcat('SHG efficiency vs sideband frequency for \Omega = ',num2str(OMEGA./1e9), ' GHz'), 'Fontsize',14)
grid