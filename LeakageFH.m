%%RUN SECTIONS SEPERATELY FOR SEPERATRE PLOTS
%FH leakage to SHG

%% LO phase and coherent leakage vs measured variance
close all;
%Convert from/to dB
dB_to_pow = @(d) 10.^(d./10);
pow_to_dB = @(p) 10.*log10(p);

%Input parameters
num_points = 100000;
phi_LO = linspace(0,2.*pi,num_points);
eps = [0 0.1 0.5 0.9];

init_sqz = 3; %Intitial squeezing/anti-squeezing in dB
%Variance in squeezer path
X = dB_to_pow(init_sqz);
Y = dB_to_pow(-init_sqz);

%Variance in LO path
X_LO = dB_to_pow(0);
Y_LO = dB_to_pow(0);

for i = 1:numel(eps)
    phi_2 = -acos(2.*sqrt(eps(i)).*cos(phi_LO)./(sqrt(4.*eps(i).*cos(phi_LO).^2 + eps(i).^2 - 2.*eps(i) +1)));

    p_var = @(V) 4.*(cos(phi_2).^2.*(cos(phi_LO).^2.*X_LO + sin(phi_LO).^2.*Y_LO +4.*eps(i).*V(1)) + ...  
    sin(phi_2).^2.*(cos(phi_LO).^2.*V(1) + sin(phi_LO).^2.*V(2) +4.*eps(i).*X_LO) + ...
    2.*sin(2.*phi_2).*sqrt(eps(i)).*cos(phi_LO).*(X_LO-V(1)));

    %Shot noise has unity variance
    shot_noise = p_var([1,1]);

    hold on

    plot(phi_LO,p_var([X,Y])./shot_noise,'DisplayName', strcat('\epsilon = ', num2str(eps(i))))

end

%Plotting
plot(phi_LO,ones(1,num_points),'k','DisplayName', 'Shot Noise')

lgd = legend;
lgd.FontSize = 10;
lgd.Location = 'best';

ylabel('Noise level (normalized)','FontSize',16); 
xlabel('LO phase','Fontsize',16);

title(strcat('Noise vs LO phase for different coherent leakage levels, initial squeezing = ',num2str(init_sqz), ' dB'))

ax=gca;
ax.FontSize = 12;
set(ax, 'xlim', [0 2*pi]);
set(ax,'XTick',0:pi/2:2*pi) 
set(ax,'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})

%% Squeezing vs leakage POWER, assuming an intial squeeze/anti-squeeze of X dB
close all;

dB_to_pow = @(d) 10.^(d./10);
pow_to_dB = @(p) 10.*log10(p);

eps = linspace(0,0.5,100000);

LO_P = 10; %Power of the local oscillator in mW

sqz = 10; %Initial squeezing in dB

%Variance in squeezer path
X = dB_to_pow(sqz);
Y = dB_to_pow(-sqz);

%Variance in LO path
X_LO = dB_to_pow(0);
Y_LO = dB_to_pow(0);

%Antisqueezing
phi_LO = 0;
phi_2 = -acos(2.*sqrt(eps).*cos(phi_LO)./(sqrt(4.*eps.*cos(phi_LO).^2 + eps.^2 - 2.*eps +1)));

p_var = @(V) 4.*(cos(phi_2).^2.*(cos(phi_LO).^2.*X_LO + sin(phi_LO).^2.*Y_LO +4.*eps.*V(1)) + ...  
sin(phi_2).^2.*(cos(phi_LO).^2.*V(1) + sin(phi_LO).^2.*V(2) +4.*eps.*X_LO) + ...
2.*sin(2.*phi_2).*sqrt(eps).*cos(phi_LO).*(X_LO-V(1)));

%Shot noise has unity variance
shot_noise = p_var([1,1]);

plot(eps.*LO_P,pow_to_dB(p_var([X,Y])./shot_noise),'b', 'DisplayName','Antisqueezing')

%Squeezing
phi_LO = pi/2;
phi_2 = -acos(2.*sqrt(eps).*cos(phi_LO)./(sqrt(4.*eps.*cos(phi_LO).^2 + eps.^2 - 2.*eps +1)));

p_var = @(V) 4.*(cos(phi_2).^2.*(cos(phi_LO).^2.*X_LO + sin(phi_LO).^2.*Y_LO +4.*eps.*V(1)) + ...  
sin(phi_2).^2.*(cos(phi_LO).^2.*V(1) + sin(phi_LO).^2.*V(2) +4.*eps.*X_LO) + ...
2.*sin(2.*phi_2).*sqrt(eps).*cos(phi_LO).*(X_LO-V(1)));

%Shot noise has unity variance
shot_noise = p_var([1,1]);

hold on

plot(eps.*LO_P,pow_to_dB(p_var([X,Y])./shot_noise),'b', 'DisplayName','Antisqueezing')
plot(eps.*LO_P,eps.*0,'k','Displayname','Shot noise')

xlabel('Leakage power (mW)','FontSize',16);
ylabel('Variance (dB)','Fontsize',16);

title(strcat('Squeezing vs power leakage. Initial squeezing=  ',num2str(sqz) ,' dB') ,'Fontsize',14)

ax=gca;
