%FH leakage to SHG
% LO phase and coherent leakage vs measured variance
close all;
%Convert from/to dB
dB_to_pow = @(d) 10.^(d./10);
pow_to_dB = @(p) 10.*log10(p);

%Input parameters
num_points = 100000;
phi_LO = linspace(0,2.*pi,num_points);
eps = [0.01 0.1 0.5];

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

ax=gca;
ax.FontSize = 12;
set(ax, 'xlim', [0 2*pi]);
set(ax,'XTick',0:pi/2:2*pi) 
set(ax,'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})
