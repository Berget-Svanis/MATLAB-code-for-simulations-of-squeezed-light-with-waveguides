% Squeezing vs power MRAD

dB_to_pow = @(d) 10.^(d./10);
pow_to_dB = @(p) 10.*log10(p);

P_in = linspace(0,600e-3,1000);


eta_norm = 0.4e4; %Normalized efficiency (1/(W m^2)
L = 45e-3; %Length

%Defining squeezing and anti-squeezing
%How to add eta?
sqz = @(P,n) n.*exp(-2.*L.*sqrt(eta_norm.*P).*tanh(L.*sqrt(eta_norm.*P))) + (1-n);
asqz = @(P,n) n.*exp(2.*L.*sqrt(eta_norm.*P).*tanh(L.*sqrt(eta_norm.*P))) + (1-n);
sqz_rad = @(P,n,rad) sqz(P,n).*cos(rad).^2 + asqz(P,n).*sin(rad).^2;
asqz_rad = @(P,n,rad) asqz(P,n).*cos(rad).^2 + sqz(P,n).*sin(rad).^2;

col = ['r','b', 'g']; %Defines colours for the graph
eta = 1;
phi_ls = [0 40 150];
for i = 1:3
    phi = phi_ls(i)./1e3;
    sqz_dB = pow_to_dB(sqz_rad(P_in,eta,phi));
    asqz_dB = pow_to_dB(asqz_rad(P_in,eta,phi));

    %Plotting
    hold on
    plot(P_in.*1e3,sqz_dB,col(i),'DisplayName',strcat('Squeezing, mrad = ', num2str(phi_ls(i))))
    plot(P_in.*1e3,asqz_dB,col(i),'DisplayName',strcat('Antisqueezing, mrad = ', num2str(phi_ls(i))))
end
plot(P_in.*1e3,P_in.*0,'k','Handlevisibility','off')

lgd = legend;
lgd.FontSize = 10;
lgd.Location = 'southwest';

ylabel('Variance (dB)','FontSize',16); 
xlabel('Input power (mW)','Fontsize',16);