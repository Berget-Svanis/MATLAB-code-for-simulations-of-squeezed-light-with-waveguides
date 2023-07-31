%SHG output vs input power as a function of length 

P_in = 150e-3; %Input power
eta_norm = 10e4; %Normalized efficiency, 1/(W*m^2)
L = linspace(0,30e-3,1000); %Length of waveguide 

%SH and FH power. See report
P_SH = P_in.*tanh(sqrt(eta_norm.*P_in).*L).^2;
P_FH = P_in.*1./cosh(sqrt(eta_norm.*P_in).*L).^2;

%Plottinng
plot(L.*1e3,P_SH.*1e3,'Displayname','Second harmonic power')
hold on 
plot(L.*1e3,P_FH.*1e3,'Displayname','Fundamental harmonic power')

xlabel('Length of waveguide (mm)','FontSize',16);
ylabel('Power (mW)','Fontsize',16);

title(strcat('SH and FH power vs waveguide length. P_{in} = ', num2str(P_in.*1e3), ' mW'),'Fontsize',14)

lgd = legend;
lgd.FontSize = 10;
lgd.Location = 'best';