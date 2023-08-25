%Cascaded squeezer waveguides
%Run sections separately for individual plots 
%% esc2 and eta_2 vs eta_eff

eff_ls = [0.2 0.4 0.6 0.8];

cmap = flip(autumn(4),1); % yellow -> red, with  4 colors (for 4 lines)
set(gca(),'ColorOrder',cmap)
num = numel(eff_ls);
c = parula(num); %Sets color gradient
figure(1) %Creates the first plot

n = 100; %Amount of data points
R2 = 1.5; %Anti-squeezing parameter

hold on
for j=1:num
    disp(j) %Shows how far the code has run
    eta_2 = linspace(0.001,1,n);
    
    sqz = zeros(1,n);
    for i=1:n
        %Solves for eta_sqz2 given a certain eta_2 and eta_eff
        syms x
        S = vpasolve(eff_ls(j) == eta_2(i).*x./((1-eta_2(i).*x).*exp(-2.*R2) + eta_2(i).*x), x);
        sqz(i) = S;
    end
    %Plotting
    plot(eta_2,sqz,'Color',c(j,:))
    text(eta_2(n-1),sqz(n-1),num2str(eff_ls(j)))
end

xlabel('\eta_2','FontSize',16);
ylabel('\eta_{sqz2}','Fontsize',16);
title('\eta_{eff} vs \eta_2 and \eta_{sqz2}','Fontsize',16);

ax=gca;
ax.FontSize = 12;
set(ax, 'xlim', [0 1]);
set(ax, 'ylim', [0 1]);
grid
%% Amplified vacuum and squuezed state comparison

%dB conversions
dB_to_pow = @(d) 10.^(d./10);
pow_to_dB = @(p) 10.*log10(p);

%Defining efficiencies
eta_1 = 1;
eta_2 = 0.3;
sqz1 = 1;
sqz2 = 1;

R1 = 1.38; %Generates 12 dB squeezing
n_R2 = 100; %How many values of R_2 we consider in the plotting
R2 = linspace(0,3,n_R2);

V_vac_eff = zeros(1,n_R2);
V_sqz_eff = zeros(1,n_R2);
for i=1:n_R2
    %Defining M_WG and M_WGL matrices. The first is squeezing and the
    %second is anti-squeezing in the X-quadrature. 
    M_WG1 = @(R_1) [sqrt(sqz1).*exp(-R_1), 0; 0, sqrt(sqz1).*exp(R_1)];
    M_L1 = @(R_1) [sqrt(1-sqz1), 0; 0, sqrt(1-sqz1)];

    M_WG2 = [sqrt(sqz2).*exp(R2(i)), 0; 0, sqrt(sqz2).*exp(-R2(i))];
    M_L2 = [sqrt(1-sqz2), 0; 0, sqrt(1-sqz2)];

    %Defining transfer matrices
    TF_in = @(R_1) sqrt(eta_1.*eta_2).*M_WG2*M_WG1(R_1);
    TF_WG1L = @(R_1) sqrt(eta_1.*eta_2).*M_WG2*M_L1(R_1);
    TF_1L = @(R_1) sqrt((1-eta_1).*eta_2).*M_WG2;
    TF_WG2L = @(R_1) sqrt(eta_2).*M_L2;
    TF_2L = @(R_1) sqrt(1-eta_2).*eye(2);

    %Amplified squeezing/anti-squeezing
    V_amp_sqz = TF_in(R1).^2 + TF_WG1L(R1).^2 + ...
        TF_1L(R1).^2 + TF_WG2L(R1).^2 + TF_2L(R1).^2; 

    %Amplified vacuum (R_1=0)
    V_amp_vac = TF_in(0).^2 + TF_WG1L(0).^2 + ...
        TF_1L(0).^2 + TF_WG2L(0).^2 + TF_2L(0).^2; 
    
    %Amplified squeezed state and amplified vacuum
    V_sqz_eff(i) = V_amp_sqz(1,1) + V_amp_sqz(1,2);
    V_vac_eff(i) = V_amp_vac(1,1) + V_amp_vac(1,2);
end

%Plotting
figure(2) %Creates the second plot
%First subplot, amplified squeezed state and amplified vacuum
subplot(3,1,[1 2])

hold on
plot(R2,pow_to_dB(V_sqz_eff),'Displayname','Amplified squeezed state')
plot(R2,pow_to_dB(V_vac_eff),'Displayname','Amplified vacuum')
plot(R2,R2.*0,'k','Displayname','Shot noise')

ylabel('Noise relative vacuum level (dB)','Fontsize',16);

title(strcat(num2str((1-eta_2).*1e2),' % detection loss, 12 dB initial squeezing'),'Fontsize',20)

ax=gca;
ax.FontSize = 12;
lgd = legend;
lgd.FontSize = 10;
lgd.Location = 'best';

%Second subplot, effective noise reduction
subplot(3,1,3)
plot(R2,pow_to_dB(V_sqz_eff./V_vac_eff),'Displayname','Effective noise reduction');

xlabel('Squeezer 2 parameter, R_2','FontSize',16);
ylabel('Noise level (dB)','Fontsize',16);
grid

%% Variance vs eta_2 loss for different squeezing parameters R2

%Cascaded squeezing

%dB conversions
dB_to_pow = @(d) 10.^(d./10);
pow_to_dB = @(p) 10.*log10(p);

n = 1000; 
%Defining efficiencies
eta_1 = 1;
eta_2 = linspace(0,1,n);
sqz1 = 1;
sqz2 = 1;

R2 = [0,0.5,1,2,5,10];

figure(3) %Creates the third plot

for j=1:numel(R2)
    V_eff = zeros(1,n);
    for i=1:n
        %Defining M_WG and M_WGL matrices. The first is squeezing and the
        %second is anti-squeezing in the X-quadrature. 
        M_WG1 = @(R_1) [sqrt(sqz1).*exp(-R_1), 0; 0, sqrt(sqz1).*exp(R_1)];
        M_L1 = @(R_1) [sqrt(1-sqz1), 0; 0, sqrt(1-sqz1)];

        M_WG2 =@(R_2) [sqrt(sqz2).*exp(R_2), 0; 0, sqrt(sqz2).*exp(-R_2)];
        M_L2 =@(R_2) [sqrt(1-sqz2), 0; 0, sqrt(1-sqz2)];

        %Defining transfer matrices
        TF_in = @(R_1,R_2) sqrt(eta_1.*eta_2(i)).*M_WG2(R_2)*M_WG1(R_1);
        TF_WG1L = @(R_1,R_2) sqrt(eta_1.*eta_2(i)).*M_WG2(R_2)*M_L1(R_1);
        TF_1L = @(R_1,R_2) sqrt((1-eta_1).*eta_2(i)).*M_WG2(R_2);
        TF_WG2L = @(R_1,R_2) sqrt(eta_2(i)).*M_L2(R_2);
        TF_2L = @(R_1,R_2) sqrt(1-eta_2(i)).*eye(2);

        %Squeezing parameter
        R1 = 0.345; %0.345 for 3dB initial squeezing, 1.15 for 10 dB, 1.727 for 15 dB 

       %Amplified squeezing/anti-squeezing
        V_amp_sqz = TF_in(R1,R2(j)).^2 + TF_WG1L(R1,R2(j)).^2 + ...
            TF_1L(R1,R2(j)).^2 + TF_WG2L(R1,R2(j)).^2 + TF_2L(R1,R2(j)).^2; 

        %Amplified vacuum (R1=0)
        V_amp_vac = TF_in(0,R2(j)).^2 + TF_WG1L(0,R2(j)).^2 + ...
            TF_1L(0,R2(j)).^2 + TF_WG2L(0,R2(j)).^2 + TF_2L(0,R2(j)).^2;  
        
        V_eff(i) = (V_amp_sqz(1,1) + V_amp_sqz(1,2))./(V_amp_vac(1,1) + V_amp_vac(1,2));
    end
    %Plotting
    plot(100.*(1-eta_2),pow_to_dB(V_eff),'Displayname',strcat('R_2 = ', num2str(R2(j))))
    hold on
end

xlabel('Loss, 1-\eta_2 (%)','FontSize',16);
ylabel('Variance (dB)','Fontsize',16);

title(strcat('Squeezing vs detection loss, initial squeezing = ',num2str(round(pow_to_dB(exp(2*R1)))),' dB' ), 'Fontsize', 14)

lgd = legend;
lgd.FontSize = 10;
lgd.Location = 'best';
grid