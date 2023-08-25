%Cascaded squeezers in WG with phase noise
%Run sections separately for separate plots! 
%% Detection loss and phase noise of 1st squeezer vs squeezing level contour plot

%Cascaded squeezing
%dB conversions
dB_to_pow = @(d) 10.^(d./10);
pow_to_dB = @(p) 10.*log10(p);

%Defining efficiencies
eta_1 = 1;
sqz1 = 1;
sqz2 = 1;

n = 100;
theta_1 = linspace(0.01,100,n)./1e3; %rad
theta_2 = linspace(0,0,n)./1e3; %rad

%Defining M_WG and M_WGL matrices. The first is squeezing and the
%second is anti-squeezing in the X-quadrature. 
M_WG1 = @(R_1) [sqrt(sqz1).*exp(-R_1), 0; 0, sqrt(sqz1).*exp(R_1)];
M_L1 = @(R_1) [sqrt(1-sqz1), 0; 0, sqrt(1-sqz1)];

M_WG2 =@(R_2) [sqrt(sqz2).*exp(R_2), 0; 0, sqrt(sqz2).*exp(-R_2)];
M_L2 =@(R_2) [sqrt(1-sqz2), 0; 0, sqrt(1-sqz2)];

V_ls_dB = [-1,-2,-3,-4,-5,-6,-7,-8,-9,-10]; %Array of squeezing levels
V_eff_ls = dB_to_pow(V_ls_dB);

num = numel(V_eff_ls);
cmap = flip(autumn(num),1); % yellow -> red
set(gca(),'ColorOrder',cmap)
c = parula(num); %Sets color gradient

figure(1) %Creates the first plot

for j=1:num
    disp(j) %Checks how far the code has run
    eta_2_ls = zeros(1,n);
    for i=1:n
        %Rotation matrices
        Rot_1 = [cos(theta_1(i)), -sin(theta_1(i)); sin(theta_1(i)), cos(theta_1(i))];
        Rot_1_inv = [cos(theta_1(i)), sin(theta_1(i)); -sin(theta_1(i)), cos(theta_1(i))];
        Rot_2 = [cos(theta_2(i)), -sin(theta_2(i)); sin(theta_2(i)), cos(theta_2(i))];
        Rot_2_inv = [cos(theta_2(i)), sin(theta_2(i)); -sin(theta_2(i)), cos(theta_2(i))];

        %Defining transfer matrices
        TF_in = @(R1,R2,eta_2) sqrt(eta_1.*eta_2).*Rot_2*M_WG2(R2)*Rot_1*M_WG1(R1)*Rot_1_inv*Rot_2_inv;
        TF_WG1L = @(R1,R2,eta_2) sqrt(eta_1.*eta_2).*Rot_2*M_WG2(R2)*Rot_1*M_L1(R1)*Rot_1_inv*Rot_2_inv;
        TF_1L = @(R1,R2,eta_2) sqrt((1-eta_1).*eta_2).*Rot_2*M_WG2(R2)*Rot_2_inv;
        TF_WG2L = @(R1,R2,eta_2) sqrt(eta_2).*Rot_2*M_L2(R2)*Rot_2_inv;
        TF_2L = @(R1,R2,eta_2) sqrt(1-eta_2).*eye(2);

        %Pump parameters
        R1 = 1.727;
        R2 = 2;
           
        %Amplified squeezing/anti-squeezing
        V_amp_sqz = @(eta_2) TF_in(R1,R2,eta_2).^2 + TF_WG1L(R1,R2,eta_2).^2 + ...
            TF_1L(R1,R2,eta_2).^2 + TF_WG2L(R1,R2,eta_2).^2 + TF_2L(R1,R2,eta_2).^2; 

        %Amplified vacuum (R1=0)
        V_amp_vac = @(eta_2) TF_in(0,R2,eta_2).^2 + TF_WG1L(0,R2,eta_2).^2 + ...
            TF_1L(0,R2,eta_2).^2 + TF_WG2L(0,R2,eta_2).^2 + TF_2L(0,R2,eta_2).^2; 

        %Takes out the (r,c) element of V
        select = @(V,r,c) V(r,c);
        V_sqz_11 = @(eta_2) select(V_amp_sqz(eta_2),1,1);
        V_vac_11 = @(eta_2) select(V_amp_vac(eta_2),1,1);
        V_sqz_12 = @(eta_2) select(V_amp_sqz(eta_2),1,2);
        V_vac_12 = @(eta_2) select(V_amp_vac(eta_2),1,2);

        find_eta2 = @(eta_2) (V_sqz_11(eta_2) + V_sqz_12(eta_2))./(V_vac_11(eta_2) + V_vac_12(eta_2)) - V_eff_ls(j); 

        %Solving for eta_2
        options = optimset('Display','off');
        eta_2_ls(i) = fsolve(find_eta2,0.5,options);
    end
    %Plotting
    plot(theta_1.*1e3,100.*(1-eta_2_ls),'Color',c(j,:))
    hold on
    text(theta_1(10).*1e3,100.*(1-eta_2_ls(10)),num2str(V_ls_dB(j)))
end

ylabel('Detection loss, 1-\eta_2 (%)','FontSize',16);
xlabel('Phase noise of WG 1, \theta_1 (mrad)','Fontsize',16);

title(strcat('Initial squeezing = ',num2str(round(pow_to_dB(eta_1.*sqz1.*exp(-2.*R1) + 1-eta_1.*sqz1),1)), ' dB. R_2 = ',num2str(R2))); 

ax=gca;
ax.FontSize = 12;
set(ax, 'xlim', [0 100]);
set(ax, 'ylim', [0 100]);
grid
%% Detection loss and phase noise of 2nd squeezer vs squeezing level contour plot

%Cascaded squeezing

%dB conversions
dB_to_pow = @(d) 10.^(d./10);
pow_to_dB = @(p) 10.*log10(p);

%Defining efficiencies
eta_1 = 1;
sqz1 = 1;
sqz2 = 1;

n = 100;
theta_1 = linspace(0,0,n)./1e3; %rad
theta_2 = linspace(0.01,100,n)./1e3; %rad

%Defining M_WG and M_WGL matrices. The first is squeezing and the
%second is anti-squeezing in the X-quadrature. 
M_WG1 = @(R_1) [sqrt(sqz1).*exp(-R_1), 0; 0, sqrt(sqz1).*exp(R_1)];
M_L1 = @(R_1) [sqrt(1-sqz1), 0; 0, sqrt(1-sqz1)];

M_WG2 =@(R_2) [sqrt(sqz2).*exp(R_2), 0; 0, sqrt(sqz2).*exp(-R_2)];
M_L2 =@(R_2) [sqrt(1-sqz2), 0; 0, sqrt(1-sqz2)];

V_ls_dB = [-1,-2,-3,-4,-5,-6,-7,-8,-9,-10]; %Array of squeezing levels
V_eff_ls = dB_to_pow(V_ls_dB);

num = numel(V_eff_ls);
cmap = flip(autumn(num),1); % yellow -> red
set(gca(),'ColorOrder',cmap)
c = parula(num); %Sets color gradient
figure(2) %Creates the second plot

for j=1:num
    disp(j) %Checks how far the code has run
    eta_2_ls = zeros(1,n);
    for i=1:n
        %Rotation matrices
        Rot_1 = [cos(theta_1(i)), -sin(theta_1(i)); sin(theta_1(i)), cos(theta_1(i))];
        Rot_1_inv = [cos(theta_1(i)), sin(theta_1(i)); -sin(theta_1(i)), cos(theta_1(i))];
        Rot_2 = [cos(theta_2(i)), -sin(theta_2(i)); sin(theta_2(i)), cos(theta_2(i))];
        Rot_2_inv = [cos(theta_2(i)), sin(theta_2(i)); -sin(theta_2(i)), cos(theta_2(i))];

        %Defining transfer matrices
        TF_in = @(R1,R2,eta_2) sqrt(eta_1.*eta_2).*Rot_2*M_WG2(R2)*Rot_1*M_WG1(R1)*Rot_1_inv*Rot_2_inv;
        TF_WG1L = @(R1,R2,eta_2) sqrt(eta_1.*eta_2).*Rot_2*M_WG2(R2)*Rot_1*M_L1(R1)*Rot_1_inv*Rot_2_inv;
        TF_1L = @(R1,R2,eta_2) sqrt((1-eta_1).*eta_2).*Rot_2*M_WG2(R2)*Rot_2_inv;
        TF_WG2L = @(R1,R2,eta_2) sqrt(eta_2).*Rot_2*M_L2(R2)*Rot_2_inv;
        TF_2L = @(R1,R2,eta_2) sqrt(1-eta_2).*eye(2);

        %Pump parameters
        R1 = 1.727;
        R2 = 2;
           
        %Amplified squeezing/anti-squeezing
        V_amp_sqz = @(eta_2) TF_in(R1,R2,eta_2).^2 + TF_WG1L(R1,R2,eta_2).^2 + ...
            TF_1L(R1,R2,eta_2).^2 + TF_WG2L(R1,R2,eta_2).^2 + TF_2L(R1,R2,eta_2).^2; 

        %Amplified vacuum (R1=0)
        V_amp_vac = @(eta_2) TF_in(0,R2,eta_2).^2 + TF_WG1L(0,R2,eta_2).^2 + ...
            TF_1L(0,R2,eta_2).^2 + TF_WG2L(0,R2,eta_2).^2 + TF_2L(0,R2,eta_2).^2; 

        %Takes out the (r,c) element of V
        select = @(V,r,c) V(r,c);
        V_sqz_11 = @(eta_2) select(V_amp_sqz(eta_2),1,1);
        V_vac_11 = @(eta_2) select(V_amp_vac(eta_2),1,1);
        V_sqz_12 = @(eta_2) select(V_amp_sqz(eta_2),1,2);
        V_vac_12 = @(eta_2) select(V_amp_vac(eta_2),1,2);

        find_eta2 = @(eta_2) (V_sqz_11(eta_2) + V_sqz_12(eta_2))./(V_vac_11(eta_2) + V_vac_12(eta_2)) - V_eff_ls(j); 

        %Solving for eta_2
        options = optimset('Display','off');
        eta_2_ls(i) = fsolve(find_eta2,0.5,options);
    end
    %Plotting
    plot(theta_2.*1e3,100.*(1-eta_2_ls),'Color',c(j,:))
    hold on
    text(theta_2(10).*1e3,100.*(1-eta_2_ls(10)),num2str(V_ls_dB(j)))
end

ylabel('Detection loss, 1-\eta_2 (%)','FontSize',16);
xlabel('Phase noise of WG 2, \theta_2 (mrad)','Fontsize',16);

title(strcat('Initial squeezing = ',num2str(round(pow_to_dB(eta_1.*sqz1.*exp(-2.*R1) + 1-eta_1.*sqz1),1)), ' dB. R_2 = ',num2str(R2))); 

ax=gca;
ax.FontSize = 12;
set(ax, 'xlim', [0 100]);
set(ax, 'ylim', [0 100]);
grid