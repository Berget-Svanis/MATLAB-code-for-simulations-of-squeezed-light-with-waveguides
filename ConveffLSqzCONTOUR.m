%Conversion efficiency vs length, squeezing contour plot

dB_to_pow = @(d) 10.^(d./10);
pow_to_dB = @(p) 10.*log10(p);

%The squeezing levels we are looking at
dB = linspace(-3,-15,13);
var = dB_to_pow(dB);

P = 100e-3; %Power

num_points = 1000;
L = linspace(0.5,5,num_points); %Length in cm

cmap = flip(autumn(16),1); % yellow -> red, with  16 colors (for 16 lines)
set(gca(),'ColorOrder',cmap)

num = numel(var);
c = parula(num); %Sets color gradient

syms x
for i = 1:num
    disp(i)
    %Solving for eta_norm
    S = vpasolve(-log(var(i)) == 2*x*tanh(x),x);
    eta_norm = ((S./L).^2)./P;
    
    %Plotting
    semilogy(L.*10,eta_norm,'Color',c(i,:))
    hold on 
    if mod(i,3) == 0
       text(L(20).*10,eta_norm(20),strcat(num2str(10.*log10(var(i))), 'dB'))
    end
end

xlabel('Length (mm)','FontSize',16);
ylabel('Normalized efficiency (1/(cm^2W)','Fontsize',16);


title('Contour plot: Normalized efficiency and Length vs Squeezing. Semilogy', 'Fontsize',12)

%% Conversion efficiency and length vs POWER. For a specific squeezing

%Conversion efficiency vs length, squeezing contour plot

dB_to_pow = @(d) 10.^(d./10);
pow_to_dB = @(p) 10.*log10(p);

%The squeezing level we are looking at
dB = -10;
var = dB_to_pow(dB);

P = linspace(10e-3,150e-3,15); %Power

num_points = 1000;
L = linspace(0.5,5,num_points); %Length in cm

cmap = flip(autumn(16),1); % yellow -> red, with  16 colors (for 16 lines)
set(gca(),'ColorOrder',cmap)

num = numel(P);
c = parula(num); %Sets color gradient

syms x
for i = 1:num
    disp(i)
    %Solving for eta_norm
    S = vpasolve(-log(var) == 2*x*tanh(x),x);
    eta_norm = ((S./L).^2)./P(i);
    
    %Plotting
    semilogy(L.*10,eta_norm,'Color',c(i,:))
    hold on 
    if i == 1 || i == 2
       text(L(20).*10,eta_norm(20),strcat(num2str(1000.*P(i)), 'mW'))
    end
end

xlabel('Length (mm)','FontSize',16);
ylabel('Normalized efficiency (1/(cm^2W)','Fontsize',16);


title('Contour plot: Normalized efficiency and Length vs Power needed for 10 dB squeezing', 'Fontsize',12)