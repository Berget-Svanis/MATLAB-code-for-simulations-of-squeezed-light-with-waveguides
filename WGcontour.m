%Contour waveguide, loss and phasenoise vs squeezing

%For ppLN
eta_norm = 0.4e4;
L = 45e-3;

%Defining squeezing, however we since we only care about the minimum of the
%function we can remove constants (eta)
sqz = @(A_0) exp(-2.*L.*sqrt(eta_norm.*A_0).*tanh(L.*sqrt(eta_norm.*A_0)));
asqz = @(A_0) exp(2.*L.*sqrt(eta_norm.*A_0).*tanh(L.*sqrt(eta_norm.*A_0)));

%Pow to dB and vice versa
dB_to_pow = @(d) 10.^(d./10);
pow_to_dB = @(p) 10.*log10(p);

%The squeezing levels we are looking at
dB = linspace(-3,-18,16);
var = dB_to_pow(dB);

num_points = 1000;
millirad = linspace(0.01,40,num_points);

cmap = flip(autumn(19),1); % yellow -> red, with  19 colors (for 19 lines)
set(gca(),'ColorOrder',cmap)

num = numel(var);
c = parula(num); %Sets color gradient

for i =1:num
    rad = millirad./1000;
    %We need to find the optimal power A_0, we do this by finding the
    %power that maximizes squeezing for a given phase noise. I.e the
    %power that the variance has a minimum at. 
    A_0 = zeros(1,num_points);
    for j=1:num_points
        sqz_min = @(A_0) sqz(A_0).*cos(rad(j)).^2 + asqz(A_0).*sin(rad(j)).^2;
        A_0(j) = fminsearch(sqz_min,0.1);
    end
    
    %Eta as a function of variance and phase noise angle
    eta = @(v_out,o) (v_out-1)./(exp(-2.*L.*sqrt(eta_norm.*A_0).*tanh(L.*sqrt(eta_norm.*A_0))).*cos(o).^2 + exp(2.*L.*sqrt(eta_norm.*A_0).*tanh(L.*sqrt(eta_norm.*A_0))).*sin(o).^2 - 1);
    
    %Plotting
    loss = 100.*(1-eta(var(i),rad));
    semilogy(millirad,loss,'Color',c(i,:))
    if mod(i,2) == 0
       text(millirad(20),loss(20),strcat(num2str(10.*log10(var(i))), 'dB'))
    end
    hold on
end

xlabel('RMS (mrad)','FontSize',16);
ylabel('Loss (%)','Fontsize',16);
xlim([0 40])
set(gca, 'ylim', [1 60]);
set(gca, 'ytick', [1 5 10 20 30 40 50 60]);

title('Contour, loss and phase noise vs squeezing','Fontsize',16)
hold off