function [s, as]=sqzWG(L ,P_in, eta_norm, eta_in, eta_out, eta_read, dB, rad)
%This function gives you squeezing and antisqueezing (dB) in a waveguide 
%for a certain efficiency.
%The efficiency can be incoupling, outcoupling,readout or propagation loss (in dB/m)
%Normalized efficiency (1/(W*m^2) and phase noise (rad) are also required. 
%Length (m) and input power (W) are also required input parameters.

a = dB.*log(10)./10; %Converting to propagation loss in 1/m

P_eff = P_in.*eta_in; %The effective power in

eta_tot = eta_out.*exp(-a.*L).*eta_read; %The total efficiency 

%Defining squeezing and anti-squeezing
sqz = @(P,n) n.*exp(-2.*L.*sqrt(eta_norm.*P).*tanh(L.*sqrt(eta_norm.*P))) + (1-n);
asqz = @(P,n) n.*exp(2.*L.*sqrt(eta_norm.*P).*tanh(L.*sqrt(eta_norm.*P))) + (1-n);
sqz_rad = @(P,n,rad) sqz(P,n).*cos(rad).^2 + asqz(P,n).*sin(rad).^2;
asqz_rad = @(P,n,rad) asqz(P,n).*cos(rad).^2 + sqz(P,n).*sin(rad).^2;

s = 10.*log10(sqz_rad(P_eff,eta_tot,rad));
as = 10.*log10(asqz_rad(P_eff,eta_tot,rad)); 
end