% MAE 158 Project   Bryce Mankovsky 28714007
clear all
clc
%airplane charachteristic
b = 93;      %ft
cr = 20;     %ft
taper = 0.27;
sweep = 33;  %degrees from 1/4 chord line
Lf = 96;     %ft
Df = 10;     %ft
W = 98000;   %lb
V = linspace(200, 900, 100);     %ft/s

%atmospheric conditions
rho = 0.0008754; %slugs / ft^3
T = 400;         % Degrees R
mu = 3.025*10^-7; % lb*s/ft^2

crex = 18.43;   %exposed root chord from solidworks
mac = (crex + cr * taper) / 2;
Swet_wing = 2.0654e+03; %from solidworks
Sref = 1181.1; %Measured in solidworks   ft^2
AR = b^2 / Sref;
t2c = .1339;    %thickness to chord ratio from solidworks
Swetf = 0.8 * pi * Df * Lf;

ReOverD = rho * V / mu;
ReWing = ReOverD * mac;
ReFus = ReOverD * Lf;
ReHtail = ReOverD * 8.07;
ReVtail = ReOverD * 14.01;
RePy = ReOverD * 16.2;
ReN = ReOverD * 16.8;

function Cdp = Cdp_1(Re, t2c, sweep, Swet, Sref, V, T)    %Cdp using t/c
%load('TC_vs_K.mat')

Cf = 0.2289 * Re.^-0.2878 + 0.0011;
Ma = V / sqrt(1.4 * 1716 * T);
Z = ((2-Ma.^2)*cosd(sweep)) / sqrt(1- Ma .^ 2 * cosd(sweep)^2);
K = (1 + Z*t2c + 100 * t2c ^ 4);
%K = TC_vs_K(t2c, sweep)
Cdp = Cf * K * Swet / Sref;
end

function Cdp = Cdp_2(Re, FR, Swet, Sref)    %Cdp using Fineness ratio
Cf = 0.2289 * Re .^ -0.2878 + 0.0011;
K = 0.0003 * FR ^ 4 - 0.0113 * FR ^ 3 + 0.1388 * FR ^ 2 - 0.7908 * FR + 2.9669; % Eqn from curve fit
Cdp = Cf * K * Swet / Sref;
end

Cdp_wing = Cdp_1(ReWing, t2c, sweep,  Swet_wing, Sref, V, T);
Cdp_f = Cdp_2(ReFus, Lf / Df, Swetf, Sref);
Cdp_ht = Cdp_1(ReHtail, 0.09, 31.6, 261 * 2 * 1.02, Sref, V, T);
Cdp_vt = Cdp_1(ReVtail, 0.09, 43.5, 161 * 2 * 1.02, Sref, V, T);
Cdp_p = Cdp_1(RePy, 0.06, 10, 117 * 2 * 1.02, Sref, V, T);   %using sweep = 10 because it is assumed equal to 0 on the chart
Cdp_n = Cdp_2(ReN, 5, 455, Sref);

Cdp = Cdp_wing + Cdp_f + Cdp_ht + Cdp_vt + Cdp_p + Cdp_n;


%Cd induced
load('AR_vs_e.mat')

e = NaN(size(V));
Cdi = NaN(size(V));
for i = 1:length(V)
    e(i) = AR_vs_e(AR, Cdp(i));  % Store e for each velocity V
end
CL = W ./ (0.5 * rho * V .^2 * Sref);
for i = 1:length(V)
    Cdi(i) = CL(i)^2 / (pi * AR * e(i));
end

Dp = 0.5 * rho * V .^ 2 .* Cdp * Sref;
Di = 0.5 * rho * V .^ 2 .* Cdi * Sref;
D = Dp + Di;

% Plot
figure;
plot(V, Dp, 'r-', 'LineWidth', 2); 
hold on;
plot(V, Di, 'b-', 'LineWidth', 2); 
plot(V, D, 'k-', 'LineWidth', 2); 

legend('Parasite Drag (Dp)', 'Induced Drag (Di)', 'Total Drag (D)', 'Location', 'northwest');

title('Drag vs Velocity');
xlabel('Velocity (ft/s)');
ylabel('Drag (lb)');

hold off;

figure;
plot(V, W./D);
title('Lift to Drag vs Velocity');
xlabel('Velocity (ft/s)');
ylabel('L/D Ratio');

%at V = 765
ReOverD765 = rho * 765 / mu;
ReWing765 = ReOverD765 * mac;
ReFus765 = ReOverD765 * Lf;
ReHtail765 = ReOverD765 * 8.07;
ReVtail765 = ReOverD765 * 14.01;
RePy765 = ReOverD765 * 16.2;
ReN765 = ReOverD765 * 16.8;
Cdp_wing765 = Cdp_1(ReWing765, t2c, sweep,  Swet_wing, Sref, 765, T);
Cdp_f765 = Cdp_2(ReFus765, Lf / Df, 0.8 * pi * Df * Lf, Sref);
Cdp_ht765 = Cdp_1(ReHtail765, 0.09, 31.6, 261 * 2 * 1.02, Sref, 765, T);
Cdp_vt765 = Cdp_1(ReVtail765, 0.09, 43.5, 161 * 2 * 1.02, Sref, 765, T);
Cdp_p765 = Cdp_1(RePy765, 0.06, 10, 117 * 2 * 1.02, Sref, 765, T);   %using sweep = 10 because it is assumed equal to 0 on the chart
Cdp_n765 = Cdp_2(ReN765, 5, 455, Sref);
Cdp765 = Cdp_wing765 + Cdp_f765 + Cdp_ht765 + Cdp_vt765 + Cdp_p765 + Cdp_n765;

e765 = AR_vs_e(AR, Cdp765);
CL765 = W / (0.5 * rho * 765^2 * Sref);
Cdi765 = CL765^2 / (pi * AR * e765);
Dp765 = 0.5 * rho * 765^2 * Cdp765 * Sref;
Di765 = 0.5 * rho * 765^2 * Cdi765 * Sref;
D765 = Dp765 + Di765
LiftToDrag765 = W / D765