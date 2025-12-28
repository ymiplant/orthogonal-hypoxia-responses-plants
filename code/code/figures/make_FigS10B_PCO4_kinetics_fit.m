clear; close all; clc

% =========================
% Data vectors for ERFVII (RAP2.12) / PCO4
% =========================
ERFVII = [15.625 31.25 62.5 125 250 500 1000 2000 4000 1000*ones(1,8)];
Ox     = [20*ones(1,9), 0 1.5 3 5.5 11 20 40 60];
vo     = [1.470724 2.709730 4.716129 9.452704 15.811300 24.016820 25.741070 21.844040 16.023140 ...
          1.6363 3.2727 6.5454 10.5273 19.6363 25.0909 31.6363 36.3273];

% xdata as 2-by-N for nlinfit
xdata = [ERFVII; Ox];

% =========================
% Nonlinear fit using nlinfit
% Model: v = Vmax * ERFVII/(KmRAP+ERFVII) * Ox/(KmO2+Ox)
% =========================
opts = statset('nlinfit');
opts.RobustWgtFun = 'Welsch';

model = @(b,x) b(1) .* (x(1,:)./(b(2)+x(1,:))) .* (x(2,:)./(b(3)+x(2,:)));
initialguess = [1 1 1];

[beta,~,~,~,MSE] = nlinfit(xdata, vo, model, initialguess, opts);

MSE
Kcat  = beta(1)/60/(1000/27311.7)  % /s, MW PCO4 = 27311.7
Kmrap = beta(2)
KmO2  = beta(3)

% define explicitly
Vmax  = beta(1);
KmRAP = beta(2);
KmO2  = beta(3);

% =========================
% Create fitted surface
% =========================
ERF = linspace(0,4500);
OX  = linspace(0,60);

MODEL = zeros(length(ERF), length(OX));
for i = 1:length(ERF)
    for j = 1:length(OX)
        MODEL(i,j) = Vmax * ERF(i)/(KmRAP+ERF(i)) * OX(j)/(KmO2+OX(j));
    end
end

figure('Position', [100, 100, 800, 600]);
surf(ERF, OX, MODEL, 'EdgeColor','none', 'DisplayName','Fitting')
alpha 0.5
xlabel('RAP2.12 Concentration (µM)','FontSize',14)
ylabel('Oxygen (%)','FontSize',14)
zlabel('PCO4 activity (µmoles/min/mg)','FontSize',14)
title('PCO4','FontSize',16);
set(gca, 'FontSize', 14);
colormap(gray)
shading flat
grid on
hold on

% =========================
% Overlay original slice data 
% =========================
ERFVII1 = [15.625 31.25 62.5 125 250 500 1000 2000 4000];
vo1     = [1.470724 2.709730 4.716129 9.452704 15.811300 24.016820 25.741070 21.844040 16.023140];

Ox1 = [0 1.5 3 5.5 11 20 40 60];
vo2 = [1.6363 3.2727 6.5454 10.5273 19.6363 25.0909 31.6363 36.3273];

line1 = 20*ones(size(ERFVII1));     % O2 fixed at 20%
line2 = 1250*ones(size(Ox1));       % RAP fixed at 1250 uM 

plot3(ERFVII1, line1, vo1, 'r:o', 'LineWidth', 1.5, 'DisplayName','ERFVII data point');
plot3(line2, Ox1, vo2, 'b:o', 'LineWidth', 1.5, 'DisplayName','Oxygen data point');

% Error bars
errlRAP = [.1 .1 .1 1.3228 2.067 3.307 3.2244 1.86 2.067]/2;
errhRAP = errlRAP;

errlOx  = [.2 .2 .2 .2 5.244 .1 4.748 3.189]/2;
errhOx  = errlOx;

plot3([ERFVII1',ERFVII1']', [line1',line1']', [-errlRAP', errhRAP']'+vo1, '-r','LineWidth', 1.5)
plot3([line2',line2']',     [Ox1',Ox1']',     [-errlOx',  errhOx']'+vo2,  '-b','LineWidth', 1.5)
