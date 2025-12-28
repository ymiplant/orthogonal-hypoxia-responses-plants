clear; close all; clc

% =========================
% Data vectors for CODD / PHD2
% =========================
CODD = [0 11.18 20.76 40.72 81.24 100 150 200 100*ones(1,9)];
Ox   = [20*ones(1,8), 0 10 20 30 40 50 60 70 80];
vo   = [0 .09 .1596 .18582 .20786 .21774 .20786 .17898 ...
        0 .0966 .1725 .21252 .23828 .27508 .29808 .32752 .3335];

% xdata as 2-by-N for nlinfit
xdata = [CODD; Ox];

% =========================
% Nonlinear fit using nlinfit
% Model: v = Vmax * CODD/(KmCODD+CODD) * Ox/(KmO2+Ox)
% =========================
opts = statset('nlinfit');
opts.RobustWgtFun = 'fair';

model = @(b,x) b(1) .* (x(1,:)./(b(2)+x(1,:))) .* (x(2,:)./(b(3)+x(2,:)));
initialguess = [1 1 1];

[beta,~,~,~,MSE] = nlinfit(xdata, vo, model, initialguess, opts);

MSE
Kcat   = beta(1)/4      % /s https://doi.org/10.1042/BJ20140779
KmCODD = beta(2)
KmO2   = beta(3)

% define explicitly
Vmax  = beta(1);
KmS   = beta(2);
KmO2  = beta(3);

% =========================
% Create fitted surface
% =========================
ERF = linspace(0,200);
OX  = linspace(0,80);

MODEL = zeros(length(ERF), length(OX));
for i = 1:length(ERF)
    for j = 1:length(OX)
        MODEL(i,j) = Vmax * ERF(i)/(KmS+ERF(i)) * OX(j)/(KmO2+OX(j));
    end
end

figure('Position', [100, 100, 800, 600]);
surf(ERF, OX, MODEL, 'EdgeColor','none', 'DisplayName', 'Fitting')
alpha 0.5
xlabel('CODD Concentration (µM)','FontSize',14)
ylabel('Oxygen (%)','FontSize',14)
zlabel('PHD2 activity (µM/s)','FontSize',14)
title('PHD2','FontSize',16);
set(gca, 'FontSize', 14);
colormap(gray)
shading flat
grid on
hold on

% =========================
% Overlay original slice data
% =========================
CODD1 = [0 11.18 20.76 40.72 81.24 100 150 200];
vo1   = [0 .09 .1596 .18582 .20786 .21774 .20786 .17898];

Ox1 = [0 10 20 30 40 50 60 70 80];
vo2 = [0 .0966 .1725 .21252 .23828 .27508 .29808 .32752 .3335];

line1 = 20*ones(size(CODD1));   % O2 fixed at 20%
line2 = 100*ones(size(Ox1));    % CODD fixed at 100 uM

plot3(CODD1, line1, vo1, 'r:o', 'LineWidth', 1.5, 'DisplayName', 'CODD data point');
plot3(line2, Ox1, vo2, 'b:o', 'LineWidth', 1.5, 'DisplayName', 'Oxygen data point');

% Error bars
errlCODD = [0 .07 .025 .025 .035 .08 .08 0.02]/2;
errhCODD = errlCODD;

errlOx = [0 .01 .03 .025 .025 .05 .07 .07 0.03]/2;
errhOx = errlOx;

plot3([CODD1',CODD1']', [line1',line1']', [-errlCODD', errhCODD']'+vo1, '-r', 'LineWidth', 1.5)
plot3([line2',line2']', [Ox1',Ox1']',     [-errlOx',   errhOx']'+vo2,  '-b', 'LineWidth', 1.5)

