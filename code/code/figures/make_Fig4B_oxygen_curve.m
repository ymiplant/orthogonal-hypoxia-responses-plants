% make_Fig4B_oxygen_curve.m
%
% Steady-state induction of a generic hypoxia-responsive gene (HRG)
% across a range of oxygen concentrations.
%
% This script reproduces the trends shown in Figure 4B.

clear; close all; clc

%% -------------------------
% Shared parameters
%% -------------------------
k1  = 0.005/1000;
k2  = 0.0002;
k4  = 0.0016/1000;
k5  = 0.00038;
KM3 = 0.05;

Enz = 0.1;   % enzyme amount (PCO4 or PHD2)

opts = odeset('RelTol',1e-8,'AbsTol',1e-10);

%% -------------------------
% Enzyme-specific parameters
%% -------------------------
% PCO4
PCO4.k3  = 26.5878;
PCO4.KM1 = 16.4217;
PCO4.KM2 = 270.2505;

% PHD2
PHD2.k3  = 0.1185;
PHD2.KM1 = 27.3944;
PHD2.KM2 = 7.8463;

%% -------------------------
% Oxygen range
%% -------------------------
O2_max  = 21;
O2_min  = 0.1;
nSteps  = 210;

O2_grid = linspace(O2_max, O2_min, nSteps);

tSS = 6 * 3600;   % 6 hours to reach steady state

mRNA_PCO4 = zeros(size(O2_grid));
mRNA_PHD2 = zeros(size(O2_grid));

%% -------------------------
% Start from steady state at 21% O2
%% -------------------------
a_PCO4 = steady_state(O2_max,[0;0],PCO4,k1,k2,k4,k5,KM3,Enz,tSS,opts);
a_PHD2 = steady_state(O2_max,[0;0],PHD2,k1,k2,k4,k5,KM3,Enz,tSS,opts);

%% -------------------------
% Oxygen sweep
%% -------------------------
for i = 1:numel(O2_grid)
    O2 = O2_grid(i);

    a_PCO4 = steady_state(O2,a_PCO4,PCO4,k1,k2,k4,k5,KM3,Enz,tSS,opts);
    mRNA_PCO4(i) = a_PCO4(2);

    a_PHD2 = steady_state(O2,a_PHD2,PHD2,k1,k2,k4,k5,KM3,Enz,tSS,opts);
    mRNA_PHD2(i) = a_PHD2(2);
end

%% -------------------------
% Plot
%% -------------------------
figure;
plot(O2_grid, mRNA_PCO4, 'b-', 'LineWidth', 2); hold on
plot(O2_grid, mRNA_PHD2, 'r-', 'LineWidth', 2);

set(gca,'XDir','reverse');
xlim([O2_min O2_max]);
set(gca,'XTick',[0.1 1 2 5 10 15 21]);

xlabel('O_2 (%)','FontSize',14);
ylabel('Simulated HRG expression (µM)','FontSize',14);
legend({'PCO4','PHD2'},'Location','northwest');
set(gca,'FontSize',14);
grid on;

%% -------------------------
% Local functions
%% -------------------------
function a_ss = steady_state(O2,a0,enz,k1,k2,k4,k5,KM3,Enz,tSS,opts)
    ode = @(t,a) rap_mrna_ode(a,O2,enz,k1,k2,k4,k5,KM3,Enz);
    [~,A] = ode45(ode,[0 tSS],a0,opts);
    a_ss = A(end,:).';
end

function da = rap_mrna_ode(a,O2,enz,k1,k2,k4,k5,KM3,Enz)
    RAP  = a(1);
    mRNA = a(2);

    v_ox = enz.k3 * Enz * (O2/(enz.KM1+O2)) * (RAP/(enz.KM2+RAP));

    dRAP  = -v_ox + k1 - k2*RAP;
    dmRNA =  k4*RAP/(KM3+RAP) - k5*mRNA;

    da = [dRAP; dmRNA];
end
