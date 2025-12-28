% make_Fig4A_timecourse.m
%
% Simulated induction of a generic hypoxia-responsive gene (HRG)
% by RAP2.12 under control of either the plant PCO oxygen-sensing
% pathway or a synthetic PHD-based sensing pathway.
%
% This script reproduces the trends shown in Figure 4A.

clear; close all; clc

%% -------------------------
% Shared parameters
%% -------------------------
k1  = 0.005/1000;
k2  = 0.0002;
k4  = 0.0016/1000;
k5  = 0.00038;
KM3 = 0.05;

Enz = 0.1;        % Enzyme amount (PCO4 or PHD2)

%% -------------------------
% Pre-equilibration at normoxia (21% O2)
%% -------------------------
O2 = 21;
t_pre = 6 * 3600;     % 6 hours

% ---- PCO4 parameters ----
k3_PCO4  = 26.5878;
KM1_PCO4 = 16.4217;
KM2_PCO4 = 270.2505;

ode_PCO4_pre = @(t,a) rap_mrna_ode(a,O2,k1,k2,k3_PCO4,Enz,KM1_PCO4,KM2_PCO4,k4,k5,KM3);
[~,Apre_PCO4] = ode45(ode_PCO4_pre,[0 t_pre],[0 0]);
a0_PCO4 = Apre_PCO4(end,:);

% ---- PHD2 parameters ----
k3_PHD2  = 0.1185;
KM1_PHD2 = 27.3944;
KM2_PHD2 = 7.8463;

ode_PHD2_pre = @(t,a) rap_mrna_ode(a,O2,k1,k2,k3_PHD2,Enz,KM1_PHD2,KM2_PHD2,k4,k5,KM3);
[~,Apre_PHD2] = ode45(ode_PHD2_pre,[0 t_pre],[0 0]);
a0_PHD2 = Apre_PHD2(end,:);

%% -------------------------
% Hypoxia step (1% O2)
%% -------------------------
O2 = 1;
t_step = 4 * 3600;    % 4 hours

% PCO4 hypoxia
ode_PCO4_step = @(t,a) rap_mrna_ode(a,O2,k1,k2,k3_PCO4,Enz,KM1_PCO4,KM2_PCO4,k4,k5,KM3);
[t_PCO4, A_PCO4] = ode45(ode_PCO4_step,[0 t_step],a0_PCO4);

% PHD2 hypoxia
ode_PHD2_step = @(t,a) rap_mrna_ode(a,O2,k1,k2,k3_PHD2,Enz,KM1_PHD2,KM2_PHD2,k4,k5,KM3);
[t_PHD2, A_PHD2] = ode45(ode_PHD2_step,[0 t_step],a0_PHD2);

%% -------------------------
% Plot
%% -------------------------
figure;
plot(t_PCO4/3600, A_PCO4(:,2), 'b-', 'LineWidth', 2); hold on
plot(t_PHD2/3600, A_PHD2(:,2), 'r-', 'LineWidth', 2);

xlabel('Time (h)','FontSize',14);
ylabel('Simulated HRG expression (µM)','FontSize',14);
legend({'PCO4','PHD2'},'Location','southeast');
set(gca,'FontSize',14);
grid on;

%% -------------------------
% Local function
%% -------------------------
function da = rap_mrna_ode(a,O2,k1,k2,k3,Enz,KM1,KM2,k4,k5,KM3)
RAP  = a(1);
mRNA = a(2);

v_ox = k3 * Enz * (O2/(KM1+O2)) * (RAP/(KM2+RAP));

dRAP  = -v_ox + k1 - k2*RAP;
dmRNA =  k4*RAP/(KM3+RAP) - k5*mRNA;

da = [dRAP; dmRNA];
end
