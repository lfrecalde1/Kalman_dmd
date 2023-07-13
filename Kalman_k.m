%% Kalman Filter matlab
clc, clear all, close all;

%% Load Data system
load('lines_features_kalman_3.mat');
load('Matrcies.mat');

Cd = 1*eye(4);
%% Samples 
n = 100;

%% World velocities
vel_w = odom(8:13, :);
quaternions = odom(4:7, :);

%% Body Velocities
for k = 1:length(vel_w)
    R_wb = quat2rotm(quaternions(:, k)');
    R_total = [R_wb, zeros(3,3);...
               zeros(3,3), R_wb];  
    vel_b(:, k) = R_total'*vel_w(:, k);
end

%% Control Values
U = vel_b(:, 1:end-n);
U_r = vel_b(:, 1:end-n);

%% Real states of the system
X = [f1_lines_inter(:, 1:end-n)];
X(3, :) = X(3, :)/313.56;
X(4, :) = X(4, :)/332.69;

X_r = [f1_lines_inter(:, 1:end-n)];
X_r(3, :) = X_r(3, :)/313.56;
X_r(4, :) = X_r(4, :)/332.69;
%% Time Definition
t(1) = 0;
ts = 0.05; 
for k = 1:length(X_r)-1
   t(k+1) = t(k) + ts; 
   
end

%% Initial Conditions
v_estimate_1 = [X_r(:, 1)];
%% Eps valu
eps_n = 0.5;
eps_d = 0.001;%% Simulation

%% Initial States Kalman
x_estimate_1 =  v_estimate_1;
P_k = zeros(4, 4, length(t));
P_k_1 = 1*eye(4);


Q = 0.01*eye(4);
R = 30*eye(4);
for k= 1:length(t)
   
    %% System Values
    v_estimate(:, k) = Ad*v_estimate_1 + Bd*U_r(:, k) + eps_d*randn(4, 1);
    y_estimate(:, k) = Cd*v_estimate(:, k)+ eps_n*randn(4, 1);
    
    %% Kalman Filter
    x_es = Ad*x_estimate_1 + Bd* U_r(:,k);
    P = Ad*P_k_1*Ad' + Q;
    
    %% Estimation
    Kk = P*Cd'*inv(Cd*P*Cd' + R);
    
    x_estimate(:, k) = x_es + Kk*(y_estimate(:, k)-Cd*x_es);
    P_k(:, :, k) = (eye(4) - Kk*Cd)*P;
    P_d(:, k) = norm(P_k(:, :, k));
    
    v_estimate_1 = v_estimate(:, k);
    x_estimate_1 = x_estimate(:, k);
    P_k_1 = P_k(:, :, k);
end

%% Parameters fancy plots
lw = 2; 
lwV = 2; 
fontsizeLabel = 9;
fontsizeLegend = 9;
fontsizeTicks = 9;
fontsizeTitel = 9;
sizeX = 900;
sizeY = 300;

%% color propreties
C1 = [246 170 141]/255;
C2 = [51 187 238]/255;
C3 = [0 153 136]/255;
C4 = [238 119 51]/255;
C5 = [204 51 17]/255;
C6 = [238 51 119]/255;
C7 = [187 187 187]/255;
C8 = [80 80 80]/255;
C9 = [140 140 140]/255;
C10 = [0 128 255]/255;
C11 = [234 52 89]/255;
C12 = [39 124 252]/255;
C13 = [40 122 125]/255;
C14 = [252 94 158]/255;
C15 = [244 171 39]/255;
C16 = [100 121 162]/255;
C17 = [255 0 0]/255;

figure('Position', [10 10 sizeX sizeY])
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8.5 11]);
set(gcf, 'PaperPositionMode', 'manual');
plot(v_estimate(1, :)' ,'-','LineWidth',lw); hold on;
plot(v_estimate(2, :)' ,'-','LineWidth',lw); hold on;
plot(y_estimate(1, :)' ,'--','LineWidth',lw); hold on;
plot(y_estimate(2, :)' ,'--','LineWidth',lw); hold on;
plot(x_estimate(1, :)' ,'--','LineWidth',lw); hold on;
plot(x_estimate(2, :)' ,'--','LineWidth',lw); hold on;

grid minor;
grid minor;
set(gca,'ticklabelinterpreter','latex',...
        'fontsize',fontsizeTicks)
%xlabel('$\textrm{Time}[s]$','interpreter','latex','fontsize',fontsizeLabel)
% ylabel('$[m/s]$','interpreter','latex','fontsize',fontsizeLabel)
% set(gca,'Xticklabel',[])
legend({'$l_1$','$l_e$'},'interpreter','latex','fontsize',fontsizeLegend)
%plot(t,ul,'-','Color',C2,'LineWidth',lw);
grid minor;
figure('Position', [10 10 sizeX sizeY])
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8.5 11]);
set(gcf, 'PaperPositionMode', 'manual');
plot(P_d(1, :) ,'-','LineWidth',lw); hold on;


grid minor;
grid minor;
set(gca,'ticklabelinterpreter','latex',...
        'fontsize',fontsizeTicks)
%xlabel('$\textrm{Time}[s]$','interpreter','latex','fontsize',fontsizeLabel)
% ylabel('$[m/s]$','interpreter','latex','fontsize',fontsizeLabel)
% set(gca,'Xticklabel',[])
legend({'$P_d$'},'interpreter','latex','fontsize',fontsizeLegend)
%plot(t,ul,'-','Color',C2,'LineWidth',lw);
grid minor;

figure('Position', [10 10 sizeX sizeY])
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8.5 11]);
set(gcf, 'PaperPositionMode', 'manual');
plot(v_estimate(3, :)' ,'-','LineWidth',lw); hold on;
plot(v_estimate(4, :)' ,'-','LineWidth',lw); hold on;
plot(y_estimate(3, :)' ,'--','LineWidth',lw); hold on;
plot(y_estimate(4, :)' ,'--','LineWidth',lw); hold on;
plot(x_estimate(3, :)' ,'--','LineWidth',lw); hold on;
plot(x_estimate(4, :)' ,'--','LineWidth',lw); hold on;

grid minor;
grid minor;
set(gca,'ticklabelinterpreter','latex',...
        'fontsize',fontsizeTicks)
%xlabel('$\textrm{Time}[s]$','interpreter','latex','fontsize',fontsizeLabel)
% ylabel('$[m/s]$','interpreter','latex','fontsize',fontsizeLabel)
% set(gca,'Xticklabel',[])
legend({'$l_1$','$l_e$'},'interpreter','latex','fontsize',fontsizeLegend)
%plot(t,ul,'-','Color',C2,'LineWidth',lw);
grid minor;