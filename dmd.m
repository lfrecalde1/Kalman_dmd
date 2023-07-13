%% Identificacion
clc, clear all, close all;

%% Lectura de datos
load('lines_features_kalman_3.mat');

% Sample time
ts = 0.03;

% End sample
n = 100;

% World velocities
vel_w = odom(8:13, :);
quaternions = odom(4:7, :);
% Body Velocities
for k = 1:length(vel_w)
    R_wb = quat2rotm(quaternions(:, k)');
    R_total = [R_wb, zeros(3,3);...
               zeros(3,3), R_wb];
           
    vel_b(:, k) = R_total'*vel_w(:, k);
end

% Control values
U = vel_b(:, 1:end-n);
U_r = vel_b(:, 1:end-n);
% Features lines
X = [f1_lines_inter(:, 1:end-n)];
X(3, :) = X(3, :)/313.56;
X(4, :) = X(4, :)/332.69;

X_r = [f1_lines_inter(:, 1:end-n)];
X_r(3, :) = X_r(3, :)/313.56;
X_r(4, :) = X_r(4, :)/332.69;
% Time generation
t(1) = 0;
for k = 1:length(X)-1
   t(k+1) = t(k) + ts; 
end

landa = 20;%
F1 = tf(landa,[1 landa]);

x1f = lsim(F1,X(1, :),t)';
x2f = lsim(F1,X(2, :),t)';
x3f = lsim(F1,X(3, :),t)';
x4f = lsim(F1,X(4, :),t)';

u1f = lsim(F1,U(1, :),t)';
u2f = lsim(F1,U(2, :),t)';
u3f = lsim(F1,U(3, :),t)';
u4f = lsim(F1,U(4, :),t)';
u5f = lsim(F1,U(5, :),t)';
u6f = lsim(F1,U(6, :),t)';



X = [x1f;...
     x2f;...
     x3f;...
     x4f];

U = [u1f; u2f; u3f; u4f; u5f; 1*u6f];
% % Solver
tic 
[sysmodel_DMDc] = Ident(X,U,ts);
toc
% 
%Initial Conditions
v_estimate(:, 1) = [X_r(:,1)];
S = 1200;
for k= 1:S
    %     f = (A*v_estimate(:,k)+B*vref(:,k));
    %     v_estimate(:, k+1) = v_estimate(:, k) + f*ts;
    v_estimate(:, k+1) = sysmodel_DMDc.A*v_estimate(:, k)+sysmodel_DMDc.B*U_r(:,k);
end
Ad = sysmodel_DMDc.A;
Bd = sysmodel_DMDc.B;
save("Matrcies.mat", 'Ad', 'Bd');
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
plot(vel_w(6, :),'-','Color',C11,'LineWidth',lw); hold on;
plot(vel_w(1, :),'-','Color',C12,'LineWidth',lw); hold on;
grid minor;
grid minor;
set(gca,'ticklabelinterpreter','latex',...
        'fontsize',fontsizeTicks)
%xlabel('$\textrm{Time}[s]$','interpreter','latex','fontsize',fontsizeLabel)
% ylabel('$[m/s]$','interpreter','latex','fontsize',fontsizeLabel)
% set(gca,'Xticklabel',[])
legend({'$\omega_w$','$\mu_{lw}$'},'interpreter','latex','fontsize',fontsizeLegend)
%plot(t,ul,'-','Color',C2,'LineWidth',lw);
grid minor;

figure('Position', [10 10 sizeX sizeY])
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8.5 11]);
set(gcf, 'PaperPositionMode', 'manual');
plot(vel_b(6, :),'-','Color',C11,'LineWidth',lw); hold on;
plot(vel_b(1, :),'-','Color',C12,'LineWidth',lw); hold on;

grid minor;
grid minor;
set(gca,'ticklabelinterpreter','latex',...
        'fontsize',fontsizeTicks)
%xlabel('$\textrm{Time}[s]$','interpreter','latex','fontsize',fontsizeLabel)
% ylabel('$[m/s]$','interpreter','latex','fontsize',fontsizeLabel)
% set(gca,'Xticklabel',[])
legend({'$\omega_b$','$\mu_{lb}$'},'interpreter','latex','fontsize',fontsizeLegend)
%plot(t,ul,'-','Color',C2,'LineWidth',lw);
grid minor;


figure('Position', [10 10 sizeX sizeY])
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8.5 11]);
set(gcf, 'PaperPositionMode', 'manual');
plot(X_r(3,1:S)','-','LineWidth',lw); hold on;



plot(v_estimate(3, 1:S)' ,'--','LineWidth',lw); hold on;


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
plot(X_r(4,1:S)','-','LineWidth',lw); hold on;



plot(v_estimate(4, 1:S)' ,'--','LineWidth',lw); hold on;


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
plot(X_r(1,1:S)','-','LineWidth',lw); hold on;
plot(X_r(2,1:S)','-','LineWidth',lw); hold on;



plot(v_estimate(1, 1:S)' ,'--','LineWidth',lw); hold on;
plot(v_estimate(2, 1:S)' ,'--','LineWidth',lw); hold on;

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
