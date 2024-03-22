%% -----------------------------------------------------
% Name    : MicroStripAntcal.m
% Purpose : 计算侧馈微带贴片天线相关的参数
% Model: 侧馈微带贴片天线
% Author: Yongwei Li / liyongwei_101@qq.com
% Date:   March 2024
%% -----------------------------------------------------

clear;
clear all;    %清空工作区与全局变量
%clc;          %清空命令区域

f=2.4*10^9; c=3*10^8; Omega=2*pi*f;
%**************自由空间中参数***************%
lambda_free_24=c/f; % 自由空间中的波长
epsilon_0=1/(pi*36e9); miu_0=pi*4e-7; % 自由空间中的 ε&μ

%**************介质板中的参数***************%
epsilon_fr4=4.4; efr4=epsilon_fr4*epsilon_0; % fr4的相对介电常数
miu_fr4=1;       mfr4=miu_fr4*miu_0; %fr4的相对磁导率
lambda_fr4_24=c/(f*sqrt(epsilon_fr4*miu_fr4)); %Fr4中的导波波长
k_fr4=Omega*sqrt(mfr4*efr4); % fr4中的波数
beta_fr4=Omega*sqrt(mfr4*efr4); %介质中的相位常数 

%% -----------------------------------------------------
%**************辐射贴片的参数***************%

h=1.6*10^-3;%基板厚度unit：mm

Z_0=6.72;%贴片天线作为传输线时的特性阻抗
W_e=w(c,f,epsilon_fr4); %辐射单元宽度w
z=W_e/2; %z为馈电中心距离贴片边缘的距离
L_e=L(c,f,h,epsilon_fr4,miu_fr4); %有效单元的长度

delta_L_value=delta_L(c,f,h,epsilon_fr4);
epsilon_e_value=epsilon_e(c,f,h,epsilon_fr4);
Y_in_z=Y_in(c,f,epsilon_fr4,k_fr4,beta_fr4,h,Z_0,z); %输入导纳

Z_in_z=1/Y_in_z;

Dir_f=(w(c,f,epsilon_fr4)*pi/lambda_fr4_24)^2*8/I(c,f,epsilon_fr4,k_fr4);

%%------------------------------输入导纳------------------------------
function Yin_fun=Y_in(c,f,epsilon_fr4,k,beta,h,Z_0,z)
Yin_fun=2*G(c,f,epsilon_fr4,k)*(cos(beta*z).^2+Z_0^2*(G(c,f,epsilon_fr4,k).^2+B(c,f,h,epsilon_fr4,k,Z_0).^2)*sin(beta*z).^2-Z_0*B(c,f,h,epsilon_fr4,k,Z_0)*sin(2*beta*z)).^-1;
end

%%------------------------------辐射电导------------------------------
function G_fun=G(c,f,epsilon_fr4,k)
G_fun=I(c,f,epsilon_fr4,k)/(120*pi^2);
end

%%------------------------------I------------------------------
function I_fun=I(c,f,epsilon_fr4,k)
I_fun=integral(@(theta)sin(0.5*k.*w(c,f,epsilon_fr4).*cos(theta)).^2.*tan(theta).^2.*sin(theta), 0, pi);
end
%%------------------------------等效电纳B------------------------------
function B_fun=B(c,f,h,epsilon_fr4,k,Z_0)
B_fun=k*delta_L(c,f,h,epsilon_fr4)*sqrt(epsilon_e(c,f,h,epsilon_fr4))/Z_0;
end

%%------------------------------实际辐射单元长度L------------------------------
function L_fun=L(c,f,h,epsilon_fr4,miu_fr4)
L_fun=0.5*c/(f*sqrt(epsilon_e(c,f,h,epsilon_fr4)*miu_fr4))-2*delta_L(c,f,h,epsilon_fr4);
end
%%------------------------------delta_L，等效辐射缝隙长度------------------------------
function deltaL_fun=delta_L(c,f,h,epsilon_fr4)
deltaL_fun=0.412*h*(epsilon_e(c,f,h,epsilon_fr4)+0.3)*(w(c,f,epsilon_fr4)/h+0.264)/((epsilon_e(c,f,h,epsilon_fr4)-0.258)*(w(c,f,epsilon_fr4)/h+0.8));
end
%%------------------------------天线基板有效相对介电常数------------------------------
function epsilonE_fun=epsilon_e(c,f,h,epsilon_fr4)
epsilonE_fun=0.5*(epsilon_fr4+1)+0.5*(epsilon_fr4-1)*(1+12*h/w(c,f,epsilon_fr4))^-0.5;
end
%%------------------------------实际辐射单元宽度w,半波长------------------------------
function w_fun=w(c,f,epsilon_fr4)
w_fun=(((epsilon_fr4+1)/2)^-0.5)*c/(2*f);
end
