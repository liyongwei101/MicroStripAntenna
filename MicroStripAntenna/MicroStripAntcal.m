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
free_lambda_24=c/f; % 自由空间中的波长
free_epsilon_0=1/(pi*36e9); free_miu_0=pi*4e-7; % 自由空间中的 ε&μ

%**************介质板中的参数 FR4***************%
Drelative_epsilon_obj=4.4; Dobj_e=Drelative_epsilon_obj*free_epsilon_0;                  % 介质板 的相对介电常数
Drelative_miu_obj=1;       Dobj_m=Drelative_miu_obj*free_miu_0; %介质板 的相对磁导率
Dobj_lambda_24=c/(f*sqrt(Drelative_epsilon_obj*Drelative_miu_obj)); %介质板中的导波波长unit:m
Dobj_k=Omega*sqrt(Dobj_m*Dobj_e); % 介质板中的波数
Dobj_beta=Omega*sqrt(Dobj_m*Dobj_e); %介质中的相位常数

%% -----------------------------------------------------
%**************辐射贴片的参数***************%

h=1.6*10^-3;                                                 %基板厚度unit：m
PatchZ_0=6.72;                                               %贴片天线作为传输线时的特性阻抗
PatchR_y0=3*10^-3;                                             %凹槽长度 unit:m

PatchW_e=w(c,f,Drelative_epsilon_obj);%辐射单元宽度w unit:m
PatchL_e=L(c,f,h,Drelative_epsilon_obj,Drelative_miu_obj); %有效单元的长度L unit:m


Patch_feed_z=PatchW_e/(2*1000); %z为馈电中心距离贴片边缘的距离 unit:m
value_delta_L=delta_L(c,f,h,Drelative_epsilon_obj); %边缩效应 unit:m
value_epsilon_e=epsilon_e(c,f,h,Drelative_epsilon_obj); %等效介电常数
value_G=G(c,f,Drelative_epsilon_obj,Dobj_k);%辐射电导 G

%Output_Y_in_z=Y_in(c,f,fr4_epsilon,fr4_k,fr4_beta,h,PatchZ_0,z); %%%%输入导纳 有问题
Output_Y_in_z=2*value_G;%输入导纳
Output_Z_in_z=1/Output_Y_in_z;%输入阻抗
Output_Imp_transformer=sqrt(Output_Z_in_z*50); %1/4阻抗变换线的阻抗
Output_Imp_transformer_recessed=sqrt(Output_Z_in_z*cos(pi*(PatchR_y0/PatchL_e))^2); %%%有问题
Output_Dir_f=(w(c,f,Drelative_epsilon_obj)*pi/Dobj_lambda_24)^2*8/I(c,f,Drelative_epsilon_obj,Dobj_k);%方向性
%%--------------------------------------------------------------

%%------------------------------输入导纳------------------------------
function Yin_fun=Y_in(c,f,relative_epsilon_obj,k,beta,h,Z_0,z)
Yin_fun=2*G(c,f,relative_epsilon_obj,k)*(cos(beta*z).^2+Z_0^2*(G(c,f,relative_epsilon_obj,k).^2+B(c,f,h,relative_epsilon_obj,k,Z_0).^2)*sin(beta*z).^2-Z_0*B(c,f,h,relative_epsilon_obj,k,Z_0)*sin(2*beta*z)).^-1;
end

%%------------------------------辐射电导 G------------------------------
function G_fun=G(c,f,relative_epsilon_obj,k)
G_fun=I(c,f,relative_epsilon_obj,k)/(120*pi^2);
end

%%------------------------------I------------------------------
function I_fun=I(c,f,relative_epsilon_obj,k)
I_fun=integral(@(theta)sin(0.5*k.*w(c,f,relative_epsilon_obj).*cos(theta)).^2.*tan(theta).^2.*sin(theta), 0, pi);
end
%%------------------------------等效电纳B------------------------------
function B_fun=B(c,f,h,relative_epsilon_obj,k,Z_0)
B_fun=k*delta_L(c,f,h,relative_epsilon_obj)*sqrt(epsilon_e(c,f,h,relative_epsilon_obj))/Z_0;
end

%%------------------------------实际辐射单元长度L------------------------------
function L_fun=L(c,f,h,relative_epsilon_obj,relative_miu_obj)
L_fun=0.5*c/(f*sqrt(epsilon_e(c,f,h,relative_epsilon_obj)*relative_miu_obj))-2*delta_L(c,f,h,relative_epsilon_obj);
end
%%------------------------------delta_L，等效辐射缝隙长度------------------------------
function deltaL_fun=delta_L(c,f,h,relative_epsilon_obj)
deltaL_fun=0.412*h*(epsilon_e(c,f,h,relative_epsilon_obj)+0.3)*(w(c,f,relative_epsilon_obj)/h+0.264)/((epsilon_e(c,f,h,relative_epsilon_obj)-0.258)*(w(c,f,relative_epsilon_obj)/h+0.8));
end
%%------------------------------天线基板有效相对介电常数 (w/h>1) ------------------------------
function epsilonE_fun=epsilon_e(c,f,h,relative_epsilon_obj)
epsilonE_fun=0.5*(relative_epsilon_obj+1)+0.5*(relative_epsilon_obj-1)*(1+12*h/w(c,f,relative_epsilon_obj))^-0.5;
end
%%------------------------------实际辐射单元宽度w,半波长------------------------------
function w_fun=w(c,f,relative_epsilon_obj)
w_fun=(((relative_epsilon_obj+1)/2)^-0.5)*c/(2*f);
end
