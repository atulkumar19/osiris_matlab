clc;
clear all;
close all;

addpath('UtilFun');

ne=1.1e22; % cm^-3


qe=-4.8032e-10; %electron charge , statcoulomb
me=9.1094e-28; %electron mass , g
mi=25;
qbymi=25;
vel_c=2.9979e10; %velocity of light ,  cm/sec
vel_b=0.9;
omp_e=sqrt((4*pi*ne*qe^2)/me) 
omp_i=sqrt((4*pi*ne*qe^2)/mi)
x_nor=vel_c/omp_e
t_nor=1/omp_e
B_norm= me*omp_e*vel_c/qe
b3=0.1;
E_norm= me*omp_e*vel_c/qe
d_e=vel_c/omp_e
d_i=vel_c/omp_i
r_l=(vel_b)/(qbymi*b3)

d_i0=d_i/d_e
