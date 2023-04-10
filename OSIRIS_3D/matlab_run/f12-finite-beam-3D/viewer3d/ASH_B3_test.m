
clc;
clear all;
close all;

ne=1.1e+22; % cm^-3 58.0e18, 6.12e+18

minch=-4;
maxch=0;

% mycmap=colormap(cmap('Jet',1024));

qe=-4.8032e-10; %electron charge , statcoulomb
me=9.1094e-28; %electron mass , g
vel_c=2.9979e10; %velocity of light ,  cm/sec
omp_e=sqrt((4*pi*ne*qe^2)/me); x_nor=vel_c/omp_e; t_nor=1/omp_e;

% len_fr=length(frm);
% scrsz = get(0,'ScreenSize');
% hfig=figure('Position',[50 10 1200 650]);
% set(gcf,'Renderer','OpenGL');

fl_nm1='b3-000250.h5'; %charge density
[xg,yg,zg,dset_b3,x1lt,x2lt,x3lt,time]=AshReadHDF5DenDat3D(fl_nm1);
% 
% pcolor3(dset_b3(1:20,1:20,1:20));

viewer3d(dset_b3);