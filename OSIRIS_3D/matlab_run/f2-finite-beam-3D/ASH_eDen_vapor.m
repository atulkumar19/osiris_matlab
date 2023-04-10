clc;
clear all;
close all;

addpath('UtilFun');

ne=1.1e+22; % cm^-3 58.0e18, 6.12e+18
frm=250;

minch=-4;
maxch=0;

pth0='/Users/akash/Desktop/OSIRIS_3D/atul_osiris';
pth1='/MS/DENSITY/Elec1/charge/';
pth2='/MS/DENSITY/Elec2/charge/';
pth3='/MS/DENSITY/Elec3/charge/';

str_h1='charge-Elec1-';
str_h2='charge-Elec2-';
str_h3='charge-Elec3-';


str_ext='.h5';
str_num='000000';
len_str_num=length(str_num);


qe=-4.8032e-10; %electron charge , statcoulomb
me=9.1094e-28; %electron mass , g
vel_c=2.9979e10; %velocity of light ,  cm/sec
omp_e=sqrt((4*pi*ne*qe^2)/me); x_nor=vel_c/omp_e; t_nor=1/omp_e;

% mycmap=colormap(cmap('blue',1024));
% mycmap=colormap(bone(2048));

len_fr=length(frm);
% scrsz = get(0,'ScreenSize');
% hfig=figure('Position',[50 10 1200 650]);
% set(gcf,'Renderer','OpenGL');
for ii=1:length(frm)    
    
    frno=frm(ii);
    str_frno=num2str(frno);
    len_frno=length(str_frno);
    str_num((len_str_num-len_frno+1):end)=str_frno;

    fl_nm1=strcat(pth0,pth1,str_h1,str_num,str_ext); %charge density
    fl_nm2=strcat(pth0,pth2,str_h2,str_num,str_ext);
    fl_nm3=strcat(pth0,pth3,str_h3,str_num,str_ext);
    
    [xg,yg,zg,dset_qe1,x1lt,x2lt,x3lt,time]=AshReadHDF5DenDat3D(fl_nm1);
    [~,~,~,dset_qe2,~,~,~,~]=AshReadHDF5DenDat3D(fl_nm2);
    [~,~,~,dset_qe3,~,~,~,~]=AshReadHDF5DenDat3D(fl_nm3);
    
    
    [xmat,ymat,zmat]=meshgrid(xg,yg,zg);
    
    dset_qe=dset_qe1+dset_qe2+dset_qe3;
    
    outpute = fopen('fnm_eden_250', 'wb'); 
     fwrite(outpute, dset_qe,'float'); 
     fclose(outpute); 


    
end