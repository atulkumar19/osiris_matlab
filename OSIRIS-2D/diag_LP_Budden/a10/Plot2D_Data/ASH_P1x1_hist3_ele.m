clc;
clear all;
close all;

addpath('UtilFun');

ne=1.1e22; % cm^-3
frm=100:100;

minch=-4;
maxch=0;

 %pth0='/home/testuser/OSIRIS/oblique_laser_nob3';
% pth0='/home/testuser/OSIRIS/Budden_m1600';
% pth0='/home/testuser/OSIRIS/LaserPlasma';
% pth0='/home/testuser/OSIRIS/Budden_m25_3pulse';
 pth0='/home/atul/OSIRIS/Budden_plane_test_nob3';
%pth0='/home/testuser/OSIRIS/oblique_laser_LH';
%  pth0='/home/testuser/OSIRIS/Weibel_symm';
pth1='/MS/DENSITY/Elec1/charge/';
pth2='/MS/DENSITY/ion/charge/';
% pth3='/MS/DENSITY/Elec3/charge/';
pth4='/MS/RAW/Elec1/';

str_h1='charge-Elec1-';
str_h2='charge-ion-';
% str_h3='charge-Elec3-';
str_h4='RAW-Elec1-';

str_ext='.h5';
str_num='000000';
len_str_num=length(str_num);


qe=-4.8032e-10; %electron charge , statcoulomb
me=9.1094e-28; %electron mass , g
vel_c=2.9979e10; %velocity of light ,  cm/sec
omp_e=sqrt((4*pi*ne*qe^2)/me); x_nor=vel_c/omp_e; t_nor=1/omp_e;

len_fr=length(frm);
scrsz = get(0,'ScreenSize');
hfig=figure('Position',[80 50 800 625]);
% Get the width and height of the figure
%============AVI-FILE-NAME==========================
writerObj = VideoWriter ('p1p2_ele.avi');
writerObj.FrameRate=5;
open(writerObj);
% aviobj = avifile ( 'eden-b-field', 'fps',5); 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for ii=1:length(frm)    
    
    frno=frm(ii);
    str_frno=num2str(frno);
    len_frno=length(str_frno);
    str_num((len_str_num-len_frno+1):end)=str_frno;
    fl_nm1=strcat(pth0,pth1,str_h1,str_num,str_ext);
    fl_nm2=strcat(pth0,pth2,str_h2,str_num,str_ext);
%     fl_nm3=strcat(pth0,pth3,str_h3,str_num,str_ext);

    fl_nm4=strcat(pth0,pth4,str_h4,str_num,str_ext); %raw data
    
    [xg_qe,yg_qe,dset_qe1,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm1);
%     [xg_qe,yg_qe,dset_qe2,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm2);
%     [xg_qe,yg_qe,dset_qe3,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm3);

    [ene_raw,x1_raw,p1_raw,x2_raw,p2_raw,tag_raw,xlt_in,xlt_fi]=...
        AshReadRaw(fl_nm4);
%      nBins_x = 10500;
%  nBins_y = 5000;
%   [counts, bin_centers] = hist3([p1_raw(1:20000:length(p1_raw)) p2_raw(1:20000:length(p2_raw))], [nBins_x nBins_y]);
%   x_bin_centers = bin_centers{1};
%   y_bin_centers = bin_centers{2};
%   pcolor( x_bin_centers,y_bin_centers,  log(counts'));caxis([0 7]);colormap('jet(128)');shading('interp');colorbar;hold on;
   plot(x1_raw(1:50:length(x1_raw)),ene_raw(1:50:length(ene_raw)),'.b');
   

    set(gca,'FontSize',24,'LineWidth',2);
%      
   
    xlim([600.0 1200.0 ]);
%      ylim([-0.0 0.01]);
    xlabel('p1(mc)'); ylabel('p2(mc)');
%     axis('square');
    title(strcat(' t=',num2str(time_qe)));
    
  drawnow
    saveas(gcf,'CurrentFig.png');
    img4=imread('CurrentFig.png');
    writeVideo(writerObj,img4);
     L = graphicsversion('handlegraphics'); 
     currFrame = getframe;
%frame = getframe ( gcf ); % aviobj = addframe ( aviobj, frame );
%      Frame = getframe;
%      writeVideo(writerObj,Frame);

 end
% aviobj = close ( aviobj );
  close(writerObj);
