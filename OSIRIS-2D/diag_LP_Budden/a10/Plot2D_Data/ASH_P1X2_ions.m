clc;
clear all;
close all;

addpath('UtilFun');

ne=1.1e22; % cm^-3
frm=100:100;

minch=-4;
maxch=0;

% pth0='/home/testuser/OSIRIS/Budden_plane_test';
% pth0='/home/testuser/OSIRIS/Budden_m1600';
% pth0='/home/testuser/OSIRIS/LaserPlasma';
% pth0='/home/testuser/OSIRIS/Budden_m25_3pulse';
% pth0='/home/testuser/OSIRIS/laser_test_tilt';
% pth0='/home/testuser/OSIRIS/Budden_m25_4l-density_ramp';
pth0='/home/testuser/OSIRIS/Weibel_assym';
pth1='/MS/DENSITY/Elec1/charge/';
pth2='/MS/DENSITY/Elec2/charge/';
% pth3='/MS/DENSITY/Elec3/charge/';
pth4='/MS/RAW/Elec1/';

str_h1='charge-Elec1-';
str_h2='charge-Elec2-';
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
hfig=figure('Position',[50 10 1200 650]);
% Get the width and height of the figure
%============AVI-FILE-NAME==========================
writerObj = VideoWriter ('phase_space_ions_m25_zoom.avi');
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
    
%     subplot(1,2,1)
    plot(x1_raw(1:400:length(x1_raw)),ene_raw(1:400:length(ene_raw)),'.b');
    set(gca,'FontSize',12,'FontWeight','Bold');
%     xlim([150 400]);
    
%     ylim([-0.1 0.1])
   
    xlabel('x1(c/\omega_{pe})'); ylabel('p1(mc)');
    title(strcat('Phase-Space (Ions), t=',num2str(time_qe)));
    
%     subplot(1,2,2)
%      plot(x1_raw(1:50:length(x1_raw)),ene_raw(1:50:length(ene_raw)),'.b');
%     set(gca,'FontSize',12,'FontWeight','Bold');
% %    xlim([150 400]);
% %     
%     ylim([0  0.01])
%    
%     xlabel('x1(c/\omega_{pe})'); ylabel('Energy(mc^2)');
%     title(strcat('Energy (Ions), t=',num2str(time_qe)));
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