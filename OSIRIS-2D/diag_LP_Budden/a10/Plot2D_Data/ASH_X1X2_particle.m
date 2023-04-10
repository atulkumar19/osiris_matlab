clc;
clear all;
close all;

addpath('UtilFun');

ne=1.1e22; % cm^-3
frm=350:350;

minch=-4;
maxch=0;

pth0='/home/testuser/OSIRIS/f12-0.9c';
pth1='/MS/DENSITY/Elec1/charge/';
pth2='/MS/DENSITY/Elec2/charge/';
pth3='/MS/DENSITY/Elec3/charge/';
pth4='/MS/RAW/Elec1/';
pth5='/MS/RAW/Elec2/';
pth6='/MS/RAW/Elec3/';

str_h1='charge-Elec1-';
str_h2='charge-Elec2-';
str_h3='charge-Elec3-';
str_h4='RAW-Elec1-';
str_h5='RAW-Elec2-';
str_h6='RAW-Elec3-';


str_ext='.h5';
str_num='000000';
len_str_num=length(str_num);


qe=-4.8032e-10; %electron charge , statcoulomb
me=9.1094e-28; %electron mass , g
vel_c=2.9979e10; %velocity of light ,  cm/sec
omp_e=sqrt((4*pi*ne*qe^2)/me); x_nor=vel_c/omp_e; t_nor=1/omp_e;

len_fr=length(frm);
scrsz = get(0,'ScreenSize');
hfig=figure('Position',[50 10 750 650]);
% Get the width and height of the figure
%============AVI-FILE-NAME==========================
writerObj = VideoWriter ('phase_space.avi');
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
    fl_nm3=strcat(pth0,pth3,str_h3,str_num,str_ext);

    fl_nm4=strcat(pth0,pth4,str_h4,str_num,str_ext); %raw data
    fl_nm5=strcat(pth0,pth5,str_h5,str_num,str_ext); %raw data
    fl_nm6=strcat(pth0,pth6,str_h6,str_num,str_ext); %raw data
    
    [xg_qe,yg_qe,dset_qe1,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm1);
    [xg_qe,yg_qe,dset_qe2,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm2);
    [xg_qe,yg_qe,dset_qe3,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm3);

    [ene_raw,x1_raw1,p1_raw1,x2_raw1,p2_raw1,tag_raw1,xlt_in,xlt_fi]=...
        AshReadRaw(fl_nm4);
     [ene_raw,x1_raw2,p1_raw2,x2_raw2,p2_raw2,tag_raw2,xlt_in,xlt_fi]=...
        AshReadRaw(fl_nm5);
     [ene_raw,x1_raw3,p1_raw3,x2_raw3,p2_raw3,tag_raw3,xlt_in,xlt_fi]=...
        AshReadRaw(fl_nm6);
    
%    hp= plot(x1_raw1(1:100:length(x1_raw1)),x2_raw1(1:100:length(x1_raw1)),...
%         x1_raw2(1:3000:length(x1_raw2)),x2_raw2(1:3000:length(x1_raw2)),...
%         x1_raw3(1:5000:length(x1_raw3)),x2_raw3(1:5000:length(x1_raw3)));
%      set(hp(1),'linestyle','none','markeredgecolor','none',...
%     'color',[1 0 0],'markersize',2,'markerfacecolor','b','marker','o');
%     set(gca,'fontsize',20,'fontweight','b');
%     xlabel('x1(c/\omega_{pe})'); ylabel('x2(c/\omega_{pe})');
% 
% set(hp(2),'linestyle','none','markeredgecolor','none',...
%      'color',[0 0 1],'markersize',2,'markerfacecolor','g','marker','o');
%  
% set(hp(3),'linestyle','none','markeredgecolor','none',...
%     'markersize',2,'markerfacecolor','g','marker','o');
 nBins_x = 1500;
 nBins_y = 1500;
  [counts, bin_centers] = hist3([x1_raw1(1:length(x1_raw1)) x2_raw1(1:length(x2_raw1))], [nBins_x nBins_y]);
  x_bin_centers = bin_centers{1};
  y_bin_centers = bin_centers{2};
  pcolor( x_bin_centers,y_bin_centers,  log(counts'));caxis([0 2.5]);colormap('jet(128)');shading('interp');colorbar;hold on;
%     
%     plot(x1_raw1(1:100:length(x1_raw1)),x2_raw1(1:100:length(x1_raw1)),'.b');hold on;
%     plot(x1_raw2(1:3000:length(x1_raw2)),x2_raw2(1:3000:length(x1_raw2)),'.g');hold on;
%     plot(x1_raw3(1:5000:length(x1_raw3)),x2_raw3(1:5000:length(x1_raw3)),'.g');hold off;
    set(gca,'FontSize',28,'FontWeight','normal','linewidth',2,'xtick',[],'ytick',[]);
  
   xlim([6 30]);
   ylim([5 25]);
   
%     xlabel('X(c/\omega_{pe})'); ylabel('Y(c/\omega_{pe})');
%     title(strcat(' t=',num2str(time_qe)));
    
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