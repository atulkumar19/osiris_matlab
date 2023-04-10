clc;
clear all;
close all;

addpath('UtilFun');

ne=1.1e22; % cm^-3
frm=100:100;

minch=-4;
maxch=0;

pth0='/home/testuser/OSIRIS/Budden_plane_test';
pth1='/MS/DENSITY/Elec1/charge/';
pth2='/MS/DENSITY/ion/charge/';
% pth3='/MS/DENSITY/Elec3/charge/';
pth4='/MS/RAW/ion/';

str_h1='charge-Elec1-';
str_h2='charge-ion-';
% str_h3='charge-Elec3-';
str_h4='RAW-ion-';

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
writerObj = VideoWriter ('phase_space.avi');
writerObj.FrameRate=5;
open(writerObj);
% aviobj = avifile ( 'eden-b-field', 'fps',5); 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for ii=1:100:length(frm)    
    
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
    
%     ind_ylim1=round(ylim1/dy)+1; 
%     ind_ylim2=round(ylim2/dy)-1;
%     dset_qe_plt=dset_qe(:,ind_ylim1:ind_ylim2);
%     yg_cha_plt=yg(ind_ylim1:ind_ylim2);
%     dset_qe_plt=dset_qe;
%     yg_cha_plt=yg;
%     LT_a=(maxene-minene)/( ylim2-ylim1);
%     LT_b=minene-LT_a*ylim1;
%     yg_cha_plt=LT_a*yg_cha_plt+LT_b;
%     ylim1cha=LT_a*ylim1+LT_b;
%     ylim2cha=LT_a*ylim2+LT_b;
    
%     subtightplot(1,2,1)
%     [n,nout]=hist(ene_raw,128);
%     plot(nout,n,'linewidth', 2)
%     xlim([10 110])
%     xlabel('Ene(mc^2)'); ylabel('counts(arb units)');
%     title(strcat('frame=',num2str(frno),'  Hist Ene'));
%     grid on;
 plot(x1_raw(1:50:length(x1_raw)),ene_raw(1:50:length(ene_raw)),'.b');
    set(gca,'FontSize',12,'FontWeight','Bold');
  
    
%     ylim([-0.5 0.5])
   
    xlabel('x1(c/\omega_{pe})'); ylabel('p1(mc)');
    title(strcat('Phase-Space (Electrons), t=',num2str(time_qe)));
    
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