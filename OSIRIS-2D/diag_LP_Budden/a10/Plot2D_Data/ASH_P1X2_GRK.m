clc;
clear all;
close all;

addpath('UtilFun');

ne=1.1e22; % cm^-3
frm=100:100;

pth0='/home/atul/OSIRIS/oblique_laser_LH';
pth1='/MS/RAW/Elec1/';
pth2='/MS/RAW/ion/';



% pth19='/MS/FLD/b3/';

str_h1='RAW-Elec1-';
str_h2='RAW-ion-';



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

for ii=1:5:length(frm)    
    
    frno=frm(ii);
    str_frno=num2str(frno);
    len_frno=length(str_frno);
    str_num((len_str_num-len_frno+1):end)=str_frno;
   fl_nm1=strcat(pth0,pth1,str_h1,str_num,str_ext);
    fl_nm2=strcat(pth0,pth2,str_h2,str_num,str_ext);
    
    
%     fl_nm19=strcat(pth0,pth19,str_h19,str_num,str_ext);
    
%     [xg_qe,yg_qe,dset_qe1,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm1);
%     [xg_qe,yg_qe,dset_qe2,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm2);
%     [xg_qe,yg_qe,dset_qe3,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm3);

    [ene_raw,x1_raw,p1_raw1,x2_raw,p2_raw1,tag_raw1,xlt_in,xlt_fi]=...
        AshReadRaw(fl_nm1);
    [ene_raw,x1_raw,p1_raw2,x2_raw,p2_raw2,tag_raw2,xlt_in,xlt_fi]=...
        AshReadRaw(fl_nm2);
   
    
%     [xg_b3,yg_b3,dset_b3,x1lt_b3,x2lt_b3,time_b3]=...
%         ReadB313Aug2014(fl_nm19);
    
%     p1_raw=p1_raw1'+p1_raw2';
       
%     p2_raw=p2_raw1+p2_raw2+p2_raw3+p2_raw4+p2_raw5+p2_raw6+p2_raw7+p2_raw8+p2_raw9+p2_raw10+p2_raw11+p2_raw12+p2_raw13+p2_raw14+p2_raw15+p2_raw16+p2_raw17+p2_raw18;
    
    p_raw_sq1=p1_raw1.*p1_raw1+p2_raw1.*p2_raw1; p_raw_sq2=p1_raw2.*p1_raw2+p2_raw2.*p2_raw2;
    
    p_raw_sq= [p_raw_sq1; p_raw_sq2];
    KE_raw=sqrt(p_raw_sq2+1)-1;
     [counts,centers] = hist(0.5e3*KE_raw,500000);
      plot(centers,counts,'*--');hold on;
%      f=(gamma(5/2))^(3/2)*(gamma(3/2))^(-5/2)*exp(-KE_raw*gamma(5/2)*(gamma(3/2))^(-1));
%     
%      plot( KE_raw(1:1:length(KE_raw)),'.');hold on;
% hist(KE_raw,500);
% [f,x]=hist(KE_raw,500); %use hist function and get unnormalized values
% figure; plot(x,f/trapz(x,f),'b-*');
% [bincounts,ind] = histc(KE_raw,500);

    set(gca,'FontSize',20,'FontWeight','Bold');
%      xlim([0 0.05])
   
    xlabel('x2(c/\omega_{pe})'); ylabel('p1(mc)');
%     title(strcat('Phase-Space Plot, t=',num2str(time_b3)));
    
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