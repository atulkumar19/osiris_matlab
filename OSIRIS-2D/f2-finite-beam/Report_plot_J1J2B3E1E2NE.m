clc;
clear all;
close all;

addpath('UtilFun');

ne=1.1e22; % cm^-3
frm=1:500;

pth0='/Volumes/Macintosh_HD_4/infinite_counter_beam';
pth1='/MS/DENSITY/Elec1/charge/';
pth2='/MS/DENSITY/Elec2/charge/';

pth3='/MS/DENSITY/Elec1/j1/';
pth4='/MS/DENSITY/Elec2/j1/';

pth5='/MS/DENSITY/Elec1/j2/';
pth6='/MS/DENSITY/Elec2/j2/';

pth7='/MS/FLD/e1/';
pth8='/MS/FLD/e2/';
pth9='/MS/FLD/e3/';

pth10='/MS/FLD/b3/';

str_h1='charge-Elec1-';
str_h2='charge-Elec2-';

str_h3='j1-Elec1-';
str_h4='j1-Elec2-';

str_h5='j2-Elec1-';
str_h6='j2-Elec2-';

str_h7='e1-';
str_h8='e2-';
str_h9='e3-';

str_h10='b3-';

str_ext='.h5';
str_num='000000';
len_str_num=length(str_num);

qe=-4.8032e-10; %electron charge , statcoulomb
me=9.1094e-28; %electron mass , g
vel_c=2.9979e10; %velocity of light ,  cm/sec
omp_e=sqrt((4*pi*ne*qe^2)/me); x_nor=vel_c/omp_e; t_nor=1/omp_e;

minch=-5;
maxch=0;
scrsz=get(0,'ScreenSize');
figure('Position',[50 10 1200 650]);
% Get the width and height of the figure
%============AVI-FILE-NAME==========================
 aviobj = avifile ( 'mean-J1J1B3E1E2-eden', 'fps',5); 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for ii=1:length(frm)   
    
    frno=frm(ii);
    
    str_frno=num2str(frno);
    len_frno=length(str_frno);
    str_num((len_str_num-len_frno+1):end)=str_frno;
    fl_nm1=strcat(pth0,pth1,str_h1,str_num,str_ext);
    fl_nm2=strcat(pth0,pth2,str_h2,str_num,str_ext);
    fl_nm3=strcat(pth0,pth3,str_h3,str_num,str_ext);
    fl_nm4=strcat(pth0,pth4,str_h4,str_num,str_ext);
    fl_nm5=strcat(pth0,pth5,str_h5,str_num,str_ext);
    fl_nm6=strcat(pth0,pth6,str_h6,str_num,str_ext);
    fl_nm7=strcat(pth0,pth7,str_h7,str_num,str_ext);
    fl_nm8=strcat(pth0,pth8,str_h8,str_num,str_ext);
    fl_nm9=strcat(pth0,pth9,str_h9,str_num,str_ext);
    fl_nm10=strcat(pth0,pth10,str_h10,str_num,str_ext);
   
    
    
    
    [xg_qe,yg_qe,dset_qe1,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm1);
    [xg_qe,yg_qe,dset_qe2,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm2);
    
    [xg_qe,yg_qe,dset_J1,x1lt_J1,x2lt_J1,time_J1]=...
        ReadECharge13Aug2014(fl_nm3);
    [xg_qe,yg_qe,dset_J2,x1lt_J2,x2lt_J2,time_J2]=...
        ReadECharge13Aug2014(fl_nm4);
    
    [xg_qe,yg_qe,dset_J4,x1lt_J4,x2lt_J4,time_J4]=...
        ReadECharge13Aug2014(fl_nm5);
    [xg_qe,yg_qe,dset_J5,x1lt_J5,x2lt_J5,time_J5]=...
        ReadECharge13Aug2014(fl_nm6);
    
    
    [xg_e1,yg_e1,dset_e1,x1lt_e1,x2lt_e1,time_e1]=...
        ReadE113Aug2014(fl_nm7);
    [xg_e2,yg_e2,dset_e2,x1lt_e2,x2lt_e2,time_e2]=...
        ReadE213Aug2014(fl_nm8);
    [xg_e3,yg_e3,dset_e3,x1lt_e3,x2lt_e3,time_e3]=...
        ReadE213Aug2014(fl_nm9);
    
    [xg_b3,yg_b3,dset_b3,x1lt_b3,x2lt_b3,time_b3]=...
        ReadB313Aug2014(fl_nm10); 
    
    dset_qe=dset_qe1+dset_qe2;
    dset_J1=dset_J1+dset_J2;
    dset_J2=dset_J4+dset_J5;
    
    mean1_qe=mean(dset_qe);
    mean1_J1=mean(dset_J1,1);
    mean1_J2=mean(dset_J2,1);
    mean1_E1=mean(dset_e1,1);
    mean1_E2=mean(dset_e2,1);
    mean1_b3=mean(dset_b3,1);
    
    
    subplot(2,3,1)
    plot(yg_qe,mean1_J1,'b','LineWidth',2);
    set(gca,'linewidth',2.0,'fontsize',20,'fontweight','b');
    xlim([0 25]);
    title(strcat(' t= ',num2str(time_qe)));
    xlabel('x2(c/\omega_{pe})');
    
    ylabel('J_x(n_0ec)');
    subplot(2,3,2)
    plot(yg_qe,mean1_J2,'b','LineWidth',2);
    xlim([0 25]);
     set(gca,'linewidth',2.0,'fontsize',20,'fontweight','b');
    title(strcat(' t= ',num2str(time_qe)));
      xlabel('x2(c/\omega_{pe})');
    
   ylabel('J_y(n_0ec)');
    subplot(2,3,3)
    plot(yg_b3,mean1_b3,'b','LineWidth',2);
    xlim([0 25]);
     set(gca,'linewidth',2.0,'fontsize',20,'fontweight','b');
     title(strcat(' t= ',num2str(time_qe)));
      xlabel('x2(c/\omega_{pe})');
    
   ylabel('B_z(mc\omega_{pe}/e)');
   
   subplot(2,3,4)
    plot(yg_qe,mean1_E1,'b','LineWidth',2);
     set(gca,'linewidth',2.0,'fontsize',20,'fontweight','b');
     xlim([0 25]);
      title(strcat(' t= ',num2str(time_qe)));
      xlabel('x2(c/\omega_{pe})');
    
    ylabel('E_x(m\omega_{pe}/e)');
    subplot(2,3,5)
    plot(yg_e2,mean1_E2,'b','LineWidth',2);
     set(gca,'linewidth',2.0,'fontsize',20,'fontweight','b');
     xlim([0 25]);
    title(strcat(' t= ',num2str(time_qe)));
      xlabel('x2(c/\omega_{pe})');
    
    ylabel('E_y(m\omega_{pe}/e)');
    subplot(2,3,6)
    plot(yg_e3,mean1_qe,'b','LineWidth',2);
     set(gca,'linewidth',2.0,'fontsize',20,'fontweight','b');
     xlim([0 25]);
      ylim([-1.2479 -0.8330]);
      title(strcat(' t= ',num2str(time_qe)));
      xlabel('x2(c/\omega_{pe})')
       ylabel('e-den(n_0)');
    drawnow;
    frame = getframe ( gcf ); aviobj = addframe ( aviobj, frame );

    
    
    
end
 aviobj = close ( aviobj );
% close(writerObj);