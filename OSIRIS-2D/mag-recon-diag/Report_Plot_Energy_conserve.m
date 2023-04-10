clc;
clear all;
close all;

addpath('UtilFun');

ne=1.1e22; % cm^-3
frm=1:500;

pth0='/Volumes/Macintosh_HD_3/osiris_run/counter_finite_tag';
pth1='/MS/DENSITY/Elec1/charge/';
pth2='/MS/DENSITY/Elec2/charge/';
pth3='/MS/DENSITY/Elec3/charge/';
pth4='/MS/FLD/b3/';
pth5='/MS/FLD/e1/';
pth6='/MS/FLD/e2/';
pth7='/MS/PHA/p1x1/Elec1/';
pth8='/MS/PHA/p1x1/Elec2/';
pth9='/MS/PHA/p1x1/Elec3/';
pth10='/MS/PHA/p2x2/Elec1/';
pth11='/MS/PHA/p2x2/Elec2/';
pth12='/MS/PHA/p2x2/Elec3/';
pth13='/MS/PHA/p3x1/Elec1/';
pth14='/MS/PHA/p3x1/Elec2/';
pth15='/MS/PHA/p3x1/Elec3/';





str_h1='charge-Elec1-';
str_h2='charge-Elec2-';
str_h3='charge-Elec3-';

str_h4='b3-';
str_h5='e1-';
str_h6='e2-';
str_h7='p1x1-Elec1-';
str_h8='p1x1-Elec2-';
str_h9='p1x1-Elec3-';
str_h10='p2x2-Elec1-';
str_h11='p2x2-Elec2-';
str_h12='p2x2-Elec3-';
str_h13='p3x1-Elec1-';
str_h14='p3x1-Elec2-';
str_h15='p3x1-Elec3-';
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
MagE=zeros(length(frm),1);
ElecE=zeros(length(frm),1);
KinE=zeros(length(frm),1);
timE=zeros(length(frm),1);
TotE=zeros(length(frm),1);
KinE0=zeros(length(frm),1);
KinE_per=zeros(length(frm),1);
cnt=1;
% %============AVI-FILE-NAME==========================
%  aviobj = avifile ( 'squeeze-eden-b-cur', 'fps',5); 
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
    fl_nm11=strcat(pth0,pth11,str_h11,str_num,str_ext);
    fl_nm12=strcat(pth0,pth12,str_h12,str_num,str_ext);
    fl_nm13=strcat(pth0,pth13,str_h13,str_num,str_ext);
    fl_nm14=strcat(pth0,pth14,str_h14,str_num,str_ext);
    fl_nm15=strcat(pth0,pth15,str_h15,str_num,str_ext);
    
    
    
    [xg_qe,yg_qe,dset_qe1,x1lt_qe,x2lt_qe,time_qe]=...
        ReadECharge13Aug2014(fl_nm1);
    [xg_qe,yg_qe,dset_qe2,x1lt_qe,x2lt_qe,time_qe]=...
        ReadECharge13Aug2014(fl_nm2);
    [xg_qe,yg_qe,dset_qe3,x1lt_qe,x2lt_qe,time_qe]=...
        ReadECharge13Aug2014(fl_nm3);
   
    [xg_b3,yg_b3,dset_b3,x1lt_b3,x2lt_b3,time_b3]=...
        ReadB313Aug2014(fl_nm4);    
    [xg_e1,yg_e1,dset_e1,x1lt_e1,x2lt_e1,time_e1]=...
        ReadE113Aug2014(fl_nm5);
     [xg_e2,yg_e2,dset_e2,x1lt_e2,x2lt_e2,time_e2]=...
        ReadE213Aug2014(fl_nm6);  
    
    
    
    [xg_P1X1,yg_P1X1,dset_P1X1,x1lt_P1X1,x2lt_P1X1,time_P1X1]=...
        ReadE213Aug2014(fl_nm7);
    [xg_P2X2,yg_P2X2,dset_P2X2,x1lt_P2X2,x2lt_P2X2,time_P2X2]=...
        ReadE213Aug2014(fl_nm8);
    [xg_P3X3,yg_P3X3,dset_P3X3,x1lt_P3X3,x2lt_P3X3,time_P3X3]=...
        ReadE213Aug2014(fl_nm9);
    [xg_P4X4,yg_P4X4,dset_P4X4,x1lt_P4X4,x2lt_P4X4,time_P4X4]=...
        ReadE213Aug2014(fl_nm10);
    [xg_P5X5,yg_P5X5,dset_P5X5,x1lt_P5X5,x2lt_P5X5,time_P5X5]=...
        ReadE213Aug2014(fl_nm11);
    [xg_P6X6,yg_P6X6,dset_P6X6,x1lt_P6X6,x2lt_P6X6,time_P6X6]=...
        ReadE213Aug2014(fl_nm12);
    [xg_P7X7,yg_P7X7,dset_P7X7,x1lt_P7X7,x2lt_P7X7,time_P7X7]=...
        ReadE213Aug2014(fl_nm13);
    [xg_P8X8,yg_P8X8,dset_P8X8,x1lt_P8X8,x2lt_P8X8,time_P8X8]=...
        ReadE213Aug2014(fl_nm14);
    [xg_P9X9,yg_P9X9,dset_P9X9,x1lt_P9X9,x2lt_P9X9,time_P9X9]=...
        ReadE213Aug2014(fl_nm15);
    
    dset_qe=dset_qe1+dset_qe2+dset_qe3;
   
    dset_p1=dset_P1X1+dset_P2X2+dset_P3X3;
    dset_p2=dset_P4X4+dset_P5X5+dset_P6X6;
    dset_p3=dset_P7X7+dset_P8X8+dset_P9X9;
    dset_psq=dset_p1.*conj(dset_p1)+dset_p2.*conj(dset_p2)+dset_p3.*conj(dset_p3);
    dset_Totsq=dset_psq+1;
    dset_KinE=sqrt(dset_Totsq)-1;
    
    KinE0(1)=sum(sum(dset_KinE(1)));   
    KinE(cnt)=sum(sum(dset_KinE));
    MagE(cnt)=sum(sum(dset_b3.*conj(dset_b3)));
    ElecE(cnt)=sum(sum(dset_e1.*conj(dset_e1)+dset_e2.*conj(dset_e2)));
    TotE(cnt)= KinE(cnt)+MagE(cnt)+ ElecE(cnt);
    timE(cnt)=time_b3;
   
    cnt=cnt+1;
   
   
%     drawnow;
%     frame = getframe ( gcf ); aviobj = addframe ( aviobj, frame );

    %{

%     
    x1minb3=x1lt_b3(1); x1maxb3=x1lt_b3(2);
    x2minb3=x2lt_b3(1); x2maxb3=x2lt_b3(2);
    
    dx_qe=abs(xg_qe(2)-xg_qe(1));
    x1minqe=x1lt_qe(1); x1maxqe=x1lt_qe(2);
    x2minqe=x2lt_qe(1); x2maxqe=x2lt_qe(2);
    [Ngx_LOqe,Ngy_LOqe]=size(dset_qe);
    Ngyby2qe=round(Ngy_LOqe/2);
    LOqe_dat=dset_qe(:,Ngyby2qe);
    
    %{
    fs=1/dx_qe; %---->sampling frequency or sampling rate
    nyq=fs/2; %---->Nyquist frequency
    fft_LOqe=fft(LOqe_dat);
    Nfft=length(fft_LOqe);
    df=fs/Nfft;
    fvec=(0:(Nfft-1))*df;
    fvec(fvec>nyq) = fvec(fvec>nyq)-fs;
    freqs=fvec(fvec>=0);
    ampSpec=abs(fft_LOqe(fvec>=0));
    ampSpec=ampSpec/Nfft;
    ampSpec(2:ceil(Nfft/2))=2*ampSpec(2:ceil(Nfft/2));
    %}
    ax1= subplot(121);
    imagesc(dset_qe','XData',[x1minqe x1maxqe],...
        'YData',[x2minqe x2maxqe]);
    set(gca,'YDir','normal');
     set(gca,'linewidth',2.0,'fontsize',20,'fontweight','b');
%     print('autoExample', '-depsc2', '-r300');
%     z = get(gca,'position');
%     set(gca,'position',[z(1)-0.08 z(2)-0.01 z(3)+0.1 z(4)+0.01],'Xtick',[],'Ytick',[]);
    axis square;
   title('t=1.2')
    colormap('bone');
    shading(ax1,'interp');%caxis(gca,[-5 0]);
%     colormap(bone(2048))
%     colormap(cmap('Copper',1024)); 
    xlabel('x1(c/\omega_{pe})'); ylabel('x2(c/\omega_{pe})');
    xlim([0 25]); ylim([0 25]);
   
    freezeColors;
    cbfreeze(colorbar('south'));
     % colorbar('south');
     title(strcat('e-Density:',num2str(ii)))
    
    
    ax2=subplot(122);
    imagesc(dset_J1','XData',[x1minb3 x1maxb3],...
        'YData',[x2minb3 x2maxb3]); 
    set(gca,'YDir','normal');
    set(gca,'linewidth',2.0,'fontsize',20,'fontweight','b');
%     print('autoExample', '-depsc2', '-r300');
%     z = get(gca,'position');
%   set(gca,'position',[z(1)-0.08 z(2)-0.01 z(3)+0.1 z(4)+0.01],'Xtick',[],'Ytick',[]); 
    axis square;
    colormap(ax2,'jet');
    shading('interp');
    xlabel('x1(c/\omega_{pe})'); %ylabel('x2(c/\omega_{pe})');
    xlim([0 25]); ylim([0 25]);
      
    freezeColors;
    cbfreeze(colorbar('south'));
%     colorbar('south');
    title(strcat('B-field:',num2str(ii)))
    
    drawnow
%     saveas(gcf,'CurrentFig.png');
%     img4=imread('CurrentFig.png');
%     writeVideo(writerObj,img4);
% frame = getframe ( gcf ); 
% aviobj = addframe ( aviobj, frame );
    
    %}
    
    
end
figure(1);


plot(timE,MagE,'b','LineWidth',2);
set(gca,'fontsize',20,'fontweight','b');

xlabel('time(\omega_{pe})');
ylabel('Mag. Energy');




figure(2)

plot(timE,ElecE,'r','LineWidth',2);
set(gca,'fontsize',20,'fontweight','b');


xlabel('time(\omega_{pe})');
ylabel('Elect. Energy');


figure(3)

plot(timE,KinE,'m','LineWidth',2);
set(gca,'fontsize',20,'fontweight','b');


xlabel('time(\omega_{pe})');
ylabel('Kin. Energy');

figure(4)

plot(timE,TotE,'g','LineWidth',2);
set(gca,'fontsize',20,'fontweight','b');


xlabel('time(\omega_{pe})');
ylabel('Tot. Energy');


       


       
       

   
%  aviobj = close ( aviobj );
% close(writerObj);