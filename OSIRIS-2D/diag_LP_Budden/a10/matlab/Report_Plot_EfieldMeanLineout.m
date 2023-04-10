clc;
clear all;
close all;

addpath('UtilFun');

ne=1.1e22; % cm^-3
frm=1:500;

pth0='/Volumes/Macintosh_HD_4/infinite_counter_beam';

pth1='/MS/DENSITY/Elec1/charge/';
pth2='/MS/DENSITY/Elec2/charge/';

pth3='/MS/FLD/e1/';
pth4='/MS/FLD/e2/';
pth5='/MS/FLD/e3/';


str_h1='charge-Elec1-';
str_h2='charge-Elec2-';

str_h3='e1-';
str_h4='e2-';
str_h5='e3-';


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
 aviobj = avifile ( 'mean-E-field', 'fps',5); 
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
    
    
    
    [xg_qe,yg_qe,dset_qe1,x1lt_qe,x2lt_qe,time_qe]=...
        ReadECharge13Aug2014(fl_nm1);
    [xg_qe,yg_qe,dset_qe2,x1lt_qe,x2lt_qe,time_qe]=...
        ReadECharge13Aug2014(fl_nm2);
    
    [xg_e1,yg_e1,dset_e1,x1lt_e1,x2lt_e1,time_e1]=...
        ReadE113Aug2014(fl_nm3);
    [xg_e2,yg_e2,dset_e2,x1lt_e2,x2lt_e2,time_e2]=...
        ReadE213Aug2014(fl_nm4);
    [xg_e3,yg_e3,dset_e3,x1lt_e3,x2lt_e3,time_e3]=...
        ReadE213Aug2014(fl_nm5);
    
    dset_qe=dset_qe1+dset_qe2;
%     dset_J=dset_J1+dset_J2+dset_J3;
    mean1_qe=mean(dset_qe,1);
    mean1_e1=mean(dset_e1,1);
    mean1_e2=mean(dset_e2,1);
    mean1_e3=mean(dset_e3,1);
    
    subplot(2,2,1)
    plot(yg_e1,mean1_e1,'b','LineWidth',2);
    set(gca,'linewidth',2.0,'fontsize',20,'fontweight','b');
     title(strcat('mean-E_x-field, t= ',num2str(time_qe)));
      xlabel('x2(c/\omega_{pe})');
    
    ylabel('E_x(m\omega_{pe}/e)');
    subplot(2,2,2)
    plot(yg_e2,mean1_e2,'b','LineWidth',2);
     set(gca,'linewidth',2.0,'fontsize',20,'fontweight','b');
    title(strcat('mean-E_y-field, t= ',num2str(time_qe)));
      xlabel('x2(c/\omega_{pe})');
    
    ylabel('E_y(m\omega_{pe}/e)');
    subplot(2,2,3)
    plot(yg_e3,mean1_e3,'b','LineWidth',2);
     set(gca,'linewidth',2.0,'fontsize',20,'fontweight','b');
      title(strcat('mean-E_z-field, t= ',num2str(time_qe)));
      xlabel('x2(c/\omega_{pe})');
    
    ylabel('E_z(m\omega_{pe}/e)');
    subplot(2,2,4)
    plot(yg_qe,mean1_qe,'b','LineWidth',2);
     ylim([-1.2479 -0.8330]);
     set(gca,'linewidth',2.0,'fontsize',20,'fontweight','b');
      title(strcat('mean-e-density, t= ',num2str(time_qe)));
      xlabel('x2(c/\omega_{pe})');
    
    ylabel('e-den(n_0)');
    drawnow;
    frame = getframe ( gcf ); aviobj = addframe ( aviobj, frame );

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
 aviobj = close ( aviobj );
% close(writerObj);