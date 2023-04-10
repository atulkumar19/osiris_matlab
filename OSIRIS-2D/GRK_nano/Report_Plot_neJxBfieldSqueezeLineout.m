clc;
clear all;
close all;

addpath('UtilFun');

ne=1.1e22; % cm^-3
frm=1:500;

pth0='/Users/akash/Desktop/Users/atul/infinite_counter_beam';
pth1='/MS/DENSITY/Elec1/charge/';
pth2='/MS/DENSITY/Elec2/charge/';

pth3='/MS/DENSITY/Elec1/j1/';
pth4='/MS/DENSITY/Elec2/j1/';

pth5='/MS/FLD/b3/';

str_h1='charge-Elec1-';
str_h2='charge-Elec2-';
str_h3='j1-Elec1-';
str_h4='j1-Elec2-';

str_h5='b3-';

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
 aviobj = avifile ( 'squeeze-eden-b-cur', 'fps',5); 
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
   
    [xg_J1,yg_J1,dset_J1,x1lt_J1,x2lt_J1,time_J1]=...
        ReadECharge13Aug2014(fl_nm3);
    [xg_J2,yg_J2,dset_J2,x1lt_J2,x2lt_J2,time_J2]=...
        ReadECharge13Aug2014(fl_nm4);
    
    [xg_b3,yg_b3,dset_b3,x1lt_b3,x2lt_b3,time_b3]=...
        ReadB313Aug2014(fl_nm5);    
    dset_qe=dset_qe1+dset_qe2;
    dset_J=dset_J1+dset_J2;
    
    dx_b3=abs(xg_b3(2)-xg_b3(1));
    x1minb3=x1lt_b3(1); x1maxb3=x1lt_b3(2);
    x2minb3=x2lt_b3(1); x2maxb3=x2lt_b3(2);
    [Ngx_LOb3,Ngy_LOb3]=size(dset_b3);
    Ngyby2b3=round(Ngy_LOb3/2);
    LOb3_dat=dset_b3(:,Ngyby2b3);
   J=zeros(Ngx_LOb3,1);
   B3=zeros(Ngx_LOb3,1);
   qe=zeros(Ngx_LOb3,1);
   J(:)=squeeze(dset_J(150,:));
   B3(:)=squeeze(dset_b3(150,:));
   qe(:)=squeeze(dset_qe(150,:));
    
%     mean1_J1=mean(dset_J,2);
%     mean1_b3=mean(dset_b3,1);
%     mean1_qe=mean(dset_qe,1);
    
    subplot(2,2,1)
    plot(yg_qe,J(:),'b','LineWidth',2);
    set(gca,'linewidth',2.0,'fontsize',20,'fontweight','b');
    title(strcat('current, t= ',num2str(time_qe)));
    xlabel('x2(c/\omega_{pe})');
    
    ylabel('J_x(n_0ec)');
    subplot(2,2,2)
    plot(yg_b3,B3(:),'b','LineWidth',2);
     set(gca,'linewidth',2.0,'fontsize',20,'fontweight','b');
    title(strcat('B-field, t= ',num2str(time_qe)));
      xlabel('x2(c/\omega_{pe})');
    
    ylabel('B_z(mc\omega_{pe}/e)');
    subplot(2,2,[3 4])
    plot(yg_qe,qe(:),'b','LineWidth',2);
     set(gca,'linewidth',2.0,'fontsize',20,'fontweight','b');
     title(strcat('e-density, t= ',num2str(time_qe)));
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