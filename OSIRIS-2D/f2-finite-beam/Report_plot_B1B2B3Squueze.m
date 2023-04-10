clc;
clear all;
close all;

addpath('UtilFun');

ne=1.1e22; % cm^-3
frm=1:500;

pth0='/Volumes/Macintosh_HD_4/counter_finite_tag';

pth1='/MS/FLD/e1/';
pth2='/MS/FLD/e2/';
pth3='/MS/FLD/e3/';


str_h1='e1-';
str_h2='e2-';
str_h3='e3-';


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
 aviobj = avifile ( 'Squeeze-E-field', 'fps',5); 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for ii=1:length(frm)   
    
    frno=frm(ii);
    
    str_frno=num2str(frno);
    len_frno=length(str_frno);
    str_num((len_str_num-len_frno+1):end)=str_frno;
    fl_nm1=strcat(pth0,pth1,str_h1,str_num,str_ext);
    fl_nm2=strcat(pth0,pth2,str_h2,str_num,str_ext);
    fl_nm3=strcat(pth0,pth3,str_h3,str_num,str_ext);
    
    
    [xg_e1,yg_e1,dset_e1,x1lt_e1,x2lt_e1,time_e1]=...
        ReadB313Aug2014(fl_nm1);    
    
    [xg_e2,yg_e2,dset_e2,x1lt_e2,x2lt_e2,time_e2]=...
        ReadB313Aug2014(fl_nm2);    
    
    [xg_e3,yg_e3,dset_e3,x1lt_e3,x2lt_e3,time_e3]=...
        ReadB313Aug2014(fl_nm3); 
    
    dx_e3=abs(xg_e3(2)-xg_e3(1));
    x1minb3=x1lt_e3(1); x1maxe3=x1lt_e3(2);
    x2minb3=x2lt_e3(1); x2maxe3=x2lt_e3(2);
    [Ngx_LOe3,Ngy_LOe3]=size(dset_e3);
    Ngyby2e3=round(Ngy_LOe3/2);
    LOe3_dat=dset_e3(:,Ngyby2e3);
   E1=zeros(Ngx_LOe3,1);
   E2=zeros(Ngx_LOe3,1);
   E3=zeros(Ngx_LOe3,1);
   E1(:)=squeeze(dset_e1(150,:));
   E2(:)=squeeze(dset_e2(150,:));
   E3(:)=squeeze(dset_e3(150,:));
    
%     mean1_e1=mean(dset_e1,1);
%     mean1_e2=mean(dset_e2,1);
%     mean1_e3=mean(dset_e3,1);
    
    
    subplot(2,2,1)
    plot(yg_e1,E1(:),'b','LineWidth',2);
    set(gca,'linewidth',2.0,'fontsize',20,'fontweight','b');
    title(strcat('E_x-field, t= ',num2str(time_e1)));
    xlabel('x2(c/\omega_{pe})');
    
    ylabel('E_x(m\omega_{pe}/e)');
    subplot(2,2,2)
    plot(yg_e2,E2(:),'b','LineWidth',2);
     set(gca,'linewidth',2.0,'fontsize',20,'fontweight','b');
    title(strcat('E_y-field, t= ',num2str(time_e2)));
      xlabel('x2(c/\omega_{pe})');
    
    ylabel('E_y(m\omega_{pe}/e)');
    subplot(2,2,[3 4])
    plot(yg_e3,E3(:),'b','LineWidth',2);
     set(gca,'linewidth',2.0,'fontsize',20,'fontweight','b');
     title(strcat('E_z-field, t= ',num2str(time_e3)));
      xlabel('x2(c/\omega_{pe})');
    
   ylabel('E_z(m\omega_{pe}/e)');
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