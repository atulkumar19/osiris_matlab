clc;
clear all;
close all;

addpath('UtilFun');

ne=1.1e2; % cm^-3
frm=1:500;

pth0='/Users/akash/Desktop/Users/atul/counter_finite_tag';
pth1='/MS/FLD/b3/';

str_h1='b3-';
str_ext='.h5';
str_num='000000';
len_str_num=length(str_num);



qe=-4.8032e-10; %electron charge , statcoulomb
me=9.1094e-28; %electron mass , g
vel_c=2.9979e10; %velocity of light ,  cm/sec
omp_e=sqrt((4*pi*ne*qe^2)/me); x_nor=vel_c/omp_e; t_nor=1/omp_e;

len_fr=length(frm);
scrsz=get(0,'ScreenSize');
figure('Position',[50 10 1200 650]);
% Get the width and height of the figure
lbwh=get(1, 'position');
figw=lbwh(3);
figh=lbwh(4);
% Number of rows and columns of axes
ncols=len_fr;
nrows=3;
% w and h of each axis in normalized units
axisw=(1/ncols)*0.95;
axish=(1/nrows)*0.95;
% figure('Position',[50 10 1200 650]);
% Get the width and height of the figure
%============AVI-FILE-NAME==========================
 aviobj = avifile ( 'fft.avi', 'fps',5); 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for ii=1:length(frm) 
    frno=frm(ii);
    
    str_frno=num2str(frno);
    len_frno=length(str_frno);
    str_num((len_str_num-len_frno+1):end)=str_frno;
    fl_nm=strcat(pth0,pth1,str_h1,str_num,str_ext);

    [xg_b3,yg_b3,dset_b3,x1lt_b3,x2lt_b3,time_b3]=...
        ReadB313Aug2014(fl_nm);
    dx_b3=abs(xg_b3(2)-xg_b3(1));
    x1minb3=x1lt_b3(1); x1maxb3=x1lt_b3(2);
    x2minb3=x2lt_b3(1); x2maxb3=x2lt_b3(2);
    [Ngx_LOb3,Ngy_LOb3]=size(dset_b3);
    Ngyby2b3=round(Ngy_LOb3/2);
    LOb3_dat=dset_b3(:,Ngyby2b3);
    
    %-------fft----------------%
    nfft=1024;
    
    E=zeros(Ngx_LOb3,1);
   E(:)=squeeze(dset_b3(150,:));
  F3=fftshift(fft( E ,nfft));
%    kx = (-nfft/2:nfft/2-1)/(dx_b3*nfft);%wave Vector
  Pows=F3.*conj(F3)/(nfft*nfft); %computing power with proper scaling
  %gama1=smoothn(abs(F3)/(Ngx_LOb3*dx_b3),'0.1');

%      E(:,1)=mean(gama,2);%mean along y axix
%      plot((kx),abs(F3)/(Ngx_LOb3*dx_b3),'k','LineWidth',2);hold on;
   %----1sided------%
  kx=(0:nfft/2-1)/(dx_b3*nfft);	
 
%    fs=1/dx_b3; %---->sampling frequency or sampling rate
%     nyq=fs/2; %---->Nyquist frequency
%     fft_LOb3=fft(LOb3_dat);
%     Nfft=length(fft_LOb3);
%     df=fs/Nfft;
%     fvec=(0:(Nfft-1))*df;
%     fvec(fvec>nyq) = fvec(fvec>nyq)-fs;
%     freqs=fvec(fvec>=0);
%     ampSpec=abs(fft_LOb3(fvec>=0));
%     ampSpec=ampSpec/Nfft;
%     ampSpec(2:ceil(Nfft/2))=2*ampSpec(2:ceil(Nfft/2));
    %-----------end fft------------------%
    
        subplot(121)
 loglog((kx),abs( F3(nfft/2+1:nfft))/(Ngx_LOb3*dx_b3),'k','LineWidth',2);
 
 
 
%  print('autoExample', '-depsc2', '-r300');
 set(gca,'YDir','normal');
  set(gca,'linewidth',2.0,'fontsize',20,'fontweight','b');
   ylabel('|B_k|^2');
   xlabel('x2(c/\omega_{pe})');
  title(strcat('FFT Spectra, t= ',num2str(time_b3)));
     
    
subplot(122)
    
    imagesc(dset_b3','XData',[x1minb3 x1maxb3],...
        'YData',[x2minb3 x2maxb3]);
    colormap('jet'); colorbar('south')
    shading('interp');
     set(gca,'YDir','normal');
    set(gca,'linewidth',2.0,'fontsize',20,'fontweight','b');
%  print('autoExample', '-depsc2', '-r300');
    xlabel('x1(c/\omega_{pe})');
    
    ylabel('x2(c/\omega_{pe})');
   
    xlim([0 25]); ylim([0 25]);
    
   title(strcat('B-field: t= ',num2str(time_b3)));
     
    drawnow;
    
    saveas(gcf,'CurrentFig.png');
     img4=imread('CurrentFig.png');
  frame = getframe ( gcf ); aviobj = addframe ( aviobj, frame );

% gif_add_frame(gcf,'finite_beam_spectra.gif',2);

end
aviobj = close ( aviobj );
