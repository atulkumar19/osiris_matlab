

addpath('UtilFun');

ne=1.1e2; % cm^-3
frm=1000:1000;

pth0='/home/testuser/OSIRIS/Weibel_assym';
pth1='/MS/DENSITY/Elec1/charge/';
pth2='/MS/DENSITY/Elec2/charge/';
% pth3='/MS/DENSITY/Elec3/charge/';
pth4='/MS/DENSITY/Elec1/ene/';
pth5='/MS/DENSITY/Elec2/ene/';
% pth6='/MS/DENSITY/Elec3/ene/';
pth7='/MS/FLD/ene_emf/';
pth8='/MS/FLD/ene_e/';
pth9='/MS/FLD/ene_b/';
pth10='/MS/FLD/b3/';

str_h1='charge-Elec1-';
str_h2='charge-Elec2-';
% str_h3='charge-Elec3-';

str_h4='ene-Elec1-';
str_h5='ene-Elec2-';
% str_h6='ene-Elec3-';

str_h7='ene_emf-';
str_h8='ene_e-';
str_h9='ene_b-';

str_h10='b3-';
str_ext='.h5';
str_num='000000';
len_str_num=length(str_num);



qe=-4.8032e-10; %electron charge , statcoulomb
me=9.1094e-28; %electron mass , g
vel_c=2.9979e10; %velocity of light ,  cm/sec
omp_e=sqrt((4*pi*ne*qe^2)/me); x_nor=vel_c/omp_e; t_nor=1/omp_e;

len_fr=length(frm);
% scrsz=get(0,'ScreenSize');
% figure('Position',[50 10 1200 650]);
 
% Get the width and height of the figure
% lbwh=get(1, 'position');
% figw=lbwh(3);
% figh=lbwh(4);
% Number of rows and columns of axes
ncols=len_fr;
nrows=3;
% w and h of each axis in normalized units
axisw=(1/ncols)*0.95;
axish=(1/nrows)*0.95;
MagE=zeros(length(frm),1);
ElecE=zeros(length(frm),1);
EME=zeros(length(frm),1);
KinE=zeros(length(frm),1);
timE=zeros(length(frm),1);
TotE=zeros(length(frm),1);
KinE0=zeros(length(frm),1);
KinE_per=zeros(length(frm),1);
cnt=1;

% Get the width and height of the figure
%============AVI-FILE-NAME==========================
writerObj = VideoWriter ('eden-b-field.avi');
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
    fl_nm4=strcat(pth0,pth4,str_h4,str_num,str_ext);
    fl_nm5=strcat(pth0,pth5,str_h5,str_num,str_ext);
%     fl_nm6=strcat(pth0,pth6,str_h6,str_num,str_ext);
    fl_nm7=strcat(pth0,pth7,str_h7,str_num,str_ext);
    fl_nm8=strcat(pth0,pth8,str_h8,str_num,str_ext);
    fl_nm9=strcat(pth0,pth9,str_h9,str_num,str_ext);
    fl_nm10=strcat(pth0,pth10,str_h10,str_num,str_ext);
    
    
    
    [xg_qe,yg_qe,dset_qe1,x1lt_qe,x2lt_qe,time_qe]=...
        ReadECharge13Aug2014(fl_nm1);
    [xg_qe,yg_qe,dset_qe2,x1lt_qe,x2lt_qe,time_qe]=...
        ReadECharge13Aug2014(fl_nm2);
%     [xg_qe,yg_qe,dset_qe3,x1lt_qe,x2lt_qe,time_qe]=...
%         ReadECharge13Aug2014(fl_nm3);
%    
    [xg_ene1,yg_ene1,dset_ene1,x1lt_ene1,x2lt_ene1,time_ene1]=...
        ReadECharge13Aug2014(fl_nm4);
    [xg_ene2,yg_ene2,dset_ene2,x1lt_ene2,x2lt_ene2,time_ene2]=...
        ReadECharge13Aug2014(fl_nm5);
%     [xg_ene3,yg_ene3,dset_ene3,x1lt_ene3,x2lt_ene3,time_ene3]=...
%         ReadECharge13Aug2014(fl_nm6);
%     
    [xg_ene_EMF,yg_ene_EMF,dset_ene_EMF,x1lt_ene_EMF,x2lt_ene_EMF,time_ene_EMF]=...
        ReadB313Aug2014(fl_nm7);
    
    [xg_ene_e,yg_ene_e,dset_ene_e,x1lt_ene_e,x2lt_ene_e,time_ene_e]=...
        ReadE113Aug2014(fl_nm8);
    
    [xg_ene_b,yg_ene_b,dset_ene_b,x1lt_ene_b,x2lt_ene_b,time_ene_b]=...
        ReadB313Aug2014(fl_nm9);
    
    [xg_b3,yg_b3,dset_b3,x1lt_b3,x2lt_b3,time_b3]=...
        ReadB313Aug2014(fl_nm10);
    
    dx_b3=abs(xg_b3(2)-xg_b3(1));
    x1minb3=x1lt_b3(1); x1maxb3=x1lt_b3(2);
    x2minb3=x2lt_b3(1); x2maxb3=x2lt_b3(2);
    [Ngx_LOb3,Ngy_LOb3]=size(dset_b3);
    Ngyby2b3=round(Ngy_LOb3/2);
    LOb3_dat=dset_b3(:,Ngyby2b3);
    
    dset_qe=dset_qe1+dset_qe2;
    dset_ene=dset_ene1+dset_ene2;
    dset_tot=dset_ene+dset_ene_EMF;
   
   
    
    
%     KinE(cnt)=sum(sum(dset_ene));
%     MagE(cnt)=sum(sum(dset_ene_b)).*(dx_b3)^2;
%     ElecE(cnt)=sum(sum(dset_ene_e));
%     
%     TotE(cnt)=sum(sum(dset_ene_EMF+dset_ene));
%     timE(cnt)=time_ene1;
%    
%     cnt=cnt+1;
    
    %-------fft----------------%
    nfft=1250;
   
%     subplot(121)
    E=zeros(Ngx_LOb3,1);
     E1=zeros(Ngx_LOb3,1);
   E(:)=squeeze(dset_b3(45,:));
   E1(:)=squeeze(dset_ene(45,:));
  F3=fftshift(fft( E ,nfft));
  F31=fftshift(fft( E1 ,nfft));
%    kx = (-nfft/2:nfft/2-1)/(dx_b3*nfft);%wave Vector
  Pows=F3.*conj(F3)/(nfft*nfft);
  Pows1=F31.*conj(F31)/(nfft*nfft);
  %computing power with proper scaling
  %gama1=smoothn(abs(F3)/(Ngx_LOb3*dx_b3),'0.1');

%      E(:,1)=mean(gama,2);%mean along y axix
%      plot((kx),abs(F3)/(Ngx_LOb3*dx_b3),'k','LineWidth',2);hold on;
   %----1sided------%
    kx=(0:nfft/2-1)/(dx_b3*nfft);	 	 
 loglog((kx),smooth(abs( F3(nfft/2+1:nfft))/(Ngx_LOb3*dx_b3)),'b','LineWidth',1);hold on;
%    loglog((kx),abs( F31(nfft/2+1:nfft))/(Ngx_LOb3*dx_b3),'b','LineWidth',2);hold on;
%  print('autoExample', '-depsc2', '-r300');
  set(gca,'linewidth',2.0,'fontsize',20,'fontweight','b');
 title(strcat('FFT Spectra (B3), t= ',num2str(time_b3)));
 xlim([0.01 25])
%   fs=1/dx_b3; %---->sampling frequency or sampling rate
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
% subplot(122)
    
%     imagesc(dset_b3','XData',[x1minb3 x1maxb3],...
%         'YData',[x2minb3 x2maxb3]);
%     set(gca,'YDir','normal');
%     colormap('jet'); colorbar('south')
%     shading('interp');
%     set(gca,'linewidth',2.0,'fontsize',20,'fontweight','b');
%     print('autoExample', '-depsc2', '-r300');
%     xlabel('x1(c/\omega_{pe})');
%     
%     ylabel('x2(c/\omega_{pe})');
%    
%     xlim([0 50]); ylim([0 50]);
%     
% %     title(num2str(ii))
     
    drawnow;
    
%     saveas(gcf,'CurrentFig.png');
%     img4=imread('CurrentFig.png');
%     writeVideo(writerObj,img4);
% gif_add_frame(gcf,'finite_beam_spectra.gif',2);

end

% close(writerObj);
