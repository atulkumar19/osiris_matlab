% clc;
% clear all;
% close all;

addpath('UtilFun');

ne=1.1e22; % cm^-3
frm=200:200;

pth0='/home/atul/OSIRIS/Asym_Ion_Recon_RR';
% pth0='/home/atul/OSIRIS/EM_PRE_POS-0.7-bxy';
pth1='/MS2/DENSITY/Elec1/charge/';
pth2='/MS2/DENSITY/Elec2/charge/';
% pth5='/MS2/DENSITY/Beam_ion/charge/';
pth3='/MS2/FLD/b1/';
pth4='/MS2/FLD/b2/';

str_h1='charge-Elec1-';
str_h2='charge-Elec2-';
% str_h5='charge-Beam_ion-';
str_h3='b1-';
str_h4='b2-';
% str_h4='charge-Elec4-';
% str_h5='charge-Elec5-';
% str_h4='b2-';
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
writerObj = VideoWriter ('Weinel_inhomo.avi');
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
     fl_nm4=strcat(pth0,pth4,str_h4,str_num,str_ext);
%       fl_nm5=strcat(pth0,pth5,str_h5,str_num,str_ext);
   
   
    
    [xg_qe,yg_qe,dset_qe1,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm1);
    [xg_qe,yg_qe,dset_qe2,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm2);
%     [xg_qe,yg_qe,dset_qe3,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm3);
%    
     [xg_b3,yg_b3,dset_b1,x1lt_b3,x2lt_b3,time_b3]=...
         ReadB313Aug2014(fl_nm3);
      [xg_b3,yg_b3,dset_b2,x1lt_b3,x2lt_b3,time_b3]=...
         ReadB313Aug2014(fl_nm4);
%      [xg_qe,yg_qe,dset_qe3,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm5);
    
    [Br,B_theta]=cart2pol(dset_b1,dset_b2 );
    dset_qe=dset_qe1+dset_qe2;
    
    
     x1minb3=x1lt_b3(1); x1maxb3=x1lt_b3(2);
     x2minb3=x2lt_b3(1); x2maxb3=x2lt_b3(2);
    
    dx_qe=abs(xg_qe(2)-xg_qe(1));
    x1minqe=x1lt_qe(1); x1maxqe=x1lt_qe(2);
    x2minqe=x2lt_qe(1); x2maxqe=x2lt_qe(2);
    [Ngx_LOqe,Ngy_LOqe]=size(dset_qe);
    Ngyby2qe=round(Ngy_LOqe/2);
    LOqe_dat=dset_qe(:,Ngyby2qe);
    dset_qe=dset_qe1+dset_qe2;
%     a1 = reshape(B_theta,409600,1); a1(409600:819200)=0;
%     dset_qe10 = reshape(a1,640,640,2);
%     outpute = fopen('fnm', 'wb'); 
%      fwrite(outpute, dset_qe10 ,'float'); 
%      fclose(outpute);
    
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
   ax1=subplot(121);
   imagesc(B_theta','XData',[x1minqe x1maxqe],...
        'YData',[x2minqe x2maxqe]);
%       (colormap(flipud(jet(4096))));hold on;
    set(gca,'YDir','normal');

axis square;
    set(gca,'FontSize',12,'FontWeight','Bold');
    shading('interp');
    colormap(ax1,'jet(4092)');
    xlabel('x1(c/\omega_{pe})'); ylabel('x2(c/\omega_{pe})');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [x,y]=meshgrid(x1minqe:dx_qe:x1maxqe,x2minqe:dx_qe:x2maxqe);
%    quiver(x,y,dset_b1',dset_b2',20);
% %    contour(x,y,B_theta',15,'k');

% xlim([30 50]);ylim([30 50]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cb = colorbar('south'); 
    set(cb,'position',[.13 .35 .335 .02]);
    title(strcat('B-field:',num2str(time_qe)));  hold off;
    
    
    ax2= subplot(122);
    imagesc(dset_qe2','XData',[x1minqe x1maxqe],...
        'YData',[x2minqe x2maxqe]);
    set(gca,'YDir','normal');
%    caxis([-1.2 0]);
    axis square;
    set(gca,'FontSize',12,'FontWeight','Bold');
    colormap(ax2,'hot(4092)'); shading('interp');
    xlabel('x1(c/\omega_{pe})'); %ylabel('x2(c/\omega_{pe})');

    cb = colorbar('south'); 
    set(cb,'position',[.57 .35 .335 .02]);
    title(strcat('e-current :',num2str(time_qe)));
    
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
