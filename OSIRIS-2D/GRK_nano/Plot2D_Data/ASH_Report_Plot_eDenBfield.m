clc;
clear all;
close all;

addpath('UtilFun');

ne=5.80000d+19; % cm^-3
frm=0:200;

pth0='/media/Elements/a_Jan1_2016/OsirisRao/bin/PICSim';
pth1='/MS/DENSITY/He_electrons/charge/';
pth2='/MS/FLD/e2/';

str_h1='charge-He_electrons-';
str_h2='e2-';

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
%  aviobj = avifile ( 'eden-b-field', 'fps',5); 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for ii=1:length(frm)   
    
    frno=frm(ii);
    
    str_frno=num2str(frno);
    len_frno=length(str_frno);
    str_num((len_str_num-len_frno+1):end)=str_frno;
    
    fl_nm1=strcat(pth0,pth1,str_h1,str_num,str_ext);
    fl_nm2=strcat(pth0,pth2,str_h2,str_num,str_ext);
    
    
    [xg_qe,yg_qe,dset_qe,x1lt_qe,x2lt_qe,time_qe,]=...
        AshReadHDF5DenDat(fl_nm1);
    [xg_e2,yg_e2,dset_e2,x1lt_e2,x2lt_e2,time_e2]=...
        AshReadHDF5DenDat(fl_nm2);
    

    x1mine2=x1lt_e2(1); x1maxe2=x1lt_e2(2);
    x2mine2=x2lt_e2(1); x2maxe2=x2lt_e2(2);
    
    dx_qe=abs(xg_qe(2)-xg_qe(1));
    x1minqe=x1lt_qe(1); x1maxqe=x1lt_qe(2);
    x2minqe=x2lt_qe(1); x2maxqe=x2lt_qe(2);
    [Ngx_LOqe,Ngy_LOqe]=size(dset_qe);
    Ngyby2qe=round(Ngy_LOqe/2);
    LOqe_dat=dset_qe(:,Ngyby2qe);
    
    LOe2_dat=dset_e2(:,Ngyby2qe);
    
    
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
    Himg=imagesc(dset_qe','XData',[x1minqe x1maxqe],...
        'YData',[x2minqe x2maxqe]); hold on;
    set(Himg,'AlphaData',0.2);
    set(gca,'YDir','normal');
     set(gca,'linewidth',2.0,'fontsize',20,'fontweight','b');
%     print('autoExample', '-depsc2', '-r300');
%     z = get(gca,'position');
%     set(gca,'position',[z(1)-0.08 z(2)-0.01 z(3)+0.1 z(4)+0.01],'Xtick',[],'Ytick',[]);
    axis square;
    colormap(cmap('blue',1024));
%     shading(ax1,'interp');%caxis(gca,[-5 0]);
%     colormap(bone(2048))
%     colormap(cmap('Copper',1024)); 
    xlabel('x1(c/\omega_{pe})'); ylabel('x2(c/\omega_{pe})');
%     xlim([0 25]); ylim([0 25]);
   
    freezeColors;
    cbfreeze(colorbar('east'));
%     colorbar('east');
    title(strcat('e-Density:','Frame=',num2str(frno),...
        ' t=',num2str(time_qe)))
    
    plot(xg_e2,LOe2_dat);
    hold off;
    
    
    ax2=subplot(122);
    imagesc(dset_e2','XData',[x1mine2 x1maxe2],...
        'YData',[x2mine2 x2maxe2]); 
    set(gca,'YDir','normal');
    set(gca,'linewidth',2.0,'fontsize',20,'fontweight','b');
%     print('autoExample', '-depsc2', '-r300');
%     z = get(gca,'position');
%   set(gca,'position',[z(1)-0.08 z(2)-0.01 z(3)+0.1 z(4)+0.01],'Xtick',[],'Ytick',[]); 
    axis square;
    colormap(ax2,'jet');
    shading('interp');
   xlabel('x1(c/\omega_{pe})');
%     xlim([0 25]); ylim([0 25]);
    
   title(strcat('B-field: t= ',num2str(time_e2)));
%     xlim([0 25]); ylim([0 25]);
      
    freezeColors;
    cbfreeze(colorbar('east'));
%     colorbar('south');
%     title(strcat('B-field:',num2str(ii)))

  
    drawnow

% frame = getframe ( gcf ); 
% aviobj = addframe ( aviobj, frame );
end
% aviobj = close ( aviobj );
