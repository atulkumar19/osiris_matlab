clc;
clear all;
close all;

addpath('UtilFun');

ne=4.0e21; % cm^-3
frm=0:200;

pth0='/Users/akash/Desktop/Users/Osiris_run/GRK_2s_5Aug15_0.98c_periodic';
pth1='/MS/DENSITY/Elec1/charge/';
pth2='/MS/DENSITY/Elec2/charge/';


str_h1='charge-Elec1-';
str_h2='charge-Elec2-';

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

for ii=1:length(frm)   
    
    frno=frm(ii);
    
    str_frno=num2str(frno);
    len_frno=length(str_frno);
    str_num((len_str_num-len_frno+1):end)=str_frno;
    fl_nm1=strcat(pth0,pth1,str_h1,str_num,str_ext);
    fl_nm2=strcat(pth0,pth2,str_h2,str_num,str_ext);
    
    [xg_qe,yg_qe,dset_qe1,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm1);
    [xg_qe,yg_qe,dset_qe2,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm2);
    
    dset_qe=dset_qe1+dset_qe2;
    
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

    imagesc(dset_qe','XData',[x1minqe x1maxqe],...
        'YData',[x2minqe x2maxqe]);colorbar;
    shading('interp');%colorbar;caxis(gca,[-5 0]);
%     colormap(cmap('blue',1024)); colorbar;
%     caxis(gca,[minch maxch]); 
%     colormap('jet'); shading('interp');
    xlabel('x1(c/\omega_{pe})'); ylabel('x2(c/\omega_{pe})');
    xlim([0 12]); ylim([4 22]);
    title(num2str(ii))
    drawnow;
    
    
%     saveas(gcf,'CurrentFig.png');
%     img4=imread('CurrentFig.png');
%     writeVideo(writerObj,img4);
 gif_add_frame(gcf,'kh.gif',2);
end

% close(writerObj);
