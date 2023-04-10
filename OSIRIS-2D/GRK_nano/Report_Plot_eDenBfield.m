clc;
clear all;
close all;

addpath('UtilFun');

ne=1.1e22; % cm^-3
frm=0:190;

pth0='/home/testuser/GRK_oblique';
pth1='/MS/DENSITY/Elec1/charge/';
pth2='/MS/DENSITY/Elec2/charge/';
pth3='/MS/DENSITY/Elec3/charge/';
pth4='/MS/DENSITY/Elec4/charge/';
pth5='/MS/DENSITY/Elec5/charge/';
pth6='/MS/DENSITY/Elec6/charge/';
pth7='/MS/DENSITY/Elec7/charge/';
pth8='/MS/DENSITY/Elec8/charge/';
pth9='/MS/DENSITY/Elec9/charge/';
pth10='/MS/DENSITY/Elec10/charge/';
pth11='/MS/DENSITY/Elec11/charge/';
pth12='/MS/DENSITY/Elec12/charge/';
pth13='/MS/DENSITY/Elec13/charge/';
pth14='/MS/DENSITY/Elec14/charge/';
pth15='/MS/DENSITY/Elec15/charge/';
pth16='/MS/DENSITY/Elec16/charge/';
pth17='/MS/DENSITY/Elec17/charge/';
pth18='/MS/DENSITY/Elec18/charge/';

pth19='/MS/FLD/b3/';

str_h1='charge-Elec1-';
str_h2='charge-Elec2-';
str_h3='charge-Elec3-';
str_h4='charge-Elec4-';
str_h5='charge-Elec5-';
str_h6='charge-Elec6-';
str_h7='charge-Elec7-';
str_h8='charge-Elec8-';
str_h9='charge-Elec9-';
str_h10='charge-Elec10-';
str_h11='charge-Elec11-';
str_h12='charge-Elec12-';
str_h13='charge-Elec13-';
str_h14='charge-Elec14-';
str_h15='charge-Elec15-';
str_h16='charge-Elec16-';
str_h17='charge-Elec17-';
str_h18='charge-Elec18-';

str_h19='b3-';
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
writerObj = VideoWriter ('eden-b-field_oblique_36nm.avi');
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
    fl_nm16=strcat(pth0,pth16,str_h16,str_num,str_ext);
    fl_nm17=strcat(pth0,pth17,str_h17,str_num,str_ext);
    fl_nm18=strcat(pth0,pth18,str_h18,str_num,str_ext);
    
    fl_nm19=strcat(pth0,pth19,str_h19,str_num,str_ext);
    
    [xg_qe,yg_qe,dset_qe1,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm1);
    [xg_qe,yg_qe,dset_qe2,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm2);
    [xg_qe,yg_qe,dset_qe3,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm3);
    [xg_qe,yg_qe,dset_qe4,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm4);
    [xg_qe,yg_qe,dset_qe5,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm5);
     [xg_qe,yg_qe,dset_qe6,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm6);
    [xg_qe,yg_qe,dset_qe7,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm7);
    [xg_qe,yg_qe,dset_qe8,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm8);
    [xg_qe,yg_qe,dset_qe9,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm9);
    [xg_qe,yg_qe,dset_qe10,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm10);
     [xg_qe,yg_qe,dset_qe11,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm11);
    [xg_qe,yg_qe,dset_qe12,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm12);
    [xg_qe,yg_qe,dset_qe13,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm13);
    [xg_qe,yg_qe,dset_qe14,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm14);
    [xg_qe,yg_qe,dset_qe15,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm15);
    [xg_qe,yg_qe,dset_qe16,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm16);
    [xg_qe,yg_qe,dset_qe17,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm17);
    [xg_qe,yg_qe,dset_qe18,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm18);
    
    [xg_b3,yg_b3,dset_b3,x1lt_b3,x2lt_b3,time_b3]=...
        ReadB313Aug2014(fl_nm19);
    
    
    dset_qe=dset_qe1+dset_qe2+dset_qe3+dset_qe4+dset_qe5+dset_qe6+dset_qe7+dset_qe8+dset_qe9+...
    dset_qe10+dset_qe11+dset_qe12+dset_qe13+dset_qe14+dset_qe15+dset_qe16+dset_qe17+dset_qe18;
%    dset_qe=dset_qe1;
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
    ax1=subplot(121);
    imagesc(dset_qe','XData',[x1minqe x1maxqe],...
        'YData',[x2minqe x2maxqe]);
    set(gca,'YDir','normal');

    axis square;
    set(gca,'FontSize',20,'FontWeight','Bold');
    shading('interp');
    colormap(ax1,'parula');
    xlabel('x1(c/\omega_{pe})'); ylabel('x2(c/\omega_{pe})');
   xlim([199 203]);
    
    cb = colorbar('south'); 
    set(cb,'position',[.13 .35 .335 .02]);
    title(strcat('e-Density:',num2str(time_qe)));
    
    
    ax2= subplot(122);
    imagesc(dset_b3','XData',[x1minb3 x1maxb3],...
        'YData',[x2minb3 x2maxb3]);
    set(gca,'YDir','normal');

    axis square;
    set(gca,'FontSize',20,'FontWeight','Bold');
    colormap(ax2,'jet'); shading('interp');
    xlabel('x1(c/\omega_{pe})'); %ylabel('x2(c/\omega_{pe})');
%     xlim([0 25])
    cb = colorbar('south'); 
    set(cb,'position',[.57 .35 .335 .02]);
    title(strcat('B-field:',num2str(time_qe)));
    
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