clc;
clear all;
close all;

addpath('UtilFun');

ne=1.0e22; % cm^-3
frm=126:126;

pth0='/home/atul/OSIRIS/Sym_Ion_Recon';
pth1='/MS/DENSITY/Elec1/charge/';
pth2='/MS/DENSITY/Elec2/charge/';
pth4='/MS/FLD/b3/';

str_h1='charge-Elec1-';
str_h2='charge-Elec2-';


str_h4='b3-';
str_ext='.h5';
str_num='000000';
len_str_num=length(str_num);

qe=-4.8032e-10; %electron charge , statcoulomb
me=9.1094e-28; %electron mass , g
vel_c=2.9979e10; %velocity of light ,  cm/sec
omp_e=sqrt((4*pi*ne*qe^2)/me); x_nor=vel_c/omp_e; t_nor=1/omp_e;

b_ext=(me*vel_c*omp_e)./(qe);

minch=-5;
maxch=0;
scrsz=get(0,'ScreenSize');
figure('Position',[50 10 1200 650]);
% Get the width and height of the figure
% %============AVI-FILE-NAME==========================
writerObj = VideoWriter ('eden-b-field_symm.avi');
writerObj.FrameRate=5;
open(writerObj);
% aviobj = avifile ( 'eden-b-field', 'fps',5); 
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for ii=1:length(frm)   
    
    frno=frm(ii);
    
    str_frno=num2str(frno);
    len_frno=length(str_frno);
    str_num((len_str_num-len_frno+1):end)=str_frno;
    fl_nm1=strcat(pth0,pth1,str_h1,str_num,str_ext);
    fl_nm2=strcat(pth0,pth2,str_h2,str_num,str_ext);
    
    fl_nm4=strcat(pth0,pth4,str_h4,str_num,str_ext);
   
    
    [xg_qe,yg_qe,dset_qe1,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm1);
    [xg_qe,yg_qe,dset_qe2,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm2);
    
   
    [xg_b3,yg_b3,dset_b3,x1lt_b3,x2lt_b3,time_b3]=...
        ReadB313Aug2014(fl_nm4);
    
    
    dset_qe=dset_qe1+dset_qe2;
     %%%%%add in matlab%%%%%%%%%%%%%%%%%
%      a1 = reshape(dset_qe,1562500,1); a1(1562500:7812500)=0; 
%      dset_qe10 = reshape(a1,1250,1250,5); 
     outpute = fopen('fnm_ab800', 'w'); 
     fwrite(outpute, dset_b3,'float'); 
     fclose(outpute); 
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    x1minb3=x1lt_b3(1); x1maxb3=x1lt_b3(2);
    x2minb3=x2lt_b3(1); x2maxb3=x2lt_b3(2);
    
    dx_qe=abs(xg_qe(2)-xg_qe(1));
    dy_qe=abs(yg_qe(2)-yg_qe(1));
    x1minqe=x1lt_qe(1); x1maxqe=x1lt_qe(2);
    x2minqe=x2lt_qe(1); x2maxqe=x2lt_qe(2);
    [Ngx_LOqe,Ngy_LOqe]=size(dset_qe);
    Ngyby2qe=round(Ngy_LOqe/2);
    LOqe_dat=dset_qe(:,Ngyby2qe);
    
  
    ax1=subplot(121);
%     a1=dset_b3(:,150);
%     plot(a1);
    imagesc(dset_qe','XData',[x1minqe x1maxqe],...
        'YData',[x2minqe x2maxqe]);
    set(gca,'YDir','normal');

    axis square;
    set(gca,'FontSize',12,'FontWeight','Bold');
    shading('interp');
    colormap(ax1,'jet(4098)');
    xlabel('x1(c/\omega_{pe})'); ylabel('x2(c/\omega_{pe})');
    xlim([0 25])
    cb = colorbar('south'); 
    set(cb,'position',[.13 .35 .335 .02]);
    title(strcat('e-Density:',num2str(time_qe)));
    
    
    ax2= subplot(122);
    imagesc(dset_b3','XData',[x1minb3 x1maxb3],...
        'YData',[x2minb3 x2maxb3]);
    set(gca,'YDir','normal');

    axis square;
    set(gca,'FontSize',12,'FontWeight','Bold');
    colormap(ax2,'jet(4098)'); shading('interp');
    xlabel('x1(c/\omega_{pe})'); %ylabel('x2(c/\omega_{pe})');
    xlim([0 25]);
    cb = colorbar('south'); 
    set(cb,'position',[.57 .35 .335 .02]);
    title(strcat('B-field:',num2str(time_b3)));
    
    drawnow
    saveas(gcf,'CurrentFig.png');
    img4=imread('CurrentFig.png');
    writeVideo(writerObj,img4);
     L = graphicsversion('handlegraphics'); 
     currFrame = getframe;

 end

  close(writerObj);