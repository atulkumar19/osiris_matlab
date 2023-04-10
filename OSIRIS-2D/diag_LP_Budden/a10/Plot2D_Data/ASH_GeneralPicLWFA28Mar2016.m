clc;
clear all;
close all;

addpath('UtilFun');

ne=5.8e+19; % cm^-3 58.0e18, 6.12e+18
frm=105;%35 110  250
% Get the width and height of the figure
%============AVI-FILE-NAME==========================
% movObj = QTWriter('GeneralPicture.mov');
% movObj.FrameRate = 10;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

minch=-2.5;
maxch=0;

pth0='/home/bhavesh/Desktop/OSIRIS/OsirisRao/EffectRamp/Micron150';
pth1='/MS/DENSITY/He_electrons/charge/';
pth2='/MS/FLD/e2/';
pth3='/MS/FLD/e1/';
pth4='/MS/RAW/He_electrons/';

str_h1='charge-He_electrons-';
str_h2='e2-';
str_h3='e1-';
str_h4='RAW-He_electrons-';

str_ext='.h5';
str_num='000000';
len_str_num=length(str_num);


qe=-4.8032e-10; %electron charge , statcoulomb
me=9.1094e-28; %electron mass , g
vel_c=2.9979e10; %velocity of light ,  cm/sec
omp_e=sqrt((4*pi*ne*qe^2)/me); x_nor=vel_c/omp_e; t_nor=1/omp_e;



xlab='time' ;              % x-axis title
pas=5;           % Number of ticks per axe
tail=2;             % font size
mark1='d';c1='g';    
mark2='<';c2='b';
mark3='s';c3='r';
mark4='p';c4='k';
mark5='o';c5='m';
mark6='*';c6='y';
fontname='courrier';
fontsize=9.2;
fontsize_ax=8;
len_fr=length(frm);
scrsz = get(0,'ScreenSize');
hfig=figure('Position',[50 10 1500 900]);
for ii=1:length(frm)    
    
    frno=frm(ii);
    str_frno=num2str(frno);
    len_frno=length(str_frno);
    str_num((len_str_num-len_frno+1):end)=str_frno;

    fl_nm1=strcat(pth0,pth1,str_h1,str_num,str_ext); %charge density
    fl_nm2=strcat(pth0,pth2,str_h2,str_num,str_ext); %raw data
    fl_nm3=strcat(pth0,pth3,str_h3,str_num,str_ext); %raw data
    fl_nm4=strcat(pth0,pth4,str_h4,str_num,str_ext); %raw data

    [xg,yg,dset_qe,x1lt,x2lt,time]=AshReadHDF5DenDat(fl_nm1);
    [~,~,dset_e2,~,~,~]=AshReadHDF5DenDat(fl_nm2);
    [~,~,dset_e1,~,~,~]=AshReadHDF5DenDat(fl_nm3);
    [ene_raw,x1_raw,p1_raw,x2_raw,~,~,~,~]=AshReadRaw(fl_nm4);
    
    d_mm=time*t_nor*vel_c*10;

    x1min=x1lt(1); x1max=x1lt(2);
    x2min=x2lt(1); x2max=x2lt(2);
    [Ngx,Ngy]=size(dset_qe);
    Ngyby2=round(Ngy/2);
    dy=yg(2)-yg(1);
    
    LOqe=dset_qe(:,Ngyby2);
    maxqe=max(LOqe); minqe=min(LOqe);
    
    [xIndMax_e2,yIndMax_e2]=ind2sub(size(dset_e2),...
        find(dset_e2==max(max(dset_e2,[],1),[],2),1));
    Ngyby2e2=yIndMax_e2;
    LOe2=dset_e2(:,Ngyby2e2);
    maxe2=max(LOe2); mine2=min(LOe2);
    
    [xIndMax_e1,yIndMax_e1]=ind2sub(size(dset_e1),...
        find(dset_e1==max(max(dset_e1,[],1),[],2),1));
    Ngyby2e1=yIndMax_e1;
    LOe1=dset_e1(:,Ngyby2e1);
    maxe1=max(LOe1); mine1=min(LOe1);
    
    
    maxabse=max(abs([maxe1;mine1;maxe2;mine2;maxqe;minqe])); 
    maxe=maxabse;  mine=-maxabse;  
    etick=linspace(mine,maxe,5);
    
    ylim1=x2min; ylim2=x2max;
    ind_ylim1=round(x2min/dy)+1; 
    ind_ylim2=round(x2max/dy)-1;
    dset_qe_plt=dset_qe(:,ind_ylim1:ind_ylim2);
    yg_cha_plt=yg(ind_ylim1:ind_ylim2);
    LT_a=(maxe-mine)/( x2max-x2min);
    LT_b=mine-LT_a*x2min;
    yg_cha_plt=LT_a*yg_cha_plt+LT_b;
    ylim1cha=LT_a*x2min+LT_b;
    ylim2cha=LT_a*x2max+LT_b;
       
    
    ene_th_loc=ene_raw>10;
    max_p1=max(p1_raw); 
%%%%%% this section plots E1X1 and P1X1 %%%%%%%%%%%%%%%%%%%%%%%    

    if ~any(ene_th_loc)
        [AX,H1,H2]=plotyy(xg,LOe1,...
        x1_raw(1),ene_raw(1)); hold on;
    else
        [AX,H1,H2]=plotyy(xg,LOe1,...
            x1_raw(ene_th_loc),ene_raw(ene_th_loc)); hold on;
    end
    
    hch=plot(AX(1),xg,LOqe); hold on;
    plot(AX(1),xg,LOe2);   hold off;
    
    set(AX(1),'xlim',[x1min x1max]);
    set(AX(2),'xlim',[x1min x1max]);
    set(AX(1),'ylim',[ylim1cha ylim2cha]);
    set(AX(1),'ytick',etick);
    set(AX(2),'ylim',[10 110]);
    set(AX(2),'ytick',linspace(10,110,10));
    set(get(AX(2),'YLabel'),'String','Energy(mc^2)');
    set(get(AX(1),'YLabel'),'String','E-Field(mc\omega / e)');
    set(get(AX(1),'YLabel'),'color','k');
    set(AX(1),'GridLineStyle','-');
    set(H1,'linewidth',2);
    set(H1,'color','r');
    set(H2,'marker','o','linestyle','none',...
        'MarkerFaceColor','g','MarkerSize',4);
    set(hch,'linewidth',2);
    set(hch,'color','k');
    grid on; hold on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    Himg=imagesc(flipdim((dset_qe_plt),2)',...
        'XData',[xg xg],...
        'YData',[yg_cha_plt yg_cha_plt]); 
    set(Himg,'AlphaData',0.6);
    shading('interp');
    colormap(cmap('orange',1024));
    caxis(gca,[minch maxch]); 
    set(gca,'ydir','normal');
    grid on;
    freezeColors;
    cbfreeze(colorbar('north'));
    str_tmp1=strcat('Frame:',num2str(frno),...
    '     : ',...
    strcat('time=',num2str(time),'(1/\omega_{pe})'),...
    '     : ',...
    strcat('Dist(mm) =',num2str(d_mm),'mm'));
    title(str_tmp1);
    xlim([x1min x1max]); ylim([ylim1cha ylim2cha]);
    xlabel('x1(c/\omega_{pe})'); 
%         ylabel('x2(c/\omega_{pe})');
    
    hold off
    
    drawnow;
    

% frame=getframe(gcf); 
% writeMovie(movObj,frame);

end


% close(movObj);
