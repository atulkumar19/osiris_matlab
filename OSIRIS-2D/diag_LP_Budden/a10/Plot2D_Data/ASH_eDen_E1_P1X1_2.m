clc;
clear all;
close all;

addpath('UtilFun');

ne=5.8e+19; % cm^-3 58.0e18, 6.12e+18
frm=200:600;

minch=-4;
maxch=0;

pth0='/home/bhavesh/Desktop/OSIRIS/OsirisRao/EffectRamp/Micron800';
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

len_fr=length(frm);
scrsz = get(0,'ScreenSize');
hfig=figure('Position',[50 10 1200 650]);
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

    x1min=x1lt(1); x1max=x1lt(2);
    x2min=x2lt(1); x2max=x2lt(2);
    [Ngx,Ngy]=size(dset_qe);
    Ngyby2=round(Ngy/2);
    dy=yg(2)-yg(1);
    
    [xIndMax_e2,yIndMax_e2]=ind2sub(size(dset_e2),...
        find(dset_e2==max(max(dset_e2,[],1),[],2),1));
    Ngyby2e2=yIndMax_e2;
    LOe2_dat=dset_e2(:,Ngyby2e2);
    maxe2=max(LOe2_dat); mine2=min(LOe2_dat);
    
    [xIndMax_e1,yIndMax_e1]=ind2sub(size(dset_e1),...
        find(dset_e1==max(max(dset_e1,[],1),[],2),1));
    Ngyby2e1=yIndMax_e1;
    LOe1=dset_e1(:,Ngyby2e1);
    maxe1=max(LOe1); mine1=min(LOe1);
    maxabse=max(abs(maxe1 ),abs(mine1));
    maxe=maxabse; mine=-maxabse;
    etick=linspace(mine,maxe,8);
    
    ene_th_loc=ene_raw>10;
    max_p1=max(p1_raw);
    
    xlim1=x1min; xlim2=x1max;
    ylim1=x2min; ylim2=x2max;

    ind_ylim1=round(ylim1/dy)+1; 
    ind_ylim2=round(ylim2/dy)-1;
    dset_qe_plt=dset_qe(:,ind_ylim1:ind_ylim2);
    yg_cha_plt=yg(ind_ylim1:ind_ylim2);
    LT_a=(maxe-mine)/( ylim2-ylim1);
    LT_b=mine-LT_a*ylim1;
    yg_cha_plt=LT_a*yg_cha_plt+LT_b;
    ylim1cha=LT_a*ylim1+LT_b
    ylim2cha=LT_a*ylim2+LT_b

    ylim1=x2min; ylim2=x2max;
    ind_ylim1=round(ylim1/dy)+1; 
    ind_ylim2=round(ylim2/dy)-1;
    dset_e2_plt=dset_e2(:,ind_ylim1:ind_ylim2);
    yg_e2_plt=yg(ind_ylim1:ind_ylim2);
    LT_a=(maxe-mine)/( ylim2-ylim1);
    LT_b=mine-LT_a*ylim1;
    yg_e2_plt=LT_a*yg_e2_plt+LT_b;

    
    if ~any(ene_th_loc)
        [AX,H1,H2]=plotyy(xg,LOe1,...
        x1_raw(1),ene_raw(1));
    else
        [AX,H1,H2]=plotyy(xg,LOe1,...
            x1_raw(ene_th_loc),ene_raw(ene_th_loc));
    end
    set(AX(1),'xlim',[xlim1 xlim2]);
    set(AX(2),'xlim',[xlim1 xlim2]);
    set(AX(1),'ylim',[ylim1cha ylim2cha]);
    set(AX(1),'ytick',etick);
    set(AX(2),'ylim',[10 80]);
    set(AX(2),'ytick',linspace(10,80,8));
    set(get(AX(2),'YLabel'),'String','Energy(mc^2)');
    set(get(AX(1),'YLabel'),'String','WakeField(mc\omega / e)');
    set(AX(1),'GridLineStyle','-');
    set(H1,'linewidth',2);
    set(H1,'color','r');
    set(H2,'marker','o','linestyle','none',...
        'MarkerFaceColor','g','MarkerSize',4);
    grid on; hold on

    Himg=imagesc(flipdim((dset_qe_plt),2)',...
        'XData',[xg xg],...
        'YData',[yg_cha_plt yg_cha_plt]); 
    set(Himg,'AlphaData',0.9);
    shading('interp');
    colormap(cmap('orange',1024));
    caxis(gca,[minch maxch]); 
    set(gca,'ydir','normal');
    grid on;
    freezeColors;
%         cbfreeze(colorbar('north'));
    str_tmp1=strcat('Frame:',num2str(frno),...
        '     : ',...
        strcat('time=',num2str(time),'(1/\omega_{pe})'));
    title(str_tmp1);
    xlim([xlim1 xlim2]); ylim([ylim1cha ylim2cha]);
    xlabel('x1(c/\omega_{pe})'); 
%         ylabel('x2(c/\omega_{pe})');
    hold on;

    Himg=imagesc(flipdim((dset_e2_plt),2)',...
        'XData',[xg xg],...
        'YData',[yg_e2_plt yg_e2_plt]); 
    set(Himg,'AlphaData',0.35);
    shading('interp');
    newmap=b2r(mine2,maxe2);
    colormap(newmap);
    set(gca,'ydir','normal');
    grid on;
    freezeColors;
%         cbfreeze(colorbar('west'));
    xlim([xlim1 xlim2]); ylim([ylim1cha ylim2cha]);
    hold off;
    
    drawnow;
    
%     gif_add_frame(gcf,'eDenE1E2P1X1_WithoutDLA_28May2015.gif');

end