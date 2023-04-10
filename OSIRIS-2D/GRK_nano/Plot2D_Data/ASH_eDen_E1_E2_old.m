clc;
clear all;
close all;

addpath('UtilFun');

% writerObj = VideoWriter('eDenE1E2P1X1_Ionization250_20Nov2014_2.avi');
% writerObj.FrameRate=2;
% open(writerObj);

DataFolder={'WithDLA_PreformedPlasma_3June2015';...
    'NairGupta_Ionization250_20Nov2014';...
    'DLA_PreformedPlasma_24Nov2014';...
    'NairGupta_Preformed250_14Oct2014';...
    'NairGupta_Ionization250_14Oct2014';...
    'NairGupta_28Mar2014';...
    'NairGupta_Ionization250_13Oct2014';...
    'NairGupta_Preformed250FixedIon_30Sep2014';...
    'NairGupta_Preformed250_25Sep2014'};

dat_pth0='/media/SIMULATION/BhaveshData/';

ne=6.12e18; % cm^-3 58.0e18
frm=190;
DirNo=1;

minch=-4;
maxch=0;

pth0='/media/SIMULATION/BhaveshData/';
pth1='/MS/DENSITY/He_electrons/charge/';
pth2='/MS/FLD/e2/';
pth3='/MS/FLD/e1/';
pth4='/MS/RAW/He_electrons/';

qe=-4.8032e-10; %electron charge , statcoulomb
me=9.1094e-28; %electron mass , g
vel_c=2.9979e10; %velocity of light ,  cm/sec
omp_e=sqrt((4*pi*ne*qe^2)/me); x_nor=vel_c/omp_e; t_nor=1/omp_e;

len_fr=length(frm);
scrsz = get(0,'ScreenSize');
hfig=figure('Position',[50 10 1200 650]);
for ii=1:length(frm)    
    frno=frm(ii);
    for jj=DirNo:DirNo
        
        fold_pth1=strcat(pth0,DataFolder{jj},pth1);
        fold_pth2=strcat(pth0,DataFolder{jj},pth2);
        fold_pth3=strcat(pth0,DataFolder{jj},pth3);
        fold_pth4=strcat(pth0,DataFolder{jj},pth4);
        
        [xg_cha,yg_cha,dset_qe,x1lt_qe,x2lt_qe,time_qe]=...
            ReadECharge13Aug2014(frno,fold_pth1);
        x1min_qe=x1lt_qe(1); x1max_qe=x1lt_qe(2);
        x2min_qe=x2lt_qe(1); x2max_qe=x2lt_qe(2);
        [Ngx_cha,Ngy_cha]=size(dset_qe);
        Ngyby2=round(Ngy_cha/2);
        dy_cha=yg_cha(2)-yg_cha(1);
        
        [xg_e2,yg_e2,time_e2,dset_e2,x1lt_e2,x2lt_e2]=...
            ReadE213Aug2014(frno,fold_pth2);
        x1mine2=x1lt_e2(1); x1maxe2=x1lt_e2(2);
        x2mine2=x2lt_e2(1); x2maxe2=x2lt_e2(2);
        [Ngx_LOe2,Ngy_LOe2]=size(dset_e2);
        [xIndMax_e2,yIndMax_e2]=ind2sub(size(dset_e2),...
            find(dset_e2==max(max(dset_e2,[],1),[],2),1));
        Ngyby2e2=yIndMax_e2;
        LOe2_dat=dset_e2(:,Ngyby2e2);
        maxe2=max(LOe2_dat); mine2=min(LOe2_dat);
        
        [xg_e1,yg_e1,time_e1,dset_e1,x1lt_e1,x2lt_e1]=...
            ReadE113Aug2014(frno,fold_pth3);
        x1mine1=x1lt_e1(1); x1maxe1=x1lt_e1(2);
        x2mine1=x2lt_e1(1); x2maxe1=x2lt_e1(2);
        [Ngx_LOe1,Ngy_LOe1]=size(dset_e1);
        [xIndMax_e1,yIndMax_e1]=ind2sub(size(dset_e1),...
            find(dset_e1==max(max(dset_e1,[],1),[],2),1));
        Ngyby2e1=yIndMax_e1;
        LOe1=dset_e1(:,Ngyby2e1);
        maxe1=max(LOe1); mine1=min(LOe1);
        
        maxe=maxe1; mine=mine1;
%         maxe=max(maxe2,maxe1); mine=min(mine2,mine1);
        
        [ene_raw,x1_raw,p1_raw,x2_raw,~,~,xlt_in,xlt_fi]=...
            ReadRaw13Aug2014(frno,fold_pth4);
        ene_th_loc=ene_raw>20;
 
        max_p1=max(p1_raw);
        
        xlim1=x1min_qe; xlim2=x1max_qe;
        ylim1=x2min_qe; ylim2=x2max_qe;
        
        ind_ylim1=round(ylim1/dy_cha)+1; 
        ind_ylim2=round(ylim2/dy_cha)-1;
        dset_qe_plt=dset_qe(:,ind_ylim1:ind_ylim2);
        yg_cha_plt=yg_cha(ind_ylim1:ind_ylim2);
        LT_a=(maxe-mine)/( ylim2-ylim1);
        LT_b=mine-LT_a*ylim1;
        yg_cha_plt=LT_a*yg_cha_plt+LT_b;
        ylim1cha=LT_a*ylim1+LT_b;
        ylim2cha=LT_a*ylim2+LT_b;
        
%         xlim1=x1max_qe-60; xlim2=x1max_qe-0;
        ylim1=x2mine2; ylim2=x2maxe2;
        ind_ylim1=round(ylim1/dy_cha)+1; 
        ind_ylim2=round(ylim2/dy_cha)-1;
        dset_e2_plt=dset_e2(:,ind_ylim1:ind_ylim2);
        yg_e2_plt=yg_e2(ind_ylim1:ind_ylim2);
        LT_a=(maxe-mine)/( ylim2-ylim1);
        LT_b=mine-LT_a*ylim1;
        yg_e2_plt=LT_a*yg_e2_plt+LT_b;
%         ylim1=LT_a*ylim1+LT_b;
%         ylim2=LT_a*ylim2+LT_b;

        subtightplot(1,1,1,[0,0.05],0.1,0.05) 
        [AX,H1,H2]=plotyy(xg_e1,LOe1,xg_e2,LOe2_dat);hold on;
        set(AX(1),'xlim',[xlim1 xlim2]);
        set(AX(2),'xlim',[xlim1 xlim2]);
        set(AX(1),'ylim',[ylim1cha ylim2cha]);
%         set(AX(2),'ylim',[0 80]);
        set(get(AX(2),'YLabel'),'String','Laser(mc\omega / e)');
        set(get(AX(1),'YLabel'),'String','WakeField(mc\omega / e)');
        set(AX(1),'GridLineStyle','-');
        set(H1,'linewidth',2);
        set(H1,'color','b');
        set(H2,'linewidth',2);
        set(H2,'color','g');
        grid on; hold on
        
      
        Himg=imagesc(flipdim((dset_qe_plt),2)',...
            'XData',[xg_cha xg_cha],...
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
            strcat('time=',num2str(time_qe),'(1/\omega_{pe})'));
        title(str_tmp1);
        xlim([xlim1 xlim2]); ylim([ylim1cha ylim2cha]);
        xlabel('x1(c/\omega_{pe})'); 
%         ylabel('x2(c/\omega_{pe})');
        hold on;
         
        Himg=imagesc(flipdim((dset_e2_plt),2)',...
            'XData',[xg_e2 xg_e2],...
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
 
    end
    
    drawnow;
    
%     gif_add_frame(gcf,'eDenE1E2P1X1_Ionization250_20Nov2014.gif');

%     saveas(gcf,'CurrentFig.png');
%     img4=imread('CurrentFig.png');
%     writeVideo(writerObj,img4);

end

% close(writerObj);