clc
clear all
close all

qe=-4.8032e-10; %electron charge , statcoulomb
me=9.1094e-28; %electron mass , g
vel_c=2.9979e10; %velocity of light ,  cm/sec

addpath('UtilFun');
pth0='/media/SIMULATION/BhaveshData/';
pth1='/MS/DENSITY/He_electrons/charge/';
pth2='/MS/FLD/e2/';
pth3='/MS/FLD/e1/';

%%%%%%% Input section starts here %%%%%%%%%%%%%%%%%%%%
DataFolder={'NairGupta_Ionization250_20Nov2014';...
    'NairGupta_Ionization250_14Oct2014';...
    'NairGupta_Preformed250_14Oct2014';...
    'NairGupta_28Mar2014';...
    'NairGupta_Ionization250_13Oct2014';...
    'NairGupta_Preformed250FixedIon_30Sep2014';...
    'NairGupta_Preformed250_25Sep2014'};


DirNo=1;
ne=58.0e18; %58.0e18 cm^-3 
frm=0;

minch=-5;
maxch=0;

% writerObj = VideoWriter('Wig_Ionization250_20Nov2014.avi');
% writerObj.FrameRate=2;
% open(writerObj);

%%%%%%% Input section end here %%%%%%%%%%%%%%%%%%%%

omp_e=sqrt((4*pi*ne*qe^2)/me); x_nor=vel_c/omp_e; 
t_nor=1/omp_e;
len_fr=length(frm);
scrsz = get(0,'ScreenSize');
figure('Position',[50 10 1200 650]);

for ii=1:length(frm)    
    frno=frm(ii);
    for jj=DirNo:DirNo
        fold_pth1=strcat(pth0,DataFolder{jj},pth1);
        fold_pth2=strcat(pth0,DataFolder{jj},pth2);
        fold_pth3=strcat(pth0,DataFolder{jj},pth3);

        [xg_e2,yg_e2,time_e2,dset_e2,x1lt_e2,x2lt_e2]=...
            ReadE213Aug2014(frno,fold_pth2);
        dx=xg_e2(2)-xg_e2(1);sps=1/dx;
        x1mine2=x1lt_e2(1); x1maxe2=x1lt_e2(2);
        x2mine2=x2lt_e2(1); x2maxe2=x2lt_e2(2);
        [Ngx_LOe2,Ngy_LOe2]=size(dset_e2);
        [xIndMax_e2,yIndMax_e2]=ind2sub(size(dset_e2),...
            find(dset_e2==max(max(dset_e2,[],1),[],2),1));
        Ngyby2e2=yIndMax_e2;
        LOe2=dset_e2(:,Ngyby2e2);
        maxe2=max(LOe2); mine2=min(LOe2);
%         y_w0=dset_e2(xIndMax_e2,:);
        [wv,wvfreq,wvtime]=wvdc(LOe2,1,length(LOe2)/2,sps);
        wvfreq=2*pi*wvfreq;
        
        
        [xg_e1,yg_e1,time_e1,dset_e1,x1lt_e1,x2lt_e1]=...
            ReadE113Aug2014(frno,fold_pth3);
        dx=xg_e1(2)-xg_e1(1);sps=1/dx;
        x1mine1=x1lt_e1(1); x1maxe1=x1lt_e1(2);
        x2mine1=x2lt_e1(1); x2maxe1=x2lt_e1(2);
        LOe1=dset_e1(:,Ngyby2e2);
        maxe1=max(LOe1); mine1=min(LOe1);
        
        maxe=max(maxe1,maxe2); mine=min(mine1,mine2);

        [xg_cha,yg_cha,dset_qe,x1lt_qe,x2lt_qe,time_qe]=...
            ReadECharge13Aug2014(frno,fold_pth1);
        dy_cha=yg_cha(2)-yg_cha(1);
        x1min_qe=x1lt_qe(1); x1max_qe=x1lt_qe(2);
        x2min_qe=x2lt_qe(1); x2max_qe=x2lt_qe(2);
        [Ngx_cha,Ngy_cha]=size(dset_qe);
        LOqe=dset_qe(:,Ngyby2e2);
        
        
        xlim1=x1max_qe-60; xlim2=x1max_qe-10;
        ylim1=44; ylim2=54;
        
        ind_ylim1=round(ylim1/dy_cha)+1; 
        ind_ylim2=round(ylim2/dy_cha)-1;
        dset_qe_plt=dset_qe(:,ind_ylim1:ind_ylim2);
        yg_cha_plt=yg_cha(ind_ylim1:ind_ylim2);
        LT_a=(maxe-mine)/( ylim2-ylim1);
        LT_b=mine-LT_a*ylim1;
        yg_cha_plt=LT_a*yg_cha_plt+LT_b;
        ylim1=LT_a*ylim1+LT_b;
        ylim2=LT_a*ylim2+LT_b;
        
%         xlim1=x1min_qe; xlim2=x1max_qe;
%         ylim1=x2min_qe; ylim2=x2max_qe;

        subtightplot(1,2,1,[0,0.05],0.1,0.05);
        plot(xg_e2,LOe2);
        xlim([xlim1 xlim2]); ylim([ylim1 ylim2]);
        grid on; hold on
        
        Himg=imagesc(flipdim((dset_qe_plt),2)',...
            'XData',[xg_cha xg_cha],...
            'YData',[yg_cha_plt yg_cha_plt]); 
        set(Himg,'AlphaData',0.5);
        shading('interp');
        colormap(cmap('orange',256));
        caxis(gca,[minch maxch]); 
        set(gca,'ydir','normal');
        grid on;
        freezeColors;
%         cbfreeze(colorbar('north'));
        str_tmp1=strcat('frame=',num2str(frno),...
            '     : ',...
            strcat(num2str(time_qe),'1/\omega_{pe}'),...
            '     : ',...
            strcat(num2str(time_qe*t_nor*vel_c*1e4,'%10.3e'),'\mu m'));
        title(str_tmp1);
        xlim([xlim1 xlim2]); 
        ylim([ylim1 ylim2]);
        xlabel('x1(c/\omega_{pe})'); ylabel('x2(c/\omega_{pe})');
        hold off
      
        
        subtightplot(1,2,2,[0.1,0.05],0.1,0.05);
        imagesc(xg_e2,wvfreq,wv); shading('interp');
        axis xy;colormap('jet'); colorbar('north');
        xlim([xlim1 xlim2]); ylim([2.5 9.5]);
        xlabel('x1(c/\omega_{pe})'); ylabel('k_x(\omega_{pe}/c)');
        grid on;
        
        drawnow;
        
    end
    
%     gif_add_frame(gcf,'Wig_Ionization250_20Nov2014.gif');
    
%     saveas(gcf,'CurrentFig.png');
%     img4=imread('CurrentFig.png');
%     writeVideo(writerObj,img4);

end

% close(writerObj);
