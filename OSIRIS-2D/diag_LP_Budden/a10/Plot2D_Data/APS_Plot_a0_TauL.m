clc;
clear all;
close all;

addpath('UtilFun');

writerObj = VideoWriter('a0TauL_Ionization250_20Nov2014_2.avi');
writerObj.FrameRate=2;
open(writerObj);

DataFolder={'NairGupta_Ionization250_20Nov2014';...
    'NairGupta_Preformed250_14Oct2014';...
    'NairGupta_Ionization250_14Oct2014';...
    'NairGupta_28Mar2014';...
    'NairGupta_Ionization250_13Oct2014';...
    'NairGupta_Preformed250FixedIon_30Sep2014';...
    'NairGupta_Preformed250_25Sep2014'};

dat_pth0='/media/SIMULATION/BhaveshData/';

ne=58.0e18; % cm^-3
LamLas=800*1e-9*1e2;
frm=0:100;
DirNo=1;

minch=-5;
maxch=0;

pth0='/media/SIMULATION/BhaveshData/';
pth1='/MS/DENSITY/He_electrons/charge/';
pth2='/MS/FLD/e2/';

qe=-4.8032e-10; %electron charge , statcoulomb
me=9.1094e-28; %electron mass , g
vel_c=2.9979e10; %velocity of light ,  cm/sec
omp_e=sqrt((4*pi*ne*qe^2)/me) 
x_nor=vel_c/omp_e; t_nor=1/omp_e;
lam_plasma=2*pi*x_nor;
nor_lam_plasma=lam_plasma/x_nor;
omLas=2*pi*vel_c/LamLas

len_fr=length(frm);
avec=zeros(len_fr,1); tvec=zeros(len_fr,1); tauvec=zeros(len_fr,1);
scrsz = get(0,'ScreenSize');
figure('Position',[50 10 1200 650]);
for ii=1:length(frm)    
    frno=frm(ii);
    for jj=DirNo:DirNo
        
        fold_pth1=strcat(pth0,DataFolder{jj},pth1);
        fold_pth2=strcat(pth0,DataFolder{jj},pth2);
        
        
        [xg_cha,yg_cha,dset_qe,x1lt_qe,x2lt_qe,time_qe]=...
            ReadECharge13Aug2014(frno,fold_pth1);
        x1min_qe=x1lt_qe(1); x1max_qe=x1lt_qe(2);
        x2min_qe=x2lt_qe(1); x2max_qe=x2lt_qe(2);
        [Ngx_cha,Ngy_cha]=size(dset_qe);
        Ngyby2=round(Ngy_cha/2);
        tvec(ii)=time_qe;
        xlim1=x1max_qe-60; xlim2=x1max_qe-20;
        ylim1=42.5; ylim2=55;
        
        
        [xg_e2,yg_e2,time_e2,dset_e2,x1lt_e2,x2lt_e2]=...
            ReadE213Aug2014(frno,fold_pth2);
        x1mine2=x1lt_e2(1); x1maxe2=x1lt_e2(2);
        x2mine2=x2lt_e2(1); x2maxe2=x2lt_e2(2);
        [Ngx_LOe2,Ngy_LOe2]=size(dset_e2);
%         Ngyby2e2=round(Ngy_LOe2/2)
%         LOe2_dat=dset_e2(:,Ngyby2e2);
        dx_e2=abs(xg_e2(2)-xg_e2(1));
        [xIndMax_e2,yIndMax_e2]=ind2sub(size(dset_e2),...
            find(dset_e2==max(max(dset_e2,[],1),[],2),1));
        Ngyby2e2=yIndMax_e2;
        LOe2_dat=dset_e2(:,Ngyby2e2);
        y_w0=dset_e2(xIndMax_e2,:);
%         plot(yg_e2,y_w0); drawnow;
        
        %%%%%% calculation of FFT of LO of E2-field%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        fs=1/dx_e2; %---->sampling frequency or sampling rate
        nyq=fs/2; %---->Nyquist frequency
        fft_LOe2=fft(LOe2_dat);
        Nfft=length(fft_LOe2);
        df=fs/Nfft;
        fvec=(0:(Nfft-1))*df;
        fvec(fvec>nyq) = fvec(fvec>nyq)-fs;
        freqs=fvec(fvec>=0);
        ampSpec=abs(fft_LOe2(fvec>=0));
        ampSpec=ampSpec/Nfft;
        ampSpec(2:ceil(Nfft/2))=2*ampSpec(2:ceil(Nfft/2));
        %%%%%% FFT of LO of E2-field%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%% calculation of a0%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tvec=xg_e2-time_e2;
        int_vec=3:length(tvec); szintvec=length(int_vec);
        a0tmp=zeros(szintvec,1); timetmp=zeros(szintvec,1);
        for kk=1:szintvec
            timetmp(kk)=tvec(int_vec(kk));
            a0tmp(kk)=trapz(tvec(1:int_vec(kk)),...
                LOe2_dat(1:int_vec(kk)),1);
        end
        avec(ii)=max(a0tmp);

%         a0=trapz(yg_cha,trapz(xg_e2,dset_e2,1),2);
%         a0=trapz(xg_e2,LOe2_dat,1); %incorrect
%         a0=trapz(xg_e2,mean(dset_e2,2),1); %incorrect
        %%%%%% calculation of a0%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%% calculation of \tau_L%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        SlopeThreshold=1e-3; AmpThreshold=(1/100)*max(LOe2_dat); 
        SmoothWidth=7e1*dx_e2; PeakGroup=3; Smoothtype=3;

        [pks] = findpeaks(xg_e2,LOe2_dat,SlopeThreshold,AmpThreshold,...
            SmoothWidth,PeakGroup,Smoothtype);
        x_gfit=pks(:,2); f_gfit=pks(:,3);hf_max=0.5*max(f_gfit);
        pp=spline(x_gfit,f_gfit-hf_max);
        f_sp=ppval(pp,x_gfit);
        
        z=fnzeros(pp,[x_gfit(1) x_gfit(end)]);
        tauvec(ii)=0.5*(z(1,2)+z(2,2))-0.5*(z(1,1)+z(2,1));
        
        
        x_gfit=pks(:,2); 
        f_gfit=pks(:,3);
        am_g=max(f_gfit); si_g=0.25*abs(x_gfit(end)-x_gfit(1)); 
        xc_g=0.5*(x_gfit(end)+x_gfit(1));
        par_gg=[am_g;si_g;xc_g];
        par_gg_lb=[1;2;10];
        par_gg_ub=[20;10;300];
        option=optimset('maxiter',1500,'FunValCheck','off');
        % [c,r,j]=nlinfit(x_gfit,f_gfit,'multi_gauss',par_g,option);
        par_go = lsqcurvefit('multi_gauss',par_gg,...
            x_gfit,f_gfit,par_gg_lb,par_gg_ub);
        f_fit=multi_gauss(par_go,xg_e2');
        PulseDuration=2*sqrt(2*log(2))*par_go(2);
        tau_L=PulseDuration*x_nor/vel_c

        
        %%%%%% calculation of \tau_L%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        subplot(2,3,1);
        imagesc(flipdim((dset_qe),2)',...
            'XData',[xg_cha xg_cha],...
            'YData',[yg_cha yg_cha]);     
        shading('interp');
        colormap(cmap('blue',256));
        caxis(gca,[minch maxch]); 
        set(gca,'ydir','normal');
        grid on;
        freezeColors;
        cbfreeze(colorbar('north'));
        str_tmp1=strcat('frame=',num2str(frno),...
            '     : ',...
            strcat(num2str(time_qe*t_nor/1e-15),'fs'),...
            '     : ',...
            strcat(num2str(time_qe*t_nor*vel_c*1e4,'%10.3e'),'\mu m'));
        %{
        str_tmp2=strcat('dist(\mu m)=',...
            num2str(time_qe*t_nor*vel_c*1e4,'%10.3e'));
        str_tmp3=strcat('xnor=', num2str(x_nor,'%10.3e'));
        str_tmp4=strcat('tnor=', num2str(t_nor,'%10.3e'));
        
        str_ti2=char(str_tmp2,str_tmp3,str_tmp4);
        text(xlim1,0.5*(ylim1+ylim2),str_ti2);
        %}
        title(str_tmp1);
        ylim([ylim1 ylim2]); 
        xlim([xlim1 xlim2]);
        xlabel('x1(c/\omega_{pe})'); ylabel('x2(c/\omega_{pe})');
        
        subplot(2,3,4);
        plot(x_gfit,f_sp+hf_max); hold on;
        plot(xg_e2,LOe2_dat,'r'); hold off;
        xlim([xlim1 xlim2]);
        xlabel('x1(c/\omega_{pe})'); 
        ylabel('E2(mc\omega_{pe}/e)');
        grid on;
                
        subplot(2,3,5);
        plot(tvec(1:ii),avec(1:ii));
        grid on;
        xlabel('t(1/\omega_{pe})'); ylabel('a0');
        
        subplot(2,3,6);
        if ii<=35
            tmpii=ii;
        else
            tmpii=35;
        end
        plot(tvec(1:tmpii),tauvec(1:tmpii)/(2*pi));
        grid on;
        xlabel('t(1/\omega_{pe})'); ylabel('\tau/\lambda_{pe}');
        
        
        subplot(2,3,[2 3]);
        plot(yg_e2,y_w0);grid on; 
        
        
        
        drawnow;
        
    end
    
    
    
    saveas(gcf,'CurrentFig.png');
    img4=imread('CurrentFig.png');
    writeVideo(writerObj,img4);

end

close(writerObj);