clc;
clear all;
close all;

addpath('UtilFun');


ne=5.8e+19; % cm^-3 58.0e18, 6.12e+18
frm=0:95;

minch=-4;
maxch=0;

pth0='/home/bhavesh/Desktop/OSIRIS/OsirisRao/bin_data/PICSim_200MicronRamp';
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
LamLas=800*1e-9*1e2;


omp_e=sqrt((4*pi*ne*qe^2)/me) 
x_nor=vel_c/omp_e; t_nor=1/omp_e;
lam_plasma=2*pi*x_nor;
nor_lam_plasma=lam_plasma/x_nor;
omLas=2*pi*vel_c/LamLas

cind=cmap('orange',1024);
len_fr=length(frm);
avec=zeros(len_fr,1); timvec=zeros(len_fr,1); tauvec=zeros(len_fr,1);
SpotSz=zeros(len_fr,1);
scrsz = get(0,'ScreenSize');
figure('Position',[50 10 1200 650]);
for ii=1:length(frm)    
    frno=frm(ii);
    str_frno=num2str(frno);
    len_frno=length(str_frno);
    str_num((len_str_num-len_frno+1):end)=str_frno;

    fl_nm1=strcat(pth0,pth1,str_h1,str_num,str_ext); %charge density
    fl_nm2=strcat(pth0,pth2,str_h2,str_num,str_ext); %raw data
    
    [xg_cha,yg_cha,dset_qe,x1lt_qe,x2lt_qe,time_qe]=...
        AshReadHDF5DenDat(fl_nm1);
    [xg_e2,yg_e2,dset_e2,x1lt_e2,x2lt_e2,time_e2]=...
        AshReadHDF5DenDat(fl_nm2);
    
    dy_cha=yg_cha(2)-yg_cha(1);
    x1min_qe=x1lt_qe(1); x1max_qe=x1lt_qe(2);
    x2min_qe=x2lt_qe(1); x2max_qe=x2lt_qe(2);
    [Ngx_cha,Ngy_cha]=size(dset_qe);
    Ngyby2=round(Ngy_cha/2);
    timvec(ii)=time_qe;
    dx_cha=abs(xg_cha(2)-xg_cha(1));
    Ngyby2cha=Ngyby2;
    LOcha_x1=dset_qe(:,Ngyby2cha);
    
    x1mine2=x1lt_e2(1); x1maxe2=x1lt_e2(2);
    x2mine2=x2lt_e2(1); x2maxe2=x2lt_e2(2);
    [Ngx_LOe2,Ngy_LOe2]=size(dset_e2);
    dx_e2=abs(xg_e2(2)-xg_e2(1));
    [xIndMax_e2,yIndMax_e2]=ind2sub(size(dset_e2),...
        find(dset_e2==max(max(dset_e2,[],1),[],2),1));
    Ngyby2e2=yIndMax_e2;

%     Ngyby2e2=round(Ngy_LOe2/2);
    LOe2_dat=dset_e2(:,Ngyby2e2);
    maxe2=max(LOe2_dat); mine2=min(LOe2_dat);
    y_w0=dset_e2(xIndMax_e2,:);
    
    
    xlim1=x1min_qe; xlim2=x1max_qe;
    ylim1=x2min_qe; ylim2=x2max_qe;
    ind_ylim1=round(ylim1/dy_cha)+1; 
    ind_ylim2=round(ylim2/dy_cha)-1;
    dset_qe_plt=dset_qe(:,ind_ylim1:ind_ylim2);
    yg_cha_plt=yg_cha(ind_ylim1:ind_ylim2);
    LT_a=(maxe2-mine2)/( ylim2-ylim1);
    LT_b=mine2-LT_a*ylim1;
    yg_cha_plt=LT_a*yg_cha_plt+LT_b;
    ylim1cha=LT_a*ylim1+LT_b;
    ylim2cha=LT_a*ylim2+LT_b;

    ylim1=x2mine2; ylim2=x2maxe2;
    ind_ylim1=round(ylim1/dy_cha)+1; 
    ind_ylim2=round(ylim2/dy_cha)-1;
    dset_e2_plt=dset_e2(:,ind_ylim1:ind_ylim2);
    yg_e2_plt=yg_e2(ind_ylim1:ind_ylim2);
    LT_a=(maxe2-mine2)/( ylim2-ylim1);
    LT_b=mine2-LT_a*ylim1;
    yg_e2_plt=LT_a*yg_e2_plt+LT_b;
    
        %%%%%% calculation of a0%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tvec=xg_e2-time_e2;
    %{
    int_vec=2:length(tvec); szintvec=length(int_vec);
    a0tmp=zeros(szintvec,1); timetmp=zeros(szintvec,1);
    for kk=1:szintvec
        timetmp(kk)=tvec(int_vec(kk));
        a0tmp(kk)=trapz(tvec(1:int_vec(kk)),LOe2_dat(1:int_vec(kk)));
    end
    avec(ii)=max(a0tmp);
    %}
    avec(ii)=max(cumtrapz(tvec,LOe2_dat));

%{ 
    % //////////// wrong ways to calculate a0 ////////////////        
    a0=trapz(yg_cha,trapz(xg_e2,dset_e2,1),2);
    a0=trapz(xg_e2,LOe2_dat,1); %incorrect
    a0=trapz(xg_e2,mean(dset_e2,2),1); %incorrect
    %%%%% calculation of a0%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}        

    %%%%%% calculation of \tau_L%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SlopeThreshold=1e-4; AmpThreshold=(1/100)*max(LOe2_dat); 
    SmoothWidth=7e1*dx_e2; PeakGroup=3; Smoothtype=3;

    [pks] = findpeaks(xg_e2,LOe2_dat,SlopeThreshold,AmpThreshold,...
        SmoothWidth,PeakGroup,Smoothtype);
    x_gfit=pks(:,2); f_gfit=pks(:,3);hf_max=0.5*max(f_gfit);
    pp=spline(x_gfit,f_gfit-hf_max);
    f_sp=ppval(pp,x_gfit);

    tauvec(ii)= fwhm(x_gfit,f_gfit);

    %{
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
    %}
    %%%%%% calculation of \tau_L%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%% calculation of spot size%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SpotSz(ii)= fwhm(yg_e2,y_w0);
    %%%%%% calculation of spot size%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subplot(2,3,1);
    plot(timvec(1:ii),SpotSz(1:ii));grid on; 
    xlabel('t(1/\omega_{pe})'); ylabel('SpotSize(w,c/\omega_{pe}');


    subplot(2,3,2);
    plot(timvec(1:ii),avec(1:ii));
    grid on;
    xlabel('t(1/\omega_{pe})'); ylabel('a0');

    subplot(2,3,3);
    if ii<=35
        tmpii=ii;
    else
        tmpii=35;
    end
    plot(timvec(1:tmpii),tauvec(1:tmpii));
    grid on;
    xlabel('t(1/\omega_{pe})'); 
    ylabel('PulseLength-fwhm (c/\omega_{pe})');
    
    subplot(2,3,[4 5 6]);
    
    H2=plot(xg_cha,LOe2_dat); 
    set(gca,'xlim',[xlim1 xlim2]);
    set(gca,'ylim',[ylim1cha ylim2cha]);
    set(get(gca,'YLabel'),'String','E-Field(mc\omega / e)');
    set(get(gca,'YLabel'),'color','k');
    set(gca,'GridLineStyle','-');
    set(H2,'linewidth',1);
    set(H2,'color','k');
    grid on; hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    Himg=imagesc(flipdim((dset_qe_plt),2)',...
        'XData',[xg_cha xg_cha],...
        'YData',[yg_cha_plt yg_cha_plt]); 
    set(Himg,'AlphaData',0.8);
    shading('interp');
    colormap(cind);
    caxis(gca,[minch maxch]); 
    set(gca,'ydir','normal');
    grid on;
    freezeColors;
%         cbfreeze(colorbar('north'));
    str_tmp1=strcat('Frame:',num2str(frno),'     : ',...
        strcat('time=',num2str(time_qe),'(','1/{\omega_{pe}}',')'),...
        strcat('Dist=',num2str(time_qe*t_nor*vel_c*1e4),'(',...
        '\mu m',')'));
    title(str_tmp1);
    xlim([xlim1 xlim2]); ylim([ylim1cha ylim2cha]);
    xlabel('x1(c/\omega_{pe})'); 
%         ylabel('x2(c/\omega_{pe})');
    hold off;
    

   
    drawnow;

end

% close(writerObj);