clc;
clear all;
close all;

addpath('UtilFun');

%%%%%%%%%%% All constants go here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pth0='/home/bhavesh/Desktop/OSIRIS/OsirisRao/bin_data/PICSim_200MicronRamp';
pth1='/MS/DENSITY/He_electrons/charge/';
pth2='/MS/RAW/He_electrons/';

str_h1='charge-He_electrons-';
str_h2='RAW-He_electrons-';

str_ext='.h5';
str_num='000000';
len_str_num=length(str_num);


qe=-4.8032e-10; %electron charge , statcoulomb
me=9.1094e-28; %electron mass , g
vel_c=2.9979e10; %velocity of light ,  cm/sec
RestMassEnergy_ele=8.19e-14; % Joules or 511kev or 0.511 Mev
cmapval=cmap('orange',1024);
%######################################################################

%%%%%%%%%%% Input from user go here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ne=5.8e19; % cm^-3
frm=0:4:200; %%%%%%% 177
% Get the width and height of the figure
%============AVI-FILE-NAME==========================
%  aviobj=avifile('EnergySpectrum','fps',10); 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
movObj = QTWriter('EnergySpectrum.mov');
movObj.FrameRate = 10;

%######################################################################



omp_e=sqrt((4*pi*ne*qe^2)/me); x_nor=vel_c/omp_e; t_nor=1/omp_e;
minch=-3;
maxch=0;

scrsz = get(0,'ScreenSize');
figure('Position',[50 10 1500 900]);
for ii=1:length(frm)    
    frno=frm(ii);
        
    str_frno=num2str(frno);
    len_frno=length(str_frno);
    str_num((len_str_num-len_frno+1):end)=str_frno;

    fl_nm1=strcat(pth0,pth1,str_h1,str_num,str_ext); %charge density
    fl_nm2=strcat(pth0,pth2,str_h2,str_num,str_ext); %raw data

    [xg,yg,dset_qe,x1lt,x2lt,time]=...
        AshReadHDF5DenDat(fl_nm1);
    [ene_raw,x1_raw,~,x2_raw,~,tag_raw,~,~]=...
        AshReadRaw(fl_nm2);
    
    dx=xg(2)-xg(1); dy=yg(2)-yg(1);
    x1min=x1lt(1); x1max=x1lt(2);
    x2min=x2lt(1); x2max=x2lt(2);
    xlim1=x1min; xlim2=x1max;
    ylim1=x2min; ylim2=x2max;
    maxene=110; minene=10;
    
    d_mm=time*t_nor*vel_c*10;
    
%     ind_ylim1=round(ylim1/dy)+1; 
%     ind_ylim2=round(ylim2/dy)-1;
%     dset_qe_plt=dset_qe(:,ind_ylim1:ind_ylim2);
%     yg_cha_plt=yg(ind_ylim1:ind_ylim2);
    dset_qe_plt=dset_qe;
    yg_cha_plt=yg;
    LT_a=(maxene-minene)/( ylim2-ylim1);
    LT_b=minene-LT_a*ylim1;
    yg_cha_plt=LT_a*yg_cha_plt+LT_b;
    ylim1cha=LT_a*ylim1+LT_b;
    ylim2cha=LT_a*ylim2+LT_b;
    
    subtightplot(1,2,1)
    [n,nout]=hist(ene_raw,128);
    plot(nout,n,'linewidth', 2)
    xlim([10 110])
    xlabel('Ene(mc^2)'); ylabel('counts(arb units)');
    title(strcat('frame=',num2str(frno),'  Hist Ene'));
    grid on;
%}    
    subtightplot(1,2,2)   
    Himg=imagesc(flipdim((dset_qe),2)',...
        'XData',[xg xg],...
        'YData',[yg_cha_plt yg_cha_plt]); 
%     set(Himg,'AlphaData',0.5);
    shading('interp');
    colormap(cmapval);
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
    xlim([xlim1 xlim2]); 
    ylim([ylim1cha ylim2cha]);
    xlabel('x1(c/\omega_{pe})'); 
%         ylabel('x2(c/\omega_{pe})');
    hold on;
    
    ene_th_loc=ene_raw>10;
    if ~any(ene_th_loc)
        H1=plot(x1_raw(1),ene_raw(1)); hold on;
    else
        H1=plot(x1_raw(ene_th_loc),ene_raw(ene_th_loc)); 
    end
    set(gca,'xlim',[xlim1 xlim2]);
    set(gca,'YAxisLocation','right');
    ylabel(gca,'Energy(mc^2)');
    set(H1,'marker','o','linestyle','none',...
        'MarkerFaceColor','g','MarkerSize',4);
    hold off;
    
    drawnow;
    
%     frame = getframe ( gcf ); 
%     aviobj = addframe ( aviobj, frame );

    frame=getframe(gcf); 
    writeMovie(movObj,frame);

end

% aviobj = close ( aviobj );

close(movObj);
