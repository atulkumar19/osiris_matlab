clc;
clear all;
close all;

addpath('UtilFun','UtilFun/QTWriter');

ne=5.8e+19; % cm^-3 58.0e18, 6.12e+18
frm=0:70;

%============AVI-FILE-NAME==========================
% movObj = QTWriter('LaserElectricSpectrum.mov');
% movObj.FrameRate = 10;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 

minch=-2.5;
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
lamL=800*1e-9*1e2; % wavelength of laser in cm


omp_e=sqrt((4*pi*ne*qe^2)/me); x_nor=vel_c/omp_e; t_nor=1/omp_e;

omL_nor=2*pi/(lamL/x_nor);
omP_nor=2*pi/(
omMax=omL_nor+3;
omMin=omL_nor-3;
omTick=[omL_nor-3 omL_nor-2 omL_nor-1 omL_nor omL_nor+1 omL_nor+2 omL_nor+3];

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
    

    [xg,yg,dset_qe,x1lt,x2lt,time]=AshReadHDF5DenDat(fl_nm1);
    [~,~,dset_e2,~,~,~]=AshReadHDF5DenDat(fl_nm2);
    [Ngx,Ngy]=size(dset_qe);
    Nk=2^nextpow2(Ngx);
    
    d_mm=time*t_nor*vel_c*10;
        
    Lx=xg(end)-xg(1);
    dx=xg(2)-xg(1);
    minKx=2*pi/Lx;
    maxKx=2*pi/dx; % this corresponds to sampling frequency (samples per unit length)
    dkx=maxKx/Nk;
    kk=dkx*(0:(Nk/2))';
                
    x1min=x1lt(1); x1max=x1lt(2);
    x2min=x2lt(1); x2max=x2lt(2);
    mine2=min(min(dset_e2));  maxe2=max(max(dset_e2));
    
    Ngyby2=round(Ngy/2);
    dy=yg(2)-yg(1);
    
%     [xIndMax_e2,yIndMax_e2]=ind2sub(size(dset_e2),...
%         find(dset_e2==max(max(dset_e2,[],1),[],2),1));
%     Ngyby2e2=yIndMax_e2;
    Ngyby2e2=round(Ngy/2);
    LOe2=dset_e2(:,Ngyby2e2);
       
    FFT_LOe2=fft(LOe2,Nk);
    powFFT=abs(FFT_LOe2(1:Nk/2+1));
      
    subplot(2,2,1)
    plot(xg,LOe2);
    title('Lineout of LaserField'); grid on;
    ylabel('E-Field(mc\omega / e)');
    xlabel('x1(c/\omega_{pe})'); 
    
    subplot(2,2,2)
    Himg=imagesc(flipdim((dset_e2),2)',...
        'XData',[xg xg],...
        'YData',[yg yg]); 
%     set(Himg,'AlphaData',0.35);
    shading('interp');  
    newmap=b2r(mine2,maxe2);
    colormap(newmap);
    set(gca,'ydir','normal');
    str_tmp1=strcat('Frame:',num2str(frno),...
    '     : ',...
    strcat('time=',num2str(time),'(1/\omega_{pe})'),...
    '     : ',...
    strcat('Dist(mm) =',num2str(d_mm),'mm'));
    title(str_tmp1);
    xlabel('x1(c/\omega_{pe})'); 
    ylabel('x2(c/\omega_{pe})');
    grid on;
    colorbar('west')
%     freezeColors;
%         cbfreeze(colorbar('west'));
%     xlim([xlim1 xlim2]); ylim([ylim1cha ylim2cha]);
    
    subplot(2,2,[3 4]);
    
    plot(kk,powFFT);
    xlim([omMin omMax]); 
    set(gca,'xtick',omTick);
    set(gca, 'XTickLabel', num2str(get(gca, 'XTick')', '%.2f'));
    grid on;
    xlabel('\omega (\omega_p)'); ylabel('Power(arb units)');
   
    
    drawnow;
    

%     frame=getframe(gcf); 
%     writeMovie(movObj,frame);

end


% close(movObj);
