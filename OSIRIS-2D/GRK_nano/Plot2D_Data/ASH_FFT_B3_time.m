clc;
clear all;
close all;

addpath('UtilFun','UtilFun/QTWriter');

ne=1.2e+22; % cm^-3 58.0e18, 6.12e+18
frm=0:500;

%============AVI-FILE-NAME==========================
% movObj = QTWriter('LaserElectricSpectrum.mov');
% movObj.FrameRate = 10;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 

minch=-2.5;
maxch=0;

pth0='/Users/repulsion/Desktop/OSIRIS/OSIRIS_2D/atul/counter_finite_0.9c';
pth1='/MS/DENSITY/Elec1/charge/';
pth2='/MS/DENSITY/Elec2/charge/';
pth3='/MS/DENSITY/Elec3/charge/';
pth4='/MS/FLD/b3/';

str_h1='charge-Elec1-';
str_h2='charge-Elec2-';
str_h3='charge-Elec3-';
% str_h4='charge-Elec4-';
% str_h5='charge-Elec5-';
str_h4='b3-';
str_ext='.h5';
str_num='000000';
len_str_num=length(str_num);
qe=-4.8032e-10; %electron charge , statcoulomb
me=9.1094e-28; %electron mass , g
vel_c=2.9979e10; %velocity of light ,  cm/sec
lamL=800*1e-9*1e2; % wavelength of laser in cm


omp_e=sqrt((4*pi*ne*qe^2)/me); x_nor=vel_c/omp_e; t_nor=1/omp_e;

omL_nor=2*pi/(lamL/x_nor);
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
    fl_nm2=strcat(pth0,pth2,str_h2,str_num,str_ext);
    fl_nm3=strcat(pth0,pth3,str_h3,str_num,str_ext);
    fl_nm4=strcat(pth0,pth4,str_h4,str_num,str_ext);
    

    [xg,yg,dset_qe1,x1lt,x2lt,time]=AshReadHDF5DenDat(fl_nm1);
    [xg,yg,dset_qe2,x1lt,x2lt,time]=AshReadHDF5DenDat(fl_nm2);
    [xg,yg,dset_qe3,x1lt,x2lt,time]=AshReadHDF5DenDat(fl_nm3);
    [~,~,dset_b3,~,~,~]=AshReadHDF5DenDat(fl_nm4);
    [Ngx,Ngy]=size(dset_qe1);
    
    d_mm=time*t_nor*vel_c*10;
    
    Lx=xg(end)-xg(1);
    dx=xg(2)-xg(1)
    tmax=max(time);
    tmin=dx./(sqrt(2));
    minOmega=1/tmax;
    maxOmega=1/tmin; % this corresponds to sampling frequency (samples per unit length)
    dt=tmax/Ngx;
    kk=dkx*(0:(Ngx-1))';
    kkplot=kk(kk<(maxKx/2));
    omplot=2*pi*kkplot;
            
    x1min=x1lt(1); x1max=x1lt(2);
    x2min=x2lt(1); x2max=x2lt(2);
    minb3=min(min(dset_b3));  maxb3=max(max(dset_b3));
    
    Ngyby2=round(Ngy/2);
    dy=yg(2)-yg(1);
    
    Ngyby2b3=round(Ngy/2);
    LOb3=dset_b3(:,Ngyby2b3);
       
    FFT_LOb3=fft(LOb3,Ngx);
    powFFT=abs(FFT_LOb3(1:round(Ngx/2)));
      
   
   
    
    plot(omplot,powFFT);
    xlim([omMin omMax]); 
    set(gca,'xtick',omTick);
    set(gca, 'XTickLabel', num2str(get(gca, 'XTick')', '%.2f'));
    grid on;
    xlabel('\omega (\omega_p)'); ylabel('Power(arb units)');
    set(gca,'YDir','normal');
    set(gca,'FontSize',16,'FontWeight','Bold');
   
    title(strcat('FFT,t= ',num2str(time)))
   
    
    drawnow;
    

%     frame=getframe(gcf); 
%     writeMovie(movObj,frame);

end


% close(movObj);
