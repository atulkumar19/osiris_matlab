clc;
clear all;
close all;

%%%%The physical constants go here. CGS only!!!!!%%%%%%%%%%%%%%%%%%%%%%
amu=1.660538782*10^-27; %1 amu = 1.660538782*10^-27 kg
kb=1.3807e-16; %Boltzmann Constant,  erg/deg(K)
qe=-4.8032e-10; %electron charge , statcoulomb
me=9.1094e-28; %electron mass , g
mp=1.6726e-24; %proton mass , g
vel_c=2.9979e10; %velocity of light ,  cm/sec
h=6.6261e-27;%Planck constant , erg-sec
dim=2;
PushTime=1e-6; % time taken to push(move) single particle in seconds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%All input parameters go here%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wav_L=800*1e-9*1e2; % wavelength of laser 800nm
tau_L=45*1e-15; % pulse length of laser 30fs
a0=1.0;
w0=10.0*1e-4; % spotsize of the focussed Gaussian Beam
% p=3*1e12; % power in Watts (1Terawatt=1e12 watts)

ne=5.8e+19; % cm^-3
zch_i=2; % charge state of ion 
mi=4.00; % mass of ion in amu

node_number=[5;4]; %# of processors in x-dir and y-dir
NumParPerCell_x=4;
NumParPerCell_y=4;
Ncx=2100;
Ncy=300;
MovWin_x=60.0;
MovWin_y=40.0;
TotSimTim=1710;

NumFiles=855; % number of data files I want for each physical quantity

Pre_rmp=0.0;
rmpLx=800*1e-6*1e2;
LenPlaCol=200*1e-6*1e2; %length(cm) of preformed plasma column
rmpRx=200*1e-6*1e2;
Pos_rmp=10*1e-6*1e2;
rmpTy=2*1e-6*1e2;
rmpBy=2*1e-6*1e2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tau_L=sqrt(2)*tau_L;
lon_rise=tau_L; lon_fall=tau_L;
per_focus=-(Pre_rmp+0.5*rmpLx);

%%%%Laser Parameters calculated here. CGS only!!!!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%
om_L=2*pi*vel_c/wav_L; % frequency of laser
T_L=(wav_L/vel_c); % time period of laser wavelength
k_L=2*pi/wav_L; %  Laser wavenumber
% I0=p/((pi/2)*w0^2); % Intensity of the focussed Gaussian Beam (watt/cm^2)
% a0_sq=7.3e-19*((wav_L*1e4)^2)*I0;
% a0=sqrt(a0_sq);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%Plasma parameters calculated here%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ni=ne;
mi=mi*amu*1e3; % mass of ion in grams
qi=zch_i*abs(qe);
rqm_e=me/qe;
rqm_i=mi/qi;
omp_e=sqrt((4*pi*ne*qe^2)/me); fp_e=omp_e/(2*pi);
omp_i=sqrt((4*pi*ni*qi^2)/mi); fp_i=omp_i/(2*pi);
kp_e=omp_e/vel_c; % plasma wave number;
wav_P=2*pi*vel_c/omp_e; %plasma wavelength

%%%%Normalization constants calculated here %%%%%%%%%%%%%%%%%%%%%%%%%%%
%     x_nor=1; t_nor=1; om_nor=1; rqm_nor=1;
    x_nor=vel_c/omp_e; t_nor=1/omp_e; om_nor=omp_e; rqm_nor=abs(rqm_e);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%Plasma parameters calculated here. CGS only!!!!!!!%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Some simulation parameters go here !!!!!!!%%%%%%%%%%%%%%%%%%%%%%%%%
TotNumParPerCell=NumParPerCell_x*NumParPerCell_y;
Tot_proc=node_number(1)*node_number(2); % total # of processors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



omPebyomL=omp_e/om_L;


%%%%Laser Power corresponding to given a0%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
powc=17*(om_L/omp_e)^2; % critical power is in GigaWatt
% p=powc*((0.5*a0)^3); % power to have given a0
% a0=2*(p/powc)^(1/3);
if (~exist('a0','var'))
    powc=powc*1e9; %critical power in watts
    a0=2*(p/powc)^(1/3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%Spot size of focused Laser pulse for stable self guiding %%%%%%%%%%
if (~exist('w0','var'))
    w0=2*sqrt(a0)/kp_e;
%     w0=1.1236*w0;
end
R=w0; %radius of blow out region
if (~exist('tau_L','var'))
    tau_L=0.5*R/vel_c;
end
lpl=vel_c*tau_L; % laser pulse length
fprintf('PlasmaFrequecyByLaserFrequency)=%f\n',omPebyomL);
fprintf('Laser Spot Size(w0(micron))=%f\n',w0*1e4);
fprintf('Laser Pulse duration(tau(fs))=%f\n',tau_L/1e-15);
fprintf('Laser Pulse Length(lpl)=%f\n',lpl);
fprintf('Radius of the bubble(R)=%f\n',R);
fprintf('For spherical bubble we need LaserPulseLength<RadiusOfBubble\n');
fprintf('Wakefield generation is efficient iff LaserPulseLength=PlasmaWaveLengthBy2\n');
fprintf('Plasma wavelength by 2 (wav_Pby2)=%f\n',wav_P/2);
fprintf('Ideal Laser Pulse Duration (tau(fs)=%f\n',...
    (pi/omp_e)/1e-15);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%Pump Depletion Length%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lpd=((om_L/omp_e)^2)*vel_c*tau_L;
display(Lpd)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%Dephasing length%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ld=(2/3)*((om_L/omp_e)^2)*R;
display(Ld)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%Energy gain%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delE=(2/3)*me*(vel_c^2)*((om_L/omp_e)^2)*a0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%Discretization of space and time domain%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%x-size of moving window%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (~exist('MovWin_x','var'))
    OneBucket=2*w0;
    MovWin_x1=2.6*OneBucket;
    MovWin_x2=4*lpl;
    % MovWin_x1=3.5*OneBucket;
    % MovWin_x2=11.32*lpl; % size of moving window in x direction. DEFINED ARB..
    MovWin_x=max(MovWin_x1,MovWin_x2);
else
     MovWin_x= MovWin_x*x_nor; % user gives MovWin_x in simulation units
end

if (~exist('MovWin_y','var'))
    %%%%%%y-size of moving window%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MovWin_y=3*OneBucket; % size of moving window in y direction. DEFINED ARB..
    % MovWin_y=2.0*OneBucket; % size of moving window in y direction. DEFINED ARB..
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    MovWin_y= MovWin_y*x_nor-(2*rmpTy+2*rmpBy); % user gives MovWin_y in simulation units
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zr=pi*w0^2/wav_L; %rayleigh length for circularly polarized light
% LenPlaCol=5*zr;
if (~exist('LenPlaCol','var'))
    LenPlaCol=max([5*zr;Lpd;Ld]); % length of plasma column through Laser pluse propagates (cm)
end
tim_Lpd_cst=Lpd/vel_c;
tim_Ld_cst=Ld/vel_c;
tim_resol_cst=abs(tim_Lpd_cst-tim_Ld_cst);

%%%%remember laser propagates in x-direction%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (exist('Ncx','var')&&exist('Ncy','var'))
    dx=MovWin_x/Ncx; dy=MovWin_y/Ncy;
else    
    dx1=wav_L/40; %ARB criterion: 20 to 40 points in laser wavelength.
    dx2=0.2/k_L; % Laser prpagates in x-dir
    dx=min(dx1,dx2);
    dy1=w0/40; %ARB criterion: 40 points in laser spot size.
    dy2=0.116/kp_e; % kp_e is plasma wave  number
    dy=min(dy1,dy2);
    Ncx=ceil(MovWin_x/dx); dx=MovWin_x/Ncx; MovWin_x=dx*Ncx;
    Ncy=ceil(MovWin_y/dy); dy=MovWin_y/Ncy; MovWin_y=dy*Ncy;

end

ncell_rmpLx=ceil(rmpLx/dx); rmpLx=ncell_rmpLx*dx;
ncell_rmpRx=ceil(rmpRx/dx); rmpRx=ncell_rmpRx*dx;
ncell_rmpTy=ceil(rmpTy/dy); rmpTy=ncell_rmpTy*dy;
ncell_rmpBy=ceil(rmpBy/dy); rmpBy=ncell_rmpBy*dy;

% LenPlaCol=LenPlaCol+0.1*LenPlaCol;
TotLenY=MovWin_y+2*rmpTy+2*rmpBy; % length of plasma column through Laser pluse propagates

TotLenX=Pre_rmp+rmpLx+LenPlaCol+rmpRx+Pos_rmp;
if (~exist('TotSimTim','var'))
    TotSimTim=TotLenX/vel_c; % total simulation time
else
    TotSimTim=TotSimTim*t_nor; %user gives time in normalised units
end

TotNumPar=TotNumParPerCell*Ncx*Ncy;
num_par_max= TotNumPar/Tot_proc; %number of particles per node
num_par_max=num_par_max+(5/100)*num_par_max; %number of particles per node
num_par_max=2*num_par_max; % helium ionizes to give 2 electrons

dt0=(0.99)/(vel_c*sqrt((1/dx)^2+(1/dy)^2));%dt calculated from courant condition
dt1=0.2/omp_e;
dt=min(dt0,dt1);
TotalTimeSteps=round(TotSimTim/dt);

ComSimTim_hr=(TotalTimeSteps*(TotNumPar*PushTime)/Tot_proc)/3600;

ndump=round(TotalTimeSteps/NumFiles);
del_t_img=TotSimTim/NumFiles;
NoFiles_Bet_Lpd_Ld=tim_resol_cst/del_t_img;

%ratio_dt must be greater than 30
ratio_dt=T_L/dt; % ARB criterion: 30-40 pts in time period
fprintf('ratio_dt=%f\n\n\n\n',ratio_dt);
fprintf('Num of files betn Lpd and Ld=%f\n\n\n\n',NoFiles_Bet_Lpd_Ld);

if (ratio_dt<30)
    fprintf('warning: consider readjusting dx and dy. poor resolution\n');
    fprintf('warning: No calculations done further\n');   
else
    cgs_par=[ne;omp_e;MovWin_x;Ncx;dx;MovWin_y;Ncy;dy;...
        Lpd;Ld;LenPlaCol;...
        dt;TotalTimeSteps;TotSimTim;...
        a0;wav_L;om_L;w0;tau_L;...
        om_L/omp_e;];
    dt=dt/t_nor;
    TotalTime=TotSimTim/t_nor;
    w0=w0/x_nor; lpl=lpl/x_nor;
    tau_L=tau_L/t_nor; lon_rise=lon_rise/t_nor; lon_fall=lon_fall/t_nor;
    om_L=om_L/om_nor;
     
    MovWin_x=MovWin_x/x_nor; dx=dx/x_nor; Ncx=round(MovWin_x/dx);
    MovWin_y=MovWin_y/x_nor; dy=dy/x_nor; Ncy=round(MovWin_y/dy);
    TotLenX=TotLenX/x_nor; TotLenY=TotLenY/x_nor;
    Pre_rmp=Pre_rmp/x_nor; rmpLx=rmpLx/x_nor; LenPlaCol=LenPlaCol/x_nor;
    rmpRx=rmpRx/x_nor; Pos_rmp=Pos_rmp/x_nor;
    rmpTy=rmpTy/x_nor; rmpBy=rmpBy/x_nor;
    per_focus=per_focus/x_nor;
   
    rqm_i=rqm_i/rqm_nor; rqm_e=rqm_e/rqm_nor;
    
    nor_par=[ne;omp_e;MovWin_x;Ncx;dx;MovWin_y;Ncy;dy;...
        Lpd/x_nor;Ld/x_nor;TotLenX;...
        dt;TotalTimeSteps;TotalTime;...
        a0;wav_L/x_nor;om_L;w0;tau_L;...
        om_L/omp_e;x_nor;t_nor;om_nor];
    
    %%%%%profile calculations%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PrX1=zeros(6,1); 
    PrX1(1)=0.0; PrX1(2)=PrX1(1)+MovWin_x+Pre_rmp; 
    PrX1(3)=PrX1(2)+rmpLx; PrX1(4)=PrX1(3)+LenPlaCol; 
    PrX1(5)=PrX1(4)+rmpRx; PrX1(6)=PrX1(5)+Pos_rmp;
    PrF1=zeros(6,1);
    PrF1(1)=0.0; PrF1(2)=0.0; PrF1(3)=1.0;
    PrF1(4)=1.0; PrF1(5)=0.0; PrF1(6)=0.0;
    
    PrX2=zeros(6,1); 
    PrX2(1)=0.0; PrX2(2)=PrX2(1)+rmpBy; 
    PrX2(3)=PrX2(2)+rmpBy; PrX2(4)=PrX2(3)+MovWin_y; 
    PrX2(5)=PrX2(4)+rmpTy; PrX2(6)=PrX2(5)+rmpTy;
    PrF2=zeros(6,1);
    PrF2(1)=0.0; PrF2(2)=0.0; PrF2(3)=1.0;
    PrF2(4)=1.0; PrF2(5)=0.0; PrF2(6)=0.0;
    
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %now calculate number of nodes required in each direction
    % ideally there should be 40cells per processor.
    % this means that---> Ncx/40~   and Ncy/40~
    %     Nnodex=round(Ncx/40); Nnodey=round(Ncy/40);
    %%%%%%%%%%%os-stdin starts here%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('x_nor=%10.4e, t_nor=%10.4e, omega_nor=%10.4e\n',...
        x_nor,t_nor,om_nor);
    fprintf('MovWin_x(cm)=%10.4e, MovWin_y(cm)=%10.4e\n',...
        MovWin_x*x_nor,MovWin_y*x_nor);
    fprintf('Length of Plasma Column(cm)=%10.4e\n',...
        LenPlaCol*x_nor);
    fprintf('dx(cm)=%10.4e, dy(cm)=%10.4e, dt(s)=%10.4e\n\n\n',...
        dx*x_nor,dy*x_nor,dt*t_nor);
    
    fprintf('!!!!!!!!!!!!!!!!!!!!!!os-stdin startss here!!!!!!!!!!!!!!!!!!!!!! \n');
    fprintf('plasma density(n0)(cm^{-3})=%10.5e \n',ne);
    fprintf('node_number(1:2)=%f,%f \n',node_number(1),node_number(2));
    fprintf('nx_p(1:2)=%f, %f \n',Ncx,(Ncy+ncell_rmpTy+ncell_rmpBy));
    fprintf('dt=%10.5e\n',dt);
    fprintf('ndump=%10.5e\n',ndump);
    fprintf('xmax(1:2)=%10.5e ,%10.5e \n',MovWin_x,TotLenY);
    fprintf('tmax=%10.5e\n',TotalTime);
    fprintf('num_par_max=%f\n',num_par_max);
    fprintf('rqm_e=%f\n',rqm_e);
    fprintf('num_par_x(1:2)=%f,%f\n',NumParPerCell_x,NumParPerCell_y);
    fprintf('num_x=%f\n',6);
    fprintf('x(1:6,1)=%8.3f,%8.3f,%8.3f,%8.3f,%8.3f,%8.3f\n',...
        PrX1(1),PrX1(2),PrX1(3),PrX1(4),PrX1(5),PrX1(6));
    fprintf('fx(1:6,1)=%8.3f,%8.3f,%8.3f,%8.3f,%8.3f,%8.3f\n',...
        PrF1(1),PrF1(2),PrF1(3),PrF1(4),PrF1(5),PrF1(6));
    fprintf('x(1:6,2)=%8.3f,%8.3f,%8.3f,%8.3f,%8.3f,%8.3f\n',...
        PrX2(1),PrX2(2),PrX2(3),PrX2(4),PrX2(5),PrX2(6));
    fprintf('fx(1:6,2)=%8.3f,%8.3f,%8.3f,%8.3f,%8.3f,%8.3f\n',...
        PrF2(1),PrF2(2),PrF2(3),PrF2(4),PrF2(5),PrF2(6));
    fprintf('ps_xmax(1:2)=%10.5e ,%10.5e \n',MovWin_x,PrX2(6));
    
    fprintf('num_par_max=%f\n',num_par_max);
    fprintf('rqm_i=%f\n',rqm_i);
    fprintf('num_par_x(1:2)=%f,%f\n',NumParPerCell_x,NumParPerCell_y);
    fprintf('num_x=%f\n',6);
    fprintf('x(1:6,1)=%8.3f,%8.3f,%8.3f,%8.3f,%8.3f,%8.3f\n',...
        PrX1(1),PrX1(2),PrX1(3),PrX1(4),PrX1(5),PrX1(6));
    fprintf('fx(1:6,1)=%8.3f,%8.3f,%8.3f,%8.3f,%8.3f,%8.3f\n',...
        PrF1(1),PrF1(2),PrF1(3),PrF1(4),PrF1(5),PrF1(6));
    fprintf('x(1:6,2)=%8.3f,%8.3f,%8.3f,%8.3f,%8.3f,%8.3f\n',...
        PrX2(1),PrX2(2),PrX2(3),PrX2(4),PrX2(5),PrX2(6));
    fprintf('fx(1:6,2)=%8.3f,%8.3f,%8.3f,%8.3f,%8.3f,%8.3f\n',...
        PrF2(1),PrF2(2),PrF2(3),PrF2(4),PrF2(5),PrF2(6));
    fprintf('ps_xmax(1:2)=%10.5e ,%10.5e \n',MovWin_x,PrX2(6));
    
    fprintf('\n\nzpulse parameters\n');
    fprintf('a0=%10.5e \n',a0);
    fprintf('Omega0=%10.5e\n',om_L);
    fprintf('lon_rise=%10.5e\tlon_fall=%10.5e\n',lon_rise,lon_fall);
    fprintf('spot size of focussed Laser(per_w0)=%10.5e\n',w0);
    fprintf('per_focus=%10.5e\n',per_focus);
    fprintf('zpulse parameters ends here\n');
    fprintf('!!!!!!!!!!!!!!!!!!!!!!os-stdin ends here!!!!!!!!!!!!!!!!!!!!!! \n');
    %%%%%%%%%%%os-stdin ends here%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fprintf('\n\n\nsome other relevant parameters\n');
    fprintf('TotalSimulationTimeInHours=%10.5e\n',ComSimTim_hr);
    fprintf('TotalTimeSteps=%10.5e\n',TotalTimeSteps);
    fprintf('tau_L=%10.5e\n',tau_L);
    fprintf('TotalNumberOfParticles=%f\n',TotNumPar);    
    fprintf('Nnodex=%f \t Nnodey=%f \n',node_number(1),node_number(2));
    fprintf('dx=%10.5e \t dy=%10.5e \n',dx,dy);
    fprintf('TotalLengthofSimBox=%10.5e \n\n\n',TotLenX);
   
end
