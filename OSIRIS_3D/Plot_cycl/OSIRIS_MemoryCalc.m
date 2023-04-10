clc
clear all
close all

%{
OpPerCycle=4; % op in one secons
cpuClockSpeed=2.9*1e9; % CyclePerSecond
OpPerPICcyclePerPart=300; % number of op in one PIC cycle
TimeOnePICcyclePerPart1=OpPerPICcyclePerPart/(cpuClockSpeed*OpPerCycle);
%}


BytePerScalar=8; % one scalar occupies 8 bytes of memory
TimeOnePICcyclePerPart=1e-6; % Time take by a cpu to complete one PIC cycle for single particle in seconds


%%%%%%%%%%%%% Input goes here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nx = 2100; % 2550.
ny = 360; % 652.
nz = 360; % 652.

ppc = 64; %  8  particle per cell

cpuNumber = 1580; % 2 cpuPernode * 790 nodes= 1580 cpu (1PetaFlop system)

Nscalars = 12; % number of scalars saved in 3D = (x,y,z,px,py,pz,Ex,Ey,Ez,Bx,By,Bz)

Nspecies = 2; %  3

tmax = 1.0e+03; %  Simulation time in normalised units
dt = 2.5e-02; % dt in normalised units

Nfiles=200;
rawfrac=10; % percentage of total number of particles for which info is stored
%%%%%%%%%%%%%%%Input ends here   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% calculations done here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tsteps=round(tmax/dt);

Noc = nx*ny*nz; %(* number of occupied cells *)
TNc = nx*ny*nz; %(* total number of cells *)
TotalParticles=Nspecies*Noc*ppc; % total number of particles 
rawParticles=(rawfrac/100)*TotalParticles;


MemOccPart = ...
    (((TotalParticles*Nscalars*BytePerScalar)/1024.0)/1024.0)/1024.0; %(* Gb *)

MemOccFields = ...
    ((TNc*9*BytePerScalar/1024.0)/1024.0)/1024.0;   %(* Gb *)

MemRawFiles=...
    Nfiles*((14*rawParticles*BytePerScalar/1024.0)/1024.0)/1024.0; %(* Gb *)
MemChargeDenFiles =...
    Nfiles*((Nspecies*TNc*BytePerScalar/1024.0)/1024.0)/1024.0;   %(* Gb *)
MemFieldFiles =...
    Nfiles*((6*TNc*BytePerScalar/1024.0)/1024.0)/1024.0;   %(* Gb *)

MemTot=MemRawFiles+MemChargeDenFiles+MemFieldFiles;

TimeOnePICcycle = TotalParticles*TimeOnePICcyclePerPart;    %(* s *)

NumberPICcycles = tmax/dt;

tsim = TimeOnePICcycle*NumberPICcycles/3600;
tsim =tsim/cpuNumber;




%%%%%%%%%%%%%%%%% we now print here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('tsim  =%f hours \n',tsim);
fprintf('tsteps  =%i  \n',tsteps);

fprintf('Total number of cells=%e \n', Noc);
fprintf('Total number of particles=%e \n', TotalParticles);
fprintf('Total number of raw particles=%e \n', rawParticles);
fprintf('Mem needed for particles =%f Gb \n', MemOccPart);
fprintf('Mem needed for fields =%f Gb \n', MemOccFields);
fprintf('total mem =%f Gb \n', MemOccFields+MemOccPart);

fprintf('mem needed for (fields and particles) / cpu  =%f Gb per cpu \n',...
    (MemOccFields + MemOccPart)/cpuNumber);

fprintf('Memory required to store %i raw files =%f Gb \n',Nfiles,...
    MemRawFiles);
fprintf('Memory required to store %i charge density files =%f Gb \n',Nfiles,...
    MemChargeDenFiles);
fprintf('Memory required to store %i field files =%f Gb \n',Nfiles,...
    MemFieldFiles);
fprintf('Total Memory required to store data =%f Gb \n',...
    MemRawFiles+MemChargeDenFiles+MemFieldFiles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%