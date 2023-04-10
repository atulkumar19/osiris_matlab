clc;
clear all;
close all;

addpath('UtilFun');

ne=1.1e+22; % cm^-3 58.0e18, 6.12e+18

minch=-4;
maxch=0;

% mycmap=colormap(cmap('Jet',1024));

qe=-4.8032e-10; %electron charge , statcoulomb
me=9.1094e-28; %electron mass , g
vel_c=2.9979e10; %velocity of light ,  cm/sec
omp_e=sqrt((4*pi*ne*qe^2)/me); x_nor=vel_c/omp_e; t_nor=1/omp_e;

% len_fr=length(frm);
% scrsz = get(0,'ScreenSize');
% hfig=figure('Position',[50 10 1200 650]);
% set(gcf,'Renderer','OpenGL');

% mycmap=colormap(cmap('blue',1024));


fl_nm1='b3-000150.h5'; %charge density
[xg,yg,zg,dsetb3,x1lt,x2lt,x3lt,time]=AshReadHDF5DenDat3D(fl_nm1);


% Mat3DtoVTK(dsetb3,'./B250.vtk');


%{
BB3=dsetb3(1:5:end,1:5:end,1:5:end);



maxB3=max(max(max(BB3)));
minB3=min(min(min(BB3)));

isovalues = linspace(minB3,maxB3,10);
nisoval= numel(isovalues);
hfig=figure('Position',[50 10 1800 1200]);
set(gcf,'Renderer','OpenGL');
p = zeros(nisoval,1);
for i=1:nisoval
    p(i) = patch(isosurface(BB3,isovalues(i)));
    isonormals(BB3,p(i))
    set(p(i), 'CData',i);
end
set(p, 'CDataMapping','direct', 'FaceColor','flat', 'EdgeColor','none')

mycmap=colormap(jet(nisoval));
colormap(mycmap)
caxis([0 nisoval])
colorbar('YTick',(1:nisoval)-0.5, 'YTickLabel',num2str(isovalues(:)))
box on; grid on; axis tight; %daspect([1 1 1])
view(3); camproj perspective
camlight; lighting gouraud; 
alpha(0.65);
rotate3d on

%}

%{
%# volumetric data, and iso-levels we want to visualize
[x,y,z,v] = flow(25);
isovalues = linspace(-2.5,1.5,6);
num = numel(isovalues);

%# plot isosurfaces at each level, using direct color mapping
figure('Renderer','opengl')
p = zeros(num,1);
for i=1:num
    p(i) = patch( isosurface(x,y,z,v,isovalues(i)) );
    isonormals(x,y,z,v,p(i))
    set(p(i), 'CData',i);
end
set(p, 'CDataMapping','direct', 'FaceColor','flat', 'EdgeColor','none')

%# define the colormap
clr = hsv(num);
colormap(clr)

%# legend of the isolevels
%#legend(p, num2str(isovalues(:)), ...
%#  'Location','North', 'Orientation','horizontal')

%# fix the colorbar to show iso-levels and their corresponding color
caxis([0 num])
colorbar('YTick',(1:num)-0.5, 'YTickLabel',num2str(isovalues(:)))

%# tweak the plot and view
box on; grid on; axis tight; daspect([1 1 1])
view(3); camproj perspective
camlight; lighting gouraud; alpha(0.75);
rotate3d on
%}


%{
szB3=size(dsetb3);
xlen=szB3(1); ylen=szB3(2); zlen=szB3(3);
% xlen=100; ylen=szB3(2); zlen=szB3(3);
tlen=xlen*ylen*zlen;
b3vec=zeros(tlen,1);
xvec=zeros(tlen,1);
yvec=zeros(tlen,1);
zvec=zeros(tlen,1);

for dd=1:tlen
    [ii,jj,kk] = fast_ind2sub(szB3,dd);
    b3vec(dd)=dsetb3(ii,jj,kk);
    xvec(dd)=xg(ii);
    yvec(dd)=yg(jj);
    zvec(dd)=zg(kk);
end

fscatter3(xvec,yvec,zvec,b3vec);
%}