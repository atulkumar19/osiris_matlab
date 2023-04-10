
clc;
clear all;
close all;

addpath('UtilFun');

ne=1.1e+22; % cm^-3 58.0e18, 6.12e+18
frm=50;

minch=-4;
maxch=0;

pth0='/Users/akash/Desktop/OSIRIS_3D/atul_osiris';
pth1='/MS/FLD/b3/';

str_h1='b3-';

str_ext='.h5';
str_num='000000';
len_str_num=length(str_num);

% mycmap=colormap(cmap('Jet',1024));

qe=-4.8032e-10; %electron charge , statcoulomb
me=9.1094e-28; %electron mass , g
vel_c=2.9979e10; %velocity of light ,  cm/sec
omp_e=sqrt((4*pi*ne*qe^2)/me); x_nor=vel_c/omp_e; t_nor=1/omp_e;

len_fr=length(frm);
% scrsz = get(0,'ScreenSize');
% hfig=figure('Position',[50 10 1200 650]);
% set(gcf,'Renderer','OpenGL');
for ii=1:length(frm)    
    
    frno=frm(ii);
    str_frno=num2str(frno);
    len_frno=length(str_frno);
    str_num((len_str_num-len_frno+1):end)=str_frno;

    fl_nm1=strcat(pth0,pth1,str_h1,str_num,str_ext); %charge density
    
    
    [xg,yg,zg,dset_b3,x1lt,x2lt,x3lt,time]=AshReadHDF5DenDat3D(fl_nm1);
    
    BB3=dset_b3(1:5:end,1:5:end,1:5:end);



    maxB3=max(max(max(BB3)));
    minB3=min(min(min(BB3)));

    isovalues = linspace(minB3,maxB3,100);
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

    mycmap=colormap(jet);
    colormap(mycmap)
%     caxis([0 nisoval])
    colorbar;
    box on; grid on; axis tight; %daspect([1 1 1])
    view(3); camproj perspective
    camlight; lighting gouraud; 
    alpha(0.65);
    rotate3d on
end
