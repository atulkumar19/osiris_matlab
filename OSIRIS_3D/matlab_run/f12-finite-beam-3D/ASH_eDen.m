clc;
clear all;
close all;

addpath('UtilFun');

ne=1.1e+22; % cm^-3 58.0e18, 6.12e+18
frm=150;

minch=-4;
maxch=0;

pth0='/Users/repulsion/Desktop/OSIRIS_3D/atul_osiris/f12-finite-beam';
pth1='/MS/DENSITY/Elec1/charge/';
pth2='/MS/DENSITY/Elec2/charge/';
pth3='/MS/DENSITY/Elec3/charge/';

str_h1='charge-Elec1-';
str_h2='charge-Elec2-';
str_h3='charge-Elec3-';


str_ext='.h5';
str_num='000000';
len_str_num=length(str_num);


qe=-4.8032e-10; %electron charge , statcoulomb
me=9.1094e-28; %electron mass , g
vel_c=2.9979e10; %velocity of light ,  cm/sec
omp_e=sqrt((4*pi*ne*qe^2)/me); x_nor=vel_c/omp_e; t_nor=1/omp_e;

% mycmap=colormap(cmap('blue',1024));
% mycmap=colormap(bone(2048));

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
    fl_nm2=strcat(pth0,pth2,str_h2,str_num,str_ext);
    fl_nm3=strcat(pth0,pth3,str_h3,str_num,str_ext);
    
    [xg,yg,zg,dset_qe1,x1lt,x2lt,x3lt,time]=AshReadHDF5DenDat3D(fl_nm1);
    [~,~,~,dset_qe2,~,~,~,~]=AshReadHDF5DenDat3D(fl_nm2);
    [~,~,~,dset_qe3,~,~,~,~]=AshReadHDF5DenDat3D(fl_nm3);
    
    
    [xmat,ymat,zmat]=meshgrid(xg,yg,zg);
    
    dset_qe=dset_qe1+dset_qe2+dset_qe3;
%     dset_qe=dset_qe1;
    xele=length(xg); xeleby2=round(xele/15);
    yele=length(yg);
    zele=length(zg);
    
%     Mat3DtoVTK(dset_qe,'qe250.vtk');
%{        
    qqe=dset_qe(1:xeleby2,1:5:end,1:5:end);

    maxqe=max(max(max(qqe)));
    minqe=min(min(min(qqe)));

    isovalues = linspace(minqe,maxqe,5);
    nisoval= numel(isovalues);
    hfig=figure('Position',[50 10 1800 1200]);
    set(gcf,'Renderer','OpenGL');
    p = zeros(nisoval,1);
    for i=1:nisoval
        p(i) = patch(isosurface(qqe,isovalues(i)));
        isonormals(qqe,p(i))
        set(p(i), 'CData',i);
    end
    set(p, 'CDataMapping','direct', 'FaceColor','flat', 'EdgeColor','none')

    mycmap=colormap(bone(nisoval));
    colormap(mycmap)
    caxis([0 nisoval])
    colorbar('YTick',(1:nisoval)-0.5, 'YTickLabel',num2str(isovalues(:)))
    box on; grid on; axis tight; 
%     daspect([1 1 1])
    view(3); 
    camproj perspective
    camlight; 
    lighting gouraud; 
    alpha(0.75);
    rotate3d on

    drawnow
%}    

    seg = dset_qe;
    segRange = [min(seg(:)) max(seg(:))];
    sizeSeg = size(seg); %will never need to change.
    sn = [100 300 0];
    visSeg = seg;
    visSeg = permute(visSeg, [3 1 2]);
%     visSeg = flipdim(visSeg, 2);
    xyview = seg;
    xzview = seg;
    yzview = seg;
    
%     xyview = flipdim(xyview, 1);
%     xyview = permute(xyview, [2 3 1]);
%     xzview = flipdim(xzview, 1);
    xzview = permute(xzview, [3 2 1]);
    xzview = flipdim(xzview, 2);
    yzview = permute(yzview, [1 3 2]);
    yzview = flipdim(yzview, 1);
    
    
    mycmap=colormap(jet(4098));
    subplot(2,2,1)
    imHandles1 = imagesc(squeeze(xyview(:,:,200)),...
        'XData',[yg yg],'YData',[xg xg],...
        segRange); 
    set(gca,'YDir','normal');
    set(gca,'FontSize',25,'FontWeight','Bold');
    shading('interp');
    colormap(mycmap); 
    ylabel('x1(c/\omega_{pe})'); xlabel('x2(c/\omega_{pe})');
    colorbar; 
    title(strcat('e-Density, XY-View, t=',num2str(time)));
    axis equal;
    axis image; 
    axis xy;
    
    
    subplot(2,2,2);
    imHandles2 = imagesc(squeeze(yzview(sn(1),:,:)),...
        'XData',[yg yg],'YData',[zg zg],...
        segRange);  
    set(gca,'YDir','normal');
    set(gca,'FontSize',25,'FontWeight','Bold');
    shading('interp');
    colormap(mycmap); 
    ylabel('x3(c/\omega_{pe})'); xlabel('x2(c/\omega_{pe})');
    colorbar; 
    title(strcat('e-Density, YZ-View, t=',num2str(time)));
    axis equal;
    axis image; 
    axis xy;
    
    subplot(2,2,3);
    imHandles3 = imagesc(squeeze(xzview(:,sn(2),:)),...
        'XData',[xg xg],'YData',[zg zg],...
        segRange);  
    set(gca,'YDir','normal');
    set(gca,'FontSize',25,'FontWeight','Bold');
    shading('interp');
    colormap(mycmap); 
    ylabel('x3(c/\omega_{pe})'); xlabel('x1(c/\omega_{pe})');
    colorbar; 
    title(strcat('e-Density, XZ-View, t=',num2str(time)));
    axis equal;
    axis image; 
    axis xy;
    
   subplot(2,2,4);
    hslc=slice(visSeg,sn(1),sn(2),sn(3),'cubic'); 
    set(gca,'FontSize',25,'FontWeight','Bold');
    shading('interp');
    colormap(mycmap); 
    ylabel('x2(c/20\omega_{pe})'); xlabel('x1(c/20 \omega_{pe})');zlabel('x3(c/20\omega_{pe})');
    colorbar; 
    title(strcat('e-Density, t=',num2str(time)));
    axis equal; 
    axis vis3d; 
   
    % set(hslc(1:3),'LineStyle','none');

    
%{    
    x1min=x1lt(1); x1max=x1lt(2);
    x2min=x2lt(1); x2max=x2lt(2);
    [Ngx,Ngy]=size(dset_qe);
    Ngyby2=round(Ngy/2);
    dy=yg(2)-yg(1);
    
        
    xlim1=x1min; xlim2=x1max;
    ylim1=x2min; ylim2=x2max;

    Himg=imagesc(flipdim((dset_qe_plt),2)',...
        'XData',[xg xg],...
        'YData',[yg yg]); 
%     set(Himg,'AlphaData',0.8);
    shading('interp');
    colormap(cmap('orange',1024));
    caxis(gca,[minch maxch]); 
    set(gca,'ydir','normal');
    grid on;
    freezeColors;
%         cbfreeze(colorbar('north'));
    str_tmp1=strcat('Frame:',num2str(frno),...
        '     : ',...
        strcat('time=',num2str(time),'(1/\omega_{pe})'));
    title(str_tmp1);
    xlim([xlim1 xlim2]); 
%     ylim([ylim1cha ylim2cha]);
    xlabel('x1(c/\omega_{pe})'); 
%         ylabel('x2(c/\omega_{pe})');


    drawnow;
    
%     gif_add_frame(gcf,'eDenE1E2P1X1_WithoutDLA_28May2015.gif');
%}
end