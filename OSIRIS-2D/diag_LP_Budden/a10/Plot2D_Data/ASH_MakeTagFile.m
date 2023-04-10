clc
clear all
close all

addpath('UtilFun');

%%%%%%%%%%% All constants go here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pth0='/home/bhavesh/Desktop/OSIRIS/OsirisRao/bin/PICSim';
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
%######################################################################


%%%%%%%%%%% Input from user go here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ne=5.8e19; % cm^-3
frno=177; %%%%%%% 75
NTE=50;
%######################################################################


omp_e=sqrt((4*pi*ne*qe^2)/me); x_nor=vel_c/omp_e; t_nor=1/omp_e;
minch=-5;
maxch=0;

str_frno=num2str(frno);
len_frno=length(str_frno);
str_num((len_str_num-len_frno+1):end)=str_frno;

fl_nm1=strcat(pth0,pth1,str_h1,str_num,str_ext); %charge density
fl_nm2=strcat(pth0,pth2,str_h2,str_num,str_ext); %raw data


[xg_cha,yg_cha,dset_qe,x1lt_qe,x2lt_qe,time_qe]=...
    AshReadHDF5DenDat(fl_nm1);
[~,x1_raw2,~,x2_raw2,~,tag_raw2,~,~]=...
    AshReadRaw(fl_nm2);

x1min_qe=x1lt_qe(1); x1max_qe=x1lt_qe(2);
x2min_qe=x2lt_qe(1); x2max_qe=x2lt_qe(2);
[Ngx_cha,Ngy_cha]=size(dset_qe);
Ngyby2=round(Ngy_cha/2);

imagesc(flipdim((dset_qe),2)',...
    'XData',[xg_cha xg_cha],...
    'YData',[yg_cha yg_cha]);      
shading('interp');
colormap(cmap('blue',256));
% caxis(gca,[minch maxch]); 
set(gca,'ydir','normal');
grid on;
freezeColors;
cbfreeze(colorbar('north'));
str_tmp1=strcat('frame=',num2str(frno),':',...
    'time=',num2str(time_qe));
str_tmp2=strcat('dist(cm)=',num2str(time_qe*t_nor*vel_c,'%10.3e'));
str_tmp3=strcat('xnor=', num2str(x_nor,'%10.3e'));
str_tmp4=strcat('tnor=', num2str(t_nor,'%10.3e'));
title(str_tmp1);
str_ti2=char(str_tmp2,str_tmp3,str_tmp4);
text(x1max_qe-61,0.465*x2max_qe,str_ti2);
ylim([15 25]); 
xlim([650 660]);
xlabel('x1'); ylabel('x2');

[xin,yin]=ginputc(NTE);

cnt1=0;loc_tag=zeros(NTE,1);
for ii=1:NTE
    dis=(x1_raw2-xin(ii)).^2+(x2_raw2-yin(ii)).^2;
    min_dis=min(dis);
    tmp=find(dis==min_dis,1);  
    if (any(tmp==loc_tag(1:(ii-1))))
        
    else
        cnt1=cnt1+1;
        loc_tag(cnt1)=tmp;
    end    
end

[fid,msg]=fopen('electrons.tags','w+');
fprintf(fid,'%10d\n',cnt1);
for ii=1:cnt1
    fprintf(fid,'%10d%c%10d\n',tag_raw2(loc_tag(ii),1),' ',...
        tag_raw2(loc_tag(ii),2));
end
fclose(fid);


