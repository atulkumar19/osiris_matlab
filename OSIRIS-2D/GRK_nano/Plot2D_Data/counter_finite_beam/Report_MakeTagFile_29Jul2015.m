clc
clear all
close all

addpath('UtilFun');

%%%%%%%%%%% All constants go here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pth0='/Volumes/Macintosh_HD_3/osiris_run/counter_finite_tag';
pth1='/MS/DENSITY/Elec1/charge/';
pth2='/MS/DENSITY/Elec2/charge/';
pth3='/MS/DENSITY/Elec3/charge/';
pth4='/MS/RAW/Elec1/';
pth5='/MS/RAW/Elec2/';
pth6='/MS/RAW/Elec3/';
% pth2='/MS/DENSITY/Elec2/charge/';
% pth3='/MS/DENSITY/Elec3/charge/';
% pth4='/MS/DENSITY/Elec4/charge/';

str_h1='charge-Elec1-';
str_h2='charge-Elec2-';
str_h3='charge-Elec3-';
str_h4='RAW-Elec1-';
str_h5='RAW-Elec2-';
str_h6='RAW-Elec3-';
% str_h3='charge-Elec3-';
% str_h4='charge-Elec4-';
str_ext='.h5';
str_num='000000';
len_str_num=length(str_num);


qe=-4.8032e-10; %electron charge , statcoulomb
me=9.1094e-28; %electron mass , g
vel_c=2.9979e10; %velocity of light ,  cm/sec
RestMassEnergy_ele=8.19e-14; % Joules or 511kev or 0.511 Mev
%######################################################################

%%%%%%%%%%% Input from user go here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ne=1.1e22; % cm^-3
frno=0; %%%%%%% 75
% NTE=500;
% x00=15.0;
% y00=15.0;
% dx=30.0;
% dy=30.0;

NTE4=1000;
d4a4=12.5;
d5a5=112.5;
d6a6=500.0;

NTE5=round((d5a5/d4a4)*NTE4)
NTE6=round((d6a6/d4a4)*NTE4)



%######################################################################


omp_e=sqrt((4*pi*ne*qe^2)/me); x_nor=vel_c/omp_e; t_nor=1/omp_e;
minch=-5;
maxch=0;

str_frno=num2str(frno);
len_frno=length(str_frno);
str_num((len_str_num-len_frno+1):end)=str_frno;

fl_nm1=strcat(pth0,pth1,str_h1,str_num,str_ext); %charge density
fl_nm2=strcat(pth0,pth2,str_h2,str_num,str_ext); %charge density
fl_nm3=strcat(pth0,pth3,str_h3,str_num,str_ext); %charge density
fl_nm4=strcat(pth0,pth4,str_h4,str_num,str_ext); %raw data
fl_nm5=strcat(pth0,pth5,str_h5,str_num,str_ext); %raw data
fl_nm6=strcat(pth0,pth6,str_h6,str_num,str_ext); %raw data


[xg_cha,yg_cha,dset_qe1,x1lt_qe,x2lt_qe,time_qe]=...
    ReadECharge13Aug2014(fl_nm1);
[xg_cha,yg_cha,dset_qe2,x1lt_qe,x2lt_qe,time_qe]=...
    ReadECharge13Aug2014(fl_nm2);
[xg_cha,yg_cha,dset_qe3,x1lt_qe,x2lt_qe,time_qe]=...
    ReadECharge13Aug2014(fl_nm3);



x1min_qe=x1lt_qe(1); x1max_qe=x1lt_qe(2);
x2min_qe=x2lt_qe(1); x2max_qe=x2lt_qe(2);
[Ngx_cha,Ngy_cha]=size(dset_qe3);
Ngyby2=round(Ngy_cha/2);


%{
xlim1=x1max_qe-40; xlim2=x1max_qe-25;
ylim1=42.5; ylim2=55;

imagesc(flipdim((dset_qe3),2)',...
    'XData',[xg_cha xg_cha],...
    'YData',[yg_cha yg_cha]); hold on;     
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
% ylim([ylim1 ylim2]); 
% xlim([xlim1 xlim2]);
xlabel('x1'); ylabel('x2');

%}


[~,x1_raw4,~,x2_raw4,~,tag_raw4,~,~]=ReadRaw13Aug2014(fl_nm4);
[~,x1_raw5,~,x2_raw5,~,tag_raw5,~,~]=ReadRaw13Aug2014(fl_nm5);
[~,x1_raw6,~,x2_raw6,~,tag_raw6,~,~]=ReadRaw13Aug2014(fl_nm6);


% par_Loc_xy4= (((x1_raw4>(x00-dx/2))&(x1_raw4<(x00+dx/2)))&...
%     ((x2_raw4>(y00-dy/2))&(x2_raw4<(y00+dy/2))));
% profile
% {
% 
% profile_type(1) = "math func", 
% 
% math_func_expr = "if((((x2-15.0)^2)<=49.0)&&(x1<=6),0.1,0.0)",
% 
% }
par_Loc_xy4= (((x2_raw4-12.5).^2)<=2.5^2)&(x1_raw4<=25.0);


% par_Loc_xy5= (((x1_raw5>(x00-dx/2))&(x1_raw5<(x00+dx/2)))&...
%     ((x2_raw5>(y00-dy/2))&(x2_raw5<(y00+dy/2))));
% profile
% {
% profile_type(1) = "math func", 
% math_func_expr = "if((((x2-15.0)^2)<=49.0)&&(x1<=6),0.9,0.0)",
% }
par_Loc_xy5= (((x2_raw5-12.5).^2)<=2.5^2)&(x1_raw5<=25.0);


% par_Loc_xy6= (((x1_raw6>(x00-dx/2))&(x1_raw6<(x00+dx/2)))&...
%     ((x2_raw6>(y00-dy/2))&(x2_raw6<(y00+dy/2))));
% profile
% {
% profile_type(1) = "math func", 
% math_func_expr = "if((((x2-15.0)^2)>49.0)||(x1>6),1.0,0.0)",
% }
par_Loc_xy6= (((x2_raw6-12.5).^2)>=2.5^2)&(x1_raw6<=25.0);


[fid,msg]=fopen('electrons1.tags','w+');
fprintf(fid,'%10d\n',NTE4);
MyTag=tag_raw4(par_Loc_xy4,:);
xx=x1_raw4(par_Loc_xy4,:);
yy=x2_raw4(par_Loc_xy4,:);
len_MyTag=length(MyTag);
ind=(randperm(len_MyTag,NTE4))';
xx4=xx(ind);
yy4=yy(ind);
% len_MyTag=length(MyTag);
% in=len_MyTag/NTE4;
% ind=round(1:in:len_MyTag);
for ii=1:NTE4
    fprintf(fid,'%10d%c%10d\n',MyTag(ind(ii),1),' ',MyTag(ind(ii),2));
end
fclose(fid);

[fid,msg]=fopen('electrons2.tags','w+');
fprintf(fid,'%10d\n',NTE5);
MyTag=tag_raw5(par_Loc_xy5,:);
xx=x1_raw5(par_Loc_xy5,:);
yy=x2_raw5(par_Loc_xy5,:);
len_MyTag=length(MyTag);
ind=(randperm(len_MyTag,NTE5))';
xx5=xx(ind);
yy5=yy(ind);
% len_MyTag=length(MyTag);
% in=len_MyTag/NTE5;
% ind=round(1:in:len_MyTag);
for ii=1:NTE5
    fprintf(fid,'%10d%c%10d\n',MyTag(ind(ii),1),' ',MyTag(ind(ii),2));
end
fclose(fid);


[fid,msg]=fopen('electrons3.tags','w+');
fprintf(fid,'%10d\n',NTE6);
MyTag=tag_raw6(par_Loc_xy6,:);
xx=x1_raw6(par_Loc_xy6,:);
yy=x2_raw6(par_Loc_xy6,:);
len_MyTag=length(MyTag);
ind=(randperm(len_MyTag,NTE6))';
xx6=xx(ind);
yy6=yy(ind);
% len_MyTag=length(MyTag);
% in=len_MyTag/NTE6;
% ind=round(1:in:len_MyTag);
for ii=1:NTE6
    fprintf(fid,'%10d%c%10d\n',MyTag(ind(ii),1),' ',MyTag(ind(ii),2));
end
fclose(fid);


hp=plot(xx4,yy4,xx5,yy5,xx6,yy6);
%     xlim([0 30]); ylim([0 30]);

set(hp(1),'linestyle','none','markeredgecolor','none',...
    'markersize',4,'markerfacecolor','r','marker','o');
set(hp(2),'linestyle','none','markeredgecolor','none',...
    'markersize',4,'markerfacecolor','b','marker','o');
set(hp(3),'linestyle','none','markeredgecolor','none',...
    'markersize',4,'markerfacecolor','g','marker','o');
% title(strcat('timestep=',num2str(tind)));




