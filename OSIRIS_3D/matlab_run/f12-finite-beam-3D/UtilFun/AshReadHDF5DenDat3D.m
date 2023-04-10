function [xg_cha,yg_cha,zg_cha,dset_cha,x1lt,x2lt,x3lt,time_cha]=...
    AshReadHDF5DenDat3D(fl_nm_cha)

%{
str1_cha='charge-Elec1-';
str2_cha='000000';
len_str2_cha=length(str2_cha);
str3_cha='.h5';
str_tmp_cha=num2str(frno);
len_tmp_cha=length(str_tmp_cha);
str2a_cha=str2_cha;
str2a_cha((len_str2_cha-len_tmp_cha+1):end)=str_tmp_cha;

fl_nm_cha=strcat(str_pth,str1_cha,str2a_cha,str3_cha);
%}

hinfo_cha=hdf5info(fl_nm_cha);

%{
tmp_struct_cha = hinfo_cha.GroupHierarchy.Attributes;
time_cha=tmp_struct_cha(3).Value;
iter_step_cha=tmp_struct_cha(4).Value;
xmin_cha=tmp_struct_cha(7).Value;
xmax_cha=tmp_struct_cha(8).Value;
x1min_cha=xmin_cha(1); x1max_cha=xmax_cha(1);
x2min_cha=xmin_cha(2); x2max_cha=xmax_cha(2);
%}

time_cha=hdf5read(hinfo_cha.GroupHierarchy.Attributes(3));
% iter_step_cha=hdf5read(hinfo_cha.GroupHierarchy.Attributes(4));
x1lt=hdf5read(hinfo_cha.GroupHierarchy.Groups.Datasets(1));
x2lt=hdf5read(hinfo_cha.GroupHierarchy.Groups.Datasets(2));
x3lt=hdf5read(hinfo_cha.GroupHierarchy.Groups.Datasets(3));
x1min_cha=x1lt(1); x1max_cha=x1lt(2);
x2min_cha=x2lt(1); x2max_cha=x2lt(2);
x3min_cha=x3lt(1); x3max_cha=x3lt(2);
dset_cha = double(hdf5read(hinfo_cha.GroupHierarchy.Datasets));
[Ngx_cha,Ngy_cha,Ngz_cha]=size(dset_cha);
% Ncxby2_cha=round(Ncx_cha/2); Ncyby2_cha=round(Ncy_cha/2);
dx_cha=(x1max_cha-x1min_cha)/(Ngx_cha-1);
dy_cha=(x2max_cha-x2min_cha)/(Ngy_cha-1);
dz_cha=(x3max_cha-x3min_cha)/(Ngz_cha-1);

xg_cha=x1min_cha+dx_cha*(0:(Ngx_cha-1)); 
yg_cha=x2min_cha+dy_cha*(0:(Ngy_cha-1));
zg_cha=x3min_cha+dz_cha*(0:(Ngz_cha-1));

end
