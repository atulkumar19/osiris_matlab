function [xg_b3,yg_b3,dset_b3,x1lt,x2lt,time_b3]=ReadB313Aug2014(fl_nm_b3)

%{
str1_b3='b3-';
str2_b3='000000';
len_str2_b3=length(str2_b3);
str3_b3='.h5';
str_tmp_b3=num2str(fr);
len_tmp_b3=length(str_tmp_b3);
str2a_b3=str2_b3;
str2a_b3((len_str2_b3-len_tmp_b3+1):end)=str_tmp_b3;
fl_nm_b3=strcat(str_pth,str1_b3,str2a_b3,str3_b3);
%}

hinfo_b3=hdf5info(fl_nm_b3);
%{
tmp_struct_b3 = hinfo_b3.GroupHierarchy.Attributes;
xmin_b3=tmp_struct_b3(7).Value;
xmax_b3=tmp_struct_b3(8).Value;
x1min_b3=xmin_b3(1); x1max_b3=xmax_b3(1);
x2min_b3=xmin_b3(2); x2max_b3=xmax_b3(2);
%}
time_b3=hdf5read(hinfo_b3.GroupHierarchy.Attributes(3));
% iter_step_b3=hdf5read(hinfo_b3.GroupHierarchy.Attributes(4));
x1lt=hdf5read(hinfo_b3.GroupHierarchy.Groups.Datasets(1));
x2lt=hdf5read(hinfo_b3.GroupHierarchy.Groups.Datasets(2));
x1min_b3=x1lt(1); x1max_b3=x1lt(2);
x2min_b3=x2lt(1); x2max_b3=x2lt(2);
dset_b3 = hdf5read(hinfo_b3.GroupHierarchy.Datasets);
[Ngx_b3,Ngy_b3]=size(dset_b3);

dx_b3=(x1max_b3-x1min_b3)/(Ngx_b3-1);
dy_b3=(x2max_b3-x2min_b3)/(Ngy_b3-1);

xg_b3=x1min_b3+dx_b3*(0:(Ngx_b3-1)); 
yg_b3=x2min_b3+dy_b3*(0:(Ngy_b3-1));

end