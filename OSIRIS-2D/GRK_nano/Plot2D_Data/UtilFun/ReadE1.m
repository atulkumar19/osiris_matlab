function [xg_e1,yg_e1,time_e1,dset_e1,x1lt,x2lt]=ReadE1(fr)

str1_e1='e1-';
str2_e1='000000';
len_str2_e1=length(str2_e1);
str3_e1='.h5';
str_tmp_e1=num2str(fr);
len_tmp_e1=length(str_tmp_e1);
str2a_e1=str2_e1;
str2a_e1((len_str2_e1-len_tmp_e1+1):end)=str_tmp_e1;
fl_nm_e1=strcat(str1_e1,str2a_e1,str3_e1);

hinfo_e1=hdf5info(fl_nm_e1);
%{
tmp_struct_e1 = hinfo_e1.GroupHierarchy.Attributes;
xmin_e1=tmp_struct_e1(7).Value;
xmax_e1=tmp_struct_e1(8).Value;
x1min_e1=xmin_e1(1); x1max_e1=xmax_e1(1);
x2min_e1=xmin_e1(2); x2max_e1=xmax_e1(2);
%}
time_e1=hdf5read(hinfo_e1.GroupHierarchy.Attributes(3));
% iter_step_e1=hdf5read(hinfo_e1.GroupHierarchy.Attributes(4));
x1lt=hdf5read(hinfo_e1.GroupHierarchy.Groups.Datasets(1));
x2lt=hdf5read(hinfo_e1.GroupHierarchy.Groups.Datasets(2));
x1min_e1=x1lt(1); x1max_e1=x1lt(2);
x2min_e1=x2lt(1); x2max_e1=x2lt(2);
dset_e1 = hdf5read(hinfo_e1.GroupHierarchy.Datasets);

[Ngx_e1,Ngy_e1]=size(dset_e1);
dx_e1=(x1max_e1-x1min_e1)/(Ngx_e1-1);
dy_e1=(x2max_e1-x2min_e1)/(Ngy_e1-1);
xg_e1=x1min_e1+dx_e1*(0:(Ngx_e1-1)); 
yg_e1=x2min_e1+dy_e1*(0:(Ngy_e1-1));

end