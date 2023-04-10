function [xg_e2,yg_e2,dset_e2,x1lt,x2lt,time_e2]=ReadE213Aug2014(fl_nm_e2)

%{
str1_e2='e2-';
str2_e2='000000';
len_str2_e2=length(str2_e2);
str3_e2='.h5';
str_tmp_e2=num2str(fr);
len_tmp_e2=length(str_tmp_e2);
str2a_e2=str2_e2;
str2a_e2((len_str2_e2-len_tmp_e2+1):end)=str_tmp_e2;
fl_nm_e2=strcat(str_pth,str1_e2,str2a_e2,str3_e2);
%}

hinfo_e2=hdf5info(fl_nm_e2);
%{
tmp_struct_e2 = hinfo_e2.GroupHierarchy.Attributes;
xmin_e2=tmp_struct_e2(7).Value;
xmax_e2=tmp_struct_e2(8).Value;
x1min_e2=xmin_e2(1); x1max_e2=xmax_e2(1);
x2min_e2=xmin_e2(2); x2max_e2=xmax_e2(2);
%}
time_e2=hdf5read(hinfo_e2.GroupHierarchy.Attributes(3));
% iter_step_e2=hdf5read(hinfo_e2.GroupHierarchy.Attributes(4));
x1lt=hdf5read(hinfo_e2.GroupHierarchy.Groups.Datasets(1));
x2lt=hdf5read(hinfo_e2.GroupHierarchy.Groups.Datasets(2));
x1min_e2=x1lt(1); x1max_e2=x1lt(2);
x2min_e2=x2lt(1); x2max_e2=x2lt(2);
dset_e2 = hdf5read(hinfo_e2.GroupHierarchy.Datasets);
[Ngx_e2,Ngy_e2]=size(dset_e2);

dx_e2=(x1max_e2-x1min_e2)/(Ngx_e2-1);
dy_e2=(x2max_e2-x2min_e2)/(Ngy_e2-1);

xg_e2=x1min_e2+dx_e2*(0:(Ngx_e2-1)); 
yg_e2=x2min_e2+dy_e2*(0:(Ngy_e2-1));

end