function [xcc_e3,ycc_e3,time_e3,dset_e3,x1lt,x2lt]=ReadE3(fr)

str1_e3='e3-';
str2_e3='000000';
len_str2_e3=length(str2_e3);
str3_e3='.h5';
str_tmp_e3=num2str(fr);
len_tmp_e3=length(str_tmp_e3);
str2a_e3=str2_e3;
str2a_e3((len_str2_e3-len_tmp_e3+1):end)=str_tmp_e3;
fl_nm_e3=strcat(str1_e3,str2a_e3,str3_e3);

hinfo_e3=hdf5info(fl_nm_e3);
%{
tmp_struct_e3 = hinfo_e3.GroupHierarchy.Attributes;
xmin_e3=tmp_struct_e3(7).Value;
xmax_e3=tmp_struct_e3(8).Value;
x1min_e3=xmin_e3(1); x1max_e3=xmax_e3(1);
x2min_e3=xmin_e3(2); x2max_e3=xmax_e3(2);
%}
time_e3=hdf5read(hinfo_e3.GroupHierarchy.Attributes(3));
% iter_step_e3=hdf5read(hinfo_e3.GroupHierarchy.Attributes(4));
x1lt=hdf5read(hinfo_e3.GroupHierarchy.Groups.Datasets(1));
x2lt=hdf5read(hinfo_e3.GroupHierarchy.Groups.Datasets(2));
x1min_e3=x1lt(1); x1max_e3=x1lt(2);
x2min_e3=x2lt(1); x2max_e3=x2lt(2);
dset_e3 = hdf5read(hinfo_e3.GroupHierarchy.Datasets);
[Ncx_e3,Ncy_e3]=size(dset_e3);
% Ncxby2_e3=round(Ncx_e3/2); Ncyby2_e3=round(Ncy_e3/2);
dx_e3=(x1max_e3-x1min_e3)/Ncx_e3; 
xcc_e3=x1min_e3+dx_e3*(1:Ncx_e3)-0.5*dx_e3;
dy_e3=(x2max_e3-x2min_e3)/Ncy_e3; 
ycc_e3=dy_e3*(1:Ncy_e3);

%{
[xmat_e3,ymat_e3]=meshgrid(xcc_cha,ycc_cha);

minch_e3=min(min(dset_e3));
maxch_e3=max(max(dset_e3));
%}


end