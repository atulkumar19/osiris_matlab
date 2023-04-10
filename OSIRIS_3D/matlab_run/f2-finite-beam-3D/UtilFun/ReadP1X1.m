function [ene_raw,x1_raw,p1_raw,x2_raw,p2_raw,xlt_in,xlt_fi]=ReadP1X1(frno)

% fl_nm='RAW-He_electrons-000000.h5';
str1_raw='RAW-He_electrons-';
str2_raw='000000';
len_str2_raw=length(str2_raw);
str3_raw='.h5';

str_tmp_raw=num2str(frno);
len_tmp_raw=length(str_tmp_raw);
str2a_raw=str2_raw;
str2a_raw((len_str2_raw-len_tmp_raw+1):end)=str_tmp_raw;
fl_nm_raw=strcat(str1_raw,str2a_raw,str3_raw);

hinfo_raw=hdf5info(fl_nm_raw);

xlt_in=hdf5read(hinfo_raw.GroupHierarchy.Attributes(10));
xlt_fi=hdf5read(hinfo_raw.GroupHierarchy.Attributes(11));

ene_raw=hdf5read(hinfo_raw.GroupHierarchy.Datasets(1));
p1_raw=hdf5read(hinfo_raw.GroupHierarchy.Datasets(2));
p2_raw=hdf5read(hinfo_raw.GroupHierarchy.Datasets(3));
x1_raw=hdf5read(hinfo_raw.GroupHierarchy.Datasets(7));
x2_raw=hdf5read(hinfo_raw.GroupHierarchy.Datasets(8));

end