function [ene_raw,x1_raw,p1_raw,x2_raw,p2_raw,tag_raw,xlt_in,xlt_fi]=...
    ReadRaw13Aug2014(fl_nm_raw)

%{
% fl_nm='RAW-He_electrons-000000.h5';
str1_raw='RAW-He_electrons-';
str2_raw='000000';
len_str2_raw=length(str2_raw);
str3_raw='.h5';
str_tmp_raw=num2str(frno);
len_tmp_raw=length(str_tmp_raw);
str2a_raw=str2_raw;
str2a_raw((len_str2_raw-len_tmp_raw+1):end)=str_tmp_raw;
fl_nm_raw=strcat(str_pth,str1_raw,str2a_raw,str3_raw);
%}
    
hinfo_raw=hdf5info(fl_nm_raw);

xlt_in=hdf5read(hinfo_raw.GroupHierarchy.Attributes(10));
xlt_fi=hdf5read(hinfo_raw.GroupHierarchy.Attributes(11));

ene_raw=hdf5read(hinfo_raw.GroupHierarchy.Datasets(1));
p1_raw=hdf5read(hinfo_raw.GroupHierarchy.Datasets(2));
p2_raw=hdf5read(hinfo_raw.GroupHierarchy.Datasets(3));
tag_tmp=hdf5read(hinfo_raw.GroupHierarchy.Datasets(6));
tag_raw=tag_tmp';
x1_raw=hdf5read(hinfo_raw.GroupHierarchy.Datasets(7));
x2_raw=hdf5read(hinfo_raw.GroupHierarchy.Datasets(8));

end