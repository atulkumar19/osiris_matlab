function [raw_dt,raw_dmp_fac,raw_x2min,raw_x2max]=ReadRawDT(frno)

% global hinfo_raw

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

raw_dt=hdf5read(hinfo_raw.GroupHierarchy.Attributes(8));
raw_dmp_fac=hdf5read(hinfo_raw.GroupHierarchy.Attributes(7));
xtmp=hdf5read(hinfo_raw.GroupHierarchy.Attributes(11));
raw_x2max=xtmp(2);
raw_x2min=0.0;

end