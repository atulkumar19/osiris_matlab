function [ene_ele,E1_ele,E2_ele,p1_ele,p2_ele,t_ele,x1_ele,x2_ele] =...
    ReadTrack13Aug2014(ele_tag,fl_nm)

% fl_nm='He_electrons-tracks.h5';
ele_tag=strcat('/',ele_tag);

hinfo_track=hdf5info(fl_nm);
NoElec=length(hinfo_track.GroupHierarchy.Groups);
my_elec=-1;
for ii=1:NoElec
    nam=hinfo_track.GroupHierarchy.Groups(ii).Name;
    if (strcmp(nam,ele_tag))
        my_elec=ii;
        break;
    end
end

if (ii<0)
    E1_ele=0.0;E2_ele=0.0;p1_ele=0.0;p2_ele=0.0;
    t_ele=0.0;x1_ele=0.0;x2_ele=0.0;
    return;
end

%{
B1_ele=hdf5read(hinfo_track.GroupHierarchy. ...
    Groups(my_elec).Datasets(1));
B2_ele=hdf5read(hinfo_track.GroupHierarchy. ...
    Groups(my_elec).Datasets(2));
B3_ele=hdf5read(hinfo_track.GroupHierarchy. ...
    Groups(my_elec).Datasets(3));
%}

E1_ele=hdf5read(hinfo_track.GroupHierarchy. ...
    Groups(my_elec).Datasets(4));
E2_ele=hdf5read(hinfo_track.GroupHierarchy. ...
    Groups(my_elec).Datasets(5));

%{
E3_ele=hdf5read(hinfo_track.GroupHierarchy. ...
    Groups(my_elec).Datasets(6));
%}
   
ene_ele=hdf5read(hinfo_track.GroupHierarchy. ...
    Groups(my_elec).Datasets(7));

%{
n_ele=hdf5read(hinfo_track.GroupHierarchy. ...
    Groups(my_elec).Datasets(8));
%}
    
p1_ele=hdf5read(hinfo_track.GroupHierarchy. ...
    Groups(my_elec).Datasets(9));
p2_ele=hdf5read(hinfo_track.GroupHierarchy. ...
    Groups(my_elec).Datasets(10));

%{
p3_ele=hdf5read(hinfo_track.GroupHierarchy. ...
    Groups(my_elec).Datasets(11));
q_ele=hdf5read(hinfo_track.GroupHierarchy. ...
    Groups(my_elec).Datasets(12));
%}
    
t_ele=hdf5read(hinfo_track.GroupHierarchy. ...
    Groups(my_elec).Datasets(13));
x1_ele=hdf5read(hinfo_track.GroupHierarchy. ...
    Groups(my_elec).Datasets(14));
x2_ele=hdf5read(hinfo_track.GroupHierarchy. ...
    Groups(my_elec).Datasets(15));

end