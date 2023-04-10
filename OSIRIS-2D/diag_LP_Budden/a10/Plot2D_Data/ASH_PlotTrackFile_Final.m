clc;
clear all
close all

TrackFile='/home/bhavesh/Desktop/OSIRIS/OsirisRao/bin/PTPICSim/MS/TRACKS/He_electrons-tracks.h5';
hinfo_track=hdf5info(TrackFile);

NoElec=length(hinfo_track.GroupHierarchy.Groups);

LenTT=0;
for ii =1:NoElec
    ttmp=hdf5read(hinfo_track.GroupHierarchy.Groups(ii).Datasets(13));
    if (LenTT<length(ttmp))
        LenTT=length(ttmp);
    end    
end

data=zeros(LenTT,NoElec*2);
cnt=1;
for ii =1:NoElec
%     ttmp=hdf5read(hinfo_track.GroupHierarchy. ...
%         Groups(ii).Datasets(13));
    x1tmp=hdf5read(hinfo_track.GroupHierarchy. ...
        Groups(ii).Datasets(14));
    LenXX=length(x1tmp);
    data(1:LenXX,cnt)=x1tmp;
    cnt=cnt+1;
    
    x2tmp=hdf5read(hinfo_track.GroupHierarchy. ...
        Groups(ii).Datasets(15));
    data(1:LenXX,cnt)=x2tmp;
    cnt=cnt+1;
    
end

for tind=1:1
    
    xx=data(tind,1:2:end);
    yy=data(tind,2:2:end);
    
    loc=(xx==0.0&yy==0.0);
    xx(loc)=[];
    yy(loc)=[];
    
    plot(xx,yy,'or');
    xlim([0 30]);
    ylim([0 30]);
    
    drawnow;
end
