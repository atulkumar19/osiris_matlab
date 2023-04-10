function [] =FunPlotTrackFile()

clc;
clear all
close all

TrackFile1='/Users/akash/Desktop/Users/atul/counter_finite_tag/ParTrack/MS/TRACKS/Elec1-tracks.h5';
TrackFile2='/Users/akash/Desktop/Users/atul/counter_finite_tag/ParTrack/MS/TRACKS/Elec2-tracks.h5';


data1=ReadTrack(TrackFile1);
data2=ReadTrack(TrackFile2);


for tind=1:6000
    
    xx1=data1(tind,1:2:end);
    yy1=data1(tind,2:2:end);
    
    loc=(xx1==0.0&yy1==0.0);
    xx1(loc)=[];
    yy1(loc)=[];
    
    
    xx2=data2(tind,1:2:end);
    yy2=data2(tind,2:2:end);
    
    loc=(xx2==0.0&yy2==0.0);
    xx2(loc)=[];
    yy2(loc)=[];
    
    
    
    
    plot(xx1,yy1,'or'); hold on;
    plot(xx2,yy2,'ob'); hold on;
   
    xlim([0 25]);
    ylim([0 25]);
    title(strcat('timestep=',num2str(tind)));
    
    drawnow;
end

end

function [data]=ReadTrack(TrackFile)

hinfo_track=hdf5info(TrackFile);

NoElec=length(hinfo_track.GroupHierarchy.Groups);
eleind=1:1:NoElec;
LenEleind=length(eleind);

LenTT=0;
for ii =1:LenEleind
    ttmp=hdf5read(hinfo_track.GroupHierarchy. ...
        Groups(ii).Datasets(13));
    if (LenTT<length(ttmp))
        LenTT=length(ttmp);
    end    
end

data=zeros(LenTT,NoElec*2);
cnt=1;
for ii =eleind
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



end
