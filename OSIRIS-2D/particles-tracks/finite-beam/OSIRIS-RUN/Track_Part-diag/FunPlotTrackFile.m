function [] =FunPlotTrackFile()

clc;
clear all;
close all;

addpath('UtilFun');

TrackFile1='/Volumes/test/2d/f2-l60/Elec1-tracks.h5';
TrackFile2='/Volumes/test/2d/f2-l60/Elec2-tracks.h5';
TrackFile3='/Volumes/test/2d/f2-l60/Elec3-tracks.h5';


data1=ReadTrack(TrackFile1,'electrons1.tags');
data2=ReadTrack(TrackFile2,'electrons2.tags');
data3=ReadTrack(TrackFile3,'electrons3.tags');


% dlmwrite('Ele1.txt', data1, 'delimiter', '\t', ...
%          'precision', '%+10.5e');
% dlmwrite('Ele2.txt', data2, 'delimiter', '\t', ...
%          'precision', '%+10.5e');
% dlmwrite('Ele3.txt', data3, 'delimiter', '\t', ...
%          'precision', '%+10.5e');

% data1=dlmread('Ele1.txt');
% data2=dlmread('Ele2.txt');
% data3=dlmread('Ele3.txt');


% scrsz = get(0,'ScreenSize');
figure('Position',[50 10 1200 650]);
writerObj = VideoWriter('TrackPart.avi');
writerObj.FrameRate=15;
open(writerObj);

NumDataPts1=20;
NumDataPts2=30;
NumDataPts3=60;

for tind=1:1
    
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
    
    
    xx3=data3(tind,1:2:end);
    yy3=data3(tind,2:2:end);
    
    loc=(xx3==0.0&yy3==0.0);
    xx3(loc)=[];
    yy3(loc)=[];
    
    vec1=round(linspace(1,length(xx1),NumDataPts1));
    vec2=round(linspace(1,length(xx2),NumDataPts2));
    vec3=round(linspace(1,length(xx3),NumDataPts3));
    
   hp=plot(xx1(vec1),yy1(vec1),xx2(vec2),yy2(vec2),xx3(vec3),yy3(vec3));
    xlim([0 2*pi]); ylim([0 2*pi]);

    set(hp(1),'linestyle','none','markeredgecolor','none',...
    'color',[1 0 0],'markersize',2,'markerfacecolor','r','marker','o');
    set(gca,'fontsize',20,'fontweight','b');
    xlabel('x1(c/\omega_{pe})'); ylabel('x2(c/\omega_{pe})');

set(hp(2),'linestyle','none','markeredgecolor','none',...
     'color',[0 0 1],'markersize',2,'markerfacecolor','b','marker','o');
 
set(hp(3),'linestyle','none','markeredgecolor','none',...
    'markersize',2,'markerfacecolor','g','marker','o');
    
    title(strcat('e-density: t=',num2str(0.006.*tind)));
    
    saveas(gcf,'CurrentFig.png');
    img4=imread('CurrentFig.png');
    writeVideo(writerObj,img4);

    
    
    
    
    
    plot(xx1,yy1,'or'); hold on;
    plot(xx2,yy2,'ob'); hold on;
    plot(xx3,yy3,'og'); hold off;
    xlim([0 30]);
    ylim([0 30]);
    
    
    drawnow;
end

close(writerObj);


end

function [data]=ReadTrack(TrackFile,TagFile)

TagDat=dlmread(TagFile);
% hinfo_track=hdf5info(TrackFile);

NoElec=TagDat(1,1);
% NoElec=length(hinfo_track.GroupHierarchy.Groups);

str_x1=cell(NoElec,1); str_x2=cell(NoElec,1);
LenTT=0; 

tic
for ii =1:NoElec
%     GrName=hinfo_track.GroupHierarchy.Groups(ii).Name;
    elec_nm=strcat('/',num2str(TagDat((ii+1),1)),'-',...
        num2str(TagDat((ii+1),2)));
    str_x1{ii}=strcat(elec_nm,'/x1');
    str_x2{ii}=strcat(elec_nm,'/x2');
    tt_str=strcat(elec_nm,'/t');
%     ttmp=hdf5read(hinfo_track.GroupHierarchy.Groups(ii).Datasets(13));
    ttmp=hdf5read(TrackFile,tt_str);
    Lentmp=length(ttmp);
    if (LenTT<Lentmp)
        LenTT=Lentmp;
    end    
end
toc

data=zeros(LenTT,NoElec*2);
cnt=1;

tic
for ii =1:NoElec
    
    x1tmp=hdf5read(TrackFile,str_x1{ii});
    LenXX=length(x1tmp);
    data(1:LenXX,cnt)=x1tmp;
    cnt=cnt+1;
    
    x2tmp=hdf5read(TrackFile,str_x2{ii});
    data(1:LenXX,cnt)=x2tmp;
    cnt=cnt+1;   
end
toc
end
