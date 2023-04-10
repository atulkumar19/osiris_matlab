clc;
clear all;
close all;

addpath('UtilFun')

fr_no='010';
fl_nm1=strcat('charge_',fr_no,'.png');
fl_nm2=strcat('p1x1_',fr_no,'.png');

ch=imread(fl_nm1);
px=imread(fl_nm2);

%{
iptsetpref('ImshowBorder','tight');
set(0,'DefaultFigureMenu','none');
format compact; 
set(0,'Default');
%}
%figure('Position', [100, 50, 650, 500]);
subtightplot(2,1,1,[0.000,0.000],[0,0],[0,0]);
imshow(ch)

subtightplot(2,1,2,[0.000,0.000],[0,0],[0,0]);
imshow(px)



