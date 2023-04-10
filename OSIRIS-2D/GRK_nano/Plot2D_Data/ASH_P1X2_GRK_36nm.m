% clc;
% clear all;
% close all;

addpath('UtilFun');

ne=1.1e22; % cm^-3
frm=10:10;

pth0='/home/testuser/GRK_circle';
pth1='/MS/RAW/Elec1/';
pth2='/MS/RAW/Elec2/';
pth3='/MS/RAW/Elec3/';
pth4='/MS/RAW/Elec4/';
pth5='/MS/RAW/Elec5/';
pth6='/MS/RAW/Elec6/';
pth7='/MS/RAW/Elec7/';
pth8='/MS/RAW/Elec8/';
pth9='/MS/RAW/Elec9/';
pth10='/MS/RAW/Elec10/';
pth11='/MS/RAW/Elec11/';
pth12='/MS/RAW/Elec12/';
pth13='/MS/RAW/Elec13/';
pth14='/MS/RAW/Elec14/';
pth15='/MS/RAW/Elec15/';
pth16='/MS/RAW/Elec16/';
pth17='/MS/RAW/Elec17/';
pth18='/MS/RAW/Elec18/';



str_h1='RAW-Elec1-';
str_h2='RAW-Elec2-';
str_h3='RAW-Elec3-';
str_h4='RAW-Elec4-';
str_h5='RAW-Elec5-';
str_h6='RAW-Elec6-';
str_h7='RAW-Elec7-';
str_h8='RAW-Elec8-';
str_h9='RAW-Elec9-';
str_h10='RAW-Elec10-';
str_h11='RAW-Elec11-';
str_h12='RAW-Elec12-';
str_h13='RAW-Elec13-';
str_h14='RAW-Elec14-';
str_h15='RAW-Elec15-';
str_h16='RAW-Elec16-';
str_h17='RAW-Elec17-';
str_h18='RAW-Elec18-';




str_ext='.h5';
str_num='000000';
len_str_num=length(str_num);


qe=-4.8032e-10; %electron charge , statcoulomb
me=9.1094e-28; %electron mass , g
vel_c=2.9979e10; %velocity of light ,  cm/sec
omp_e=sqrt((4*pi*ne*qe^2)/me); x_nor=vel_c/omp_e; t_nor=1/omp_e;

len_fr=length(frm);
% scrsz = get(0,'ScreenSize');
% hfig=figure('Position',[50 10 1200 650]);
% % Get the width and height of the figure
% %============AVI-FILE-NAME==========================
% writerObj = VideoWriter ('phase_space.avi');
% writerObj.FrameRate=5;
% open(writerObj);
% aviobj = avifile ( 'eden-b-field', 'fps',5); 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for ii=1:5:length(frm)    
    
    frno=frm(ii);
    str_frno=num2str(frno);
    len_frno=length(str_frno);
    str_num((len_str_num-len_frno+1):end)=str_frno;
   fl_nm1=strcat(pth0,pth1,str_h1,str_num,str_ext);
    fl_nm2=strcat(pth0,pth2,str_h2,str_num,str_ext);
    fl_nm3=strcat(pth0,pth3,str_h3,str_num,str_ext);
    fl_nm4=strcat(pth0,pth4,str_h4,str_num,str_ext);
    fl_nm5=strcat(pth0,pth5,str_h5,str_num,str_ext);
    fl_nm6=strcat(pth0,pth6,str_h6,str_num,str_ext);
    fl_nm7=strcat(pth0,pth7,str_h7,str_num,str_ext);
    fl_nm8=strcat(pth0,pth8,str_h8,str_num,str_ext);
    fl_nm9=strcat(pth0,pth9,str_h9,str_num,str_ext);
    fl_nm10=strcat(pth0,pth10,str_h10,str_num,str_ext);
    fl_nm11=strcat(pth0,pth11,str_h11,str_num,str_ext);
    fl_nm12=strcat(pth0,pth12,str_h12,str_num,str_ext);
    fl_nm13=strcat(pth0,pth13,str_h13,str_num,str_ext);
    fl_nm14=strcat(pth0,pth14,str_h14,str_num,str_ext);
    fl_nm15=strcat(pth0,pth15,str_h15,str_num,str_ext);
    fl_nm16=strcat(pth0,pth16,str_h16,str_num,str_ext);
    fl_nm17=strcat(pth0,pth17,str_h17,str_num,str_ext);
    fl_nm18=strcat(pth0,pth18,str_h18,str_num,str_ext);


    [ene_raw,x1_raw1,p1_raw1,x2_raw1,p2_raw1,tag_raw1,xlt_in,xlt_fi]=...
        AshReadRaw(fl_nm1);
    [ene_raw,x1_raw2,p1_raw2,x2_raw2,p2_raw2,tag_raw2,xlt_in,xlt_fi]=...
        AshReadRaw(fl_nm2);
    [ene_raw,x1_raw3,p1_raw3,x2_raw3,p2_raw3,tag_raw3,xlt_in,xlt_fi]=...
        AshReadRaw(fl_nm3);
    [ene_raw,x1_raw4,p1_raw4,x2_raw4,p2_raw4,tag_raw4,xlt_in,xlt_fi]=...
        AshReadRaw(fl_nm4);
    [ene_raw,x1_raw5,p1_raw5,x2_raw5,p2_raw5,tag_raw5,xlt_in,xlt_fi]=...
        AshReadRaw(fl_nm5);
    [ene_raw,x1_raw6,p1_raw6,x2_raw6,p2_raw6,tag_raw6,xlt_in,xlt_fi]=...
        AshReadRaw(fl_nm6);
    [ene_raw,x1_raw7,p1_raw7,x2_raw7,p2_raw7,tag_raw7,xlt_in,xlt_fi]=...
        AshReadRaw(fl_nm7);
    [ene_raw,x1_raw8,p1_raw8,x2_raw8,p2_raw8,tag_raw8,xlt_in,xlt_fi]=...
        AshReadRaw(fl_nm8);
    [ene_raw,x1_raw9,p1_raw9,x2_raw9,p2_raw9,tag_raw9,xlt_in,xlt_fi]=...
        AshReadRaw(fl_nm9);
    [ene_raw,x1_raw10,p1_raw10,x2_raw10,p2_raw10,tag_raw10,xlt_in,xlt_fi]=...
        AshReadRaw(fl_nm4);
    [ene_raw,x1_raw11,p1_raw11,x2_raw11,p2_raw11,tag_raw,xlt_in,xlt_fi]=...
        AshReadRaw(fl_nm10);
   
    [ene_raw,x1_raw12,p1_raw12,x2_raw12,p2_raw12,tag_raw,xlt_in,xlt_fi]=...
        AshReadRaw(fl_nm12);
    [ene_raw,x1_raw13,p1_raw13,x2_raw13,p2_raw13,tag_raw13,xlt_in,xlt_fi]=...
        AshReadRaw(fl_nm13);
    [ene_raw,x1_raw14,p1_raw14,x2_raw14,p2_raw14,tag_raw14,xlt_in,xlt_fi]=...
        AshReadRaw(fl_nm14);
    [ene_raw,x1_raw15,p1_raw15,x2_raw15,p2_raw15,tag_raw15,xlt_in,xlt_fi]=...
        AshReadRaw(fl_nm15);
    [ene_raw,x1_raw16,p1_raw16,x2_raw16,p2_raw16,tag_raw16,xlt_in,xlt_fi]=...
        AshReadRaw(fl_nm16);
    [ene_raw,x1_raw17,p1_raw17,x2_raw17,p2_raw17,tag_raw17,xlt_in,xlt_fi]=...
        AshReadRaw(fl_nm17);
    [ene_raw,x1_raw18,p1_raw18,x2_raw18,p2_raw18,tag_raw18,xlt_in,xlt_fi]=...
        AshReadRaw(fl_nm18);
   
    
    

    p_raw_sq1=p1_raw1.*p1_raw1+p2_raw1.*p2_raw1; p_raw_sq2=p1_raw2.*p1_raw2+p2_raw2.*p2_raw2;
    p_raw_sq3=p1_raw3.*p1_raw3+p2_raw3.*p2_raw3; p_raw_sq4=p1_raw4.*p1_raw4+p2_raw4.*p2_raw4;
    p_raw_sq5=p1_raw5.*p1_raw5+p2_raw5.*p2_raw5; p_raw_sq6=p1_raw6.*p1_raw6+p2_raw6.*p2_raw6;
    p_raw_sq7=p1_raw7.*p1_raw7+p2_raw7.*p2_raw7; p_raw_sq8=p1_raw8.*p1_raw8+p2_raw8.*p2_raw8;
    p_raw_sq9=p1_raw9.*p1_raw9+p2_raw9.*p2_raw9; p_raw_sq10=p1_raw10.*p1_raw10+p2_raw10.*p2_raw10;
    p_raw_sq11=p1_raw11.*p1_raw11+p2_raw11.*p2_raw11; p_raw_sq12=p1_raw12.*p1_raw12+p2_raw12.*p2_raw12;
    p_raw_sq13=p1_raw13.*p1_raw13+p2_raw13.*p2_raw13; p_raw_sq14=p1_raw14.*p1_raw14+p2_raw14.*p2_raw14;
    p_raw_sq15=p1_raw15.*p1_raw15+p2_raw15.*p2_raw15; p_raw_sq16=p1_raw16.*p1_raw16+p2_raw16.*p2_raw16;
    p_raw_sq17=p1_raw17.*p1_raw17+p2_raw17.*p2_raw17; p_raw_sq18=p1_raw18.*p1_raw18+p2_raw18.*p2_raw18;
   
%     p_raw_sq= [p_raw_sq1; p_raw_sq2; p_raw_sq3; p_raw_sq4; p_raw_sq5; p_raw_sq6; p_raw_sq7....
%     ; p_raw_sq8; p_raw_sq9; p_raw_sq10; p_raw_sq11; p_raw_sq12; p_raw_sq13; p_raw_sq14....
%     ; p_raw_sq15; p_raw_sq16; p_raw_sq17; p_raw_sq18];

    X1_RAW=[x1_raw1;x1_raw2;x1_raw3;x1_raw4;x1_raw5;x1_raw6;x1_raw7;x1_raw8;x1_raw9;x1_raw10;x1_raw11;x1_raw12;x1_raw13;x1_raw14;x1_raw15;x1_raw16;x1_raw17;x1_raw18];
    X2_RAW=[x2_raw1;x2_raw2;x2_raw3;x2_raw4;x2_raw5;x2_raw6;x2_raw7;x2_raw8;x2_raw9;x2_raw10;x2_raw11;x2_raw12;x2_raw13;x2_raw14;x2_raw15;x2_raw16;x2_raw17;x2_raw18];
%     p_raw_sq= [p_raw_sq1; p_raw_sq2; p_raw_sq3; p_raw_sq4; p_raw_sq5; p_raw_sq6; p_raw_sq7....
%     ; p_raw_sq8; p_raw_sq9; p_raw_sq10; p_raw_sq11; p_raw_sq12; p_raw_sq13; p_raw_sq14....
%     ; p_raw_sq15; p_raw_sq16; p_raw_sq17; p_raw_sq18;p_raw_sq19;p_raw_sq20;p_raw_sq21;p_raw_sq22 ...
%     ; p_raw_sq23;p_raw_sq24;p_raw_sq25;p_raw_sq26;p_raw_sq27;p_raw_sq28;p_raw_sq29];
    p_raw_sq= [p_raw_sq1; p_raw_sq2; p_raw_sq3; p_raw_sq4; p_raw_sq5; p_raw_sq6; p_raw_sq7....
    ; p_raw_sq8; p_raw_sq9; p_raw_sq10; p_raw_sq11; p_raw_sq12; p_raw_sq13; p_raw_sq14....
    ; p_raw_sq15; p_raw_sq16; p_raw_sq17; p_raw_sq18];

% p_raw_sq= [p_raw_sq1; p_raw_sq2; p_raw_sq3; p_raw_sq4; p_raw_sq5; p_raw_sq6; p_raw_sq7];
% p_raw_sq= [p_raw_sq1; p_raw_sq2; p_raw_sq3; p_raw_sq4; p_raw_sq5; p_raw_sq6; p_raw_sq7];
%     KE_raw=sqrt(p_raw_sq+1)-1; 
%     
% KE_raw=sqrt(p_raw_sq+1)-1; 
%     KE_raw=0.5e3*KE_raw;
%      [counts,centers] = hist(KE_raw,170);
%      semilogy(centers,counts,'*g');hold on;
plot(X1_RAW(1:1:length(X1_RAW)),X2_RAW(1:1:length(X2_RAW)),'.b')
%      f=(gamma(5/2))^(3/2)*(gamma(3/2))^(-5/2)*exp(-KE_raw*gamma(5/2)*(gamma(3/2))^(-1));
%     
%      plot( KE_raw(1:1:length(KE_raw)),'.');hold on;
% hist(KE_raw,500);
% [f,x]=hist(KE_raw,500); %use hist function and get unnormalized values
% figure; plot(x,f/trapz(x,f),'b-*');
% [bincounts,ind] = histc(KE_raw,500);

    set(gca,'FontSize',12,'FontWeight','Bold');
%       xlim([0 0.002]);
%        ylim([0 0.5e5]);
xlim([199 203])
   
%     xlabel('KE(KeV)'); ylabel('counts');
%     title(strcat('Phase-Space Plot, t=',num2str(time_b3)));
    
  drawnow
%     saveas(gcf,'CurrentFig.png');
%     img4=imread('CurrentFig.png');
% %     writeVideo(writerObj,img4);
%      L = graphicsversion('handlegraphics'); 
%      currFrame = getframe;
%frame = getframe ( gcf ); % aviobj = addframe ( aviobj, frame );
%      Frame = getframe;
%      writeVideo(writerObj,Frame);

 end
% aviobj = close ( aviobj );
%   close(writerObj);