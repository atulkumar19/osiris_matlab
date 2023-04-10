clc;
clear all;
close all;

addpath('UtilFun');

ne=1.1e22; % cm^-3
frm=00:100;

pth0='/home/testuser/GRK_circle_10nm';
pth1='/MS/RAW/Elec1/';
pth2='/MS/RAW/Elec2/';
pth3='/MS/RAW/Elec3/';
pth4='/MS/RAW/Elec4/';
pth5='/MS/RAW/Elec5/';
pth6='/MS/RAW/Elec6/';
pth7='/MS/RAW/Elec7/';
% pth8='/MS/RAW/Elec8/';
% pth9='/MS/RAW/Elec9/';
% pth10='/MS/RAW/Elec10/';
% pth11='/MS/RAW/Elec11/';
% pth12='/MS/RAW/Elec12/';
% pth13='/MS/RAW/Elec13/';
% pth14='/MS/RAW/Elec14/';
% pth15='/MS/RAW/Elec15/';
% pth16='/MS/RAW/Elec16/';
% pth17='/MS/RAW/Elec17/';
% pth18='/MS/RAW/Elec18/';
% pth19='/MS/RAW/Elec19/';
% pth20='/MS/RAW/Elec20/';
% pth21='/MS/RAW/Elec21/';
% pth22='/MS/RAW/Elec22/';
% pth23='/MS/RAW/Elec23/';
% pth24='/MS/RAW/Elec24/';
% pth25='/MS/RAW/Elec25/';
% pth26='/MS/RAW/Elec26/';
% pth27='/MS/RAW/Elec27/';
% pth28='/MS/RAW/Elec28/';
% pth29='/MS/RAW/Elec29/';
% pth30='/MS/RAW/Elec30/';
% pth31='/MS/RAW/Elec31/';
% pth32='/MS/RAW/Elec32/';
% pth33='/MS/RAW/Elec33/';
% pth34='/MS/RAW/Elec34/';
% pth35='/MS/RAW/Elec35/';
% pth36='/MS/RAW/Elec36/';
% pth37='/MS/RAW/Elec37/';
% pth38='/MS/RAW/Elec38/';
% pth39='/MS/RAW/Elec39/';
% pth40='/MS/RAW/Elec40/';
% pth41='/MS/RAW/Elec41/';
% pth42='/MS/RAW/Elec42/';
% pth43='/MS/RAW/Elec43/';
% pth44='/MS/RAW/Elec44/';
% pth45='/MS/RAW/Elec45/';
pth46='/MS/FLD/b3/';


str_h1='RAW-Elec1-';
str_h2='RAW-Elec2-';
str_h3='RAW-Elec3-';
str_h4='RAW-Elec4-';
str_h5='RAW-Elec5-';
str_h6='RAW-Elec6-';
str_h7='RAW-Elec7-';
% str_h8='RAW-Elec8-';
% str_h9='RAW-Elec9-';
% str_h10='RAW-Elec10-';
% str_h11='RAW-Elec11-';
% str_h12='RAW-Elec12-';
% str_h13='RAW-Elec13-';
% str_h14='RAW-Elec14-';
% str_h15='RAW-Elec15-';
% str_h16='RAW-Elec16-';
% str_h17='RAW-Elec17-';
% str_h18='RAW-Elec18-';
% str_h19='RAW-Elec19-';
% str_h20='RAW-Elec20-';
% str_h21='RAW-Elec21-';
% str_h22='RAW-Elec22-';
% str_h23='RAW-Elec23-';
% str_h24='RAW-Elec24-';
% str_h25='RAW-Elec25-';
% str_h26='RAW-Elec26-';
% str_h27='RAW-Elec27-';
% str_h28='RAW-Elec28-';
% str_h29='RAW-Elec29-';
% str_h30='RAW-Elec30-';
% str_h31='RAW-Elec31-';
% str_h32='RAW-Elec32-';
% str_h33='RAW-Elec33-';
% str_h34='RAW-Elec34-';
% str_h35='RAW-Elec35-';
% str_h36='RAW-Elec36-';
% str_h37='RAW-Elec37-';
% str_h38='RAW-Elec38-';
% str_h39='RAW-Elec39-';
% str_h40='RAW-Elec40-';
% str_h41='RAW-Elec41-';
% str_h42='RAW-Elec42-';
% str_h43='RAW-Elec43-';
% str_h44='RAW-Elec44-';
% str_h45='RAW-Elec45-';
str_h46='b3-';


str_ext='.h5';
str_num='000000';
len_str_num=length(str_num);


qe=-4.8032e-10; %electron charge , statcoulomb
me=9.1094e-28; %electron mass , g
vel_c=2.9979e10; %velocity of light ,  cm/sec
omp_e=sqrt((4*pi*ne*qe^2)/me); x_nor=vel_c/omp_e; t_nor=1/omp_e;

len_fr=length(frm);
scrsz = get(0,'ScreenSize');
hfig=figure('Position',[50 10 1200 650]);
% Get the width and height of the figure
%============AVI-FILE-NAME==========================
writerObj = VideoWriter ('e-Trajectory.avi');
writerObj.FrameRate=5;
open(writerObj);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for ii=1:length(frm)    
    
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
%     fl_nm8=strcat(pth0,pth8,str_h8,str_num,str_ext);
%     fl_nm9=strcat(pth0,pth9,str_h9,str_num,str_ext);
%     fl_nm10=strcat(pth0,pth10,str_h10,str_num,str_ext);
%     fl_nm11=strcat(pth0,pth11,str_h11,str_num,str_ext);
%     fl_nm12=strcat(pth0,pth12,str_h12,str_num,str_ext);
%     fl_nm13=strcat(pth0,pth13,str_h13,str_num,str_ext);
%     fl_nm14=strcat(pth0,pth14,str_h14,str_num,str_ext);
%     fl_nm15=strcat(pth0,pth15,str_h15,str_num,str_ext);
%     fl_nm16=strcat(pth0,pth16,str_h16,str_num,str_ext);
%     fl_nm17=strcat(pth0,pth17,str_h17,str_num,str_ext);
%     fl_nm18=strcat(pth0,pth18,str_h18,str_num,str_ext);
%     fl_nm19=strcat(pth0,pth19,str_h19,str_num,str_ext);
%     fl_nm20=strcat(pth0,pth20,str_h20,str_num,str_ext);
%     fl_nm21=strcat(pth0,pth21,str_h21,str_num,str_ext);
%     fl_nm22=strcat(pth0,pth22,str_h22,str_num,str_ext);
%     fl_nm23=strcat(pth0,pth23,str_h23,str_num,str_ext);
%     fl_nm24=strcat(pth0,pth24,str_h24,str_num,str_ext);
%     fl_nm25=strcat(pth0,pth25,str_h25,str_num,str_ext);
%     fl_nm26=strcat(pth0,pth26,str_h26,str_num,str_ext);
%     fl_nm27=strcat(pth0,pth27,str_h27,str_num,str_ext);
%     fl_nm28=strcat(pth0,pth28,str_h28,str_num,str_ext);
%     fl_nm29=strcat(pth0,pth29,str_h29,str_num,str_ext);
%     fl_nm30=strcat(pth0,pth30,str_h30,str_num,str_ext);
%     fl_nm31=strcat(pth0,pth31,str_h31,str_num,str_ext);
%     fl_nm32=strcat(pth0,pth32,str_h32,str_num,str_ext);
%     fl_nm33=strcat(pth0,pth33,str_h33,str_num,str_ext);
%     fl_nm34=strcat(pth0,pth34,str_h34,str_num,str_ext);
%     fl_nm35=strcat(pth0,pth35,str_h35,str_num,str_ext);
%     fl_nm36=strcat(pth0,pth36,str_h36,str_num,str_ext);
%     fl_nm37=strcat(pth0,pth37,str_h37,str_num,str_ext);
%     fl_nm38=strcat(pth0,pth38,str_h38,str_num,str_ext);
%     fl_nm39=strcat(pth0,pth39,str_h39,str_num,str_ext);
%     fl_nm40=strcat(pth0,pth40,str_h40,str_num,str_ext);
%     fl_nm41=strcat(pth0,pth41,str_h41,str_num,str_ext);
%     fl_nm42=strcat(pth0,pth42,str_h42,str_num,str_ext);
%     fl_nm43=strcat(pth0,pth43,str_h43,str_num,str_ext);
%     fl_nm44=strcat(pth0,pth44,str_h44,str_num,str_ext);
%     fl_nm45=strcat(pth0,pth45,str_h45,str_num,str_ext);
    fl_nm46=strcat(pth0,pth46,str_h46,str_num,str_ext);
    
%     fl_nm19=strcat(pth0,pth19,str_h19,str_num,str_ext);
    
%     [xg_qe,yg_qe,dset_qe1,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm1);
%     [xg_qe,yg_qe,dset_qe2,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm2);
%     [xg_qe,yg_qe,dset_qe3,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm3);

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
%     [ene_raw,x1_raw8,p1_raw8,x2_raw8,p2_raw8,tag_raw8,xlt_in,xlt_fi]=...
%         AshReadRaw(fl_nm8);
%     [ene_raw,x1_raw9,p1_raw9,x2_raw9,p2_raw9,tag_raw9,xlt_in,xlt_fi]=...
%         AshReadRaw(fl_nm9);
%     [ene_raw,x1_raw10,p1_raw10,x2_raw10,p2_raw10,tag_raw10,xlt_in,xlt_fi]=...
%         AshReadRaw(fl_nm10);
%     [ene_raw,x1_raw11,p1_raw11,x2_raw11,p2_raw11,tag_raw11,xlt_in,xlt_fi]=...
%         AshReadRaw(fl_nm11);
%     [ene_raw,x1_raw12,p1_raw,x2_raw12,p2_raw,tag_raw,xlt_in,xlt_fi]=...
%         AshReadRaw(fl_nm12);
%     [ene_raw,x1_raw13,p1_raw13,x2_raw13,p2_raw13,tag_raw13,xlt_in,xlt_fi]=...
%         AshReadRaw(fl_nm13);
%     [ene_raw,x1_raw14,p1_raw14,x2_raw14,p2_raw14,tag_raw14,xlt_in,xlt_fi]=...
%         AshReadRaw(fl_nm14);
%     [ene_raw,x1_raw15,p1_raw15,x2_raw15,p2_raw15,tag_raw15,xlt_in,xlt_fi]=...
%         AshReadRaw(fl_nm15);
%     [ene_raw,x1_raw16,p1_raw16,x2_raw16,p2_raw16,tag_raw16,xlt_in,xlt_fi]=...
%         AshReadRaw(fl_nm16);
%     [ene_raw,x1_raw17,p1_raw17,x2_raw17,p2_raw17,tag_raw17,xlt_in,xlt_fi]=...
%         AshReadRaw(fl_nm17);
%     [ene_raw,x1_raw18,p1_raw18,x2_raw18,p2_raw18,tag_raw18,xlt_in,xlt_fi]=...
%         AshReadRaw(fl_nm18);
%     [ene_raw,x1_raw19,p1_raw19,x2_raw19,p2_raw19,tag_raw1,xlt_in,xlt_fi]=...
%         AshReadRaw(fl_nm19);
%     [ene_raw,x1_raw20,p1_raw20,x2_raw20,p2_raw20,tag_raw2,xlt_in,xlt_fi]=...
%         AshReadRaw(fl_nm20);
%     [ene_raw,x1_raw21,p1_raw21,x2_raw21,p2_raw21,tag_raw3,xlt_in,xlt_fi]=...
%         AshReadRaw(fl_nm21);
%     [ene_raw,x1_raw22,p1_raw22,x2_raw22,p2_raw22,tag_raw4,xlt_in,xlt_fi]=...
%         AshReadRaw(fl_nm22);
%     [ene_raw,x1_raw23,p1_raw23,x2_raw23,p2_raw23,tag_raw5,xlt_in,xlt_fi]=...
%         AshReadRaw(fl_nm23);
%     [ene_raw,x1_raw24,p1_raw24,x2_raw24,p2_raw24,tag_raw6,xlt_in,xlt_fi]=...
%         AshReadRaw(fl_nm24);
%     [ene_raw,x1_raw25,p1_raw25,x2_raw25,p2_raw25,tag_raw7,xlt_in,xlt_fi]=...
%         AshReadRaw(fl_nm25);
%     [ene_raw,x1_raw26,p1_raw26,x2_raw26,p2_raw26,tag_raw8,xlt_in,xlt_fi]=...
%         AshReadRaw(fl_nm26);
%     [ene_raw,x1_raw27,p1_raw27,x2_raw27,p2_raw27,tag_raw9,xlt_in,xlt_fi]=...
%         AshReadRaw(fl_nm27);
%     [ene_raw,x1_raw28,p1_raw28,x2_raw28,p2_raw28,tag_raw10,xlt_in,xlt_fi]=...
%         AshReadRaw(fl_nm28);
%     [ene_raw,x1_raw29,p1_raw29,x2_raw29,p2_raw29,tag_raw,xlt_in,xlt_fi]=...
%         AshReadRaw(fl_nm29);
%     [ene_raw,x1_raw30,p1_raw30,x2_raw30,p2_raw30,tag_raw1,xlt_in,xlt_fi]=...
%         AshReadRaw(fl_nm30);
%     [ene_raw,x1_raw31,p1_raw31,x2_raw31,p2_raw31,tag_raw2,xlt_in,xlt_fi]=...
%         AshReadRaw(fl_nm31);
%     [ene_raw,x1_raw32,p1_raw32,x2_raw32,p2_raw32,tag_raw3,xlt_in,xlt_fi]=...
%         AshReadRaw(fl_nm32);
%     [ene_raw,x1_raw33,p1_raw33,x2_raw33,p2_raw33,tag_raw4,xlt_in,xlt_fi]=...
%         AshReadRaw(fl_nm33);
%     [ene_raw,x1_raw34,p1_raw34,x2_raw34,p2_raw34,tag_raw5,xlt_in,xlt_fi]=...
%         AshReadRaw(fl_nm34);
%     [ene_raw,x1_raw35,p1_raw35,x2_raw35,p2_raw35,tag_raw6,xlt_in,xlt_fi]=...
%         AshReadRaw(fl_nm35);
%     [ene_raw,x1_raw36,p1_raw36,x2_raw36,p2_raw36,tag_raw7,xlt_in,xlt_fi]=...
%         AshReadRaw(fl_nm36);
%     [ene_raw,x1_raw37,p1_raw37,x2_raw37,p2_raw37,tag_raw8,xlt_in,xlt_fi]=...
%         AshReadRaw(fl_nm37);
%     [ene_raw,x1_raw38,p1_raw38,x2_raw38,p2_raw38,tag_raw9,xlt_in,xlt_fi]=...
%         AshReadRaw(fl_nm38);
%     [ene_raw,x1_raw39,p1_raw39,x2_raw39,p2_raw39,tag_raw10,xlt_in,xlt_fi]=...
%         AshReadRaw(fl_nm39);
%     [ene_raw,x1_raw40,p1_raw40,x2_raw40,p2_raw40,tag_raw,xlt_in,xlt_fi]=...
%         AshReadRaw(fl_nm40);
%     [ene_raw,x1_raw41,p1_raw41,x2_raw41,p2_raw41,tag_raw11,xlt_in,xlt_fi]=...
%         AshReadRaw(fl_nm41);
%     [ene_raw,x1_raw42,p1_raw42,x2_raw42,p2_raw42,tag_raw12,xlt_in,xlt_fi]=...
%         AshReadRaw(fl_nm42);
%     [ene_raw,x1_raw43,p1_raw43,x2_raw43,p2_raw43,tag_raw,xlt_in,xlt_fi]=...
%         AshReadRaw(fl_nm43);
%     [ene_raw,x1_raw44,p1_raw44,x2_raw44,p2_raw44,tag_raw13,xlt_in,xlt_fi]=...
%         AshReadRaw(fl_nm44);
%     [ene_raw,x1_raw45,p1_raw45,x2_raw45,p2_raw45,tag_raw14,xlt_in,xlt_fi]=...
%         AshReadRaw(fl_nm45);
    
        [xg_b3,yg_b3,dset_b3,x1lt_b3,x2lt_b3,time_b3]=...
        ReadB313Aug2014(fl_nm46);
    
    X1_RAW= [x1_raw1;x1_raw2;x1_raw3;x1_raw4;x1_raw5;x1_raw6;x1_raw7];
     X2_RAW= [x2_raw1;x2_raw2;x2_raw3;x2_raw4;x2_raw5;x2_raw6;x2_raw7];
    
%     X1_RAW= [x1_raw1;x1_raw2;x1_raw3;x1_raw4;x1_raw5;x1_raw6;x1_raw7;x1_raw8;x1_raw9;x1_raw10;x1_raw11;x1_raw12;x1_raw13;x1_raw14;x1_raw15;x1_raw16;x1_raw17;x1_raw18 ...
%              ;x1_raw19;x1_raw20;x1_raw21;x1_raw22;x1_raw23;x1_raw24;x1_raw25;x1_raw26;x1_raw27;x1_raw28;x1_raw29;x1_raw30;x1_raw31;x1_raw32;x1_raw33;x1_raw34;x1_raw35;x1_raw36 ...
%              ;x1_raw37;x1_raw38;x1_raw39;x1_raw40;x1_raw41;x1_raw42;x1_raw43;x1_raw44;x1_raw45];
%     
%      X2_RAW= [x2_raw1;x2_raw2;x2_raw3;x2_raw4;x2_raw5;x2_raw6;x2_raw7;x2_raw8;x2_raw9;x2_raw10;x2_raw11;x2_raw12;x2_raw13;x2_raw14;x2_raw15;x2_raw16;x2_raw17;x2_raw18 ...
%              ;x2_raw19;x2_raw20;x2_raw21;x2_raw22;x2_raw23;x2_raw24;x2_raw25;x2_raw26;x2_raw27;x2_raw28;x2_raw29;x2_raw30;x2_raw31;x2_raw32;x2_raw33;x2_raw34;x2_raw35;x2_raw36 ...
%              ;x2_raw37;x2_raw38;x2_raw39;x2_raw40;x2_raw41;x2_raw42;x2_raw43;x2_raw44;x2_raw45];

   
plot(X1_RAW(1:length(X1_RAW)),X2_RAW(1:length(X2_RAW)),'.r');
    set(gca,'FontSize',14,'FontWeight','Bold');
%       xlim([0 0.002]);
%        ylim([0 0.5e5]);

xlim([199.5 202]);
   
    xlabel('x1'); ylabel('x2');
    title(strcat('e-Trajectory-10nm, t=',num2str(time_b3)));
    
  drawnow
    saveas(gcf,'CurrentFig.png');
    img4=imread('CurrentFig.png');
    writeVideo(writerObj,img4);
     L = graphicsversion('handlegraphics'); 
     currFrame = getframe;
%frame = getframe ( gcf ); % aviobj = addframe ( aviobj, frame );
%      Frame = getframe;
%      writeVideo(writerObj,Frame);

 end
% aviobj = close ( aviobj );
  close(writerObj);