% clc;
% clear all;
% close all;

addpath('UtilFun');

ne=1.1e22; % cm^-3
frm=1:600;



pth0='/home/testuser/GRK_circle_62nm';
pth1='/MS/DENSITY/Elec1/ene/';
pth2='/MS/DENSITY/Elec2/ene/';
pth3='/MS/DENSITY/Elec3/ene/';
pth4='/MS/DENSITY/Elec4/ene/';
pth5='/MS/DENSITY/Elec5/ene/';
pth6='/MS/DENSITY/Elec6/ene/';
pth7='/MS/DENSITY/Elec7/ene/';
pth8='/MS/DENSITY/Elec8/ene/';
pth9='/MS/DENSITY/Elec9/ene/';
pth10='/MS/DENSITY/Elec10/ene/';
pth11='/MS/DENSITY/Elec11/ene/';
pth12='/MS/DENSITY/Elec12/ene/';
pth13='/MS/DENSITY/Elec13/ene/';
pth14='/MS/DENSITY/Elec14/ene/';
pth15='/MS/DENSITY/Elec15/ene/';
pth16='/MS/DENSITY/Elec16/ene/';
pth17='/MS/DENSITY/Elec17/ene/';
pth18='/MS/DENSITY/Elec18/ene/';

pth19='/MS/DENSITY/Elec19/ene/';
pth20='/MS/DENSITY/Elec20/ene/';
pth21='/MS/DENSITY/Elec21/ene/';
pth22='/MS/DENSITY/Elec22/ene/';
pth23='/MS/DENSITY/Elec23/ene/';
pth24='/MS/DENSITY/Elec24/ene/';
pth25='/MS/DENSITY/Elec25/ene/';
pth26='/MS/DENSITY/Elec26/ene/';
pth27='/MS/DENSITY/Elec27/ene/';
pth28='/MS/DENSITY/Elec28/ene/';
pth29='/MS/DENSITY/Elec29/ene/';


pth37='/MS/FLD/ene_emf/';
pth38='/MS/FLD/ene_e/';
pth39='/MS/FLD/ene_b/';


str_h1='ene-Elec1-';
str_h2='ene-Elec2-';
str_h3='ene-Elec3-';
str_h4='ene-Elec4-';
str_h5='ene-Elec5-';
str_h6='ene-Elec6-';
str_h7='ene-Elec7-';
str_h8='ene-Elec8-';
str_h9='ene-Elec9-';
str_h10='ene-Elec10-';
str_h11='ene-Elec11-';
str_h12='ene-Elec12-';
str_h13='ene-Elec13-';
str_h14='ene-Elec14-';
str_h15='ene-Elec15-';
str_h16='ene-Elec16-';
str_h17='ene-Elec17-';
str_h18='ene-Elec18-';

str_h19='ene-Elec19-';
str_h20='ene-Elec20-';
str_h21='ene-Elec21-';
str_h22='ene-Elec22-';
str_h23='ene-Elec23-';
str_h24='ene-Elec24-';
str_h25='ene-Elec25-';
str_h26='ene-Elec26-';
str_h27='ene-Elec27-';
str_h28='ene-Elec28-';
str_h29='ene-Elec29-';

% str_h30='ene-Elec12-';
% str_h31='ene-Elec13-';
% str_h32='ene-Elec14-';
% str_h33='ene-Elec15-';
% str_h34='ene-Elec16-';
% str_h35='ene-Elec17-';
% str_h36='ene-Elec18-';



str_h37='ene_emf-';
str_h38='ene_e-';
str_h39='ene_b-';

str_ext='.h5';
str_num='000000';
len_str_num=length(str_num);

qe=-4.8032e-10; %electron charge , statcoulomb
me=9.1094e-28; %electron mass , g
vel_c=2.9979e10; %velocity of light ,  cm/sec
omp_e=sqrt((4*pi*ne*qe^2)/me); x_nor=vel_c/omp_e; t_nor=1/omp_e;

minch=-5;
maxch=0;
scrsz=get(0,'ScreenSize');
figure('Position',[50 10 1200 650]);
% Get the width and height of the figure
MagE=zeros(length(frm),1);
ElecE=zeros(length(frm),1);
EME=zeros(length(frm),1);
KinE=zeros(length(frm),1);
timE=zeros(length(frm),1);
TotE=zeros(length(frm),1);
KinE0=zeros(length(frm),1);
KinE_per=zeros(length(frm),1);
Absorption=zeros(length(frm),1);
cnt=1;
% %============AVI-FILE-NAME==========================
%  aviobj = avifile ( 'squeeze-eden-b-cur', 'fps',5); 
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
    fl_nm19=strcat(pth0,pth19,str_h19,str_num,str_ext);
    fl_nm20=strcat(pth0,pth20,str_h20,str_num,str_ext);
    fl_nm21=strcat(pth0,pth21,str_h21,str_num,str_ext);
    fl_nm22=strcat(pth0,pth22,str_h22,str_num,str_ext);
    fl_nm23=strcat(pth0,pth23,str_h23,str_num,str_ext);
    fl_nm24=strcat(pth0,pth24,str_h24,str_num,str_ext);
    fl_nm25=strcat(pth0,pth25,str_h25,str_num,str_ext);
    fl_nm26=strcat(pth0,pth26,str_h26,str_num,str_ext);
    fl_nm27=strcat(pth0,pth27,str_h27,str_num,str_ext);
    fl_nm28=strcat(pth0,pth28,str_h28,str_num,str_ext);
    fl_nm29=strcat(pth0,pth29,str_h29,str_num,str_ext);
%     fl_nm30=strcat(pth0,pth30,str_h30,str_num,str_ext);
%     fl_nm31=strcat(pth0,pth31,str_h31,str_num,str_ext);
%     fl_nm32=strcat(pth0,pth32,str_h32,str_num,str_ext);
%     fl_nm33=strcat(pth0,pth33,str_h33,str_num,str_ext);
%     fl_nm34=strcat(pth0,pth34,str_h34,str_num,str_ext);
%     fl_nm35=strcat(pth0,pth35,str_h35,str_num,str_ext);
%     fl_nm36=strcat(pth0,pth36,str_h36,str_num,str_ext);
%     
    
    fl_nm37=strcat(pth0,pth37,str_h37,str_num,str_ext);
    fl_nm38=strcat(pth0,pth38,str_h38,str_num,str_ext);
    fl_nm39=strcat(pth0,pth39,str_h39,str_num,str_ext);
    
    
    
%     [xg_qe,yg_qe,dset_qe1,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm1);
%     [xg_qe,yg_qe,dset_qe2,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm2);
%     [xg_qe,yg_qe,dset_qe3,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm3);
%     [xg_qe,yg_qe,dset_qe4,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm4);
%     [xg_qe,yg_qe,dset_qe5,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm5);
%      [xg_qe,yg_qe,dset_qe6,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm6);
%     [xg_qe,yg_qe,dset_qe7,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm7);
%     [xg_qe,yg_qe,dset_qe8,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm8);
%     [xg_qe,yg_qe,dset_qe9,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm9);
%     [xg_qe,yg_qe,dset_qe10,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm10);
%      [xg_qe,yg_qe,dset_qe11,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm11);
%     [xg_qe,yg_qe,dset_qe12,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm12);
%     [xg_qe,yg_qe,dset_qe13,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm13);
%     [xg_qe,yg_qe,dset_qe14,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm14);
%     [xg_qe,yg_qe,dset_qe15,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm15);
%     [xg_qe,yg_qe,dset_qe16,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm16);
%     [xg_qe,yg_qe,dset_qe17,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm17);
%     [xg_qe,yg_qe,dset_qe18,x1lt_qe,x2lt_qe,time_qe,]=...
%         ReadECharge13Aug2014(fl_nm18);
    
   
    [xg_ene1,yg_ene1,dset_ene1,x1lt_ene1,x2lt_ene1,time_ene1]=...
        ReadECharge13Aug2014(fl_nm1);
    [xg_ene2,yg_ene2,dset_ene2,x1lt_ene2,x2lt_ene2,time_ene2]=...
        ReadECharge13Aug2014(fl_nm2);
    [xg_ene3,yg_ene3,dset_ene3,x1lt_ene3,x2lt_ene3,time_ene3]=...
        ReadECharge13Aug2014(fl_nm3);
     [xg_ene1,yg_ene1,dset_ene4,x1lt_ene1,x2lt_ene1,time_ene1]=...
        ReadECharge13Aug2014(fl_nm4);
    [xg_ene2,yg_ene2,dset_ene5,x1lt_ene2,x2lt_ene2,time_ene2]=...
        ReadECharge13Aug2014(fl_nm5);
    [xg_ene3,yg_ene3,dset_ene6,x1lt_ene3,x2lt_ene3,time_ene3]=...
        ReadECharge13Aug2014(fl_nm6);
    [xg_ene1,yg_ene1,dset_ene7,x1lt_ene1,x2lt_ene1,time_ene1]=...
        ReadECharge13Aug2014(fl_nm7);
    [xg_ene2,yg_ene2,dset_ene8,x1lt_ene2,x2lt_ene2,time_ene2]=...
        ReadECharge13Aug2014(fl_nm8);
    [xg_ene3,yg_ene3,dset_ene9,x1lt_ene3,x2lt_ene3,time_ene3]=...
        ReadECharge13Aug2014(fl_nm9);
     [xg_ene1,yg_ene1,dset_ene10,x1lt_ene1,x2lt_ene1,time_ene1]=...
        ReadECharge13Aug2014(fl_nm10);
    [xg_ene2,yg_ene2,dset_ene11,x1lt_ene2,x2lt_ene2,time_ene2]=...
        ReadECharge13Aug2014(fl_nm11);
    [xg_ene3,yg_ene3,dset_ene12,x1lt_ene3,x2lt_ene3,time_ene3]=...
        ReadECharge13Aug2014(fl_nm12);
    [xg_ene1,yg_ene1,dset_ene13,x1lt_ene1,x2lt_ene1,time_ene1]=...
        ReadECharge13Aug2014(fl_nm13);
    [xg_ene2,yg_ene2,dset_ene14,x1lt_ene2,x2lt_ene2,time_ene2]=...
        ReadECharge13Aug2014(fl_nm14);
    [xg_ene3,yg_ene3,dset_ene15,x1lt_ene3,x2lt_ene3,time_ene3]=...
        ReadECharge13Aug2014(fl_nm15);
     [xg_ene1,yg_ene1,dset_ene16,x1lt_ene1,x2lt_ene1,time_ene1]=...
        ReadECharge13Aug2014(fl_nm16);
    [xg_ene2,yg_ene2,dset_ene17,x1lt_ene2,x2lt_ene2,time_ene2]=...
        ReadECharge13Aug2014(fl_nm17);
    [xg_ene3,yg_ene3,dset_ene18,x1lt_ene3,x2lt_ene3,time_ene3]=...
        ReadECharge13Aug2014(fl_nm18);
    [xg_ene1,yg_ene1,dset_ene19,x1lt_ene1,x2lt_ene1,time_ene1]=...
        ReadECharge13Aug2014(fl_nm19);
    [xg_ene2,yg_ene2,dset_ene20,x1lt_ene2,x2lt_ene2,time_ene2]=...
        ReadECharge13Aug2014(fl_nm20);
    [xg_ene3,yg_ene3,dset_ene21,x1lt_ene3,x2lt_ene3,time_ene3]=...
        ReadECharge13Aug2014(fl_nm21);
     [xg_ene1,yg_ene1,dset_ene22,x1lt_ene1,x2lt_ene1,time_ene1]=...
        ReadECharge13Aug2014(fl_nm22);
    [xg_ene2,yg_ene2,dset_ene23,x1lt_ene2,x2lt_ene2,time_ene2]=...
        ReadECharge13Aug2014(fl_nm23);
    [xg_ene3,yg_ene3,dset_ene24,x1lt_ene3,x2lt_ene3,time_ene3]=...
        ReadECharge13Aug2014(fl_nm24);
    [xg_ene1,yg_ene1,dset_ene25,x1lt_ene1,x2lt_ene1,time_ene1]=...
        ReadECharge13Aug2014(fl_nm25);
    [xg_ene2,yg_ene2,dset_ene26,x1lt_ene2,x2lt_ene2,time_ene2]=...
        ReadECharge13Aug2014(fl_nm26);
    [xg_ene3,yg_ene3,dset_ene27,x1lt_ene3,x2lt_ene3,time_ene3]=...
        ReadECharge13Aug2014(fl_nm27);
     [xg_ene1,yg_ene1,dset_ene28,x1lt_ene1,x2lt_ene1,time_ene1]=...
        ReadECharge13Aug2014(fl_nm28);
    [xg_ene2,yg_ene2,dset_ene29,x1lt_ene2,x2lt_ene2,time_ene2]=...
        ReadECharge13Aug2014(fl_nm29);
    
    
    [xg_ene_EMF,yg_ene_EMF,dset_ene_EMF,x1lt_ene_EMF,x2lt_ene_EMF,time_ene_EMF]=...
        ReadB313Aug2014(fl_nm37);
    
    [xg_ene_e,yg_ene_e,dset_ene_e,x1lt_ene_e,x2lt_ene_e,time_ene_e]=...
        ReadE113Aug2014(fl_nm38);
    
    [xg_ene_b,yg_ene_b,dset_ene_b,x1lt_ene_b,x2lt_ene_b,time_ene_b]=...
        ReadB313Aug2014(fl_nm39);
    
    
    dset_qe=dset_qe1+dset_qe2+dset_qe3;
    dset_ene=dset_ene1+dset_ene2+dset_ene3+dset_ene4+dset_ene5+dset_ene6+dset_ene7+dset_ene8+dset_ene9+...
             dset_ene10+dset_ene11+dset_ene12+dset_ene13+dset_ene14+dset_ene15+dset_ene16+dset_ene17+dset_ene18+...
             dset_ene19+dset_ene20+dset_ene21+dset_ene22+dset_ene23+dset_ene24+dset_ene25+dset_ene26+dset_ene27+dset_ene28+dset_ene29;
    dset_tot=dset_ene+dset_ene_EMF;
   
    dx_qe=abs(xg_ene1(2)-xg_ene1(1));
    
   
    
    KinE(cnt)=sum(sum(dset_ene))*dx_qe.^2;
%     MagE(cnt)=sum(sum(dset_ene_b))*dx_qe.^2;
%     ElecE(cnt)=sum(sum(dset_ene_e))*dx_qe.^2;
%     
%     TotE(cnt)=sum(sum((1/2)*dset_ene_EMF+dset_ene))*dx_qe.^2;
    timE(cnt)=time_ene1;
   
    cnt=cnt+1;


end
% 
% for i=1:length(frm)-1
% 
% Absorption(i)=KinE(i+1)-KinE(i);
% 
% end

figure(1)
plot(timE, KinE,'--r','LineWidth',2);
hold on;
set(gca,'fontsize',20,'fontweight','b');
xlabel('time(\omega_{pe})');
ylabel('KinEnergy');


% hold on;
% plot(timE,MagE,'r','LineWidth',2);
% set(gca,'fontsize',20,'fontweight','b');
% xlabel('time(\omega_{pe})');
% ylabel('Energy');
% hold on;
% 
% plot(timE,ElecE,'g','LineWidth',2);
% set(gca,'fontsize',20,'fontweight','b');
% xlabel('time(\omega_{pe})');
% 
% hold on;
% plot(timE,TotE,'m','LineWidth',2);
% set(gca,'fontsize',20,'fontweight','b');
% xlabel('time(\omega_{pe})');
% 
% hold off;
% 
% figure(2)
% 
% 
% plot(timE,MagE,'r','LineWidth',2);
% set(gca,'fontsize',20,'fontweight','b');
% xlabel('time(\omega_{pe})');
% ylabel('Energy');
% hold on;
% 
% plot(timE,ElecE,'g','LineWidth',2);
% set(gca,'fontsize',20,'fontweight','b');
% xlabel('time(\omega_{pe})');
% 
% 
% hold off;
% 
% 
%        
% 
% 
       
       
