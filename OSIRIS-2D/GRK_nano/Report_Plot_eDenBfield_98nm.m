clc;
clear all;
close all;

addpath('UtilFun');

ne=1.1e22; % cm^-3
frm=0:150;

pth0='/home/testuser/GRK_circle_98nm';
pth1='/MS/DENSITY/Elec1/charge/';
pth2='/MS/DENSITY/Elec2/charge/';
pth3='/MS/DENSITY/Elec3/charge/';
pth4='/MS/DENSITY/Elec4/charge/';
pth5='/MS/DENSITY/Elec5/charge/';
pth6='/MS/DENSITY/Elec6/charge/';
pth7='/MS/DENSITY/Elec7/charge/';
pth8='/MS/DENSITY/Elec8/charge/';
pth9='/MS/DENSITY/Elec9/charge/';
pth10='/MS/DENSITY/Elec10/charge/';
pth11='/MS/DENSITY/Elec11/charge/';
pth12='/MS/DENSITY/Elec12/charge/';
pth13='/MS/DENSITY/Elec13/charge/';
pth14='/MS/DENSITY/Elec14/charge/';
pth15='/MS/DENSITY/Elec15/charge/';
pth16='/MS/DENSITY/Elec16/charge/';
pth17='/MS/DENSITY/Elec17/charge/';
pth18='/MS/DENSITY/Elec18/charge/';
pth19='/MS/DENSITY/Elec19/charge/';
pth20='/MS/DENSITY/Elec20/charge/';
pth21='/MS/DENSITY/Elec21/charge/';
pth22='/MS/DENSITY/Elec22/charge/';
pth23='/MS/DENSITY/Elec23/charge/';
pth24='/MS/DENSITY/Elec24/charge/';
pth25='/MS/DENSITY/Elec25/charge/';
pth26='/MS/DENSITY/Elec26/charge/';
pth27='/MS/DENSITY/Elec27/charge/';
pth28='/MS/DENSITY/Elec28/charge/';
pth29='/MS/DENSITY/Elec29/charge/';
pth30='/MS/DENSITY/Elec30/charge/';
pth31='/MS/DENSITY/Elec31/charge/';
pth32='/MS/DENSITY/Elec32/charge/';
pth33='/MS/DENSITY/Elec33/charge/';
pth34='/MS/DENSITY/Elec34/charge/';
pth35='/MS/DENSITY/Elec35/charge/';
pth36='/MS/DENSITY/Elec36/charge/';
pth37='/MS/DENSITY/Elec37/charge/';
pth38='/MS/DENSITY/Elec38/charge/';
pth39='/MS/DENSITY/Elec39/charge/';
pth40='/MS/DENSITY/Elec40/charge/';
pth41='/MS/DENSITY/Elec41/charge/';
pth42='/MS/DENSITY/Elec42/charge/';
pth43='/MS/DENSITY/Elec43/charge/';
pth44='/MS/DENSITY/Elec44/charge/';
pth45='/MS/DENSITY/Elec45/charge/';

pth47='/MS/FLD/b3/';

str_h1='charge-Elec1-';
str_h2='charge-Elec2-';
str_h3='charge-Elec3-';
str_h4='charge-Elec4-';
str_h5='charge-Elec5-';
str_h6='charge-Elec6-';
str_h7='charge-Elec7-';
str_h8='charge-Elec8-';
str_h9='charge-Elec9-';
str_h10='charge-Elec10-';
str_h11='charge-Elec11-';
str_h12='charge-Elec12-';
str_h13='charge-Elec13-';
str_h14='charge-Elec14-';
str_h15='charge-Elec15-';
str_h16='charge-Elec16-';
str_h17='charge-Elec17-';
str_h18='charge-Elec18-';
str_h19='charge-Elec19-';
str_h20='charge-Elec20-';
str_h21='charge-Elec21-';
str_h22='charge-Elec22-';
str_h23='charge-Elec23-';
str_h24='charge-Elec24-';
str_h25='charge-Elec25-';
str_h26='charge-Elec26-';
str_h27='charge-Elec27-';
str_h28='charge-Elec28-';
str_h29='charge-Elec29-';
str_h30='charge-Elec30-';
str_h31='charge-Elec31-';
str_h32='charge-Elec32-';
str_h33='charge-Elec33-';
str_h34='charge-Elec34-';
str_h35='charge-Elec35-';
str_h36='charge-Elec36-';
str_h37='charge-Elec37-';
str_h38='charge-Elec38-';
str_h39='charge-Elec39-';
str_h40='charge-Elec40-';
str_h41='charge-Elec41-';
str_h42='charge-Elec42-';
str_h43='charge-Elec43-';
str_h44='charge-Elec44-';
str_h45='charge-Elec45-';

str_h47='b3-';
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
%============AVI-FILE-NAME==========================
writerObj = VideoWriter ('eden-b-field_circle_101nm.avi');
writerObj.FrameRate=5;
open(writerObj);
% aviobj = avifile ( 'eden-b-field', 'fps',5); 
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
    fl_nm30=strcat(pth0,pth30,str_h30,str_num,str_ext);
    fl_nm31=strcat(pth0,pth31,str_h31,str_num,str_ext);
    fl_nm32=strcat(pth0,pth32,str_h32,str_num,str_ext);
    fl_nm33=strcat(pth0,pth33,str_h33,str_num,str_ext);
    fl_nm34=strcat(pth0,pth34,str_h34,str_num,str_ext);
    fl_nm35=strcat(pth0,pth35,str_h35,str_num,str_ext);
    fl_nm36=strcat(pth0,pth36,str_h36,str_num,str_ext);
    fl_nm37=strcat(pth0,pth37,str_h37,str_num,str_ext);
    fl_nm38=strcat(pth0,pth38,str_h38,str_num,str_ext);
    fl_nm39=strcat(pth0,pth39,str_h39,str_num,str_ext);
    fl_nm40=strcat(pth0,pth40,str_h40,str_num,str_ext);
    fl_nm41=strcat(pth0,pth41,str_h41,str_num,str_ext);
    fl_nm42=strcat(pth0,pth42,str_h42,str_num,str_ext);
    fl_nm43=strcat(pth0,pth43,str_h43,str_num,str_ext);
    fl_nm44=strcat(pth0,pth44,str_h44,str_num,str_ext);
    fl_nm45=strcat(pth0,pth45,str_h45,str_num,str_ext);
   
    
    fl_nm47=strcat(pth0,pth47,str_h47,str_num,str_ext);
    
    [xg_qe,yg_qe,dset_qe1,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm1);
    [xg_qe,yg_qe,dset_qe2,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm2);
    [xg_qe,yg_qe,dset_qe3,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm3);
    [xg_qe,yg_qe,dset_qe4,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm4);
    [xg_qe,yg_qe,dset_qe5,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm5);
     [xg_qe,yg_qe,dset_qe6,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm6);
    [xg_qe,yg_qe,dset_qe7,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm7);
    [xg_qe,yg_qe,dset_qe8,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm8);
    [xg_qe,yg_qe,dset_qe9,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm9);
    [xg_qe,yg_qe,dset_qe10,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm10);
     [xg_qe,yg_qe,dset_qe11,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm11);
    [xg_qe,yg_qe,dset_qe12,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm12);
    [xg_qe,yg_qe,dset_qe13,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm13);
    [xg_qe,yg_qe,dset_qe14,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm14);
    [xg_qe,yg_qe,dset_qe15,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm15);
    [xg_qe,yg_qe,dset_qe16,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm16);
    [xg_qe,yg_qe,dset_qe17,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm17);
    [xg_qe,yg_qe,dset_qe18,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm18);
    [xg_qe,yg_qe,dset_qe19,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm19);
    [xg_qe,yg_qe,dset_qe20,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm20);
     [xg_qe,yg_qe,dset_qe21,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm21);
    [xg_qe,yg_qe,dset_qe22,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm22);
    [xg_qe,yg_qe,dset_qe23,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm23);
    [xg_qe,yg_qe,dset_qe24,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm24);
    [xg_qe,yg_qe,dset_qe25,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm25);
    [xg_qe,yg_qe,dset_qe26,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm26);
    [xg_qe,yg_qe,dset_qe27,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm27);
    [xg_qe,yg_qe,dset_qe28,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm28);
    [xg_qe,yg_qe,dset_qe29,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm29);
    [xg_qe,yg_qe,dset_qe30,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm30);
    [xg_qe,yg_qe,dset_qe31,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm31);
    [xg_qe,yg_qe,dset_qe32,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm32);
    [xg_qe,yg_qe,dset_qe33,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm33);
    [xg_qe,yg_qe,dset_qe34,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm34);
     [xg_qe,yg_qe,dset_qe35,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm35);
    [xg_qe,yg_qe,dset_qe36,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm36);
    [xg_qe,yg_qe,dset_qe37,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm37);
    [xg_qe,yg_qe,dset_qe38,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm38);
    [xg_qe,yg_qe,dset_qe39,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm39);
     [xg_qe,yg_qe,dset_qe40,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm40);
    [xg_qe,yg_qe,dset_qe41,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm41);
    [xg_qe,yg_qe,dset_qe42,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm42);
    [xg_qe,yg_qe,dset_qe43,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm43);
    [xg_qe,yg_qe,dset_qe44,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm44);
    [xg_qe,yg_qe,dset_qe45,x1lt_qe,x2lt_qe,time_qe,]=...
        ReadECharge13Aug2014(fl_nm45);
    
    
%     
    [xg_b3,yg_b3,dset_b3,x1lt_b3,x2lt_b3,time_b3]=...
        ReadB313Aug2014(fl_nm47);
    
    
    dset_qe=dset_qe1+dset_qe2+dset_qe3+dset_qe4+dset_qe5+dset_qe6+dset_qe7+dset_qe8+...
            dset_qe9+dset_qe10+dset_qe11+dset_qe12+dset_qe13+dset_qe14+dset_qe15+dset_qe16+dset_qe17+dset_qe18+...
            dset_qe19+dset_qe20+dset_qe21+dset_qe22+dset_qe23+dset_qe24+dset_qe25+dset_qe26+dset_qe27+dset_qe28+ ...
            dset_qe29+dset_qe30+dset_qe31+dset_qe32+dset_qe33+dset_qe34+dset_qe35+dset_qe36+dset_qe37+...
            dset_qe38+dset_qe39+dset_qe40+dset_qe41+dset_qe42+dset_qe43+dset_qe44+dset_qe45;
%     dset_qe=dset_qe1;
    x1minb3=x1lt_b3(1); x1maxb3=x1lt_b3(2);
    x2minb3=x2lt_b3(1); x2maxb3=x2lt_b3(2);
    
    dx_qe=abs(xg_qe(2)-xg_qe(1));
    x1minqe=x1lt_qe(1); x1maxqe=x1lt_qe(2);
    x2minqe=x2lt_qe(1); x2maxqe=x2lt_qe(2);
    [Ngx_LOqe,Ngy_LOqe]=size(dset_qe);
    Ngyby2qe=round(Ngy_LOqe/2);
    LOqe_dat=dset_qe(:,Ngyby2qe);
    
    %{
    fs=1/dx_qe; %---->sampling frequency or sampling rate
    nyq=fs/2; %---->Nyquist frequency
    fft_LOqe=fft(LOqe_dat);
    Nfft=length(fft_LOqe);
    df=fs/Nfft;
    fvec=(0:(Nfft-1))*df;
    fvec(fvec>nyq) = fvec(fvec>nyq)-fs;
    freqs=fvec(fvec>=0);
    ampSpec=abs(fft_LOqe(fvec>=0));
    ampSpec=ampSpec/Nfft;
    ampSpec(2:ceil(Nfft/2))=2*ampSpec(2:ceil(Nfft/2));
    %}
    ax1=subplot(121);
    imagesc(dset_qe','XData',[x1minqe x1maxqe],...
        'YData',[x2minqe x2maxqe]);
    set(gca,'YDir','normal');

    axis square;
    set(gca,'FontSize',20,'FontWeight','Bold');
    shading('interp');
    colormap(ax1,'parula');
    xlabel('x1(c/\omega_{pe})'); ylabel('x2(c/\omega_{pe})');
   xlim([199 204]);
    
    cb = colorbar('south'); 
    set(cb,'position',[.13 .35 .335 .02]);
    title(strcat('e-Density:',num2str(time_qe)));
    
    
    ax2= subplot(122);
    imagesc(dset_b3','XData',[x1minb3 x1maxb3],...
        'YData',[x2minb3 x2maxb3]);
    set(gca,'YDir','normal');

    axis square;
    set(gca,'FontSize',20,'FontWeight','Bold');
    colormap(ax2,'jet'); shading('interp');
    xlabel('x1(c/\omega_{pe})'); %ylabel('x2(c/\omega_{pe})');
%     xlim([0 25])
    cb = colorbar('south'); 
    set(cb,'position',[.57 .35 .335 .02]);
    title(strcat('B-field:',num2str(time_qe)));
    
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