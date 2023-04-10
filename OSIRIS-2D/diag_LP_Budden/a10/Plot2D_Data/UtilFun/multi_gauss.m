function [in_pred]=multi_gauss(par,t)

%function [int_pred]=multi_gauss(no_guass,in_parm);
% multi_guass.m function will provide no of guassian calculation to main
% program for fitting in nlinfit, nlinfit will converge according to this
% function and no of guassian
% no_guass= need from main program how many guassian needed to model data
% d_t= this is the data(x) you want to model 
%int_pred= this fuction will return data(y)after fitting data according to model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the model fitted is
%
%           y=a1*exp(-0.5*((x-c1)/b1)^2))+...+a4*exp(-0.5*((x-c4)/b4)^2))
%
% par is the parameter vector characterising the model
% t is the vector at which model values are evaluated by the nlinfit
% function and this t has no bearance with the t of the original data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rt=length(t);
rpar=length(par);
n_gauss=rpar/3;
in_pred=zeros(rt,1);

cnt=1;
for i=1:1:n_gauss
    a=par(cnt); b=par(cnt+1); c=par(cnt+2);
%     cst1=sqrt(b/pi); a=a*cst1;
    in_pred=in_pred+a*exp(-0.5*((((1/b)*(t-c)).^2)));
    cnt=cnt+3; 
end
return;