function []=WriteOsirisInput()
clc
clear all
close all

%{
simulation=struct('random_seed',0,'omega_p0',0.0,...
    'n0',0.0,'gamma',0.0);
fmt_str_sim={'%d'; '%f'; '%f' ; '%f'};
WriteStructFile('simulation',simulation,fmt_str_sim);
%}

node_conf=struct('node_number',[16;8],'if_periodic',...
    {{'.false.','.false.'}});
fmt_str_ncon={'%d'; '%s'};
WriteStructFile('node_conf',node_conf,fmt_str_ncon);

end

function []=WriteStructFile(st_nm,st,fmt_st)

fid=fopen('os-stdin','a+');
lensim=length(st);
fld_nm=fieldnames(st);
lenfnm=length(fld_nm);

for i=1:lensim
    fprintf(fid,'%s\n%s\n',st_nm,'{');
    for j=1:lenfnm
        fnm=fld_nm{j};
        fval=eval(strcat('st(',num2str(i),...
            ').',fnm));
        len_fval=length(fval);
        tmp_str=strcat(fmt_st{j},',');
        if (len_fval>1)            
            fprintf(fid,'%s=',strcat(fnm,'(',num2str(1),...
                ':',num2str(len_fval),')'));
            for k=1:len_fval                
                fprintf(fid,tmp_str,fval(k));
            end
            fprintf(fid,'\n');
        else
            fprintf(fid,'%s=',fnm);
            fprintf(fid,tmp_str,fval);
            fprintf(fid,'\n');
        end
    end
    fprintf(fid,'%s\n','}');
end
fclose(fid);

end