function [ output_args ] = infocus_func( pad,subwells)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
D = dir('*r.tif');
col='cyg';
delcol='rp'; % color which will be deleted
for i=1:length(col)
    a=strrep(D(1).name, '-r.', ['-',col(i),'.']);
    if exist(a, 'file') == 2
        delcol(length(delcol)+1)=col(i);
    end
end

v=[];
mkdir('Deleted Images');
currentfolder=pwd;
delfolder=[currentfolder,'\Deleted Images'];
fileID = fopen('names.txt','w');
re_start=(pad-1)*subwells*5+1;
for i=re_start:size(D)
    v=focus(D(i).name);
    disp(num2str(i))
    if v<90
        fprintf(fileID,'%12s\n',D(i).name);
        for j=1:length(delcol)
            file=strrep(D(i).name, '-r.', ['-',delcol(j),'.']);
            movefile(file,delfolder);
        end     
    end
end
fclose(fileID);

end

