function [ pads, namef ] = getting_prefix
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%finds log file
D=dir('*.nd');

%imports data
fileID = fopen(D(1).name);
data=textscan(fileID,'%s %s','HeaderLines',11,'Delimiter',',');

%find rows with pad names
rowname=strfind(data{1},'"Stage');

posname={};
for i=1:length(rowname)
    if rowname{i}==1
        leng=length(posname);
        posname{leng+1}=data{2}{i};
    end
end


%getting # of pads
pads=0;
for i=1:length(posname)
    var=strfind(posname{i},'_1"');
      if var>0
        pads=pads+1;
      end
end

%getting pad names
subwells=length(posname)/pads;

for i=1:pads
    r=strfind(posname{i*subwells},'_');
    namef{i}=posname{i*subwells}(2:r(length(r))-1);
end

for i=1:pads
    r=strfind(namef{i},'%');
    if length(r)>0
        namef{i}=[namef{i}(1:r(1)-1),namef{i}(r(1)+1:end)];
    end
end
disp(namef);


end

