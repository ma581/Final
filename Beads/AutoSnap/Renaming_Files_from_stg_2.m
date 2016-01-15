%finds log file
D=dir('*.stg');

%imports data
fileID = fopen(D(1).name);
data=textscan(fileID,'%s %s %s %s %s %s %s %s %s %s %s %s %s','HeaderLines',4,'Delimiter',',');
posname=data{1}; % saving pads names

%number of pads =  length posname
pads=length(posname);
subwells=1;

%getting pad names
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

col=getcolours;
D=dir('awell-t-*');
subwells=length(D)/(pads*5);

RenameMultiwellGUI(pads, subwells, namef, 'y');