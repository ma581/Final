numWells = 02; %input('Enter number of wells visited: ');
numSubwells = 15; %input('Enter number of snapshots within each well: ');
inputFile = 'awell-y-'; %input('Enter original file header, surrounded by single quotes: ');
outputFile = 'sigMsigBEtho'; %input('Enter new file name, surrounded by single quotes: ');
cpath=pwd;
foldername='Renamed Files';
cpathnew=[cpath '\' foldername];
wellnames= ['sigM','sigM_Etho'];

% Check if foldername already exists
directorynames=dir;
foldersize=size(directorynames);
for  i= 1:foldersize(1)
    test = strcmp(directorynames(i).name,foldername);
    hit=directorynames(i).name;
    if test==1
        hit;
        break
    end
end

%creats folder for new files
if foldersize(1)==i
    mkdir(foldername);
end


%testing variables

if length(wellnames) ~= numWells
    disp(['You have to name all Wells'])
    return
end


channel = 'y'; %input('Enter the channel color (eg y, c, r) surrounded by single quotes: ');

for i = 1 : numWells,
    for j = 1 : numSubwells,
        cd(cpath);
        disp(['Reading input file = ' inputFile num2str(str3((i-1)*numSubwells + (j))) '.tif'])
        im = imread([inputFile num2str(str3((i-1)*numSubwells + (j))) '.tif']);
        cd(cpathnew);
        imwrite(im,[outputFile num2str(str3(i)) '-' str2(j) '-' channel '.tif'],'tiff');
    end
end
  
cd(cpath);
inputFile = 'awell-t-'; %input('Enter original file header, surrounded by single quotes: ');

channel = 'r'; %input('Enter the channel color (eg y, c, r) surrounded by single quotes: ');

for i = 1 : numWells,
    for j = 1 : numSubwells,
        disp(['Reading input file = ' inputFile num2str(str3((i-1)*numSubwells + (j))) '.tif'])
        im = imread([inputFile num2str(str3((i-1)*numSubwells + (j))) '.tif']);
        imwrite(im,[outputFile num2str(str3(i)) '-' str2(j) '-' channel '.tif'],'tiff');
    end
end

inputFile = 'awell-p-'; %input('Enter original file header, surrounded by single quotes: ');

channel = 'p'; %input('Enter the channel color (eg y, c, r) surrounded by single quotes: ');

for i = 1 : numWells,
    for j = 1 : numSubwells,
        disp(['Reading input file = ' inputFile num2str(str3((i-1)*numSubwells + (j))) '.tif'])
        im = imread([inputFile num2str(str3((i-1)*numSubwells + (j))) '.tif']);
        imwrite(im,[outputFile num2str(str3(i)) '-' str2(j) '-' channel '.tif'],'tiff');
    end
end


inputFile = 'awell-c-'; %input('Enter original file header, surrounded by single quotes: ');

channel = 'c'; %input('Enter the channel color (eg y, c, r) surrounded by single quotes: ');

for i = 1 : numWells,
    for j = 1 : numSubwells,
        disp(['Reading input file = ' inputFile num2str(str3((i-1)*numSubwells + (j))) '.tif'])
        im = imread([inputFile num2str(str3((i-1)*numSubwells + (j))) '.tif']);
        imwrite(im,[outputFile num2str(str3(i)) '-' str2(j) '-' channel '.tif'],'tiff');
    end
end

inputFile = 'awell-g-'; 
channel = 'g'; %input('Enter the channel color (eg y, c, r) surrounded by single quotes: ');

for i = 1 : numWells,
    for j = 1 : numSubwells,
        disp(['Reading input file = ' inputFile num2str(str3((i-1)*numSubwells + (j))) '.tif'])
        im = imread([inputFile num2str(str3((i-1)*numSubwells + (j))) '.tif']);
        imwrite(im,[outputFile num2str(str3(i)) '-' str2(j) '-' channel '.tif'],'tiff');
    end
end
%     for j = 1 : numFrames
%         im = imread([inputFile num2str(j) '.tif']);
%         imwrite(im,[outputFile str3(j) '.tif'],'tiff');
%         j
%     end
