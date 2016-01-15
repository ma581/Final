function RenameMultiwellGUI_start_later( numWells, Subwells,outputFile , colours)

%First the global parameters for the renaming are defined (numWells, numSubwells, numSubwells)

%numWells = 4; %input('Enter number of wells visited: ');
%numSubwells = 30; %input('Enter number of snapshots within each well: ');
%outputFile =  'sigM_4EtoH'; %input('Enter new file name, surrounded by single quotes: ');

% Then the paramters for each for loop are defined. And each for Loop
% renames the files in such a way that it is obvious from which pad the
% image is.

D=dir('awell-p*');
start_image=str2num(D(1).name(9:end-4));
startwell=floor(start_image/30)+1;
numWells=startwell+numWells+1;




numSubwells=Subwells*5; %time 5 as for each position 5 images are taken

if size(strfind(colours,'y'))==1
    inputFile = 'awell-y-'; %input('Enter original file header, surrounded by single quotes: ');

    channel = 'y'; %input('Enter the channel color (eg y, c, r) surrounded by single quotes: ');

    for i = startwell : numWells,
        for j = 1 : numSubwells,
            disp(['Reading input file = ' inputFile num2str(str3((i-1)*numSubwells + (j))) '.tif'])
            im = imread([inputFile num2str(str3((i-1)*numSubwells + (j))) '.tif']);
            imwrite(im,[outputFile{i} '_' num2str(str3(i)) '-' str2(j) '-' channel '.tif'],'tiff');
        end
    end
end

inputFile = 'awell-t-'; %input('Enter original file header, surrounded by single quotes: ');

channel = 'r'; %input('Enter the channel color (eg y, c, r) surrounded by single quotes: ');

for i = startwell : numWells,
    for j = 1 : numSubwells,
        disp(['Reading input file = ' inputFile num2str(str3((i-1)*numSubwells + (j))) '.tif'])
        im = imread([inputFile num2str(str3((i-1)*numSubwells + (j))) '.tif']);
        imwrite(im,[outputFile{i} '_' num2str(str3(i)) '-' str2(j) '-' channel '.tif'],'tiff');
    end
end

inputFile = 'awell-p-'; %input('Enter original file header, surrounded by single quotes: ');

channel = 'p'; %input('Enter the channel color (eg y, c, r) surrounded by single quotes: ');

for i = startwell : numWells,
    for j = 1 : numSubwells,
        disp(['Reading input file = ' inputFile num2str(str3((i-1)*numSubwells + (j))) '.tif'])
        im = imread([inputFile num2str(str3((i-1)*numSubwells + (j))) '.tif']);
        imwrite(im,[outputFile{i} '_' num2str(str3(i)) '-' str2(j) '-' channel '.tif'],'tiff');
    end
end

if size(strfind(colours,'c'))==1
    inputFile = 'awell-c-'; %input('Enter original file header, surrounded by single quotes: ');

    channel = 'c'; %input('Enter the channel color (eg y, c, r) surrounded by single quotes: ');

    for i = 1 : numWells,
        for j = 1 : numSubwells,
            disp(['Reading input file = ' inputFile num2str(str3((i-1)*numSubwells + (j))) '.tif'])
            im = imread([inputFile num2str(str3((i-1)*numSubwells + (j))) '.tif']);
            imwrite(im,[outputFile{i} '_' num2str(str3(i)) '-' str2(j) '-' channel '.tif'],'tiff');
        end
    end
end

if size(strfind(colours,'g'))==1
    inputFile = 'awell-g-'; 
    channel = 'g'; %input('Enter the channel color (eg y, c, r) surrounded by single quotes: ');

    for i = startwell : numWells,
        for j = 1 : numSubwells,
            disp(['Reading input file = ' inputFile num2str(str3((i-1)*numSubwells + (j))) '.tif'])
            im = imread([inputFile num2str(str3((i-1)*numSubwells + (j))) '.tif']);
            imwrite(im,[outputFile{i} '_' num2str(str3(i)) '-' str2(j) '-' channel '.tif'],'tiff');
        end
    end
%     for j = 1 : numFrames
%         im = imread([inputFile num2str(j) '.tif']);
%         imwrite(im,[outputFile str3(j) '.tif'],'tiff');
%         j
%     end
end

end

