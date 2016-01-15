function [ output_args ] = autosnap_beta1(pads, positions, padnames, colors, outputfolder)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% inut parameters:
% pads: number of pads
% positions: numbero of positions per pads
% padnames: name of pads
% colors: in the form 'ryp'
% outputfolder: name of output folder

handles.outputfolder= outputfolder;
mkdir(handles.outputfolder);
handles.current_directory=cd;
for i=1:pads
    cd([handles.current_directory]);
    for j=1:positions
        for k=1:length(colors)
           if colors(k)=='r'
              inputFile = ['awell-','t','-']; %input('Enter original file header, surrounded by single quotes: ');
           else
              inputFile = ['awell-',colors(k),'-']; %input('Enter original file header, surrounded by single quotes: ');
           end
              
           channel = colors(k); %input('Enter the channel color (eg y, c, r) surrounded by single quotes: ');
           
           if exist([handles.current_directory,'\',inputFile num2str(str3((i-1)*positions + (j))) '.tif'],'file')==0
               timepast=0;
               while exist([handles.current_directory,'\',inputFile num2str(str3((i-1)*positions + (j))) '.tif'],'file')==0
                   pause(1);
                   timepast=timepast+1;
                   if timepast>10
                       disp('Renaming Stopped');
                       return
                   end
               end
           end
                disp(['Reading input file = ' inputFile num2str(str3((i-1)*positions + (j))) '.tif'])
                im = imread([inputFile num2str(str3((i-1)*positions + (j))) '.tif']);
                imwrite(im,[handles.current_directory, '\',handles.outputfolder,'\', padnames{i} '_' num2str(str3(i)) '-' str2(j) '-' channel '.tif'],'tiff');
        end
    end
    cd([handles.current_directory, '\',handles.outputfolder]);
    fname=[padnames{i} '_' num2str(str3(i))];
    redseg_better_del_no_cells_auto(fname, 'naturalback',1);
    num_col=length(colors);
    analysesnaphist_auto(fname,colors);
end
end


