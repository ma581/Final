function [ output_args ] = mother_auto_seg_plus_copy( subwells,num_col ,cur_folder,dir_source)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
cd(cur_folder);
[anzahl_pads, namef]=getting_prefix_func;
num_names=length(namef);
im_per_pad=subwells*5*num_col;
pad_num=0;
D_old=0;
while pad_num~=anzahl_pads
    copy_files_subfunc(anzahl_pads,subwells,num_col,dir_source,cur_folder);
    cd(cur_folder);
    D = dir('awell*');
    if floor(length(D)/im_per_pad)~=pad_num && D_old ~=length(D)&&~isempty(D);
        %gets pad number
        pad=floor(length(D)/im_per_pad);
        for i=pad_num+1:pad
            %sets doit to current pad
            doit=zeros(num_names,1);
            doit(i)=1;
            %run redseg
            autoseg_func_auto( pad,subwells, namef,doit);
        end
            %upfdat
            D_old=length(D);
            pad_num=pad;
    end
        

end

