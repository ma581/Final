function [ output_args ] = mother_auto_seg( subwells,num_col ,cur_folder)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
cd(cur_folder);
[anzahl_pads, namef]=getting_prefix;
num_names=length(namef);
im_per_pad=subwells*5*num_col;
pad_num=0;
D_old=0;
while pad_num~=anzahl_pads
    D = dir('awell*');
    if mod(length(D),im_per_pad)==0 && D_old ~=length(D)&&~isempty(D);
        %gets pad number
        pad=length(D)/im_per_pad;
        %sets doit to current pad
        doit=zeros(num_names,1);
        doit(pad)=1;
        %run redseg
        autoseg_func_auto( pad,subwells, namef,doit);
        %upfdat
        D_old=length(D);
        pad_num=pad;
    end
        

end

