function [ output_args ] = autoseg_func_auto( pad,subwells, namef,doit )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
RenameMultiwellGUI_func(pad, subwells, namef, 'y');
infocus_func(pad,subwells);
redseg_better_auto_func(doit);
analysesnaphist_modified_loop_correct_title;


end

