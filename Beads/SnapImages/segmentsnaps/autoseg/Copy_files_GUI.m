function varargout = Copy_files_GUI(varargin)
% COPY_FILES_GUI MATLAB code for Copy_files_GUI.fig
%      COPY_FILES_GUI, by itself, creates a new COPY_FILES_GUI or raises the existing
%      singleton*.
%
%      H = COPY_FILES_GUI returns the handle to a new COPY_FILES_GUI or the handle to
%      the existing singleton*.
%
%      COPY_FILES_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COPY_FILES_GUI.M with the given input arguments.
%
%      COPY_FILES_GUI('Property','Value',...) creates a new COPY_FILES_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Copy_files_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Copy_files_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Copy_files_GUI

% Last Modified by GUIDE v2.5 24-Apr-2014 11:27:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Copy_files_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @Copy_files_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Copy_files_GUI is made visible.
function Copy_files_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Copy_files_GUI (see VARARGIN)


% Choose default command line output for Copy_files_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Copy_files_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Copy_files_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function pads_Callback(hObject, eventdata, handles)
% hObject    handle to pads (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pads as text
%        str2double(get(hObject,'String')) returns contents of pads as a double


% --- Executes during object creation, after setting all properties.
function pads_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pads (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pos_Callback(hObject, eventdata, handles)
% hObject    handle to pos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pos as text
%        str2double(get(hObject,'String')) returns contents of pos as a double


% --- Executes during object creation, after setting all properties.
function pos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dir_source_Callback(hObject, eventdata, handles)
% hObject    handle to dir_source (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dir_source as text
%        str2double(get(hObject,'String')) returns contents of dir_source as a double


% --- Executes during object creation, after setting all properties.
function dir_source_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dir_source (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dir_out_Callback(hObject, eventdata, handles)
% hObject    handle to dir_out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dir_out as text
%        str2double(get(hObject,'String')) returns contents of dir_out as a double


% --- Executes during object creation, after setting all properties.
function dir_out_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dir_out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in Go.
function Go_Callback(hObject, eventdata, handles)
% hObject    handle to Go (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

copy_files( handles.pads, handles.pos,3,handles.dir_source, handles.dir_out)



function [ output_args ] = copy_files( pads, images,col,dir_source, dir_out)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if exist(dir_out,'dir')==0
    mkdir(dir_out)
end

tot_num=pads*images*col*5;
im_num=0;
cd(dir_source);
colors='ytp';
c_names_num=[0,0,0];
while im_num~=tot_num
    D = dir('awell*');
    
    
    if length(D)~=im_num
        %converting D into Cell array
        for i=1:length(D)
            names{i}=D(i).name;
        end
        
        for c=1:col;
            %getting file name of colors
            for i=1:length(D)
                if names{i}(7)==colors(c)
                    index(i)=true;
                else
                    index(i)=false;
                end
            end
            c_names=names(index);
            %copying files
            if length(c_names)~=c_names_num(c)
                for i=c_names_num(c)+1:length(c_names)
                    source=[dir_source, c_names{i}];
                    copyfile(source, dir_out);
                end
                c_names_num(c)=length(c_names);
            end
        end
        im_num=length(D);
        set(handles.copied,'string',num2str(c_names_num));
        guidata(hObject, handles);
    end
end



