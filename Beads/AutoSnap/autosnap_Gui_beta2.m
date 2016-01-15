function varargout = autosnap_Gui_beta2(varargin)
% AUTOSNAP_GUI_BETA2 MATLAB code for autosnap_Gui_beta2.fig
%      AUTOSNAP_GUI_BETA2, by itself, creates a new AUTOSNAP_GUI_BETA2 or raises the existing
%      singleton*.
%
%      H = AUTOSNAP_GUI_BETA2 returns the handle to a new AUTOSNAP_GUI_BETA2 or the handle to
%      the existing singleton*.
%
%      AUTOSNAP_GUI_BETA2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AUTOSNAP_GUI_BETA2.M with the given input arguments.
%
%      AUTOSNAP_GUI_BETA2('Property','Value',...) creates a new AUTOSNAP_GUI_BETA2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before autosnap_Gui_beta2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to autosnap_Gui_beta2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help autosnap_Gui_beta2

% Last Modified by GUIDE v2.5 23-Jul-2014 16:44:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @autosnap_Gui_beta2_OpeningFcn, ...
                   'gui_OutputFcn',  @autosnap_Gui_beta2_OutputFcn, ...
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


% --- Executes just before autosnap_Gui_beta2 is made visible.
function autosnap_Gui_beta2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to autosnap_Gui_beta2 (see VARARGIN)

% Choose default command line output for autosnap_Gui_beta2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes autosnap_Gui_beta2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = autosnap_Gui_beta2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in start_push.
function start_push_Callback(hObject, eventdata, handles)
% hObject    handle to start_push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cd(get(handles.path,'String'));
handles.number_of_pos=get(handles.Pos_num,'String');%get number of positions per pad
handles.number_of_pos=str2double(handles.number_of_pos);
%number of colors
handles.colors=1+get(handles.button_RFP, 'Value') + get(handles.button_YFP, 'Value')+ get(handles.button_CFP, 'Value')+ get(handles.button_GFP, 'Value');
[handles.posname, handles.pads, handels.posname_title]=pads_stg;
if handles.pads==0
    return;
end
handles.cstring=color_string(handles);
if handles.pads==0
    errordlg('Save stage positions','Save the stage positions!');
else
    autosnap_beta_error_send_4(handles.pads,handles.number_of_pos, handles.posname, handles.cstring, get(handles.output_folder, 'String'), handels.posname_title, get(handles.auto_print,'Value'))
end

%Copying Files
if get(handles.auto_move,'Value')==1
    handles.folder_date=get_date;
    handles.move_files=get(handles.move_path,'String');
    handles.foldername=[handles.move_files,'\',handles.folder_date];
    if exist(handles.foldername,'file')==0
        mkdir(handles.foldername);
    end
    cd(handles.foldername);
    D=dir;
    handles.foldername_take=[handles.foldername,'\Take',num2str(length(D)-1)];%error?
    mkdir(handles.foldername_take);
    cd(get(handles.path,'String'));
    D=dir('*tif');
    set(handles.Done_Copying,'String','@work');
    pause(0.1);
    for i=1:length(D)
        copyfile([get(handles.path,'String'),'\',D(i).name],handles.foldername_take);
    end
    set(handles.Done_Copying,'String','33');
    pause(0.1);
    cd(get(handles.path,'String'));
    D=dir('*stg');
    for i=1:length(D)
        copyfile([get(handles.path,'String'),'\',D(i).name],handles.foldername_take);
    end
    set(handles.Done_Copying,'String','66');
    pause(0.1);
    copyfile([get(handles.path,'String'),'\',get(handles.output_folder,'String')],[handles.foldername_take,'\',get(handles.output_folder,'String')]);
    set(handles.Done_Copying,'String','100');
end
    
        
    
    




function done_Callback(hObject, eventdata, handles)
% hObject    handle to done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of done as text
%        str2double(get(hObject,'String')) returns contents of done as a double


% --- Executes during object creation, after setting all properties.
function done_CreateFcn(hObject, eventdata, handles)
% hObject    handle to done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Pos_num_Callback(hObject, eventdata, handles)
% hObject    handle to Pos_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Pos_num as text
%        str2double(get(hObject,'String')) returns contents of Pos_num as a double


% --- Executes during object creation, after setting all properties.
function Pos_num_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pos_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function path_Callback(hObject, eventdata, handles)
% hObject    handle to path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of path as text
%        str2double(get(hObject,'String')) returns contents of path as a double


% --- Executes during object creation, after setting all properties.
function path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function send_error
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Defining Account details
myaddress = 'nikon.bacillus.1@gmail.com';
mypassword = 'P455w0rd!';

%Setting Prefences
setpref('Internet','E_mail',myaddress);
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username',myaddress);
setpref('Internet','SMTP_Password',mypassword);

%Setting Properties
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');

%sending Email
sendmail(myaddress, 'Acquisition Error', 'Come and check me. I am having problems!');


% --- Executes on button press in button_CFP.
function button_CFP_Callback(hObject, eventdata, handles)
% hObject    handle to button_CFP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of button_CFP


% --- Executes on button press in button_YFP.
function button_YFP_Callback(hObject, eventdata, handles)
% hObject    handle to button_YFP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of button_YFP


% --- Executes on button press in button_GFP.
function button_GFP_Callback(hObject, eventdata, handles)
% hObject    handle to button_GFP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of button_GFP


% --- Executes on button press in button_RFP.
function button_RFP_Callback(hObject, eventdata, handles)
% hObject    handle to button_RFP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of button_RFP

function  [namef, pads, namef_title]= pads_stg
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%finds log file
D=dir('*.stg');
s=size(D);
if s(1)==0
    pads=0;
    namef=0;
    namef_title=0;
    disp('Stopped');
    errordlg('Save Stage Positions','Save Stage Positions!');
else
    %imports data
    fileID = fopen(D(1).name);
    data=textscan(fileID,'%s %s %s %s %s %s %s %s %s %s %s %s %s','HeaderLines',4,'Delimiter',',');
    posname=data{1}; % saving pads names

    %number of pads =  length posname
    pads=length(posname);
    
    %getting pad name
    for i=1:pads
        if posname{i}(end-1)=='1';
            r=strfind(posname{i},'_');
            namef{i}=posname{i}(2:r(length(r))-1);
        else
            namef{i}=posname{i}(2:end-1);
        end
    end
    namef_title=namef;
    for i=1:pads
        r=strfind(namef{i},'%');
        if length(r)>0
            namef{i}=[namef{i}(1:r(1)-1),namef{i}(r(1)+1:end)];
        end
    end
    for i=1:pads
        r=strfind(namef{i},'\');
        if length(r)>0
            namef{i}=[namef{i}(1:r(1)-1),namef{i}(r(1)+1:end)];
        end
    end
    for i=1:pads
        r=strfind(namef{i},'/');
        if length(r)>0
            namef{i}=[namef{i}(1:r(1)-1),namef{i}(r(1)+1:end)];
        end
    end
    
    
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function cstring=color_string(handles)
cstring='p';
if get(handles.button_YFP, 'Value')==1
    cstring=[cstring, 'y'];
end
if get(handles.button_CFP, 'Value')==1
    cstring=[cstring, 'c'];
end
if get(handles.button_GFP, 'Value')==1
    cstring=[cstring,'g'];
end
if get(handles.button_RFP, 'Value')==1
    cstring=[cstring,'r'];
end



function output_folder_Callback(hObject, eventdata, handles)
% hObject    handle to output_folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of output_folder as text
%        str2double(get(hObject,'String')) returns contents of output_folder as a double


% --- Executes during object creation, after setting all properties.
function output_folder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to output_folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in auto_move.
function auto_move_Callback(hObject, eventdata, handles)
% hObject    handle to auto_move (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of auto_move



function move_path_Callback(hObject, eventdata, handles)
% hObject    handle to move_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of move_path as text
%        str2double(get(hObject,'String')) returns contents of move_path as a double


% --- Executes during object creation, after setting all properties.
function move_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to move_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function folder_name=get_date
d1=date;
d_pos=strfind(d1,'-');
tag=d1(1:d_pos(1)-1);
monat=d1(d_pos(1)+1:d_pos(2)-1);
monate_jahr={'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
for i=1:12
    if sum(monate_jahr{i}==monat)==3
        if i<10
            monate_num=['0',num2str(i)];
        else
            monate_num=num2str(i);
        end
        break
    end
end
jahr=d1(d_pos(2)+1:end);
folder_name=[jahr,'_',monate_num,'_',tag];


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.folder_date=get_date;
handles.move_files=get(handles.move_path,'String');
handles.foldername=[handles.move_files,'\',handles.folder_date];
if exist(handles.foldername,'file')==0
    mkdir(handles.foldername);
end
cd(handles.foldername);
D=dir;
handles.foldername_take=[handles.foldername,'\Take',num2str(length(D)-1)];%error?
mkdir(handles.foldername_take);
cd(get(handles.path,'String'));
D=dir('*tif');
set(handles.Done_Copying,'String','@work');
pause(0.1);
for i=1:length(D)
    copyfile([get(handles.path,'String'),'\',D(i).name],handles.foldername_take);
end
set(handles.Done_Copying,'String','33');
pause(0.1);
cd(get(handles.path,'String'));
D=dir('*stg');
for i=1:length(D)
    copyfile([get(handles.path,'String'),'\',D(i).name],handles.foldername_take);
end
set(handles.Done_Copying,'String','66');
pause(0.1);
copyfile([get(handles.path,'String'),'\',get(handles.output_folder,'String')],[handles.foldername_take,'\',get(handles.output_folder,'String')]);
set(handles.Done_Copying,'String','100');
    



function Done_Copying_Callback(hObject, eventdata, handles)
% hObject    handle to Done_Copying (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Done_Copying as text
%        str2double(get(hObject,'String')) returns contents of Done_Copying as a double


% --- Executes during object creation, after setting all properties.
function Done_Copying_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Done_Copying (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in auto_print.
function auto_print_Callback(hObject, eventdata, handles)
% hObject    handle to auto_print (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of auto_print


% --- Executes on button press in print_now.
function print_now_Callback(hObject, eventdata, handles)
% hObject    handle to print_now (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cd([get(handles.path,'String'),'\',get(handles.output_folder,'String')]);
out_names=dir('output*.fig');
for i=1:length(out_names)
    h=hgload(out_names(i).name);
    set(gcf,'Units','centimeters');
    set(gcf,'PaperOrientation','landscape');
    set(gcf,'PaperType','A4');
    set(gcf, 'PaperPosition', [0.2 0.5 29 19.5 ]);
    print(h,'-PSouth Wing Colour Laserjet')
end
