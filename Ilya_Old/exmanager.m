function varargout = exmanager(varargin)
% EXMANAGER M-file for exmanager.fig
%      EXMANAGER, by itself, creates a new EXMANAGER or raises the existing
%      singleton*.
%
%      H = EXMANAGER returns the handle to a new EXMANAGER or the handle to
%      the existing singleton*.
%
%      EXMANAGER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EXMANAGER.M with the given input arguments.
%
%      EXMANAGER('Property','Value',...) creates a new EXMANAGER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before exmanager_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to exmanager_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help exmanager

% Last Modified by GUIDE v2.5 10-Jul-2003 15:22:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @exmanager_OpeningFcn, ...
                   'gui_OutputFcn',  @exmanager_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before exmanager is made visible.
function exmanager_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to exmanager (see VARARGIN)

% Choose default command line output for exmanager
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes exmanager wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = exmanager_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% ==============================
% OBJECT: pathfield
% ==============================

% --- Executes during object creation, after setting all properties.
function pathfield_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pathfield (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function pathfield_Callback(hObject, eventdata, handles)
% hObject    handle to pathfield (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pathfield as text
%        str2double(get(hObject,'String')) returns contents of pathfield as a double

% ==============================
% OBJECT: browsebutton
% ==============================

% --- Executes on button press in browsebutton.
function browsebutton_Callback(hObject, eventdata, handles)
% hObject    handle to browsebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set_pathname(handles);

% --- Executes on button press in loadworkspacechk.
function loadworkspacechk_Callback(hObject, eventdata, handles)
% hObject    handle to loadworkspacechk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of loadworkspacechk


% --- Executes on button press in loadfigureschk.
function loadfigureschk_Callback(hObject, eventdata, handles)
% hObject    handle to loadfigureschk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of loadfigureschk


% --- Executes on button press in diarychk.
function diarychk_Callback(hObject, eventdata, handles)
% hObject    handle to diarychk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of diarychk
writediary = get(hObject,'Value');
pathname = pathname_check(handles);
if writediary & ~isempty(pathname);
    filename= [pathname 'diary.txt'];
    diary(filename);
else
    diary off;
end;
if isempty(pathname)
    set(hObject,'Value',0);
end;    

% --- Executes on button press in saveworkspacechk.
function saveworkspacechk_Callback(hObject, eventdata, handles)
% hObject    handle to saveworkspacechk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of saveworkspacechk


% --- Executes on button press in savefigureschk.
function savefigureschk_Callback(hObject, eventdata, handles)
% hObject    handle to savefigureschk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of savefigureschk


% --- Executes on button press in diarybutton.
function diarybutton_Callback(hObject, eventdata, handles)
% hObject    handle to diarybutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pathname = pathname_check(handles);
if ~isempty(pathname)
    filename = [pathname 'diary.txt'];
    if exist(filename,'file')
        edit(filename);
    else
        warndlg('There is no diary in this experiment folder yet.')
    end;
end;

% --- Executes on button press in loadbutton.
function loadbutton_Callback(hObject, eventdata, handles)
% hObject    handle to loadbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pathname = pathname_check(handles);

loadworkspace = get(handles.loadworkspacechk,'Value');
loadfigures = get(handles.loadfigureschk,'Value');
if ~isempty(pathname)
    if loadworkspace
        if exist([pathname 'experiment_workspace.mat'])
            evalin('base',['load(''' pathname 'experiment_workspace.mat' ''');']);
        else
            warndlg('No saved workspace found');
        end;
    end;
    if loadfigures
        figs = dir([pathname '*.fig']);
        if isempty(figs)
            warndlg('No figure found');
        end;
        for k = 1:length(figs)
            open([pathname figs(k).name]);
        end;
    end;
end;

% --- Executes on button press in savebutton.
function savebutton_Callback(hObject, eventdata, handles)
% hObject    handle to savebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pathname = pathname_check(handles);
if ~isempty(pathname)
    saveworkspace = get(handles.saveworkspacechk,'Value');
    savefigures = get(handles.savefigureschk,'Value');
    if saveworkspace
        answer = questdlg('Save workspace (may overwrite previously saved workspace)?');
        if strcmp(answer,'Yes');
            filename = [pathname 'experiment_workspace.mat'];
            evalin('base',['save(''' filename ''')']);
        end;
    end;
    if savefigures
        answer = questdlg('Save current figures (may delete previously saved figures)?');
        if strcmp(answer,'Yes');
            delete([pathname '*.fig']);
            figs = get(0,'Children');
            figs(figs==gcbf) = [];
            for k = 1:length(figs);
                figname = [pathname num2str(figs(k)) '.fig'];
                saveas(figs(k),figname);
            end;
        end;
    end;
end;

% --- Path name check function
function pathname = pathname_check(handles);
pathname = get(handles.pathfield,'String');
if isempty(pathname)
    pathname = set_pathname(handles);
end;
if ischar(pathname) & ~isempty(pathname) & ~strcmp(pathname(end),filesep)
    pathname = [pathname filesep];
else
    pathname = '';
end;

function pathname = set_pathname(handles);
pathname = get(handles.pathfield,'String');
getdirtitle = 'Select experiment folder';
if isempty(pathname)
    pathname = uigetdir(pwd,getdirtitle);
else
    pathname = uigetdir(pathname,getdirtitle);
end;
if ischar(pathname)
    set(handles.pathfield,'String',pathname);
end;