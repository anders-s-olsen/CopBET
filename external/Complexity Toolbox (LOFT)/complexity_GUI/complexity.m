function varargout = complexity(varargin)
% toolkit for complexity analysis of resting state fMRI data 
% ----------------------------------------------------------
% Written by AnithaPriya Krishnan
% Version 1.0
% Mail to Authors:  <a href="Dawnwei.Song@gmail.com">SONG Xiao-Wei</a>; <a href="ycg.yan@gmail.com">YAN Chao-Gan</a>; <a href="dongzy08@gmail.com">DONG Zhang-Ye</a> 
% Release = 20130402
% uses "Tools for NIfTI and ANALYZE image" toolbox by Jimmy Shen from MATLAB file
% exchange for reading and writing NIFTI/ANALYZE images

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @complexity_OpeningFcn, ...
                   'gui_OutputFcn',  @complexity_OutputFcn, ...
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


% --- Executes just before complexity is made visible.
function complexity_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);
handles.banner = imread('LOFT_logo.png');
axes(handles.ax_topLeft);
image(handles.banner);
set(handles.ax_topLeft,'Visible', 'off');
handles.banner = imread('brain_fractal.png');
axes(handles.ax_topRight);
image(handles.banner);
set(handles.ax_topRight,'Visible', 'off');
clBlue = [11/255 132/255 199/255];
setbgcolor(handles.pb_preProcessing,clBlue);
setbgcolor(handles.pb_ApEn,clBlue);
setbgcolor(handles.pb_sampEn,clBlue);
setbgcolor(handles.pb_crossApEn,clBlue);
setbgcolor(handles.pb_MSE,clBlue);
setbgcolor(handles.pb_waveletMSE,clBlue);
% UIWAIT makes complexity wait for user response (see UIRESUME)
% uiwait(handles.Complexity);


% --- Outputs from this function are returned to the command line.
function varargout = complexity_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


% --- Executes on button press in pb_preProcessing.
function pb_preProcessing_Callback(hObject, eventdata, handles)
preProcessing

% --- Executes on button press in pb_ApEn.
function pb_ApEn_Callback(hObject, eventdata, handles)
complexity_ApEn

% --- Executes on button press in pb_sampEn.
function pb_sampEn_Callback(hObject, eventdata, handles)
complexity_sampEn

% --- Executes on button press in pb_crossApEn.
function pb_crossApEn_Callback(hObject, eventdata, handles)
complexity_cApEn

% --- Executes on button press in pb_MSE.
function pb_MSE_Callback(hObject, eventdata, handles)
% hObject    handle to pb_MSE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
complexity_MSE

% --- Executes on button press in pb_waveletMSE.
function pb_waveletMSE_Callback(hObject, eventdata, handles)
% hObject    handle to pb_waveletMSE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
