function varargout = displayImage(varargin)
% module for displaying 3D images and selecting ROI voxels
% ---------------------------------------------------------- 
% Written by AnithaPriya Krishnan 
% Version 1.0 
% contact: Dr.Danny JJ Wang "jjwang@loni.ucla.edu"; 
% Release = 20130402
% 1) Default configuration: Needs a 3D image as input for displaying the
% images. 
% 2) Can be used for selecting the ROI voxels by passing a 3D image and the
% number of ROI voxels 

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @displayImage_OpeningFcn, ...
                   'gui_OutputFcn',  @displayImage_OutputFcn, ...
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


% --- Executes just before displayImage is made visible.
function displayImage_OpeningFcn(hObject, eventdata, handles, varargin)
global opRCSL;
handles.output = hObject;
img = cell2mat(varargin(1));
nSlices = size(img,3);
axes(handles.ax_display);
imshow(img(:,:,ceil(nSlices/2)),[]);
set(handles.ax_display,'Visible', 'off');
set(handles.edt_slNumber,'String', num2str(ceil(nSlices/2))); 
set(handles.sl_display,'Value',0.5);
set(handles.pb_selectROI,'visible','off'); 
opRCSL = 0;
handles.img = img;
handles.nSlices = nSlices;
if (length(varargin)>1)
    clBlue = [11/255 132/255 199/255];
    set(handles.pb_selectROI,'visible','on');
    setbgcolor(handles.pb_selectROI,clBlue);
    handles.nPoints = cell2mat(varargin(2));
end
guidata(hObject,handles);
% UIWAIT makes displayImage wait for user response (see UIRESUME)
if (length(varargin)>1)
    uiwait(handles.fig_display);
end
% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = displayImage_OutputFcn(hObject, eventdata, handles) 
global opRCSl;
%if isfield(handles,'voxROI')
%    varargout{1} = handles.voxROI;
if (opRCSl ~= 0)
    varargout{1} = opRCSl;
    guidata(hObject, handles);
    delete(handles.fig_display);
else
    varargout{1} = handles.output;
end

% --- Executes on slider movement.
function sl_display_Callback(hObject, eventdata, handles)
sliderValue = get(handles.sl_display,'Value');
%puts the slider value into the edit text component
val = ceil(sliderValue*(handles.nSlices-1))+1;
set(handles.edt_slNumber,'String', num2str(val)); 
axes(handles.ax_display);
imshow(handles.img(:,:,val),[]);
set(handles.ax_display,'Visible', 'off');

% --- Executes during object creation, after setting all properties.
function sl_display_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function edt_slNumber_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edt_slNumber_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_selectROI.
function pb_selectROI_Callback(hObject, eventdata, handles)
global opRCSl;
opRCSl = zeros(handles.nPoints,3);
for i = 1:handles.nPoints
%     set(handles.pb_selectROI,'Enable','off');
    slNum = str2num(get(handles.edt_slNumber,'string'));
    [x y] = ginput(1);
    row = int32(y);
    col = int32(x);
    opRCSl(i,:) = [row col slNum];
    disp(['voxel: ', num2str([row col slNum])]);
    pause(1)
end
handles.voxROI = opRCSl;
guidata(hObject,handles);
uiresume(handles.fig_display);
