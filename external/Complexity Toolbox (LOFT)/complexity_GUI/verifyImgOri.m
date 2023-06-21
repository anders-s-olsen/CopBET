function varargout = verifyImgOri(varargin)
% module for verifying the orientation of the fMRI images and brain mask
% ---------------------------------------------------------- 
% Written by AnithaPriya Krishnan 
% Version 1.0 
% contact: Dr.Danny JJ Wang "jjwang@loni.ucla.edu"; 
% Release = 20130402
% There might be a flip in the left right orientation of the images
% depending on the tool used to generate the brain mask. This module allows
% to modify the orientation of the mask to match that of the fMRI data

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @verifyImgOri_OpeningFcn, ...
                   'gui_OutputFcn',  @verifyImgOri_OutputFcn, ...
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


% --- Executes just before verifyImgOri is made visible.
function verifyImgOri_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
img1 = cell2mat(varargin(1));
img2 = cell2mat(varargin(2));
nSlices = size(img1,3);
axes(handles.ax_left);
imshow(img1(:,:,ceil(nSlices/2)),[]);
set(handles.ax_left,'Visible', 'off');
axes(handles.ax_right);
imshow(img2(:,:,ceil(nSlices/2)),[]);
set(handles.ax_right,'Visible', 'off');
clBlue = [11/255 132/255 199/255];
setbgcolor(handles.pb_rot90mask,clBlue);
setbgcolor(handles.pb_fliplrmask,clBlue);
setbgcolor(handles.pb_flipudmask,clBlue);
setbgcolor(handles.pb_confOrient,clBlue);
setbgcolor(handles.pb_flipz,clBlue);
handles.nSlices = nSlices;
handles.img1 = img1;
handles.img2 = img2;
guidata(hObject,handles);
% UIWAIT makes verifyImgOri wait for user response (see UIRESUME)
uiwait(handles.fig_verImgOri);
% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = verifyImgOri_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.img2;
guidata(hObject,handles);
delete(handles.fig_verImgOri);


% --- Executes on slider movement.
function sl_dispSlice_Callback(hObject, eventdata, handles)
sliderValue = get(handles.sl_dispSlice,'Value');
%puts the slider value into the edit text component
val = ceil(sliderValue*(handles.nSlices-1))+1;
set(handles.edt_dispSlice,'String', num2str(val)); 
axes(handles.ax_left);
imshow(handles.img1(:,:,val),[]);
set(handles.ax_left,'Visible', 'off');
axes(handles.ax_right);
imshow(handles.img2(:,:,val),[]);
set(handles.ax_right,'Visible', 'off');


% --- Executes during object creation, after setting all properties.
function sl_dispSlice_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function edt_dispSlice_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edt_dispSlice_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_rot90mask.
function pb_rot90mask_Callback(hObject, eventdata, handles)
for i = 1:handles.nSlices
    tmp(:,:,i) = rot90(handles.img2(:,:,i));
end
handles.img2 = tmp;
guidata(hObject,handles);
update_disp(handles);

% --- Executes on button press in pb_fliplrmask.
function pb_fliplrmask_Callback(hObject, eventdata, handles)
for i = 1:handles.nSlices
    tmp(:,:,i) = fliplr(handles.img2(:,:,i));
end
handles.img2 = tmp;
guidata(hObject,handles);
update_disp(handles);

% --- Executes on button press in pb_flipudmask.
function pb_flipudmask_Callback(hObject, eventdata, handles)
for i = 1:handles.nSlices
    tmp(:,:,i) = flipud(handles.img2(:,:,i));
end
handles.img2 = tmp;
guidata(hObject,handles);
update_disp(handles);


% --- Executes on button press in pb_flipz.
function pb_flipz_Callback(hObject, eventdata, handles)
n=0;
nSlices = handles.nSlices;
for i = 1:nSlices
    tmp(:,:,i) = handles.img2(:,:,nSlices-n);
    n = n+1;
end
handles.img2 = tmp;
guidata(hObject,handles);
update_disp(handles)

% --- Executes on button press in pb_confOrient.
function pb_confOrient_Callback(hObject, eventdata, handles)
uiresume(handles.fig_verImgOri);

function update_disp(strHandle)
axes(strHandle.ax_left);
imshow(strHandle.img1(:,:,ceil(strHandle.nSlices/2)),[]);
set(strHandle.ax_left,'Visible', 'off');
axes(strHandle.ax_right);
imshow(strHandle.img2(:,:,ceil(strHandle.nSlices/2)),[]);
set(strHandle.ax_right,'Visible', 'off');
