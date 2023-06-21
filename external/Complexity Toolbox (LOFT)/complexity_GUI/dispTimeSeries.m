function varargout = dispTimeSeries(varargin)
% module for displaying processed (filtered, detrended) time series data
% of resting state fMRI data
% ---------------------------------------------------------- 
% Written by AnithaPriya Krishnan 
% Version 1.0 
% contact: Dr.Danny JJ Wang "jjwang@loni.ucla.edu"; 
% Release = 20130402
% 1) Displays the original and processed image. Click on the image to plot
% the time series data
% 2) For saving the processed volume (4-D), the filenames for 
% a) filtered data is: baseName_fil_ordN_NCFMper.nii
% b) detrended data are: baseName_detrendMean.nii,
% basename_detrendLinear.nii
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dispTimeSeries_OpeningFcn, ...
                   'gui_OutputFcn',  @dispTimeSeries_OutputFcn, ...
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


% --- Executes just before dispTimeSeries is made visible.
function dispTimeSeries_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
origImg = cell2mat(varargin(1));
processedImg = cell2mat(varargin(2));
opFname = cell2mat(varargin(3));
imgVoxDim = cell2mat(varargin(4));
nSlices = size(origImg,3);
nTimePoints = size(origImg,4);
handles.nSlices = nSlices;
handles.origImg = origImg;
handles.processedImg = processedImg;
handles.nTimePoints = nTimePoints;
axes(handles.ax_Image);
h1 = imshow(squeeze(origImg(:,:,ceil(nSlices/2),1)),[]);
set(h1, 'ButtonDownFcn', {@ax_Image_ButtonDownFcn, handles});
axes(handles.ax_processedImg);
h2 = imshow(squeeze(processedImg(:,:,ceil(nSlices/2),1)),[]);
set(h2, 'ButtonDownFcn', {@ax_processedImg_ButtonDownFcn, handles});
set(handles.edt_dispSlice,'String',num2str(ceil(nSlices/2)));
set(handles.sl_dispSlice,'Value',ceil(nSlices/2)/nSlices);
set(handles.edt_dispVol,'String',num2str(1));
set(handles.sl_dispVol,'Value',1/nTimePoints);
clBlue = [11/255 132/255 199/255];
setbgcolor(handles.pb_saveImg,clBlue);
handles.imgVoxDim = imgVoxDim;
handles.opFname = opFname;
guidata(hObject,handles);
% UIWAIT makes dispTimeSeries wait for user response (see UIRESUME)
% uiwait(handles.fig_dispTimeSeries);

% --- Outputs from this function are returned to the command line.
function varargout = dispTimeSeries_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

% --- Executes on slider movement.
function sl_dispSlice_Callback(hObject, eventdata, handles)
tpVal = get(handles.edt_dispVol,'Value');
tp = ceil(tpVal*(handles.nTimePoints-1))+1;
sliderValue = get(handles.sl_dispSlice,'Value');
%puts the slider value into the edit text component
val = ceil(sliderValue*(handles.nSlices-1))+1;
set(handles.edt_dispSlice,'String', num2str(val)); 
axes(handles.ax_Image);
h1 = imshow(handles.origImg(:,:,val,tp),[]);
set(h1, 'ButtonDownFcn', {@ax_Image_ButtonDownFcn, handles});
axes(handles.ax_processedImg);
h2 = imshow(handles.processedImg(:,:,val,tp),[]);
set(h2, 'ButtonDownFcn', {@ax_processedImg_ButtonDownFcn, handles});
guidata(hObject,handles)

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

% --- Executes on slider movement.
function sl_dispVol_Callback(hObject, eventdata, handles)
slVal = get(handles.edt_dispSlice,'Value');
sl = ceil(slVal*(handles.nSlices-1))+1;
sliderValue = get(handles.sl_dispVol,'Value');
%puts the slider value into the edit text component
val = ceil(sliderValue*(handles.nTimePoints-1))+1;
set(handles.edt_dispVol,'String', num2str(val)); 
axes(handles.ax_Image);
h1 = imshow(handles.origImg(:,:,sl,val),[]);
set(h1, 'ButtonDownFcn', {@ax_Image_ButtonDownFcn, handles});
axes(handles.ax_processedImg);
h2 = imshow(handles.processedImg(:,:,sl,val),[]);
set(h2, 'ButtonDownFcn', {@ax_processedImg_ButtonDownFcn, handles});
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function sl_dispVol_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function edt_dispVol_Callback(hObject, eventdata, handles)

function edt_dispVol_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on mouse press over axes background.
function ax_Image_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to ax_Image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
imgSize = size(handles.origImg);
axesHandle  = get(hObject,'Parent');
coordinates = get(axesHandle,'CurrentPoint');
coordinates = coordinates(1,1:2);
tmpRow = floor(coordinates(2));
tmpCol = floor(coordinates(1));
if (tmpRow>imgSize(1))
    row = imgSize(1);
elseif (tmpRow < 1)
    row = 1;
else
    row = tmpRow;
end
if (tmpCol>imgSize(2))
    col = imgSize(2);
elseif (tmpCol < 1)
    col = 1;
else
    col = tmpCol;
end
s1 = str2num(get(handles.edt_dispSlice,'String'));
origTS = squeeze(handles.origImg(row,col,s1,:));
axes(handles.ax_TS)
plot(origTS);
procTS = squeeze(handles.processedImg(row,col,s1,:));
axes(handles.ax_processedTS)
plot(procTS);

% --- Executes on mouse press over axes background.
function ax_processedImg_ButtonDownFcn(hObject, eventdata, handles)
imgSize = size(handles.origImg);
axesHandle  = get(hObject,'Parent');
coordinates = get(axesHandle,'CurrentPoint');
coordinates = coordinates(1,1:2);
tmpRow = floor(coordinates(2));
tmpCol = floor(coordinates(1));
if (tmpRow>imgSize(1))
    row = imgSize(1);
elseif (tmpRow < 1)
    row = 1;
else
    row = tmpRow;
end
if (tmpCol>imgSize(2))
    col = imgSize(2);
elseif (tmpCol < 1)
    col = 1;
else
    col = tmpCol;
end
s1 = str2num(get(handles.edt_dispSlice,'String'));
origTS = squeeze(handles.origImg(row,col,s1,:));
axes(handles.ax_TS)
plot(origTS);
procTS = squeeze(handles.processedImg(row,col,s1,:));
axes(handles.ax_processedTS)
plot(procTS);

% --- Executes on button press in pb_saveImg.
function pb_saveImg_Callback(hObject, eventdata, handles)
niiStruct = make_nii(handles.processedImg,handles.imgVoxDim,[],64,[]);
niiStruct.hdr.hk.data_type = 'float64';
save_nii(niiStruct,handles.opFname,[]);
