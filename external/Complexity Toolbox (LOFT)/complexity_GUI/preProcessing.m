function varargout = preProcessing(varargin)
% module for filtering and detrending the time series of resting state fMRI data
% ---------------------------------------------------------- 
% Written by AnithaPriya Krishnan 
% Version 1.0 
% contact: Dr.Danny JJ Wang "jjwang@loni.ucla.edu"; 
% Release = 20130402
% 1) uses a butterworth zero-phase filter for filtering. The cut off
% frequency is normalized and takes a value between 0 and 1.
% 2) the output images are displayed in the dispTimeSeries module but are
% not saved. Please save the images in dispTimeSeries module

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @preProcessing_OpeningFcn, ...
                   'gui_OutputFcn',  @preProcessing_OutputFcn, ...
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

% --- Executes just before preProcessing is made visible.
function preProcessing_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
axes(handles.ax_topLeft);
image(imread('LOFT_logo.png'));
set(handles.ax_topLeft,'Visible', 'off');
axes(handles.ax_topRight);
image(imread('brain_fractal.png'));
set(handles.ax_topRight,'Visible', 'off');
clBlue = [11/255 132/255 199/255];
setbgcolor(handles.pb_ipDir,clBlue);
setbgcolor(handles.pb_applyFilter,clBlue);
setbgcolor(handles.pb_brainMask,clBlue);
setbgcolor(handles.pb_detrend,clBlue);
setbgcolor(handles.pb_opDir,clBlue);
handles.detrendMean = 0;
handles.detrendLinear = 0;
% Update handles structure
guidata(hObject, handles);
% UIWAIT makes preProcessing wait for user response (see UIRESUME)
% uiwait(handles.fig_preProcessing);


% --- Outputs from this function are returned to the command line.
function varargout = preProcessing_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

% --- Executes on button press in pb_ipDir.
function pb_ipDir_Callback(hObject, eventdata, handles)
ipFormat = cell2mat(inputdlg('Input 3D or 4D','Input selection'));
if isempty(ipFormat)
    disp('please choose the input format: 3D or 4D');
else
    if (strcmp(ipFormat,'3D')==1)
        dirName = uigetdir;
        if (dirName == 0)
            disp('no image input directory selected');
        else
            imgStruct = readImages4D(dirName);
            handles.img_4D = imgStruct.img_4D;
            handles.baseName = imgStruct.bName;
            handles.imgVoxDim = imgStruct.voxDim;
        end
    else
        [fname,pname] = uigetfile('*.*','select the 4D image');
        if (fname==0 & pname==0)
            disp('4D image file not selected');
        else
            imgName = [pname,fname];
            imgStruct = load_nii(imgName);
            handles.img_4D = imgStruct.img;
            [p,f,e] = fileparts(imgName);
            handles.baseName = f;
            handles.imgVoxDim = imgStruct.hdr.dime.pixdim(2:4);
        end
    end
end
guidata(hObject,handles);
disp('done reading input images');

% --- Executes on button press in pb_brainMask.
function pb_brainMask_Callback(hObject, eventdata, handles)
[fname,pname] = uigetfile('*.*','select the brain mask');
if (fname==0 & pname==0)
    disp('Brain mask not selected');
else
    mask_file = [pname,fname];
    mask = load_nii(mask_file);
    if (size(mask.img) ~= size(handles.img_4D(:,:,:,1)))
        msgbox('mask and input image dimensions do not match');
    else
        handles.brainMask = mask.img;
        guidata(hObject, handles);
    end
end

% --- Executes on button press in pb_opDir.
function pb_opDir_Callback(hObject, eventdata, handles)
opFolder = uigetdir;
handles.opFolder = opFolder;
guidata(hObject, handles);

function edt_order_Callback(hObject, eventdata, handles)
order = get(hObject,'String');
order = str2num(order);
handles.order = order;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edt_order_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edt_cutOffFreq_Callback(hObject, eventdata, handles)
cutOffFreq = get(hObject,'String');
cutOffFreq = str2double(cutOffFreq);
handles.cutOffFreq = cutOffFreq;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edt_cutOffFreq_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pb_applyFilter.
function pb_applyFilter_Callback(hObject, eventdata, handles)
A = sum(strcmp(fieldnames(handles),'img_4D'));
B = sum(strcmp(fieldnames(handles),'brainMask'));
C = sum(strcmp(fieldnames(handles),'opFolder'));
D = sum(strcmp(fieldnames(handles),'order'));
E = sum(strcmp(fieldnames(handles),'cutOffFreq'));
ipChk = [A B C D E];
clear A B C D E;
set(handles.edt_order,'Enable','off');
set(handles.edt_cutOffFreq,'Enable','off');
if (sum(ipChk)==5)
    mask = handles.brainMask;
    brainVox = find(mask == max(mask(:)));
    disp('Filtering the time series');
    imgSize = size(mask);
    filteredImg_4D = zeros(size(handles.img_4D));
    [b,a]= butter(handles.order,handles.cutOffFreq,'low'); 
    msg = ['filtering the time seres: order=',num2str(handles.order),...
        ',cut off frequency',num2str(handles.cutOffFreq)];
    h = waitbar(0,msg);
    for vox = 1:length(brainVox)
        [row,col,sl] = ind2sub(imgSize,brainVox(vox));
        TS1 = squeeze(handles.img_4D(row,col,sl,:));
        filteredTS1 = filtfilt(b,a,double(TS1));
        filteredImg_4D(row,col,sl,:) = filteredTS1;
        waitbar(vox/length(brainVox));
    end
    close(h);
    disp('Done filtering the time series');
else
    ipParam = {'Input Folder','Brain Mask','output folder','order of the filter','cut off frequency'};
    ind = find(ipChk==0)
    misParam = ipParam(ind);
    if (length(ind)==1)
        msg = ['Missing input:',misParam];
    else
        msg = ['Missing inputs:',misParam];
    end
    disp(msg);
    msgbox(msg);
end
set(handles.edt_order,'String', ' ');
set(handles.edt_cutOffFreq,'String', ' ');
set(handles.edt_order,'Enable','on');
set(handles.edt_cutOffFreq,'Enable','on');
opFname = [handles.opFolder,filesep,handles.baseName,'_fil_ord',...
    num2str(handles.order),'_NCF',num2str(handles.cutOffFreq*100),'per.nii'];
rmfield(handles,'order');
rmfield(handles,'cutOffFreq');
handles.filteredImg_4D = filteredImg_4D;
guidata(hObject,handles);
dispTimeSeries(handles.img_4D,filteredImg_4D,opFname,handles.imgVoxDim);


% --- Executes on button press in rb_detrendMean.
function rb_detrendMean_Callback(hObject, eventdata, handles)
if (get(handles.rb_detrendMean,'Value')==1)
    handles.detrendMean = 1;
    guidata(hObject,handles);
end


% --- Executes on button press in rb_detrendLinear.
function rb_detrendLinear_Callback(hObject, eventdata, handles)
if (get(handles.rb_detrendLinear,'Value')==1)
    handles.detrendLinear = 1;
    guidata(hObject,handles);
end

% --- Executes on button press in pb_detrend.
function pb_detrend_Callback(hObject, eventdata, handles)
A = sum(strcmp(fieldnames(handles),'img_4D'));
B = sum(strcmp(fieldnames(handles),'brainMask'));
C = sum(strcmp(fieldnames(handles),'opFolder'));
D = (get(handles.rb_detrendMean,'Value') | get(handles.rb_detrendLinear,'Value'));
ipChk = [A B C D];
clear A B C D;
set(handles.rb_detrendMean,'Enable','off');
set(handles.rb_detrendLinear,'Enable','off');
if (sum(ipChk)==4)
    mask = handles.brainMask;
    brainVox = find(mask == max(mask(:)));
    disp('Detrending the time series');
    imgSize = size(mask);
    if (handles.detrendMean == 1)
        detrendMean_4D = zeros(size(handles.img_4D));
    end
    if (handles.detrendLinear == 1)
        detrendLinear_4D = zeros(size(handles.img_4D));
    end
    msg = ['Detrending the time seres'];
    h = waitbar(0,msg);
    for vox = 1:length(brainVox)
        [row,col,sl] = ind2sub(imgSize,brainVox(vox));
        TS1 = double(squeeze(handles.img_4D(row,col,sl,:)));
        if (handles.detrendMean==1)
            detrendMean_4D(row,col,sl,:) = detrend(TS1,'constant');
        end
        if (handles.detrendLinear == 1)
            detrendLinear_4D(row,col,sl,:) = detrend(TS1);
        end
        waitbar(vox/length(brainVox));
    end
    close(h);
    disp('Done detrending the time series');
else
    ipParam = {'Input Image','Brain Mask','remove mean or linear trend'};
    ind = find(ipChk==0)
    misParam = ipParam(ind);
    if (length(ind)==1)
        msg = ['Missing input:',misParam];
    else
        msg = ['Missing inputs:',misParam];
    end
    disp(msg);
    msgbox(msg);
end
set(handles.rb_detrendMean,'Value',0);
set(handles.rb_detrendLinear,'Value',0);
set(handles.rb_detrendMean,'Enable','on');
set(handles.rb_detrendLinear,'Enable','on');
if (handles.detrendMean == 1)
    handles.detrendMean = 0;
    opFname = [handles.opFolder,filesep,handles.baseName,'_detrendMean.nii'];
    dispTimeSeries(handles.img_4D,detrendMean_4D,opFname,handles.imgVoxDim);
end
if (handles.detrendLinear == 1)
    handles.detrendLinear = 0;
    opFname = [handles.opFolder,filesep,handles.baseName,'_detrendLinear.nii'];
    dispTimeSeries(handles.img_4D,detrendLinear_4D,opFname,handles.imgVoxDim);
end
