function varargout = complexity_cApEn(varargin)
% module for cross approximate entropy analysis of resting state fMRI data
% ---------------------------------------------------------- 
% Written by AnithaPriya Krishnan 
% Version 1.0 
% contact: Dr.Danny JJ Wang "jjwang@loni.ucla.edu"; 
% Release = 20130402
% 1)Input: fMRI time series data, brain mask, r, m, output folder, ROI voxel
% and filter option
% 2) Please verify if the orientation of the fMRI data matches the brain mask.
% 3) r and m can take multiple values (enter the values separated by comma 0.2,0.3)
% 4) The basename of the output file is the basename of the first fMRI volume. As matlab
% does not support long filenames, please name the fMRI data accordingly. 
% 5) Please use the pre-processing module for low pass filtering and detrendign
% the time series 
% 6) If you would like to use more than one ROI voxel, please decide the ROI 
% voxels before running this module as there is a pause of only three seconds
% before selecting the next ROI. If you need more time, increase the time
% in line 142 of displayImage.m

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @complexity_cApEn_OpeningFcn, ...
                   'gui_OutputFcn',  @complexity_cApEn_OutputFcn, ...
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


% --- Executes just before complexity_cApEn is made visible.
function complexity_cApEn_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
axes(handles.ax_topLeft);
image(imread('LOFT_logo.png'));
set(handles.ax_topLeft,'Visible', 'off');
axes(handles.ax_topRight);
image(imread('brain_fractal.png'));
set(handles.ax_topRight,'Visible', 'off');
clBlue = [11/255 132/255 199/255];
setbgcolor(handles.pb_ipDir,clBlue);
setbgcolor(handles.pb_BrMask,clBlue);
setbgcolor(handles.pb_verImgOri,clBlue);
setbgcolor(handles.pb_opDir,clBlue);
setbgcolor(handles.pb_ROIVox,clBlue);
setbgcolor(handles.pb_cApEn,clBlue);
setbgcolor(handles.pb_dispCApEn,clBlue);
setbgcolor(handles.pb_help,clBlue);
% UIWAIT makes complexity_cApEn wait for user response (see UIRESUME)
% uiwait(handles.fig_crossApEn);
% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = complexity_cApEn_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

% --- Executes on button press in pb_ipDir.
function pb_ipDir_Callback(hObject, eventdata, handles)
ipFormat = cell2mat(inputdlg('Input 3D or 4D','Input selection'));
if isempty(ipFormat)
    disp('please choose the input format: 3D or 4D');
else
    if (strcmp(ipFormat,'3D')==1 | strcmp(ipFormat,'3d')==1)
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

% --- Executes on button press in pb_BrMask.
function pb_BrMask_Callback(hObject, eventdata, handles)
% hObject    handle to pb_BrMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
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

% --- Executes on button press in pb_verImgOri.
function pb_verImgOri_Callback(hObject, eventdata, handles)
% hObject    handle to pb_verImgOri (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mask = verifyImgOri(handles.img_4D(:,:,:,1),handles.brainMask);
handles.brainMask = mask;
guidata(hObject,handles);

function edtTxt_r_Callback(hObject, eventdata, handles)
% hObject    handle to edtTxt_r (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtTxt_r as text
%        str2double(get(hObject,'String')) returns contents of edtTxt_r as a double
r = get(hObject,'String');
r = str2num(r);
handles.r = r;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edtTxt_r_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtTxt_r (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edtTxt_m_Callback(hObject, eventdata, handles)
% hObject    handle to edtTxt_m (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtTxt_m as text
%        str2double(get(hObject,'String')) returns contents of edtTxt_m as a double
m = get(hObject,'String');
m = str2num(m);
handles.m = m;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edtTxt_m_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtTxt_m (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pb_opDir.
function pb_opDir_Callback(hObject, eventdata, handles)
% hObject    handle to pb_opDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
opFolder = uigetdir;
handles.opFolder = opFolder;
guidata(hObject, handles);

% --- Executes on button press in pb_ROIVox.
function pb_ROIVox_Callback(hObject, eventdata, handles)
% hObject    handle to pb_ROIVox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prompt = 'Enter the number of ROI voxels';
dlg_title = 'ROI selection';
def = {'1'};
nPoints = inputdlg(prompt,dlg_title,1,def);
if isempty(nPoints)
    disp('Specify the number of ROI voxels');
else
    nPoints = str2num(char(nPoints));
    voxROI = displayImage(handles.img_4D(:,:,:,1),nPoints);
    handles.nPoints = nPoints;
    handles.voxROI = voxROI;
    guidata(hObject,handles);
end

% --- Executes on button press in pb_cApEn.
function pb_cApEn_Callback(hObject, eventdata, handles)
% hObject    handle to pb_cApEn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
A = sum(strcmp(fieldnames(handles),'img_4D'));
B = sum(strcmp(fieldnames(handles),'brainMask'));
C = sum(strcmp(fieldnames(handles),'r'));
D = sum(strcmp(fieldnames(handles),'m'));
E = sum(strcmp(fieldnames(handles),'opFolder'));
F = sum(strcmp(fieldnames(handles),'voxROI'));
ipList = [A B C D E F];
clear A B C D E F
ipParam = {'Input Image','Brain Mask','r','m','output folder','ROI voxel'};
set(handles.edtTxt_r,'Enable','off');
set(handles.edtTxt_m,'Enable','off');
if (sum(ipList)== 6)
    mask = handles.brainMask;
    r = handles.r;
    m = handles.m;
    % calculate approximate entropy
    disp('calculating cross approximate entropy');
    imgSize = size(mask);
    cApEn = zeros(imgSize);
    for i = 1:length(m)
        for j = 1:length(r)
            disp(['calculating cross ApEn: m=',num2str(m(i)),',r=',num2str(r(j))]);
            nFail = 0;
            for p = 1:handles.nPoints
                cApEn = calcCrossApEn_noFilt(mask,handles.img_4D,m(i),r(j),handles.voxROI(p,:));
                opFname = [handles.opFolder,filesep,handles.baseName,'cApEn_m',...
                    num2str(m(i)),'_r',num2str(r(j)*100),'per',...
                    num2str(handles.voxROI(p,1)),num2str(handles.voxROI(p,2)),...
                    num2str(handles.voxROI(p,3)),'.nii'];
                niiStruct = make_nii(cApEn,handles.imgVoxDim,[],64,[]);
                niiStruct.hdr.hk.data_type = 'float64';
                save_nii(niiStruct,opFname,[]);
            end
        end
    end
    guidata(hObject,handles);
    disp('done calculating cross approximate entropy');
else
    ind = find(ipList==0)
    misParam = ipParam(ind);
    if (length(ind)==1)
        msg = ['Missing input:',misParam];
    else
        msg = ['Missing inputs:',misParam];
    end
    disp(msg);
    msgbox(msg);
end
set(handles.edtTxt_r,'String', ' ');
set(handles.edtTxt_m,'String', ' ');
set(handles.edtTxt_r,'Enable','on');
set(handles.edtTxt_m,'Enable','on');
rmfield(handles,'r');
rmfield(handles,'m');
rmfield(handles,'voxROI');
guidata(hObject,handles);


% --- Executes on button press in pb_dispCApEn.
function pb_dispCApEn_Callback(hObject, eventdata, handles)
% hObject    handle to pb_dispCApEn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fname,pname] = uigetfile('*.nii','select the cApEn image');
ApEnName = [pname,fname];
if (pname==0 & fname==0)
    disp('No image selected to display')
else
    ApEn = load_nii(ApEnName);
    ApEnImg = ApEn.img;
    displayImage(ApEnImg);
end

% --- Executes on button press in pb_help.
function pb_help_Callback(hObject, eventdata, handles)
% hObject    handle to pb_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
open('helpDocuments.pdf');