function varargout = complexity_MSE(varargin)
% COMPLEXITY_MSE MATLAB code for complexity_MSE.fig
% module for multiscale sample entropy analysis of resting state fMRI data
% ---------------------------------------------------------- 
% Written by Jun Fang
% Version 1.0 
% 1) Input: fMRI time series data, brain mask, r, m, scale, output folder
% 2) Please verify if the orientation of the fMRI data matches the brain mask.
% 3) r and m can take multiple values (enter the values separated by comma 0.2,0.3)
% 4) The basename of the output file is the basename of the brain mask. As matlab
% does not support long filenames, please name the brain mask accordingly. 

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @complexity_MSE_OpeningFcn, ...
                   'gui_OutputFcn',  @complexity_MSE_OutputFcn, ...
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


% --- Executes just before complexity_MSE is made visible.
function complexity_MSE_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to complexity_MSE (see VARARGIN)

% Choose default command line output for complexity_MSE
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
setbgcolor(handles.pb_MSE,clBlue);
setbgcolor(handles.pb_dispMSE,clBlue);
setbgcolor(handles.pb_verImgOri,clBlue);
setbgcolor(handles.pb_opDir,clBlue);
setbgcolor(handles.pb_help,clBlue);
%set(handles.pb_help,'Enable','off');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes complexity_MSE wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = complexity_MSE_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pb_ipDir.
function pb_ipDir_Callback(hObject, eventdata, handles)
% hObject    handle to pb_ipDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
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


function edtTxt_a_Callback(hObject, eventdata, handles)
% hObject    handle to edtTxt_a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtTxt_a as text
%        str2double(get(hObject,'String')) returns contents of edtTxt_a as a double
a = get(hObject,'String');
a = str2num(a);
handles.a = a;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edtTxt_a_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtTxt_a (see GCBO)
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

% --- Executes on button press in pb_MSE.
function pb_MSE_Callback(hObject, eventdata, handles)
% hObject    handle to pb_MSE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
A = sum(strcmp(fieldnames(handles),'img_4D'));
B = sum(strcmp(fieldnames(handles),'brainMask'));
C = sum(strcmp(fieldnames(handles),'r'));
D = sum(strcmp(fieldnames(handles),'m'));
E = sum(strcmp(fieldnames(handles),'a'));
F = sum(strcmp(fieldnames(handles),'opFolder'));
ipChk =[A B C D E F];
clear A B C D E F
set(handles.edtTxt_r,'Enable','off');
set(handles.edtTxt_m,'Enable','off');
set(handles.edtTxt_a,'Enable','off');
if (sum(ipChk)==6)
    mask = handles.brainMask;
    brainVox = find(mask == max(mask(:)));
    r = handles.r;
    m = handles.m;
    a = handles.a;
    % calculate sample entropy
    disp('calculating multiscale sample entropy');
    imgSize = size(mask);
    MSE = zeros(imgSize);
    for i = 1:length(m)
        for j = 1:length(r)
            msg = ['calculating MSE: m=',num2str(m(i)),',r=',num2str(r(j)), ',a=', num2str(a)];
            h = waitbar(0,msg);
            nFail = 0;
            for vox = 1:length(brainVox)
                [row,col,sl] = ind2sub(imgSize,brainVox(vox));
                TS = squeeze(handles.img_4D(row,col,sl,:));
                r_val = r(j)*std(double(TS));
                tmp = sample_entropy(m(i),r_val,TS,a);
                MSE(row,col,sl) = tmp(1);
                nFail = nFail+tmp(2);
                waitbar(vox/length(brainVox));
            end
            close(h)
            opFname = [handles.opFolder,filesep,handles.baseName,'MSE_m',...
                num2str(m(i)),'_r',num2str(r(j)*100),'_a',num2str(a),'per','.nii'];
            niiStruct = make_nii(MSE,handles.imgVoxDim,[],64,[]);
            niiStruct.hdr.hk.data_type = 'float64';
            save_nii(niiStruct,opFname,[]);
        end
    end
    handles.nSlices = imgSize(3);
    guidata(hObject,handles);
    disp('done calculating multiscale sample entropy');
else
    ipParam = {'Input Image','Brain Mask','r','m','a','output folder'};
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
set(handles.edtTxt_a,'String', ' ');
set(handles.edtTxt_r,'String', ' ');
set(handles.edtTxt_m,'String', ' ');
set(handles.edtTxt_a,'Enable','on');
set(handles.edtTxt_r,'Enable','on');
set(handles.edtTxt_m,'Enable','on');
rmfield(handles,'a');
rmfield(handles,'r');
rmfield(handles,'m');
guidata(hObject,handles);

% --- Executes on button press in pb_dispMSE.
function pb_dispMSE_Callback(hObject, eventdata, handles)
% hObject    handle to pb_dispMSE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fname,pname] = uigetfile('*.nii','select the MSE image');
MSEName = [pname,fname];
if (pname==0 & fname==0)
    disp('No image selected to display')
else
    MSE = load_nii(MSEName);
    MSEImg = MSE.img;
    displayImage(MSEImg);
end

% --- Executes on button press in pb_help.
function pb_help_Callback(hObject, eventdata, handles)
% hObject    handle to pb_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
open('helpDocuments.pdf');
