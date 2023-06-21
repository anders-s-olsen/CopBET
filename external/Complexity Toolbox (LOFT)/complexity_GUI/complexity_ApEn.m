function varargout = complexity_ApEn(varargin)
% module for approximate entropy analysis of resting state fMRI data
% ---------------------------------------------------------- 
% Written by AnithaPriya Krishnan 
% Version 1.0 
% contact: Dr.Danny JJ Wang "jjwang@loni.ucla.edu"; 
% Release = 20130402
% 1) Input: fMRI time series data, brain mask, r, m, output folder
% 2) Please verify if the orientation of the fMRI data matches the brain mask.
% 3) r and m can take multiple values (enter the values separated by comma 0.2,0.3)
% 4) The basename of the output file is the basename of the brain mask. As matlab
% does not support long filenames, please name the brain mask accordingly. 


gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @complexity_ApEn_OpeningFcn, ...
                   'gui_OutputFcn',  @complexity_ApEn_OutputFcn, ...
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


% --- Executes just before complexity_ApEn is made visible.
function complexity_ApEn_OpeningFcn(hObject, eventdata, handles, varargin)
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
setbgcolor(handles.pb_ApEn,clBlue);
setbgcolor(handles.pb_dispApEn,clBlue);
setbgcolor(handles.pb_verImgOri,clBlue);
setbgcolor(handles.pb_opDir,clBlue);
setbgcolor(handles.pb_help,clBlue);
%set(handles.pb_help,'Enable','off');
% UIWAIT makes complexity_ApEn wait for user response (see UIRESUME)
% uiwait(handles.fig_ApEn);
% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = complexity_ApEn_OutputFcn(hObject, eventdata, handles) 
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

% --- Executes on button press in pb_ApEn.
function pb_ApEn_Callback(hObject, eventdata, handles)
% hObject    handle to pb_ApEn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
A = sum(strcmp(fieldnames(handles),'img_4D'));
B = sum(strcmp(fieldnames(handles),'brainMask'));
C = sum(strcmp(fieldnames(handles),'r'));
D = sum(strcmp(fieldnames(handles),'m'));
E = sum(strcmp(fieldnames(handles),'opFolder'));
ipChk =[A B C D E];
clear A B C D E
set(handles.edtTxt_r,'Enable','off');
set(handles.edtTxt_m,'Enable','off');
if (sum(ipChk)==5)
    mask = handles.brainMask;
    brainVox = find(mask == max(mask(:)));
    r = handles.r;
    m = handles.m;
    % calculate approximate entropy
    disp('calculating self approximate entropy');
    imgSize = size(mask);
    ApEn = zeros(imgSize);
    for i = 1:length(m)
        for j = 1:length(r)
            msg = ['calculating self ApEn: m=',num2str(m(i)),',r=',num2str(r(j))];
            h = waitbar(0,msg);
            nFail = 0;
            for vox = 1:length(brainVox)
                [row,col,sl] = ind2sub(imgSize,brainVox(vox));
                TS1 = squeeze(handles.img_4D(row,col,sl,:));
                TS2 = TS1;
                r_val = r(j)*std(double(TS1));
                tmp =  cross_approx_entropy(m(i),r_val,TS1,TS2);
                ApEn(row,col,sl) = tmp(1);
                nFail = nFail+tmp(2);
                waitbar(vox/length(brainVox));
            end
            close(h)
            opFname = [handles.opFolder,filesep,handles.baseName,'ApEn_m',...
                num2str(m(i)),'_r',num2str(r(j)*100),'per','.nii'];
            niiStruct = make_nii(ApEn,handles.imgVoxDim,[],64,[]);
            niiStruct.hdr.hk.data_type = 'float64';
            save_nii(niiStruct,opFname,[]);
        end
    end
    handles.nSlices = imgSize(3);
    guidata(hObject,handles);
    disp('done calculating self approximate entropy');
else
    ipParam = {'Input Image','Brain Mask','r','m','output folder'};
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
set(handles.edtTxt_r,'String', ' ');
set(handles.edtTxt_m,'String', ' ');
set(handles.edtTxt_r,'Enable','on');
set(handles.edtTxt_m,'Enable','on');
rmfield(handles,'r');
rmfield(handles,'m');
guidata(hObject,handles);


% --- Executes on button press in pb_dispApEn.
function pb_dispApEn_Callback(hObject, eventdata, handles)
% hObject    handle to pb_dispApEn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fname,pname] = uigetfile('*.nii','select the ApEn image');
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
