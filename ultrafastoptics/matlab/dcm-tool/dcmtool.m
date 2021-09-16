function varargout = dcmtool(varargin)
% DCMTOOL M-file for dcmtool.fig
%      DCMTOOL, by itself, creates a new DCMTOOL or raises the existing
%      singleton*.
%
%      H = DCMTOOL returns the handle to a new DCMTOOL or the handle to
%      the existing singleton*.
%
%      DCMTOOL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DCMTOOL.M with the given input arguments.
%
%      DCMTOOL('Property','Value',...) creates a new DCMTOOL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dcmtool_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dcmtool_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help dcmtool

% Last Modified by GUIDE v2.5 11-May-2007 15:48:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dcmtool_OpeningFcn, ...
                   'gui_OutputFcn',  @dcmtool_OutputFcn, ...
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


% --- Executes just before dcmtool is made visible.
function dcmtool_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dcmtool (see VARARGIN)

% Choose default command line output for dcmtool
handles.output = hObject;

% Create index list.
n = 6;
m = 4;
poslist = get(handles.m1, 'Position');
nlist = length(poslist);
xlist = zeros(1,nlist);
ylist = zeros(1,nlist);
for k = 1:nlist,
  xlist(k) = poslist{k}(1);
  ylist(k) = poslist{k}(2);
end
[xs, k1] = sort(xlist);
[ys, k2] = sort(ylist(k1));
handles.m1 = handles.m1(k2(k1));
for k = 1:nlist,
  set(handles.m1(k), 'UserData', [ceil(k/m) mod(k,n)+1])
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes dcmtool wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = dcmtool_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;





% --- Executes on button press in dual.
function dual_Callback(hObject, eventdata, handles)
% hObject    handle to dual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dual



% --- Executes on button press in m1pump.
function m1pump_Callback(hObject, eventdata, handles)
% hObject    handle to m1pump (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of m1pump


% --- Executes on button press in m1pass1.
function m1pass1_Callback(hObject, eventdata, handles)
% hObject    handle to m1pass1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of m1pass1


% --- Executes on button press in m1pulse.
function m1pulse_Callback(hObject, eventdata, handles)
% hObject    handle to m1pulse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of m1pulse


% --- Executes on button press in m1pass2.
function m1pass2_Callback(hObject, eventdata, handles)
% hObject    handle to m1pass2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of m1pass2



function m1_Callback(hObject, eventdata, handles)
% hObject    handle to m1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of m1 as text
%        str2double(get(hObject,'String')) returns contents of m1 as a double



% --- Executes during object creation, after setting all properties.
function m1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to m1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function m2_Callback(hObject, eventdata, handles)
% hObject    handle to m2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of m2 as text
%        str2double(get(hObject,'String')) returns contents of m2 as a double


% --- Executes during object creation, after setting all properties.
function m2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to m2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

