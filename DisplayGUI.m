function varargout = DisplayGUI(varargin)
% DisplayGUI MATLAB code for DisplayGUI.fig
%      DisplayGUI, by itself, creates a new DisplayGUI or raises the existing
%      singleton*.
%
%      H = DisplayGUI returns the handle to a new DisplayGUI or the handle to
%      the existing singleton*.
%
%      DisplayGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DisplayGUI.M with the given input arguments.
%
%      DisplayGUI('Property','Value',...) creates a new DisplayGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DisplayGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DisplayGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DisplayGUI

% Last Modified by GUIDE v2.5 21-Feb-2022 06:01:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DisplayGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @DisplayGUI_OutputFcn, ...
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

% --- Executes just before DisplayGUI is made visible.
function DisplayGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DisplayGUI (see VARARGIN)

% Choose default command line output for DisplayGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DisplayGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = DisplayGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in RUN.
function RUN_Callback(hObject, eventdata, handles) %Values Stored when run is selected.
% hObject    handle to RUN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

LBStartTime = get(handles.LBStartTime,'string'); %Lower Bound for Start Time
UBStartTime = get(handles.UBStartTime,'string'); %Upper Bound for Start Time
NumSamples = str2num(get(handles.NumOfSamples,'string')); %Number of Samples Taken
tFinal = datenum(UBStartTime); %Converting Final Start Time to DateNum
tInit = datenum(LBStartTime); %Converting Initial Start Time to DateNum
dT = (tFinal - tInit)/NumSamples; %Calculating Time Step
SampleDates(1) = tInit; %First Date of the Vector is Lower Bound Start Time
SampleDates(str2num(get(handles.NumOfSamples,'string'))) = tFinal; %Last Date of the Vector is Upper Bound Start Time

tb = tInit; %Initializing time
WB = waitbar((0)/NumSamples,'I am Loading'); %Creating Progress Bar

for i = 2:str2num(get(handles.NumOfSamples,'string'))-1 %Stepping through and Adding Values of to Start Time Vector    
    
    if getappdata(WB,'canceling')% Check for clicked Cancel button
        break
    end
    
    waitbar((i-1)/NumSamples,WB,'Initializing Start Times');
    tb = tb + dT;
    SampleDates(i) = tb;
    
end
close(WB)

%---Converts Back to String Vector for Time in STK Format
SampleDates = datetime(SampleDates, 'Format','yyyy MM dd HH:mm:ss:SSS','ConvertFrom', 'datenum'); 

%---Displaying Saved Values Into GUI Table
DisplayTable = cellstr(SampleDates');
set(handles.DisplayTable,'Data',DisplayTable,'ColumnName',{'Start Time'});


%---Assigning Selected Outbound Traj---%
OutboundTrajs = handles.OutboundTraj.String;
OutboundTrajIndex = handles.OutboundTraj.Value;
OutboundTraj = OutboundTrajs{OutboundTrajIndex};

%---Assigning Selected Return Traj---%
ReturnTrajs = handles.ReturnTraj.String;
ReturnTrajIndex = handles.ReturnTraj.Value;
ReturnTraj = ReturnTrajs{ReturnTrajIndex};

%---Assigning Spacecraft Parameters---%
EngineIsp = get(handles.Isp,'Value');
EngineThrust = get(handles.Thrust,'Value');
DryMass = get(handles.DryMass,'Value');
DryMass = get(handles.PropMass,'Value');
MiningDur = get(handles.MiningDuration,'Value');

switch OutboundTraj %Nested Switch Case for Outbound/Return Trajs
    case 'Hayabusa Mission Profile'
        switch ReturnTraj
            case 'Direct Return to Earth'
            case 'EGA to Return to Earth'
        end
    case 'Earth Escape Direct to Asteroid Plane'
        switch ReturnTraj
            case 'Direct Return to Earth'
            case 'EGA to Return to Earth'
        end
    case 'EGA to Asteroid Plane'
        switch ReturnTraj
            case 'Direct Return to Earth'
            case 'EGA to Return to Earth'
        end
    case 'EGA in Earth Plane to Raise Orbit'
        switch ReturnTraj
            case 'Direct Return to Earth'
            case 'EGA to Return to Earth'
        end
    case 'Direct Trajectory'
        switch ReturnTraj
            case 'Direct Return to Earth'
            case 'EGA to Return to Earth'
        end
end 

% --- Executes when uipanel1 is resized.
function uipanel1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to uipanel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in OutboundTraj.
function OutboundTraj_Callback(hObject, eventdata, handles)
% hObject    handle to OutboundTraj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns OutboundTraj contents as cell array
%        contents{get(hObject,'Value')} returns selected item from OutboundTraj


% --- Executes during object creation, after setting all properties.
function OutboundTraj_CreateFcn(hObject, eventdata, handles)
% hObject    handle to OutboundTraj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function LBStartTime_Callback(hObject, eventdata, handles)
% hObject    handle to LBStartTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LBStartTime as text
%        str2double(get(hObject,'String')) returns contents of LBStartTime as a double


% --- Executes during object creation, after setting all properties.
function LBStartTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LBStartTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function UBStartTime_Callback(hObject, eventdata, handles)
% hObject    handle to UBStartTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of UBStartTime as text
%        str2double(get(hObject,'String')) returns contents of UBStartTime as a double


% --- Executes during object creation, after setting all properties.
function UBStartTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UBStartTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ReturnTraj.
function ReturnTraj_Callback(hObject, eventdata, handles)
% hObject    handle to ReturnTraj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ReturnTraj contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ReturnTraj


% --- Executes during object creation, after setting all properties.
function ReturnTraj_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ReturnTraj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MiningDur_Callback(hObject, eventdata, handles)
% hObject    handle to MiningDur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MiningDur as text
%        str2double(get(hObject,'String')) returns contents of MiningDur as a double


% --- Executes during object creation, after setting all properties.
function MiningDur_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MiningDur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Uncertainty_Callback(hObject, eventdata, handles)
% hObject    handle to Uncertainty (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Uncertainty as text
%        str2double(get(hObject,'String')) returns contents of Uncertainty as a double


% --- Executes during object creation, after setting all properties.
function Uncertainty_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Uncertainty (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Isp_Callback(hObject, eventdata, handles)
% hObject    handle to Isp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Isp as text
%        str2double(get(hObject,'String')) returns contents of Isp as a double


% --- Executes during object creation, after setting all properties.
function Isp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Isp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Thrust_Callback(hObject, eventdata, handles)
% hObject    handle to Thrust (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Thrust as text
%        str2double(get(hObject,'String')) returns contents of Thrust as a double


% --- Executes during object creation, after setting all properties.
function Thrust_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Thrust (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function SaveLocation_Callback(hObject, eventdata, handles)
% hObject    handle to SaveLocation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of SaveLocation as text
%        str2double(get(hObject,'String')) returns contents of SaveLocation as a double
OutboundTrajs = handles.OutboundTraj.String;


% --- Executes during object creation, after setting all properties.
function SaveLocation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SaveLocation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SAVE.
function SAVE_Callback(hObject, eventdata, handles)
% hObject    handle to SAVE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
LBStartTime = get(handles.LBStartTime,'string'); %Lower Bound for Start Time
UBStartTime = get(handles.UBStartTime,'string'); %Upper Bound for Start Time
NumSamples = str2num(get(handles.NumOfSamples,'string')); %Number of Samples Taken
tFinal = datenum(UBStartTime); %Converting Final Start Time to DateNum
tInit = datenum(LBStartTime); %Converting Initial Start Time to DateNum
dT = (tFinal - tInit)/NumSamples; %Calculating Time Step
SampleDates(1) = tInit; %First Date of the Vector is Lower Bound Start Time
SampleDates(str2num(get(handles.NumOfSamples,'string'))) = tFinal; %Last Date of the Vector is Upper Bound Start Time

tb = tInit; %Initializing time
for i = 2:str2num(get(handles.NumOfSamples,'string'))-1 %Stepping through and Adding Values of to Start Time Vector
    tb = tb + dT;
    SampleDates(i) = tb;
end

SampleDates = datetime(SampleDates, 'Format','yyyy MM dd HH:mm:ss:SSS','ConvertFrom', 'datenum') %Converts Back to String Vector for Time in STK Format
DateString = string(SampleDates)' 
[file,path] = uiputfile('*.xlsx');
filename = fullfile(path,file);
xlswrite(filename,DateString);


% --- Executes on button press in Clear.
function Clear_Callback(hObject, eventdata, handles)
% hObject    handle to Clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%OutboundTrajs = handles.OutboundTraj.String;
%ReturnTrajs = handles.ReturnTraj.String;
set(handles.LBStartTime,'string','07-Jun-2022 00:00:00');
set(handles.UBStartTime,'string','07-Jun-2022 00:00:00');
set(handles.NumOfSamples,'string','');
set(handles.OutboundTraj,'Value',1);
set(handles.ReturnTraj,'Value',1);
set(handles.Uncertainty,'string','');
set(handles.Isp,'string','');
set(handles.Thrust,'string','');
set(handles.DryMass,'string','');
set(handles.PropMass,'string','');
set(handles.MiningDuration,'string','');
set(handles.SaveLoc,'string','');


% --- Executes when Clear is resized.
function Clear_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to Clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Clear_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on OutboundTraj and none of its controls.
function OutboundTraj_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to OutboundTraj (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

function PropMass_Callback(hObject, eventdata, handles)
% hObject    handle to PropMass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of PropMass as text
%        str2double(get(hObject,'String')) returns contents of PropMass as a double

% --- Executes during object creation, after setting all properties.
function PropMass_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PropMass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function DryMass_Callback(hObject, eventdata, handles)
% hObject    handle to DryMass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DryMass as text
%        str2double(get(hObject,'String')) returns contents of DryMass as a double


% --- Executes during object creation, after setting all properties.
function DryMass_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DryMass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function MiningDuration_Callback(hObject, eventdata, handles)
% hObject    handle to MiningDuration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MiningDuration as text
%        str2double(get(hObject,'String')) returns contents of MiningDuration as a double


% --- Executes during object creation, after setting all properties.
function MiningDuration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MiningDuration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function NumOfSamples_Callback(hObject, eventdata, handles)
% hObject    handle to NumOfSamples (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NumOfSamples as text
%        str2double(get(hObject,'String')) returns contents of NumOfSamples as a double


% --- Executes during object creation, after setting all properties.
function NumOfSamples_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NumOfSamples (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function SaveLoc_Callback(hObject, eventdata, handles)
% hObject    handle to SaveLoc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SaveLoc as text
%        str2double(get(hObject,'String')) returns contents of SaveLoc as a double


% --- Executes during object creation, after setting all properties.
function SaveLoc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SaveLoc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object deletion, before destroying properties.
function DisplayTable_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to DisplayTable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
