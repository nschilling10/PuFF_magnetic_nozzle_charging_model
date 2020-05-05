function varargout = pulsed_power(varargin)
% PULSED_POWER MATLAB code for pulsed_power.fig
%      PULSED_POWER, by itself, creates a new PULSED_POWER or raises the existing
%      singleton*.
%
%      H = PULSED_POWER returns the handle to a new PULSED_POWER or the handle to
%      the existing singleton*.
%
%      PULSED_POWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PULSED_POWER.M with the given input arguments.
%
%      PULSED_POWER('Property','Value',...) creates a new PULSED_POWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pulsed_power_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pulsed_power_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pulsed_power

% Last Modified by GUIDE v2.5 12-Jan-2018 08:54:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pulsed_power_OpeningFcn, ...
                   'gui_OutputFcn',  @pulsed_power_OutputFcn, ...
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

% --- Executes just before pulsed_power is made visible.
function pulsed_power_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pulsed_power (see VARARGIN)

% Choose default command line output for pulsed_power
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using pulsed_power.
if strcmp(get(hObject,'Visible'),'off')
    Run_ode_solver_Callback(hObject, eventdata, handles)
end

% UIWAIT makes pulsed_power wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = pulsed_power_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in Run_ode_solver.
function Run_ode_solver_Callback(hObject, eventdata, handles)
% hObject    handle to Run_ode_solver (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);
cla;
V0 = eval(handles.Charging_voltage.String);
I0 = 0;
C = eval(handles.Capacitance.String);
R = eval(handles.Resistance_transmission_line.String) + eval(handles.Resistance_load.String);
L = eval(handles.Inductance_transmission_line.String) + eval(handles.Inductance_load.String);

%time constant, based on half period of circuit
t_half_period = pi*sqrt(L*C);
tf = t_half_period;
dt = tf/100;

%initial conditions
IV = [V0 I0];

[t,y]=ode45(@(t,y) rlcfuns(t,y,R,L,C),(0:dt:tf),IV);   
V = y(:,1);
I = y(:,2);

plot(t*1e6,V/1e3,'--'), hold on
plot(t*1e6,I/1e3,'-'), hold off
xlabel('t (\mus)');
ylabel('V (kV), I (kA)')
lh = legend('V','I','location','northeast'); legend boxoff
grid on


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', {'plot(rand(5))', 'plot(sin(1:0.01:25))', 'bar(1:.5:10)', 'plot(membrane)', 'surf(peaks)'});



function Resistance_load_Callback(hObject, eventdata, handles)
% hObject    handle to Resistance_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Resistance_load as text
%        str2double(get(hObject,'String')) returns contents of Resistance_load as a double


% --- Executes during object creation, after setting all properties.
function Resistance_load_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Resistance_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Inductance_load_Callback(hObject, eventdata, handles)
% hObject    handle to Inductance_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Inductance_load as text
%        str2double(get(hObject,'String')) returns contents of Inductance_load as a double


% --- Executes during object creation, after setting all properties.
function Inductance_load_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Inductance_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Inductance_transmission_line_Callback(hObject, eventdata, handles)
% hObject    handle to Inductance_transmission_line (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Inductance_transmission_line as text
%        str2double(get(hObject,'String')) returns contents of Inductance_transmission_line as a double


% --- Executes during object creation, after setting all properties.
function Inductance_transmission_line_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Inductance_transmission_line (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Resistance_transmission_line_Callback(hObject, eventdata, handles)
% hObject    handle to Resistance_transmission_line (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Resistance_transmission_line as text
%        str2double(get(hObject,'String')) returns contents of Resistance_transmission_line as a double


% --- Executes during object creation, after setting all properties.
function Resistance_transmission_line_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Resistance_transmission_line (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Capacitance_Callback(hObject, eventdata, handles)
% hObject    handle to Capacitance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Capacitance as text
%        str2double(get(hObject,'String')) returns contents of Capacitance as a double


% --- Executes during object creation, after setting all properties.
function Capacitance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Capacitance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Charging_voltage_Callback(hObject, eventdata, handles)
% hObject    handle to Charging_voltage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Stringkkkkkkkk') returns contents of Charging_voltage as text
%        str2double(get(hObject,'String')) returns contents of Charging_voltage as a double


% --- Executes during object creation, after setting all properties.
function Charging_voltage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Charging_voltage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%this is the ode circuit solver
    function dy = rlcfuns(t,y,R,L,C)
        %R circuit and load  resistance
        %L circuit and load inductance
        %C capacitance
        
        dy = zeros(2,1);
        
        voltage = y(1);
        current = y(2);
        
        %circuit equations
        dy(1) = -current/C;
        dy(2) = (voltage - current*R)/L;
        
