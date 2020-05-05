function varargout = flux_compression(varargin)
% FLUX_COMPRESSION MATLAB code for flux_compression.fig
%      FLUX_COMPRESSION, by itself, creates a new FLUX_COMPRESSION or raises the existing
%      singleton*.
%
%      H = FLUX_COMPRESSION returns the handle to a new FLUX_COMPRESSION or the handle to
%      the existing singleton*.
%
%      FLUX_COMPRESSION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FLUX_COMPRESSION.M with the given input arguments.
%
%      FLUX_COMPRESSION('Property','Value',...) creates a new FLUX_COMPRESSION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before flux_compression_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to flux_compression_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help flux_compression

% Last Modified by GUIDE v2.5 12-Jan-2018 16:33:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @flux_compression_OpeningFcn, ...
                   'gui_OutputFcn',  @flux_compression_OutputFcn, ...
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

% --- Executes just before flux_compression is made visible.
function flux_compression_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to flux_compression (see VARARGIN)

% Choose default command line output for flux_compression
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

%turn on the toolbar
set(handles.figure1,'Toolbar','figure');

% This sets up the initial plot - only do when we are invisible
% so window can get raised using flux_compression.
if strcmp(get(hObject,'Visible'),'off')
    Run_ode_solver_Callback(hObject, eventdata, handles)
end

% UIWAIT makes flux_compression wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = flux_compression_OutputFcn(hObject, eventdata, handles)
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
mu0 = 4*pi*1e-7;
v = eval(handles.Conductor_velocity.String);
I0 = eval(handles.I0.String);
width = eval(handles.width.String);
height0 = eval(handles.height.String);
tau = height0/v*0.95;
handles.dwell_time.String = num2str(tau*1e6,3);
depth = eval(handles.depth.String);
L_FCG = mu0 * width * height0 /depth;
handles.Inductance_FCG.String = [num2str(L_FCG*1e9,3)];
%L_FCG = eval(handles.Inductance_FCG.String);
L_load = eval(handles.Inductance_load.String);
R = eval(handles.Resistance_transmission_line.String) + eval(handles.Resistance_load.String);
L = L_FCG + L_load;
B0 = 4*pi*1e-7 * I0 /depth;
handles.B_initial.String = num2str(B0);
%eval(handles.B_initial.String);

flux = B0*width*height0;

%time constant, based on half period of circuit
%t_half_period = pi*sqrt(L*C);
tf = tau;
dt = tf/100;

%initial conditions
Iode = [I0];

[t,y]=ode45(@(t,y) fcgfuns(t,y,R,L_load,v,width,height0,depth,flux),(0:dt:tf),Iode);   
I = y(:,1);
height = height0 - v*t;  
A = width*height;
B = flux./A;
Energy = 1/2 * L_FCG * I0^2 * height0./height; %energy in the FCG

%update final values
handles.final_current.String = num2str(I(end)/1e6,3);
Bend = B(end);
if abs(Bend)<1
    handles.final_B.String = num2str(B(end),'%.2f');
elseif abs(Bend)<10
    handles.final_B.String = num2str(B(end),'%.1f');
else
    handles.final_B.String = num2str(B(end),'%.0f');
end
    
%handles.final_LFCG.String = num2str(I(end)/1e6,3);

semilogy(t*1e6,I/1e6,'-'), hold on
plot(t*1e6,B,'--'), hold on
plot(t*1e6,(4*pi*1e-7 * width.*height./depth + L_load)/1e-9,'--'), hold off
xlabel('t (\mus)');
ylabel('I (MA), B (T), L(nH)')
lh = legend('I (MA)','B (T)','L (nH)','location','northeast'); legend boxoff
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



function Inductance_FCG_Callback(hObject, eventdata, handles)
% hObject    handle to Inductance_FCG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Inductance_FCG as text
%        str2double(get(hObject,'String')) returns contents of Inductance_FCG as a double


% --- Executes during object creation, after setting all properties.
function Inductance_FCG_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Inductance_FCG (see GCBO)
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



function dwell_time_Callback(hObject, eventdata, handles)
% hObject    handle to text6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text6 as text
%        str2double(get(hObject,'String')) returns contents of text6 as a double


% --- Executes during object creation, after setting all properties.
function text6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Conductor_velocity_Callback(hObject, eventdata, handles)
% hObject    handle to Conductor_velocity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Conductor_velocity as text
%        str2double(get(hObject,'String')) returns contents of Conductor_velocity as a double


% --- Executes during object creation, after setting all properties.
function Conductor_velocity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Conductor_velocity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%this is the ode circuit solver
    function dy = fcgfuns(t,y,R,L_load,v,w,h0,d,flux)
        %R circuit and load  resistance
        %L circuit and load inductance
        %v velocity of the plate
        
        dy = zeros(1,1);
        
        current = y(1);
        
        h = h0 - v*t;
        A = w*h;
        B = flux./A;
        
        mu0 = 4*pi*1e-7;
        L_FCG = mu0 * w * h /d;

        L = L_FCG + L_load;
        
        emf_motional = -B*w*v;  %negative because work is being done on the circuit, so current should increase

        %circuit equations (it's important to get L as a function of time correctly)
        dy(1) = -(emf_motional + current*R)/L;
        


function B_initial_Callback(hObject, eventdata, handles)
% hObject    handle to B_initial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of B_initial as text
%        str2double(get(hObject,'String')) returns contents of B_initial as a double


% --- Executes during object creation, after setting all properties.
function B_initial_CreateFcn(hObject, eventdata, handles)
% hObject    handle to B_initial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function width_Callback(hObject, eventdata, handles)
% hObject    handle to width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of width as text
%        str2double(get(hObject,'String')) returns contents of width as a double


% --- Executes during object creation, after setting all properties.
function width_CreateFcn(hObject, eventdata, handles)
% hObject    handle to width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function I0_Callback(hObject, eventdata, handles)
% hObject    handle to I0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of I0 as text
%        str2double(get(hObject,'String')) returns contents of I0 as a double


% --- Executes during object creation, after setting all properties.
function I0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to I0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function depth_Callback(hObject, eventdata, handles)
% hObject    handle to depth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of depth as text
%        str2double(get(hObject,'String')) returns contents of depth as a double


% --- Executes during object creation, after setting all properties.
function depth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to depth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function height_Callback(hObject, eventdata, handles)
% hObject    handle to height (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of height as text
%        str2double(get(hObject,'String')) returns contents of height as a double


% --- Executes during object creation, after setting all properties.
function height_CreateFcn(hObject, eventdata, handles)
% hObject    handle to height (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
