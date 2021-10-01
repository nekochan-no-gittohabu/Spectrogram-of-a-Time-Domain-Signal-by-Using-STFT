%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% OUR CODE STARTS AT LINE 392 %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% AYÇA AYDOĞAN 2093375
% YİĞİT DİNCER 2188100


function varargout = myapp(varargin)
% MYAPP MATLAB code for myapp.fig
%      MYAPP, by itself, creates a new MYAPP or raises the existing
%      singleton*.
%
%      H = MYAPP returns the handle to a new MYAPP or the handle to
%      the existing singleton*.
%
%      MYAPP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MYAPP.M with the given input arguments.
%
%      MYAPP('Property','Value',...) creates a new MYAPP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before myapp_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to myapp_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help myapp

% Last Modified by GUIDE v2.5 15-Dec-2020 18:52:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @myapp_OpeningFcn, ...
                   'gui_OutputFcn',  @myapp_OutputFcn, ...
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


% --- Executes just before myapp is made visible.
function myapp_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to myapp (see VARARGIN)

% Choose default command line output for myapp
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes myapp wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = myapp_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function text2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in Hamming.
function Hamming_Callback(hObject, eventdata, handles)
% hObject    handle to Hamming (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Hamming


% --- Executes on button press in Hann.
function Hann_Callback(hObject, eventdata, handles)
% hObject    handle to Hann (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Hann


% --- Executes on button press in Tukey.
function Tukey_Callback(hObject, eventdata, handles)
% hObject    handle to Tukey (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Tukey


% --- Executes on button press in Cosine.
function Cosine_Callback(hObject, eventdata, handles)
% hObject    handle to Cosine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Cosine


% --- Executes on button press in Triangular.
function Triangular_Callback(hObject, eventdata, handles)
% hObject    handle to Triangular (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Triangular


% --- Executes on button press in Gaussian.
function Gaussian_Callback(hObject, eventdata, handles)
% hObject    handle to Gaussian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Gaussian


% --- Executes on button press in Blackman.
function Blackman_Callback(hObject, eventdata, handles)
% hObject    handle to Blackman (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Blackman


% --- Executes on button press in Kaiser.
function Kaiser_Callback(hObject, eventdata, handles)
% hObject    handle to Kaiser (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Kaiser


% --- Executes on button press in Chebwin.
function Chebwin_Callback(hObject, eventdata, handles)
% hObject    handle to Chebwin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Chebwin


% --- Executes on button press in Rectangular.
function Rectangular_Callback(hObject, eventdata, handles)
% hObject    handle to Rectangular (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Rectangular

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in Record.
function Record_Callback(hObject, eventdata, handles)
% hObject    handle to Record (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Record


% --- Executes on button press in Sinusoidal.
function Sinusoidal_Callback(hObject, eventdata, handles)
% hObject    handle to Sinusoidal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Sinusoidal


% --- Executes on button press in Square.
function Square_Callback(hObject, eventdata, handles)
% hObject    handle to Square (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Square


% --- Executes on button press in Chirp.
function Chirp_Callback(hObject, eventdata, handles)
% hObject    handle to Chirp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Chirp


% --- Executes on button press in Sawtooth.
function Sawtooth_Callback(hObject, eventdata, handles)
% hObject    handle to Sawtooth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Sawtooth


% --- Executes on button press in Audio.
function Audio_Callback(hObject, eventdata, handles)
% hObject    handle to Audio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Audio


% --- Executes during object creation, after setting all properties.
function Window_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function uibuttongroup2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uibuttongroup2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes on button press in timed.
function timed_Callback(hObject, eventdata, handles)
% hObject    handle to timed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of timed


% --- Executes on button press in freq.
function freq_Callback(hObject, eventdata, handles)
% hObject    handle to freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of freq


% --- Executes on button press in spec.
function spec_Callback(hObject, eventdata, handles)
% hObject    handle to spec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of spec


% --- Executes on button press in ChSq.
function ChSq_Callback(hObject, eventdata, handles)
% hObject    handle to ChSq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ChSq


% --- Executes on button press in SinSaw.
function SinSaw_Callback(hObject, eventdata, handles)
% hObject    handle to SinSaw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SinSaw


% --- Executes on button press in Sqsin.
function Sqsin_Callback(hObject, eventdata, handles)
% hObject    handle to Sqsin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Sqsin

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Fs = str2double(get(handles.edit1,'String'));
L = str2double(get(handles.edit2,'String'));
N = str2double(get(handles.edit3,'String'));
M = str2double(get(handles.edit4,'String'));
if M > N
    msgbox('M should be smaller than N');
end

x = zeros(1,L); %initialize x, the vector that holds sound data
time=0;

%%%%%%%%% PLOTS %%%%%%%%%%
timed = 0;
freq = 0;
spec = 0;
if get(handles.timed,'Value') == get(handles.timed,'Max')
    timed = true;
end

if get(handles.freq,'Value') == get(handles.freq,'Max')
    freq = true;
end

if get(handles.spec,'Value') == get(handles.spec,'Max')
    spec = true;
end

%%%%%%%%%% WINDOW TYPE %%%%%%%%%%
if get(handles.Hamming,'Value') == get(handles.Hamming,'Max')
        win = hamming(N).';

elseif get(handles.Hann,'Value') == get(handles.Hann,'Max')
        win = hanning(N).';  
    
elseif get(handles.Tukey,'Value') == get(handles.Tukey,'Max')
        win=tukeywin(N).'; 
    
elseif get(handles.Cosine,'Value') == get(handles.Cosine,'Max')
        % As there is no cosine filter in MATLAB, 
        % we created our own filter
        for n = 1:N
            win(n) = cos(2*pi*500*n);  
        end
    
elseif get(handles.Triangular,'Value') == get(handles.Triangular,'Max')
        win = triang(N).'; 
    
elseif get(handles.Gaussian,'Value') == get(handles.Gaussian,'Max')
        win = gausswin(N).';
    
elseif get(handles.Blackman,'Value') == get(handles.Blackman,'Max')
        win = blackman(N).';
    
elseif get(handles.Kaiser,'Value') == get(handles.Kaiser,'Max')
        beta = 5;
        win = kaiser(N,beta).';
    
elseif get(handles.Chebwin,'Value') == get(handles.Chebwin,'Max')
        win = chebwin(N).';
         
elseif get(handles.Rectangular,'Value') == get(handles.Rectangular,'Max')
        win = ones(1,N);
    
end


%%%%%%%% INPUT TYPE %%%%%%%%%%%
if get(handles.Record,'Value') == get(handles.Record,'Max')
    recObj = audiorecorder(Fs,16,1); %record the sound the Fs Hz
    disp('Start speaking.');
    recordblocking(recObj, (L/Fs));
    disp('End of Recording.');
    x = getaudiodata(recObj);
    time = 1;
    
elseif get(handles.Sinusoidal,'Value') == get(handles.Sinusoidal,'Max')
    t = 0:1/Fs:(L-1)/Fs;
    x = sin(2*pi*500*t);

elseif get(handles.Square,'Value') == get(handles.Square,'Max')
    t = linspace(0,L,L-1)';
    x = square(t);
    
elseif get(handles.Chirp,'Value') == get(handles.Chirp,'Max')
   t = 0:1/Fs:(L-1)/Fs;
   x = chirp(t,0,(L-1)/Fs,250);

elseif get(handles.Sawtooth,'Value') == get(handles.Sawtooth,'Max')
    t = 0:1/Fs:(L-1)/Fs;
    x = sawtooth(2*pi*50*t);

    
elseif get(handles.Audio,'Value') == get(handles.Audio,'Max')
    [z, F] = audioread('andean-flute.wav');  % Default Fs=44100Hz
    % Estimating the Fs to the desired value
    m = min(L-1,int32((length(z)+1)*Fs/44100));
    disp(m);
    for n=0:m-1;
    x(n+1) = z(int32(44100*n/Fs + 1));      
    end
    time = 1;
    
    
elseif get(handles.ChSq,'Value') == get(handles.ChSq,'Max')
    t = linspace(0,L,L)';
    z = square(t).';
    t = 0:1/Fs:(L-1)/Fs;
    y = chirp(t,0,(L-1)/Fs,250);
    x = y+z;
    
elseif get(handles.SinSaw,'Value') == get(handles.SinSaw,'Max')
    t = 0:1/Fs:(L-1)/Fs;
    y = sin(2*pi*1000*t);
    z = sin(2*pi*4000*t);
    x = y + z;
    
elseif get(handles.Sqsin,'Value') == get(handles.Sqsin,'Max')
    % Square followed by a sin wave
    t = 0:1/Fs:(L-2)/(2*Fs);
    z = 0.5*sin(2*pi*500*t);
    disp(length(z));
    t = linspace(0,(L/2)+1,L/2)';
    y = square(t).';
    disp(length(y));
    x = [y, z];
    
end
    

Y = myspectrogram(Fs, L, N, M, win, x, time, timed, freq, spec);
set(handles.text2,'string',Y);

function [y] = myspectrogram(Fs, L, N, M, win, x, time, timed, freq, spec);
    y='The results are shown';
    
    % Play the sampled sound signal
    if time ==1
        soundsc(x, Fs);
    elseif time == 0
        soundsc(x);
    end
    
    if timed == true
        %Time domain plot of the signal
        figure
        plot(x);    
        title('Time domain plot')
        xlabel('n, samples')
        ylabel('x[n]')
    end

    if freq == true
        [X, w] = freqz(x);
        %Frequency domain plot of the signal
        figure
        plot(abs(X));
        title('Frequency domain plot')
        xlabel('w')
        ylabel('X(w)')
    end
    
    % Normalize the signal
    x = x.'/max(abs(x));

    % Make the signal a row vector
    x = x(:).';

    % Number of segments (frames) the signal is divided to.
    K = floor((L-M)/(N-M)); 

    %STFT
    X = zeros(N,K);
    S = zeros(N,K);

    for k=1:K
        X(:,k) = x((k-1)*(N-M)+1:k*N - (k-1)*M).*win;
        S(:,k) = fft(X(:,k));
    end
    S = S(1:N/2,:); % STFT is saved in matrix S 

    %Time axis is determined differently for computer generated data
    %like sinusoidal, square wave etc. So there are 2 cases;
    if time == 1        % For recording and audio data
        t = (N-1)/Fs + (0:K-1)*(N-M)/Fs;
    elseif time == 0    % For computer generated data.
        t =(0:K-1)*L/(K*10000);
    end

    %Frequency points in Hz
    f = (0:N/2-1)*Fs/N;  

    if spec == true
        %Plot the spectogram
        h = figure('Name','STFT - Spectogram');
        colormap('jet');

        [T,F] = meshgrid(t,f/1000); % f in KHz.
        % For spectogram, square magnitude of STFT is used
        surface(T,F,10*log10(abs(S.^2) + eps),'EdgeColor','none');

        axis tight;
        grid on;
        title(['   Fs: ',num2str(Fs),', Signal Length: ',num2str(L),', Window Length: ', num2str(N),', Overlap: ', num2str(M), ' samples.']);
        xlabel('Time (sec)');
        ylabel('Frequency (KHz)');
        colorbar('Limits',[-80, 40]);
        cbar_handle = findobj(h,'tag','Colorbar');
        set(get(cbar_handle,'YLabel'),'String','(dB)','Rotation',0);
        zlim([-80 40]);
    end
