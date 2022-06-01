function varargout = GROUP_GUI(varargin)
% GROUP_GUI MATLAB code for GROUP_GUI.fig
%      GROUP_GUI, by itself, creates a new GROUP_GUI or raises the existing
%      singleton*.
%
%      H = GROUP_GUI returns the handle to a new GROUP_GUI or the handle to
%      the existing singleton*.
%
%      GROUP_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GROUP_GUI.M with the given input arguments.
%
%      GROUP_GUI('Property','Value',...) creates a new GROUP_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GROUP_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GROUP_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GROUP_GUI

% Last Modified by GUIDE v2.5 04-Feb-2021 22:39:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GROUP_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @GROUP_GUI_OutputFcn, ...
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


% --- Executes just before GROUP_GUI is made visible.
function GROUP_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GROUP_GUI (see VARARGIN)

% Choose default command line output for GROUP_GUI
handles.output = hObject;

% set(hObject,'toolbar','figure');
set(gcf,'HandleVis','off');
set(findall(handles.uipanel1, '-property', 'Enable'), 'Enable', 'off');
set(findall(handles.uipanel2, '-property', 'Enable'), 'Enable', 'off');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GROUP_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GROUP_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Read Data Hub.
function ReadXls_Callback(hObject, eventdata, handles)
% hObject    handle to ReadXls (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clear global group_result;
if ismember('control',get(gcf,'currentModifier'))  % 按下ctrl再点击，进行调试，打开Load_Mat_data.m
    edit Load_Mat_data;
    return;
end
if strcmpi(get(gcf,'CurrentCharacter'),'e') % 按下e再点击，打开excel
    if exist('Z:\Data\MOOG\Results\Result_Lwh.xlsm')
        winopen('Z:\Data\MOOG\Results\Result_Lwh.xlsm');
    elseif exist('Z:\Data\MOOG\Results\Result_LJY.xlsm')
        winopen('Z:\Data\MOOG\Results\Result_LJY.xlsm');
    end
    return;
end
    
set(handles.num_entries,'string',0); drawnow;
xlRange{1} = 'A1:DG3000';
xlRange{2} = 'A1:DG3000';
monkey_choose = str2num(get(handles.sheetN,'string')); % 第一个sheet默认ringbell，第二个sheet默认Arthas
% monkey_choose = 1;

if exist('Z:\BaiduNetdiskWorkspace\data\Result_Lwh.xlsm')
    handles.XlsData = ReadXls('Z:\BaiduNetdiskWorkspace\data\Result_Lwh.xlsm',monkey_choose,str2num(get(handles.headerN,'string')),xlRange);
elseif exist('Z:\Data\MOOG\Results\Result_Lwh.xlsm')
    handles.XlsData = ReadXls('Z:\Data\MOOG\Results\Result_Lwh.xlsm',monkey_choose,str2num(get(handles.headerN,'string')),xlRange);
elseif exist('Z:\Data\MOOG\Results\Result_LJY.xlsm')
    handles.XlsData = ReadXls('Z:\Data\MOOG\Results\Result_LJY.xlsm',monkey_choose,str2num(get(handles.headerN,'string')),xlRange);
end

handles.N = size(handles.XlsData.num,1);

% ----------------------------------------------
% data_hub = matching_file_backup202102(handles.XlsData); % old code,
% combine load file and matching file, now separate, Lwh 202001
data_hub = Load_Mat_data(handles.XlsData,monkey_choose);
handles.data_hub = data_hub;
guidata(hObject,handles);
% ----------------------------------------------

set(handles.num_entries,'string',num2str(handles.N));

set(findall(handles.uipanel2, '-property', 'Enable'), 'Enable', 'on');
guidata(hObject,handles);


% --- Executes on button press in matchingfile.
function matchingfile_Callback(hObject, eventdata, handles)
% hObject    handle to matchingfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ismember('control',get(gcf,'currentModifier'))
    edit matching_file;
    return;
end

% ----------------------------------------------
data_hub = handles.data_hub;
loadM_data = matching_file(data_hub);
handles.loadM_data = loadM_data;
guidata(hObject,handles);

cc = get(gcbf,'color');
set(gcbf,'color','y');
pause(0.1);
set(gcbf,'color',cc);
drawnow;
% ----------------------------------------------

set(handles.num_entries,'string',num2str(handles.N));

set(findall(handles.uipanel2, '-property', 'Enable'), 'Enable', 'on');
guidata(hObject,handles);


% --- Executes on button press in classification.
function classification_Callback(hObject, eventdata, handles)
% hObject    handle to classification (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ismember('control',get(gcf,'currentModifier'))
    edit Group_Classification;
    return;
end

set([handles.level1 handles.level2 handles.fig_selected_txt],'string',[]);
set([handles.num_all_units],'string',[]);
handles.level1_selected = [];
handles.level2_selected = [];

loadM_data = handles.loadM_data;
p_cri = [];
p_cri = str2num(get(handles.p_cri,'string'));
loadM_data = Group_Classification(loadM_data,p_cri);
handles.loadM_data = loadM_data;

cc = get(gcbf,'color');
set(gcbf,'color','y');
pause(0.1);
set(gcbf,'color',cc);
drawnow;

% hs = msgbox('Classification Compleeted');
% ht = findobj(hs, 'Type', 'text');
% set(ht, 'FontSize', 10, 'Unit', 'normal');
% set(hs, 'Resize', 'on');

set(findall(handles.uipanel1, '-property', 'Enable'), 'Enable', 'on');
set(findall(handles.cell_counter, '-property', 'Enable'), 'Enable', 'on');
set(handles.close_figs,'enable','on');
guidata(hObject,handles);


% --- Executes on button press in Microstimulation.
function Microstimulation_Callback(hObject, eventdata, handles)
% hObject    handle to Microstimulation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if ismember('control',get(gcf,'currentModifier')) % ����ctrlֱ�ӱ༭code
    edit Group_MicroStim_step;
    return;
end

set([handles.level1 handles.level2 handles.fig_selected_txt],'string',[]);
set([handles.num_all_units],'string',[]);
handles.protocol_selected = 'Stim';
handles.level1_selected = [];
handles.level2_selected = [];

set(findall(handles.uipanel1, '-property', 'Enable'), 'Enable', 'on');
set(findall(handles.cell_counter, '-property', 'Enable'), 'Enable', 'on');

data = handles.XlsData;
loadM_data = handles.loadM_data;

%^ auto run 'Group_SpiT'
if ~isfield(handles,'function_handles_spi')
    disp('****** Run SpiT ******');
    p_cri = [];
    p_cri = str2num(get(handles.p_cri,'string'));
    function_handles_spi = Group_SpiT(loadM_data,handles,p_cri);
    handles.function_handles_spi = function_handles_spi;
    disp('****** Done SpiT ******'); 
end

disp('****** Run Microstimulation ******');
function_handles_stim = Group_MicroStim_step(loadM_data,handles);

handles.function_handles_stim = function_handles_stim;
guidata(hObject,handles);
load_function_handles(hObject,handles,handles.function_handles_stim); % Update possible function_handles
disp('****** Done Microstimulation ******');

% --- Executes on button press in SpiT.
function SpiT_Callback(hObject, eventdata, handles)
disp('****** Run SpiT ******');
% hObject    handle to SpiT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ismember('control',get(gcf,'currentModifier'))
    edit Group_SpiT;
    return;
end

set([handles.level1 handles.level2 handles.fig_selected_txt],'string',[]);
set([handles.num_all_units],'string',[]);
handles.protocol_selected = 'SpiT';
handles.level1_selected = [];
handles.level2_selected = [];

set(handles.close_figs,'enable','on');

set(findall(handles.uipanel1, '-property', 'Enable'), 'Enable', 'on');
set(findall(handles.cell_counter, '-property', 'Enable'), 'Enable', 'on');

data = handles.XlsData;
loadM_data = handles.loadM_data;
p_cri = [];
p_cri = str2num(get(handles.p_cri,'string'));

function_handles_spi = Group_SpiT(loadM_data,handles,p_cri);
handles.function_handles_spi = function_handles_spi;

guidata(hObject,handles);
load_function_handles(hObject,handles,handles.function_handles_spi); % Update possible function_handles
disp('****** Done SpiT ******');


% --- Executes on button press in CueS.
function CueS_Callback(hObject, eventdata, handles)
% hObject    handle to LIP_HD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ismember('control',get(gcf,'currentModifier'))
    edit Group_CueS;
    return;
end

set([handles.level1 handles.level2 handles.fig_selected_txt],'string',[]);
set([handles.num_all_units],'string',[]);
handles.level1_selected = [];
handles.level2_selected = [];
handles.protocol_selected = 'CueS';

set(findall(handles.uipanel1, '-property', 'Enable'), 'Enable', 'on');
set(findall(handles.cell_counter, '-property', 'Enable'), 'Enable', 'on');

data = handles.XlsData;
loadM_data = handles.loadM_data;
function_handles_cues = Group_CueS(loadM_data,handles);

handles.function_handles_cues = function_handles_cues;
guidata(hObject,handles);
load_function_handles(hObject,handles,function_handles_cues);

% % Lwh, 20220507
% function load_function_handles(hObject, handles)
% % Clear previous values
% set([handles.level1 handles.level2],'value',[]);
% 
% % Update possible function_handles
% handles.level1_n = size(handles.function_handles,1)-1; % The last one reserved for functions that is accessible to GROUP_GUI 
%                                                        % but not to manual selection. (cell_selection undapte) HH20150723
%                                                        % for this one:  'NoShow',{@cell_selection};                                                       
% 
% handles.level1_txt = cell(handles.level1_n,1);
% 
% handles.fig_all = [];
% for i = 1:handles.level1_n
%     handles.level1_txt{i} = handles.function_handles{i,1};
%     level2_len = size(handles.function_handles{i,2},1);
%     handles.fig_all = [handles.fig_all; repmat(i,level2_len,1) (1:level2_len)'];
% end
% set(handles.level1,'string',handles.level1_txt);
% 
% % Clear fig_selected_txt cache
% handles.fig_selected = [];
% guidata(hObject,handles);
% 
% cc = get(gcbf,'color');
% set(gcbf,'color','y');
% pause(0.1);
% set(gcbf,'color',cc);
% drawnow;

% seperate function handle for different task
% Lwh 20220507
% handles_all: same as original "handles", use for level1, level2, and so on
% handles_this: selected task for specific function_handles
function load_function_handles(hObject, handles_all, handles_this)
% Clear previous values
set([handles_all.level1 handles_all.level2],'value',[]);

% Update possible function_handles
handles_all.level1_n = size(handles_this,1)-1; % The last one reserved for functions that is accessible to GROUP_GUI 
                                                       % but not to manual selection. (cell_selection undapte) HH20150723
                                                       % for this one:  'NoShow',{@cell_selection};
                                                      
% Lwh20220507 : handles_all.level1_n,level1_txt,fig_selected: 
% use the same location for saving the level1_n
% and other for differenttask
                                                       
handles_all.level1_txt = cell(handles_all.level1_n,1);

handles_all.fig_all = [];
for i = 1:handles_all.level1_n
    handles_all.level1_txt{i} = handles_this{i,1};
    level2_len = size(handles_this{i,2},1);
    handles_all.fig_all = [handles_all.fig_all; repmat(i,level2_len,1) (1:level2_len)'];
end
set(handles_all.level1,'string',handles_all.level1_txt);

% Clear fig_selected_txt cache
handles_all.fig_selected = [];
guidata(hObject,handles_all);

cc = get(gcbf,'color');
set(gcbf,'color','y');
pause(0.1);
set(gcbf,'color',cc);
drawnow;

% --- Executes on selection change in level1.
function level1_Callback(hObject, eventdata, handles)
% hObject    handle to level1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns level1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from level1


% if ~isfield(handles,'function_handles'); return; end   % Lwh, 20220507
% seperate function handle for different task
if strcmpi(handles.protocol_selected,'SpiT')
    function_handles_this = handles.function_handles_spi;
elseif strcmpi(handles.protocol_selected,'stim')
    function_handles_this = handles.function_handles_stim;
elseif strcmpi(handles.protocol_selected,'cues')
    function_handles_this = handles.function_handles_cues;
else
%     return
keyborad
end

select_item = get(hObject,'value');

handles.level1_selected = select_item;

set(handles.level2,'value',[]);
handles.level2_selected = [];
handles.level2_txt = [];

if length(select_item) == 1
    %     handles.level2_n = size(handles.function_handles{select_item,2},1); % Lwh, 20220507
    handles.level2_n = size(function_handles_this{select_item,2},1);
    for i = 1:handles.level2_n
        %         handles.level2_txt{i} = handles.function_handles{select_item,2}{i,1}; % Lwh, 20220507
        handles.level2_txt{i} = function_handles_this{select_item,2}{i,1};
    end
end

set(handles.level2,'string',handles.level2_txt);
set(handles.debug,'enable','off');

guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function level1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to level1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in level2.
function level2_Callback(hObject, eventdata, handles)
% hObject    handle to level2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns level2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from level2

handles.level2_selected = get(hObject,'value');
if length(handles.level2_selected) == 1
    set(handles.debug,'enable','on');
else
    set(handles.debug,'enable','off');
end

guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function level2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to level2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in fig_add.
function fig_add_Callback(hObject, eventdata, handles)
% hObject    handle to fig_add (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fig_add = [];

if ~isempty(handles.level2_selected) % Add selected level 2 entries
    fig_add = [repmat(handles.level1_selected,length(handles.level2_selected),1) handles.level2_selected'];    
else % Add all level 2 entries for each level 1 entry
    for j = 1:length(handles.level1_selected)
        level2_len = size(handles.function_handles{handles.level1_selected(j),2},1);
        fig_add = [fig_add ; repmat(handles.level1_selected(j),level2_len,1) (1:level2_len)'];
    end
end

handles.fig_selected = [handles.fig_selected; fig_add];
handles.fig_selected = sortrows(munique(handles.fig_selected));

guidata(hObject,handles);
update_fig_selected(handles);


function update_fig_selected(handles)
txt = [];
for i = 1:size(handles.fig_selected,1)
    l1 = handles.fig_selected(i,1);
    l2 = handles.fig_selected(i,2);
    txt{i} = sprintf('%s | %s\n',handles.function_handles{l1,1},handles.function_handles{l1,2}{l2,1});
end
set(handles.fig_selected_txt,'string',txt);
set(handles.fig_selected_txt,'value',[]);


% --- Executes on button press in fig_remove_all.
function fig_remove_all_Callback(hObject, eventdata, handles)
% hObject    handle to fig_remove_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.fig_selected = [];
guidata(hObject,handles);
update_fig_selected(handles);


% --- Executes on selection change in fig_selected_txt.
function fig_selected_txt_Callback(hObject, eventdata, handles)
% hObject    handle to fig_selected_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns fig_selected_txt contents as cell array
%        contents{get(hObject,'Value')} returns selected item from fig_selected_txt


% --- Executes during object creation, after setting all properties.
function fig_selected_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fig_selected_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sheetN_Callback(hObject, eventdata, handles)
% hObject    handle to sheetN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sheetN as text
%        str2double(get(hObject,'String')) returns contents of sheetN as a double


% --- Executes during object creation, after setting all properties.
function sheetN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sheetN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function headerN_Callback(hObject, eventdata, handles)
% hObject    handle to headerN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of headerN as text
%        str2double(get(hObject,'String')) returns contents of headerN as a double


% --- Executes during object creation, after setting all properties.
function headerN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to headerN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% remove by Lwh 202102
% % --- Executes during object creation, after setting all properties.
% function monkey_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to monkey (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: popupmenu controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end


% --- Executes on button press in close_figs.
function close_figs_Callback(hObject, eventdata, handles)
% hObject    handle to close_figs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h_figs = findall(0,'type','figure');
h_figs(h_figs == gcbf) = [];
close(h_figs);


% --- Executes on button press in fig_add_all.
function fig_add_all_Callback(hObject, eventdata, handles)
% hObject    handle to fig_add_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.fig_selected = handles.fig_all;
guidata(hObject,handles);
update_fig_selected(handles);

% --- Executes on button press in fig_reverse.
function fig_reverse_Callback(hObject, eventdata, handles)
% hObject    handle to fig_reverse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
all = handles.fig_all;
sel = handles.fig_selected;

if isempty(sel)
    handles.fig_selected = all;
else
    [~,reverse_ind] = setdiff(all , sel ,'rows');
    handles.fig_selected = all(reverse_ind,:);
end
guidata(hObject,handles);
update_fig_selected(handles);

% --- Executes on button press in fig_remove.
function fig_remove_Callback(hObject, eventdata, handles)
% hObject    handle to fig_remove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
to_remove = get(handles.fig_selected_txt,'value');
handles.fig_selected(to_remove,:) = [];
guidata(hObject,handles);
update_fig_selected(handles);

% --- Executes on button press in go.
function go_Callback(hObject, eventdata, handles)
% hObject    handle to go (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.fig_selected)
    if isempty(handles.level2_selected)
        to_do = handles.fig_all(ismember(handles.fig_all(:,1) ,handles.level1_selected),:);
    else
        to_do = [repmat(handles.level1_selected,length(handles.level2_selected),1) handles.level2_selected'];
    end
%     handles.fig_selected = to_do;
%     guidata(hObject,handles);
%     update_fig_selected(handles);
%    
else
    to_do = handles.fig_selected;
end

% Do analysis
% seperate function handle for different task
if strcmpi(handles.protocol_selected,'SpiT')
    function_handles_this = handles.function_handles_spi;
elseif strcmpi(handles.protocol_selected,'stim')
    function_handles_this = handles.function_handles_stim;
elseif strcmpi(handles.protocol_selected,'cues')
    function_handles_this = handles.function_handles_cues;
else
    keyboard
end

cell_num = str2num(get(handles.num_all_units,'String'));

if isempty(cell_num) || cell_num(end,1) > 0 % More than one cell included
    for i = 1:size(to_do,1)
        feval(function_handles_this{to_do(i,1),2}{to_do(i,2),2},false);  % Not debugging
    end
else
    beep; pause(0.1); beep;
end

figure(gcbf);


% --- Executes during object creation, after setting all properties.
function go_CreateFcn(hObject, eventdata, handles)
% hObject    handle to go (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text1.
function text1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to text1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guide GROUP_GUI;
edit GROUP_GUI;

% remove by Lwh 202102
% % --- Executes on selection change in monkey.
% function monkey_Callback(hObject, eventdata, handles)
% % hObject    handle to monkey (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: contents = cellstr(get(hObject,'String')) returns monkey contents as cell array
% %        contents{get(hObject,'Value')} returns selected item from monkey


% --- Executes on button press in debug.
function debug_Callback(hObject, eventdata, handles)
% hObject    handle to debug (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% seperate function handle for different task
if strcmpi(handles.protocol_selected,'SpiT')
    function_handles_this = handles.function_handles_spi;
elseif strcmpi(handles.protocol_selected,'stim')
    function_handles_this = handles.function_handles_stim;
elseif strcmpi(handles.protocol_selected,'cues')
    function_handles_this = handles.function_handles_cues;
else
    keyboard
end

if length(handles.level2_selected) == 1
    feval(function_handles_this{handles.level1_selected,2}{handles.level2_selected,2},true);  % Debug mode
end


% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    dbquit('all');
catch
end

% --- Executes on button press in Ringbell_data.
function Ringbell_data_Callback(hObject, eventdata, handles)
% hObject    handle to Ringbell_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% seperate function handle for different task
if strcmpi(handles.protocol_selected,'SpiT')
    function_handles_this = handles.function_handles_spi;
elseif strcmpi(handles.protocol_selected,'stim')
    function_handles_this = handles.function_handles_stim;
elseif strcmpi(handles.protocol_selected,'cues')
    function_handles_this = handles.function_handles_cues;
else
    keyboard
end

typical_cell_selection_num = get(handles.t_criterion,'value'); 
feval(function_handles_this{end,2}{1},typical_cell_selection_num);  % Update cell_selection
% ���ص�ǰ�����Group��������һ��end��'NoShow',{@cell_selection}��ִ��cell_selection

% Hint: get(hObject,'Value') returns toggle state of Ringbell_data


% --- Executes on button press in Arthas_data.
function Arthas_data_Callback(hObject, eventdata, handles)
% hObject    handle to Arthas_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% seperate function handle for different task
if strcmpi(handles.protocol_selected,'SpiT')
    function_handles_this = handles.function_handles_spi;
elseif strcmpi(handles.protocol_selected,'stim')
    function_handles_this = handles.function_handles_stim;
elseif strcmpi(handles.protocol_selected,'cues')
    function_handles_this = handles.function_handles_cues;
else
    keyboard
end

typical_cell_selection_num = get(handles.t_criterion,'value'); 
feval(function_handles_this{end,2}{1},typical_cell_selection_num);  % Update cell_selection

% Hint: get(hObject,'Value') returns toggle state of Arthas_data


% --- Executes on selection change in t_criterion.
function t_criterion_Callback(hObject, eventdata, handles)  % HH20160912
% hObject    handle to t_criterion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% seperate function handle for different task
if strcmpi(handles.protocol_selected,'SpiT')
    function_handles_this = handles.function_handles_spi;
elseif strcmpi(handles.protocol_selected,'stim')
    function_handles_this = handles.function_handles_stim;
elseif strcmpi(handles.protocol_selected,'cues')
    function_handles_this = handles.function_handles_cues;
else
    keyboard
end

typical_cell_selection_num = get(handles.t_criterion,'value'); 
feval(function_handles_this{end,2}{1},typical_cell_selection_num);  % Update cell_selection

% Hints: contents = cellstr(get(hObject,'String')) returns t_criterion contents as cell array
%        contents{get(hObject,'Value')} returns selected item from t_criterion


% --- Executes during object creation, after setting all properties.
function t_criterion_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t_criterion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% remove by Lwh 202102
% % --- Executes on button press in transparent.
% function transparent_Callback(hObject, eventdata, handles)
% % hObject    handle to transparent (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hint: get(hObject,'Value') returns toggle state of transparent
% % LIP_HD_Callback(hObject, eventdata, handles);
% SpiT_Callback(hObject, eventdata, handles);


function p_cri_Callback(hObject, eventdata, handles)
% hObject    handle to p_cri (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of p_cri as text
%        str2double(get(hObject,'String')) returns contents of p_cri as a double


% --- Executes during object creation, after setting all properties.
function p_cri_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p_cri (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in fine_task.
function fine_task_Callback(hObject, eventdata, handles)
% hObject    handle to fine_task (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% seperate function handle for different task
if strcmpi(handles.protocol_selected,'SpiT')
    function_handles_this = handles.function_handles_spi;
elseif strcmpi(handles.protocol_selected,'stim')
    function_handles_this = handles.function_handles_stim;
elseif strcmpi(handles.protocol_selected,'cues')
    function_handles_this = handles.function_handles_cues;
else
    keyboard
end

typical_cell_selection_num = get(handles.t_criterion,'value'); 
feval(function_handles_this{end,2}{1},typical_cell_selection_num);  % Update task type selection
% Hint: get(hObject,'Value') returns toggle state of fine_task


% --- Executes on button press in coarse_task.
function coarse_task_Callback(hObject, eventdata, handles)
% hObject    handle to coarse_task (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% seperate function handle for different task
if strcmpi(handles.protocol_selected,'SpiT')
    function_handles_this = handles.function_handles_spi;
elseif strcmpi(handles.protocol_selected,'stim')
    function_handles_this = handles.function_handles_stim;
elseif strcmpi(handles.protocol_selected,'cues')
    function_handles_this = handles.function_handles_cues;
else
    keyboard
end

typical_cell_selection_num = get(handles.t_criterion,'value'); 
feval(function_handles_this{end,2}{1},typical_cell_selection_num);  % Update task type selection
% Hint: get(hObject,'Value') returns toggle state of coarse_task


% --- Executes on button press in two_target.
function two_target_Callback(hObject, eventdata, handles)
% hObject    handle to two_target (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% seperate function handle for different task
if strcmpi(handles.protocol_selected,'SpiT')
    function_handles_this = handles.function_handles_spi;
elseif strcmpi(handles.protocol_selected,'stim')
    function_handles_this = handles.function_handles_stim;
elseif strcmpi(handles.protocol_selected,'cues')
    function_handles_this = handles.function_handles_cues;
else
    keyboard
end

typical_cell_selection_num = get(handles.t_criterion,'value'); 
feval(function_handles_this{end,2}{1},typical_cell_selection_num);  % Update task type selection
% Hint: get(hObject,'Value') returns toggle state of two_target


% --- Executes on button press in four_target.
function four_target_Callback(hObject, eventdata, handles)
% hObject    handle to four_target (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% seperate function handle for different task
if strcmpi(handles.protocol_selected,'SpiT')
    function_handles_this = handles.function_handles_spi;
elseif strcmpi(handles.protocol_selected,'stim')
    function_handles_this = handles.function_handles_stim;
elseif strcmpi(handles.protocol_selected,'cues')
    function_handles_this = handles.function_handles_cues;
else
    keyboard
end

typical_cell_selection_num = get(handles.t_criterion,'value'); 
feval(function_handles_this{end,2}{1},typical_cell_selection_num);  % Update task type selection
% Hint: get(hObject,'Value') returns toggle state of four_target
