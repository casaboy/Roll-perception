% matching spiral file and cueswtiching (2AFC/4AFC) and u-stim. task  Lwh 201906 
% update for two monkey Lwh20200727
% seperate Load_Mat_data and matching_file  Lwh 202102

function data_hub = Load_Mat_data(data,monkey_choose)

%% Get excel data: cell_info
num = data.num;
txt = data.txt;
raw = data.raw;
header = data.header;

data_num = size(num,1);

cell_info.Date = txt(:,header.Date);
cell_info.Area = txt(:,header.Area);
cell_info.FileNo = txt(:,header.FileNo);
cell_info.Protocol = txt(:,header.Protocol);
cell_info.Monkey = num(:,header.Monkey);
cell_info.Session = num(:,header.Session);
cell_info.Hemisphere = num(:,header.Hemisphere);
cell_info.xloc = num(:,header.Xloc);
cell_info.yloc = num(:,header.Yloc);
cell_info.Depth = num(:,header.Depth);
cell_info.Chan = num(:,header.Chan1);
cell_info.RF = num(:,header.X : header.H);
cell_info.coherence = num(:,header.coherence);
cell_info.stimAmp = num(:,header.stimAmp);
cell_info.guidetube = num(:,header.guidetube);
cell_info.offset = num(:,header.offset);

eye_x = num(:,header.Eyex);
eye_y = num(:,header.Eyey);
eye_x(isnan(eye_x)) = 0;
eye_y(isnan(eye_y)) = 0;
cell_info.eye_x = eye_x;
cell_info.eye_y = eye_y;

note_MUSU = txt(:,header.Note);
MU_mark = strfind(note_MUSU,'MU','ForceCellOutput',true);
MU_id = logical(cellfun(@length,MU_mark));

SU_mark = strfind(note_MUSU,'SU','ForceCellOutput',true);
SU_id_temp = cellfun(@length,SU_mark);
SU_id = logical(~MU_id & SU_id_temp);

% SU and MU marker; SU=1, MU=2, other=0
SUMU = zeros(data_num,1);
SUMU(SU_id) = 1;
SUMU(MU_id) = 2;

cell_info.SUMU = SUMU;

unique_monkey = unique(cell_info.Monkey);
unique_Protocol = unique(cell_info.Protocol);

% recode Protocol:用数字表示常用protocol
cell_info.Protocol_code = zeros(data_num,1);
cell_info.Protocol_code(strcmpi(cell_info.Protocol,'SpiT')) = 1;
cell_info.Protocol_code(strcmpi(cell_info.Protocol,'CueS4T')) = 2;
cell_info.Protocol_code(strcmpi(cell_info.Protocol,'Honly')) = 3; % two target，单独的H
cell_info.Protocol_code(strcmpi(cell_info.Protocol,'Ronly')) = 4; % two target，单独的R
cell_info.Protocol_code(strcmpi(cell_info.Protocol,'Coarse4T')) = 5;
cell_info.Protocol_code(strcmpi(cell_info.Protocol,'HCoarse2T')) = 6;
cell_info.Protocol_code(strcmpi(cell_info.Protocol,'RCoarse2T')) = 7;
cell_info.Protocol_code(strcmpi(cell_info.Protocol,'stimCueS')) = 8;
cell_info.Protocol_code(strcmpi(cell_info.Protocol,'stimHonly')) = 9;
cell_info.Protocol_code(strcmpi(cell_info.Protocol,'stimRonly')) = 10;
cell_info.Protocol_code(strcmpi(cell_info.Protocol,'stimCo')) = 11;
cell_info.Protocol_code(strcmpi(cell_info.Protocol,'stimCoH')) = 12;
cell_info.Protocol_code(strcmpi(cell_info.Protocol,'stimCoR')) = 13;
cell_info.Protocol_code(strcmpi(cell_info.Protocol,'afterStimSpiT')) = 14;
cell_info.Protocol_code(strcmpi(cell_info.Protocol,'RFplot')) = 15;
cell_info.Protocol_code(strcmpi(cell_info.Protocol,'MemSac')) = 16;

cell_info.Protocol_code(strcmpi(cell_info.Protocol,'stimRonly(MultiAmp)')) = 10;
cell_info.Protocol_code(strcmpi(cell_info.Protocol,'stimCoR(MultiAmp)')) = 13;

unique_Protocol_code = unique(cell_info.Protocol_code);  

if sum(unique_Protocol_code==0) > 0
    disp('CHECK EXECL PROTOCOL NAME !!!!');
    keyboard
end

% if sum(unique_Protocol_code==0) > 0 % 0 为没有protocol,例如MemSac
% find([cell_info.Protocol_code]==0)
%     keyboard
% end


% 文件名例子 m7c117r1
for i = 1:data_num
    c_loc = find(cell_info.FileNo{i}=='c'); % 找到c的位置
    r_loc = find(cell_info.FileNo{i}=='r'); % 找到r的位置
    dot_loc = find(cell_info.FileNo{i}=='.'); % 找到r的位置
    
    if cell_info.Monkey(i) == 7 % ringbell
        try
            cell_index(i,1) = str2num(cell_info.FileNo{i}(c_loc+1:r_loc-1));
            run_index(i,1) = str2num(cell_info.FileNo{i}(r_loc+1:dot_loc-1));
        catch
            cell_index(i,1) = nan; % only RF, no file
            run_index(i,1) = nan;
        end
        cell_index_all_monkey(i,1) = cell_index(i,1);
    elseif cell_info.Monkey(i) == 16 % Arthas
        cell_index(i,2) = str2num(cell_info.FileNo{i}(c_loc+1:r_loc-1));
        run_index(i,2) = str2num(cell_info.FileNo{i}(r_loc+1:dot_loc-1));
        cell_index_all_monkey(i,1) = cell_index(i,2);
    end
    
    % 补全缺省项
    % 1. coherence，monkey7：SpiT默认100，电刺激默认15，task默认15/30
    if isnan(cell_info.coherence(i))
        if cell_info.Monkey(i) == 7
            if sum(strcmp(cell_info.Protocol(i),{'Coarse','CueS2T','CueS4T','Honly','Ronly'}))>0
                cell_info.coherence(i) = 15;
            elseif strncmp(cell_info.Protocol(i),'stim',4)
                cell_info.coherence(i) = 15;
            elseif sum(strcmp(cell_info.Protocol(i),{'SpiT','afterStimSpiT','MemSac','DelSac'}))>0
                cell_info.coherence(i) = 100;
            else
                cell_info.coherence(i) = 99999;
            end
        end
    end
    
    % 2. stimAmp, 前面用的是20uA(未标注)，session 255后改为40uA（已标注）
    if isnan(cell_info.stimAmp(i))
        if cell_info.Monkey(i) == 7 % Ringbell
            if strncmp(cell_info.Protocol(i),'stim',4)
                cell_info.stimAmp(i) = 20;
            else
                cell_info.stimAmp(i) = nan;
            end
        elseif cell_info.Monkey(i) == 16 % Arthas
            if strncmp(cell_info.Protocol(i),'stim',4)
                cell_info.stimAmp(i) = 40;
            else
                cell_info.stimAmp(i) = nan;
            end
        end
    end
    
    % RF
    if isnan(cell_info.RF(i,1))
        cell_info.RF(i,:) = nan(1,4);
    end
end

%% load .mat data here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FnameCode=2; % load all mat files from interested folder
SpiT_data = popul_load_data_LWH(FnameCode,1,monkey_choose);
CueS_data = popul_load_data_LWH(FnameCode,2,monkey_choose);
Stim_data = popul_load_data_LWH(FnameCode,4,monkey_choose);
ASTun_data = popul_load_data_LWH(FnameCode,5,monkey_choose); % after stim tuning data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% add stim amplitutude to stim_data
% 202002 后续工作可能同一个unit做两种不同强度的电刺激，如何处理？
% 1. 不同强度的电刺激分block进行，有两个file  （此情况暂不考虑，默认只有一个file）
% 2. 不同强度的电刺激在同一个block中进行，只有一个file
% 3. 只做一个强度的电刺激，只有一个file

% 20210316
% 不同电刺激强度
for i = 1:length(Stim_data)
    temp_stimAmp = cell_info.stimAmp(strcmp(Stim_data(i).FILE,cell_info.FileNo));
    if temp_stimAmp == 20
        Stim_data(i).stimAmp = 20;
    elseif temp_stimAmp == 40
        Stim_data(i).stimAmp = 40;
    elseif temp_stimAmp == 100
        Stim_data(i).stimAmp = 100;
    elseif isempty(temp_stimAmp)
        Stim_data(i).stimAmp = nan;
    end
end

%% save and export
data_hub.cell_info = cell_info;
data_hub.SpiT_data = SpiT_data;
data_hub.CueS_data = CueS_data;
data_hub.Stim_data = Stim_data;
data_hub.ASTun_data = ASTun_data; % after stim tuning data

return;