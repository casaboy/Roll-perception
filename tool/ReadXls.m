% example
% xlRange{1} = 'A1:DG3000';
% xlRange{2} = 'A1:DG3000';
% handles.XlsData = ReadXls('Z:\Data\MOOG\Results\Result_Lwh.xlsm',[1 2],3,xlRange);

function XlsData = ReadXls(FILE,Sheet,HEADS_N,xlRange)

if nargin<4
    xlRange = [];
end

%% Read xls. HH20140624
% global num txt raw header;

%% 同时读取多个sheet Lwh 202102
num = []; 
txt = []; 
raw = []; 
for s = 1:length(Sheet)
    num_temp = []; txt_temp = []; raw_temp = [];
    if isempty(xlRange)
        [num_temp,txt_temp,raw_temp] = xlsread(FILE,Sheet(s));
    else
        if ~iscell(xlRange)
            [num_temp,txt_temp,raw_temp] = xlsread(FILE,Sheet(s),xlRange); % 原来的xlRange输入为char，懒得改了保留
        else
            [num_temp,txt_temp,raw_temp] = xlsread(FILE,Sheet(s),xlRange{s}); % 多个sheet时循环
        end
    end
    
    % Get Header infomation
    % HEADS_N = 3;
    if HEADS_N == 3
        header_all_this_sheet = txt_temp(HEADS_N-1,:);
        header_all_this_sheet(strcmp(header_all_this_sheet,'')) = txt_temp(HEADS_N-2,strcmp(header_all_this_sheet,''));
    else
        header_all_this_sheet = txt_temp(1,:);
    end
    
    header_this = []; types_this = [];
    for i = 1:length(header_all_this_sheet)
        if ~isempty(header_all_this_sheet{i}) % 可能有的是预留的空各自，跳过
            % 判断第一个字符是不是数字，如果是数字则不能作为变量名
            if ~isempty(str2num(header_all_this_sheet{i}(1)))
                disp('Header info error...');
                keyboard
            end
            
            if i == num_temp(1,i)
                eval(['header_this.' header_all_this_sheet{i} '=' num2str(i) ';']);
                
                if sum(~isnan(num_temp(HEADS_N+1:end,i))) > 0 % Number type
                    eval(['types_this.' header_all_this_sheet{i} '= 1;']); % num, =1
                else
                    eval(['types_this.' header_all_this_sheet{i} '= 2;']); % not num, =2
                end
            else
                disp('Header info error...');
                keyboard;
            end
        end
    end
    
    % Delete headers
    if isfield(header_this,'Monkey')
        end_line = find(~isnan(num_temp(:,header_this.Monkey)),1,'last');
    else
        end_line = find(~isnan(num_temp(:,1)),1,'last');
    end
    
    % 多个sheet统一一起
    num = [num;num_temp(HEADS_N+1:end_line,:)];
    txt = [txt;txt_temp(HEADS_N : end_line - 1,:)]; % Note here
    raw = [raw;raw_temp(HEADS_N+1:end_line,:)];
    
    % save header, type，hName
    header_this_sheet{s} = header_this; % header and code num
    %     types_this_sheet{s} = types_this; % 暂时不输出type类型
    hName_this_sheet{s} = fieldnames(header_this); % header name
end

% Save and Output
if length(Sheet)>1
    if isequal(header_this_sheet{1:length(Sheet)}) && isequal(hName_this_sheet{1:length(Sheet)}) % 不同sheet之间的heading，hName是否一致
        XlsData.header = header_this_sheet{1};
        XlsData.hName = hName_this_sheet{1};
    end
else
    XlsData.header = header_this_sheet{1};
    XlsData.hName = hName_this_sheet{1};
end

XlsData.num = num;
XlsData.txt = txt; % Note here
XlsData.raw = raw;

