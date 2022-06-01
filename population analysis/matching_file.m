% matching spiral file and cueswtiching (2AFC/4AFC) and u-stim. task
% 201906 Lwh
% update for two monkey Lwh20200727
% add some tuning properties here for classification ()

% cell_info.Protocol_code(strcmp(cell_info.Protocol,'SpiT')) = 1;
% cell_info.Protocol_code(strcmp(cell_info.Protocol,'CueS4T')) = 2;
% cell_info.Protocol_code(strcmp(cell_info.Protocol,'Honly')) = 3; % two target，单独的H
% cell_info.Protocol_code(strcmp(cell_info.Protocol,'Ronly')) = 4; % two target，单独的R
% cell_info.Protocol_code(strcmp(cell_info.Protocol,'Coarse4T')) = 5;
% cell_info.Protocol_code(strcmp(cell_info.Protocol,'HCoarse2T')) = 6;
% cell_info.Protocol_code(strcmp(cell_info.Protocol,'RCoarse2T')) = 7;
% cell_info.Protocol_code(strcmp(cell_info.Protocol,'stimCueS')) = 8;
% cell_info.Protocol_code(strcmp(cell_info.Protocol,'stimHonly')) = 9;
% cell_info.Protocol_code(strcmp(cell_info.Protocol,'stimRonly')) = 10;
% cell_info.Protocol_code(strcmp(cell_info.Protocol,'stimCo')) = 11;
% cell_info.Protocol_code(strcmp(cell_info.Protocol,'stimCoH')) = 12;
% cell_info.Protocol_code(strcmp(cell_info.Protocol,'stimCoR')) = 13;
% cell_info.Protocol_code(strcmp(cell_info.Protocol,'afterStimSpiT')) = 14;
% cell_info.Protocol_code(strcmp(cell_info.Protocol,'RFplot')) = 15;
% ============================================================
%    loadM_data.Data:
%    1. different depth: 绝对深度和depth_code，depth_code = 0: main depth
%    2. different spike channel
%    3. different eye postion: centerEye_Protocol, varyEye_Protocol
%    4. different protocol
%    Lwh 201909
% ============================================================

function loadM_data = matching_file(data_hub)

%% Get data
cell_info = data_hub.cell_info; % excel
SpiT_data = data_hub.SpiT_data;
CueS_data = data_hub.CueS_data;
Stim_data = data_hub.Stim_data;
ASTun_data = data_hub.ASTun_data; % after stim tuning data

spifile_num = length(SpiT_data);
spifile_name = {SpiT_data.FILE}';
spifile_chan = cell2mat({SpiT_data.SpikeChan})';

if isfield(CueS_data,'FILE') % 借此判断结构体CueS_data是否为空
    cuesfile_num = length(CueS_data);
    cuesfile_name = {CueS_data.FILE}';
    cuesfile_chan = cell2mat({CueS_data.SpikeChan})';
else
    cuesfile_num = 0;
    cuesfile_name = nan;
    cuesfile_chan = nan;
end

stimfile_num = length(Stim_data);
stimfile_name = {Stim_data.FILE}';
stimfile_chan = cell2mat({Stim_data.SpikeChan})';

if isfield(ASTun_data,'FILE')
    ASTunfile_num = length(ASTun_data);
    ASTunfile_name = {ASTun_data.FILE}';
    ASTunfile_chan = cell2mat({ASTun_data.SpikeChan})';
else
    ASTunfile_num = 0;
    ASTunfile_name = nan;
    ASTunfile_chan = nan;
end

excel_data_num = length(cell_info.Date);
unique_monkey = unique(cell_info.Monkey);

for i = 1:excel_data_num
    c_loc = find(cell_info.FileNo{i}=='c'); % 找到c的位置
    r_loc = find(cell_info.FileNo{i}=='r'); % 找到r的位置
    dot_loc = find(cell_info.FileNo{i}=='.'); % 找到r的位置
    if cell_info.Protocol_code(i) ~= 15 % RFplot
        cell_index(i,1) = str2num(cell_info.FileNo{i}(c_loc+1:r_loc-1));
    else
        cell_index(i,1) = nan; % only RF, no file
    end
end

unique_cell_all_monkey = [];
for m = 1:length(unique_monkey)
    select_monkey = logical([cell_info.Monkey] == unique_monkey(m));
    cell_this_monkey = cell_index(select_monkey);
    cell_this_monkey(isnan(cell_this_monkey)) = [];
    unique_cell_this_monkey{m} = unique(cell_this_monkey);
    cell_num_unique_monkey(m) = length(unique_cell_this_monkey{m});
    unique_cell_all_monkey = [unique_cell_all_monkey;unique_cell_this_monkey{m}];
end
unique_cell_num_all_monkey = sum(cell_num_unique_monkey);



%% matching file for different task
% 每个cell循环，以cellID,depth,eye_position, SpikeChan来获取不同的protocol结果

progressbar('Matching file...');

pass_index = false(excel_data_num,1); % 确保不会筛选到monkey不同但cell_index相同的细胞
cell_info.newDepth = cell_info.Depth;
cell_info.Depth_code = nan(size(cell_info.Depth));
depth_code_temp = [];

for n = 1:unique_cell_num_all_monkey
    for mk = 1:length(unique_monkey)
        %% matching file
        % 1. cell ID
        % *******************************************************************
        this_cell = logical(cell_index == unique_cell_all_monkey(n) & cell_info.Monkey == unique_monkey(mk) & pass_index == false);
        this_cell_index = find(cell_index == unique_cell_all_monkey(n) & cell_info.Monkey == unique_monkey(mk) & pass_index == false);

        if ~isempty(this_cell_index)
            this_monkey = unique_monkey(mk);
            pass_index(this_cell) = true;
            break
        end
    end

    % 2. depth
    % *******************************************************************
    this_depth = cell_info.Depth(this_cell);
    [unique_depth,~,ic] = unique(this_depth); % unique_depth(ic) = this_depth, this_depth(ia) = unique_depth
    
    % 有时候在depth1测了tuing，挪动一点电极后再测其他，此时认为在同一depth,微调允许范围depthTol = 25um
    depthTol = 25; % Depth tolerance
    merge_depth = false(length(unique_depth),1);
    new_unique_depth_temp = zeros(length(unique_depth),1);
    merge_num = 0;
    new_this_depth = [];
    new_unique_depth = [];
    if length(unique_depth)>1
        for dp = 1:length(unique_depth)-1
            if (unique_depth(dp+1) - unique_depth(dp))<depthTol % 微调电极了
                merge_num = merge_num+1;
                merge_depth(dp) = merge_num;
                merge_depth(dp+1) = merge_num;
            end
        end
        
        % adjustment the depth: 微调深度的depth，取平均值
        if sum(merge_depth)>0
            for i = 1:length(unique(nonzeros(merge_depth)))
                %                     temp = [];
                %                     temp = unique_depth(merge_depth==i);
                %                     merge_temp = mean(temp);
                
                merge_temp= mean(unique_depth(merge_depth==i));
                
                % rewrite the uniqeu_depth
                new_unique_depth_temp(merge_depth) = merge_temp;
                new_unique_depth_temp(~merge_depth) = unique_depth(~merge_depth);
            end
            new_this_depth = new_unique_depth_temp(ic);
            new_unique_depth = unique(new_this_depth);
            
            % rewrite the depth
            cell_info.newDepth(this_cell) = new_this_depth;
        else
            new_unique_depth = unique_depth;
            new_this_depth = this_depth;
        end
    else
        new_unique_depth = unique_depth;
        new_this_depth = this_depth;
        cell_info.newDepth(this_cell) = new_this_depth; % update depth after adjustment
    end
    
    % 不同深度下protocol的数目
    protocol_depth = [];
    protocol_depth_num = [];
    protocol_all = [];
    unique_protocol_depth = [];
    main_depth_index{n} = [];
    
    for dp = 1:length(new_unique_depth)
        select = [];
        select = logical(this_cell & cell_info.newDepth == new_unique_depth(dp));
        protocol_depth{dp} = cell_info.Protocol_code(select);
        protocol_depth_num(dp) = length(unique(protocol_depth{dp}));
        protocol_all = [protocol_all;protocol_depth{dp}];
    end
    
    max_protocol_num = max(protocol_depth_num);
    main_depth_index{n} = (protocol_depth_num == max_protocol_num)';
    
    if sum(main_depth_index{n})==1 % only one main depth with maximum protocol
        main_depth_loc = 1;
    elseif sum(main_depth_index{n})>1 % 存在两个深度做的task数目一样多
        index = [];
        index = 1:length(new_unique_depth);
        main_depth_loc = [];
        if sum(sum(protocol_all == [7:12]))>0 % 某个深度下做了电刺激,取该深度作为main depth
            ustim_depth_temp = cellfun(@(x) sum(sum(x == [7:12])),protocol_depth,'uniformoutput',0); % 看哪个深度下做了电刺激
            main_depth_loc = find(cell2mat(ustim_depth_temp)');
            main_depth_index{n} = false(size(new_unique_depth));
            main_depth_index{n}(main_depth_loc) = true;
            
        else % 其他情况一律取中间的深度作为main depth
            main_depth_loc = ceil(length(new_unique_depth)/2);
            main_depth_index{n} = false(size(new_unique_depth));
            main_depth_index{n}(main_depth_loc) = true;
        end
    end
    
    % 根据main depth，recode depth，main为0，main上方为-1，-2，..., 下方为+1,+2,...(Depth_code)
    % (-1和0之间的间隔理论约100um，但是不一定，大于depthTol就会认为是新的一个depth_code)
    unique_depth_code_temp = [];
    unique_depth_code = [];
    for i = 1:length(new_unique_depth)
        unique_depth_code_temp(i) = i - find(main_depth_index{n});
    end
    unique_depth_code = unique_depth_code_temp';
    
    ic = [];
    [~,~,ic] = unique(new_this_depth);

    % rewrite the depth code
    cell_info.Depth_code(this_cell) = unique_depth_code(ic);
    
    % 3. eye position
    % *******************************************************************
    %         % eye position  找到每一个unit做task时候使用的eye position参数组合
    %         % unique_eyePos 其实只有两种，(0,0) 和 根据RF调整的(x,y),默认所有细胞都做了这两种
    %         % 且每个深度最多只有两种眼位置
    %         this_eyePos = [];unique_eyePos = [];
    %         this_eyePos(:,1) = cell_info.eye_x(this_cell);
    %         this_eyePos(:,2) = cell_info.eye_y(this_cell);
    %         unique_eyePos_temp = unique(this_eyePos,'rows');
    %         eyePos_dim = size(unique_eyePos_temp);
    %         unique_eyePos_num(n) = eyePos_dim(1);
    unique_eyePos_num_new = 2; % 每个深度最多只有两种眼位置
    
    % 4. Spike channel
    % *******************************************************************
    this_cell_spichan = cell_info.Chan(this_cell);
    unique_spikechan = unique(this_cell_spichan);
    
    
    % *******************************************************************
    % 以以上4个指标，循环并找到不同potocol的.m文件（实验结果）
    for dp = 1:length(unique_depth_code) % 每个深度进行循环
        % for save
        select_depth = [];
        select_depth_code = [];
        select_depth = new_unique_depth(dp);
        select_depth_code = unique_depth_code(dp) ;
        
        for ey = 1:unique_eyePos_num_new  % 每个眼位置进行循环
            for c = 1:length(unique_spikechan) % 每个spike chan循环
                select_temp = logical(this_cell & cell_info.Depth_code == unique_depth_code(dp));
                select_find_temp  = find(this_cell & cell_info.Depth_code == unique_depth_code(dp));
                
                if ey == 1 % 永远把eye position (0.0)放第一个位置
                    select = logical(select_temp & cell_info.eye_x == 0 & cell_info.eye_y == 0);
                    select_find = find(select_temp & cell_info.eye_x == 0 & cell_info.eye_y == 0);
                    
                elseif ey == 2 % 第二个位置为vary eye position （每个细胞的rf不一样，vary也不同）
                    this_depth_eyePos = [];
                    this_depth_eyePos(:,1) = cell_info.eye_x(select_temp);
                    this_depth_eyePos(:,2) = cell_info.eye_y(select_temp);
                    
                    alteR_eyePos_temp = this_depth_eyePos(sum(this_depth_eyePos == [0 0],2) ~= 2,:); % 非（0，0）的眼位置
                    alteR_eyePos = unique(alteR_eyePos_temp,'rows');
                    
                    if ~isempty(alteR_eyePos)
                        select  = logical(select_temp & cell_info.eye_x == alteR_eyePos(1) & cell_info.eye_y == alteR_eyePos(2));
                        select_find  = find(select_temp & cell_info.eye_x == alteR_eyePos(1) & cell_info.eye_y == alteR_eyePos(2));
                    else
                        select = [];
                        select_find = [];
                    end
                end
                
                % 此时select只可能对应不同的且唯一的protocol （因为不会相同的条件下重复做两个相同的任务）
                % 所以只需要在spiralfiel.m中找到select对应的文件名和chan，就一定是该条件下的protocol result
                % Load .mat
                
                % 202002 后续工作可能同一个unit做两种不同强度的电刺激，如何处理？
                % 1. 不同强度的电刺激分block进行，有两个file  （此情况暂不考虑，默认只有一个file）
                % 2. 不同强度的电刺激在同一个block中进行，只有一个file
                
                protocol_select = cell_info.Protocol_code(select);
                for i = 1:14 % 14个常用task
                    this_name = [];
                    this_name = cell_info.FileNo(select_find(protocol_select == i));
                    this_SUMU = cell_info.SUMU(select_find(protocol_select == i));
                    
                    if ~isempty(this_name)
                        if length(this_name)>1  % 可能原因：more than one spike channel but same file name, 同一个深度做了两次相同的任务，同一个任务微调深度后再做了一次
                            % 原始深度
                            depth_ori = cell_info.Depth(select_find(protocol_select == i));
                            
                            % 电刺激protocol原始深度
                            depth_stim = []; depth_stim_temp = [];
                            for ii = 8:13 % stimCueS stimHonly stimRonly stimCo stimCoH stimCoR
                                depth_stim_temp = cell_info.Depth(select_find(protocol_select == ii));
                                depth_stim = [depth_stim;depth_stim_temp];
                            end
                            unique_depth_stim = unique(depth_stim);
                            
                            % recording protocol原始深度
                            depth_record = []; depth_record_temp = [];
                            for ii = 2:7 % CueS4T Honly Ronly Coarse4T HCoarse2T RCoarse2T
                                depth_record_temp = cell_info.Depth(select_find(protocol_select == ii));
                                depth_record = [depth_record;depth_record_temp];
                            end
                            unique_depth_record = unique(depth_record);
                            
                            if length(unique(depth_ori))==1 % more than one spike channel but same file name, 同一个深度做了两次相同的任务
                                this_name = this_name(1); % 直接选择第一个
                                
                            else % 对于同一个任务微调深度后再做了一次，直接找电刺激（如果做了）相同位置或最近位置的spit
                                % 1. 首先找到其他任务和电刺激任务的原始深度，其次找到其他任务和recording task的原始深度
                                if ~isempty(unique_depth_stim) % 做了电刺激
                                    if length(unique_depth_stim)>1 % 电刺激也调了深度。。。
                                        unique_depth_stim = unique_depth_stim(1); % 暂时选取第一个电刺激深度
                                    end
                                    % 2. 直接找电刺激相同位置或最近位置的spit
                                    [~,use_this] = min(abs(depth_ori-unique_depth_stim));
                                    if length(use_this)>1
                                        use_this = use_this(1);
                                    end
                                    this_name = this_name(use_this);
                                    
                                elseif ~isempty(unique_depth_record) % 没有做电刺激但是做了recording
                                    if length(unique_depth_record)>1 % recording也微调了深度。。。
                                        keyboard
                                        unique_depth_record = unique_depth_record(1); % 暂时选取第一个recording深度  
                                    end
                                    % 3. 找recording相同位置或最近位置的spit
                                    [~,use_this] = min(abs(depth_ori-unique_depth_record));
                                    if length(use_this)>1
                                        use_this = use_this(1);
                                    end
                                    this_name = this_name(use_this);
                                    
                                else % 电刺激和recording都没有做
                                    this_name = this_name(1);
                                end
                            end
                        end
                        
                        file_set = [];file_chan_set = [];
                        if sum(i == [1 14])>0
                            file_set = spifile_name;
                            file_chan_set = spifile_chan;
                            file_data_set = SpiT_data;
                        elseif sum(i == [2 3 4 5 6 7])>0
                            file_set = cuesfile_name;
                            file_chan_set = cuesfile_chan;
                            file_data_set = CueS_data;
                        elseif sum(i == [8 9 10 11 12 13])>0
                            file_set = stimfile_name;
                            file_chan_set = stimfile_chan;
                            file_data_set = Stim_data;
                        end
                        
                        if sum(protocol_select == i)~=0
                            loadM_data(n).Data(dp).Depth = select_depth;
                            loadM_data(n).Data(dp).Depth_code = select_depth_code;
                            if ey == 1
                                loadM_data(n).Data(dp).SpikeChan(c).centerEye_Protocol{i} = file_data_set(strcmp(this_name,file_set) & file_chan_set == unique_spikechan(c));
                                % add SU/MU marker, Lwh 20200803
                                try
                                    loadM_data(n).Data(dp).SpikeChan(c).centerEye_Protocol{i}.SUMU = this_SUMU(c);
                                end
                            else
                                loadM_data(n).Data(dp).SpikeChan(c).varyEye_Protocol{i} = file_data_set(strcmp(this_name,file_set) & file_chan_set == unique_spikechan(c));
                                % add SU/MU marker, Lwh 20200803
                                try
                                    loadM_data(n).Data(dp).SpikeChan(c).varyEye_Protocol{i}.SUMU = this_SUMU(c);
                                end
                            end
                        else
                            loadM_data(n).Data(dp).Depth = select_depth;
                            loadM_data(n).Data(dp).Depth_code = select_depth_code;
                            if ey == 1
                                loadM_data(n).Data(dp).SpikeChan(c).centerEye_Protocol{i} = struct;
                            else
                                loadM_data(n).Data(dp).SpikeChan(c).varyEye_Protocol{i} = struct;
                            end
                        end
                    else
                        loadM_data(n).Data(dp).Depth = select_depth;
                        loadM_data(n).Data(dp).Depth_code = select_depth_code;
                        if ey == 1
                            loadM_data(n).Data(dp).SpikeChan(c).centerEye_Protocol{i} = struct;
                        else
                            loadM_data(n).Data(dp).SpikeChan(c).varyEye_Protocol{i} = struct;
                        end
                    end
                end
            end
        end
    end
    
    % add afterStim tuning data, 只找centerEye_Protocol, Lwh 20211231
    if size(loadM_data(n).Data(main_depth_index{n}).SpikeChan(1).centerEye_Protocol{1},2) && ...
            isfield(loadM_data(n).Data(main_depth_index{n}).SpikeChan(1).centerEye_Protocol{1},'FILE') % 电刺激前做在centerEye_Protocol下做了spiral tuning, 有FILE字段且不为空值
        beforeStim_tun_file_temp = loadM_data(n).Data(main_depth_index{n}).SpikeChan(1).centerEye_Protocol{1}.FILE;
        r_loc_temp = find(beforeStim_tun_file_temp=='r'); % 找到r的位置;
        beforeStim_tun_file = beforeStim_tun_file_temp(1:r_loc_temp-1);
        beforeStim_tun_file_length = length(beforeStim_tun_file);
        
        select_matchAB_file = strncmp(ASTunfile_name,beforeStim_tun_file,beforeStim_tun_file_length); % use file name to find the matching ASTunfile
        find_matchAB_file = find(select_matchAB_file);
        
        if sum(select_matchAB_file) == 0
            loadM_data(n).Data(main_depth_index{n}).SpikeChan(1).AfterStim_Tuning = struct;
        elseif sum(select_matchAB_file) == 1
            ASTun_data_this = ASTun_data(select_matchAB_file);
            if ASTun_data_this.eyepos_x == 0 && ASTun_data_this.eyepos_y == 0 % centerEye_Protocol
                loadM_data(n).Data(main_depth_index{n}).SpikeChan(1).AfterStim_Tuning = ASTun_data_this;
            else
                loadM_data(n).Data(main_depth_index{n}).SpikeChan(1).AfterStim_Tuning = struct;
            end
        else % 一般是在不同眼位置同时做了after stim tuning
            have_this = 0;
            for iii = 1:sum(select_matchAB_file)
                ASTun_data_this = ASTun_data(find_matchAB_file(iii));
                if ASTun_data_this.eyepos_x == 0 && ASTun_data_this.eyepos_y == 0 % centerEye_Protocol
                    loadM_data(n).Data(main_depth_index{n}).SpikeChan(1).AfterStim_Tuning = ASTun_data_this;
                    have_this = 1;
                end
            end
            if ~have_this
                loadM_data(n).Data(main_depth_index{n}).SpikeChan(1).AfterStim_Tuning = struct;
            end
        end
    else
        loadM_data(n).Data(main_depth_index{n}).SpikeChan(1).AfterStim_Tuning = struct;
    end
    
    
    % Save result
    % one time parameter
    loadM_data(n).Date = unique(cell_info.Date(this_cell));
    loadM_data(n).monkey = this_monkey;
    loadM_data(n).cellID = unique_cell_all_monkey(n);
    loadM_data(n).GridX = unique(cell_info.xloc(this_cell));
    loadM_data(n).GridY = unique(cell_info.yloc(this_cell));
    loadM_data(n).Area = cell2mat(unique(cell_info.Area(this_cell)));
    loadM_data(n).main_depth = new_this_depth(main_depth_index{n});
    
    loadM_data(n).guidetube = unique(cell_info.guidetube(this_cell));
    loadM_data(n).offset = unique(cell_info.offset(this_cell));
    

    if length(loadM_data(n).GridX)>1 || length(loadM_data(n).GridX)>1
        keyboard
    end

    % 关于stimAmp的参数输入转到上方写进stim_data中（写进 loadM_data(mk,n).Data(dp).SpikeChan(c).centerEye_Protocol{i} 等中，不单独提取出来到loadM_data.stimAmp）
    %         % u-stim amplitude, 暂时一个细胞只有一个电刺激强度（不考虑一个细胞做了两种电刺激强度,如果做了，选40uA的）
    %         temp_stimAmp = cell_info.stimAmp(this_cell);
    %         temp_stimAmp = temp_stimAmp(~isnan(temp_stimAmp)); % 去掉不是电刺激的file （nan）
    %         temp_stimAmp = unique(temp_stimAmp);
    %
    %         if temp_stimAmp == 20
    %             loadM_data(mk,n).stimAmp = 20;
    %         elseif temp_stimAmp == 40
    %             loadM_data(mk,n).stimAmp = 40;
    %         elseif temp_stimAmp == 2040 % 同一个file同时做了20/40uA电刺激
    %             loadM_data(mk,n).stimAmp = [20 40];  % 暂定！！！
    %         end
    
    % RF (与深度，眼位置基本无关)
    this_RF_temp = cell_info.RF(this_cell,:);
    have_rf = ~isnan(this_RF_temp(:,1));
    if sum(have_rf)==1
        loadM_data(n).RF = cell_info.RF(this_cell_index(have_rf),:);
    elseif sum(have_rf)>1 % 多次测量RF（不同深度），取平均
        temp = cell_info.RF(this_cell_index(have_rf),:);
        loadM_data(n).RF = mean(temp);
    else
        loadM_data(n).RF = nan(1,4);
    end
    
    loadM_data(n).unique_SpikeChan = unique_spikechan;
    
    % 计数
    progressbar(n/unique_cell_num_all_monkey);
end

disp(' ');
disp('Matching files finish!');
disp('---------------------------------------------------------------------------------------------');

return;