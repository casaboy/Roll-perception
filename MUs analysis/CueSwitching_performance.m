% Plots performance for H/R discrimination coarse task
% % Lwh, 04,2018
% % Lwh, 07, 2018 added microstimulation
% % Lwh  09, 2018 added coarse task
% % Lwh  12, 2020 前：固定Rot_amplitude = 45, rotation定义为实际旋转为固定Rot_amplitude的百分比（≈0.01-0.2） true_rot_amplitude = Rot_amplitude*rotation，横坐标：trotation，即与Expasion偏离的polar角度
% % Lwh  12, 2020 后：rotation定义为与Expasion偏离的polar角度，true_rot_amplitude = Rot_amplitude*cos(90-rotation)，横坐标：rotation
% % Lwh  1, 2021  session 319, m7c1274(含)以后的task改用圆点视觉刺激，并且调整旋转中心为fixation point，用stars_type定义
% % Lwh  03, 2021 add difficult trial condition, modify old code
% % for positive value assign preferred direction (manmul input)

% %-------------------------------------------------------
function CueSwitching_performance(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, batch_flag);

% manmul input the preferred direction
% global_pref_code 1: 左，CW，rotation；  2: 右，CCW，heading
global_pref_code = nan(1,3);
global_pref_code = [1 2 2];

TEMPO_Defs;
Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

IF_FIGURE = 1;  % 1：画图； 0：不画图
Is_step = 0; % 是否分step画psy curve并计算

%get the column of values for azimuth and rotation and stim_type
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
temp_heading   = data.moog_params(HEADING, :, MOOG);
temp_rotation = data.moog_params(ROTATION,:,MOOG);
temp_amplitude = data.moog_params(AMPLITUDE,:,CAMERAS);
temp_switch_index = data.moog_params(SWITCH_INDEX,:,MOOG);
temp_coherence = data.moog_params(COHERENCE,:,MOOG);
temp_microstim = data.misc_params(MICROSTIM,:);
temp_microstim2 = data.misc_params(MICROSTIM2,:);

% 更换刺激从三角形变成圆点，STARS_TYPE从0变成1
temp_stars_type = data.moog_params(STARS_TYPE,:,MOOG);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for old setting
%
% 由于TEMPO代码没写好，改了rotation定义后，TEMPO输出的为实际呈现的AMPLITUDE和ROT_AMPLITUDE（也就是true_amplitude和true_rot_amplitude）
% Rotation时：RotAzimuth = 90 - ROTATION;
% true_amplitude = sin(RotAzimuth)*Amplitude
% true_rot_amplitude = cos(RotAzimuth)*Rot_Amplitude
% 后面改为设定的amplitude和rot_amplitude，实际刺激按照上式计算即可
% 只涉及大概202011-202012的file，手动设置为temp_rot_amplitude=45
try
    temp_rot_amplitude = data.moog_params(ROT_AMPLITUDE,:,CAMERAS);
catch
    temp_rot_amplitude = ones(size(temp_rotation))*45;
end
unique_amplitude = munique(temp_amplitude');
unique_rot_amplitude = munique(temp_rot_amplitude');
if length(unique_amplitude)>1
    temp_amplitude = ones(size(temp_amplitude))*0.2;
end
if length(unique_rot_amplitude)>1
    temp_rot_amplitude = ones(size(temp_rot_amplitude))*45;
end

if max(unique(temp_rotation')) < 2 % 202012前的old setting
    Rot_Amplitude = 45; % 原来默认是45
    true_rot_amplitude = temp_rotation*Rot_Amplitude; % true_rot_amplitude
    % 转化为偏离Expasion的polar角度rotation:
    % true_rot_amplitude = cos(90-rotation)*Rot_Amplitude
    % rotation = 90-acos(true_rot_amplitude/Rot_Amplitude)
    temp_rotation = 90-acosd(true_rot_amplitude/Rot_Amplitude);
end

try
    temp_eyepos_x = data.targ_params(1,:,1);
    temp_eyepos_y = data.targ_params(3,:,1);
    
    eyepos_x = munique(temp_eyepos_x');
    eyepos_y = munique(temp_eyepos_y');
catch % for rescue file
    eyepos_x = data.targ_params(1);
    eyepos_y = data.targ_params(2);
end

% mistake in file m7c960r48
if strncmp(FILE,'m7c960r48.htb',9)
    eyepos_y=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MOOG 中，self-motion为CW，stimulu为CCW时候，rotation>0，我们画图的时候self-motion CW为rotation<0
% 将rotation值调转，让rotation>0对应CCW self-motion ,rotation<0对应CW self-motion    Lwh20181214
temp_rotation = temp_rotation.*(-1);

unique_heading = munique(temp_heading');
unique_rotation = munique(temp_rotation'); % h+r task中rotation没有0度
unique_heading_exc0 = nonzeros(unique_heading);
unique_rotation_exc0 = nonzeros(unique_rotation);
unique_stim_type = munique(temp_stim_type');
unique_coherence = munique(temp_coherence');
unique_switch_index = munique(temp_switch_index');
unique_amplitude = munique(temp_amplitude');
unique_rot_amplitude = munique(temp_rot_amplitude');

% add temp_microstim2 to control interleave amplitude stimulation， Lwh 2021014
%        MICROSTIM   MICROSTIM2
% amp1      1            0
% amp2      0            1
% amp3      1            1
% amp1,2,3对应的amplitude值（uA）看ao的script
unique_microstim = munique(temp_microstim');
unique_microstim2 = munique(temp_microstim2');
temp_stimcondition = zeros(size(temp_microstim));
temp_stimcondition(temp_microstim == 1 & temp_microstim2 == 0) = 1; % amp1
temp_stimcondition(temp_microstim == 0 & temp_microstim2 == 1) = 2; % amp1
temp_stimcondition(temp_microstim == 1 & temp_microstim2 == 1) = 3; % amp1,如果有的话
unique_stimcondition = munique(temp_stimcondition');


% For coarse task, adjustment coherence from negative to positive
% use unique_coherence as condition like unique_heading or unique_rotation
if length(unique_coherence) > 1 % coarse task
    % 分两种情况，一种用正负角度来控制负值coherence(coherence都是正值)，另一种直接用负值coherence
    if sum(unique_coherence<0)==0  % 第一种情况
        % 将角度为负值的trial的coherence改为负值表示,使得coherence从负值到0到正值，heading\rotation角度全部用正值表示
        select_neg = find( temp_heading<0 | temp_rotation<0);
        select_pos = find( temp_heading>0 | temp_rotation>0);
        temp_coherence(select_neg) = -1*temp_coherence(select_neg);
        
        temp_heading(temp_heading~=0)=90; 
        temp_rotation(temp_rotation~=0)=90; 
         
        unique_coherence = munique(temp_coherence');
        unique_heading = munique(temp_heading');
        unique_rotation = munique(temp_rotation');
        unique_heading_exc0 = nonzeros(unique_heading);
        unique_rotation_exc0 = nonzeros(unique_rotation);
        unique_coherence_exc0 = nonzeros(unique_coherence);
    else % 第二种情况，不需要调整
    end
else % find task
end

%% Find repetition num, zero_cond_num
% Assume trials from first to the last one
% To find how many condition dose one repetition have (one_repetion).   Zero degree condition may be two or three or...
if length(unique_coherence) == 1  % fine task
    zero_trial_num = sum(logical(temp_heading==0 & temp_rotation==0)); % include ctrl and stim
    if ~isempty(unique_heading_exc0) % have translation part: heading only 2-target task or 4-target task
        for n = 1:length(unique_stimcondition)
            for i = 1:length(unique_heading_exc0)
                h_trial_num(n,i)=sum(logical(temp_heading==unique_heading_exc0(i) & temp_rotation==0 & temp_stimcondition == unique_stimcondition(n)));
            end
        end
    else % 单独做一种motion type
        h_trial_num = nan;
    end
    
    if ~isempty(unique_rotation_exc0) % have rotation part: rotation only 2-target task or 4-target task
        for n = 1:length(unique_stimcondition)
            for i = 1:length(unique_rotation_exc0)
                r_trial_num(n,i)=sum(logical(temp_rotation==unique_rotation_exc0(i) & temp_heading==0 & temp_stimcondition == unique_stimcondition(n)));
            end
        end
    else % 单独做一种motion type
        r_trial_num = nan;
    end

    repetition_temp = min([r_trial_num(:);h_trial_num(:)]); % ctrl重复了repetition_temp，stim也重复了repetition_temp
    zero_cond_num = round(zero_trial_num/repetition_temp / length(unique_stimcondition)); % Number of zero degree condition (0,1,2,3,...)
    one_repetition = (zero_cond_num + length(unique_heading_exc0) + length(unique_rotation_exc0))*length(unique_stim_type)*length(unique_stimcondition)*length(unique_coherence); % include ctrl and stim
    
elseif length(unique_coherence) > 1 % coarse task
    zero_trial_num = sum(temp_coherence==0);
%     if zero_trial_num > 0
        for n = 1:length(unique_stimcondition)
            for j = 1:length(unique_switch_index)
                for i = 1:length(unique_coherence_exc0)
                    c_trial_num(n,j,i)=sum(logical(temp_coherence==unique_coherence_exc0(i) & temp_stimcondition == unique_stimcondition(n)) & temp_switch_index == unique_switch_index(j));
                end
            end
        end
%     else
%         
%         
%     end
    
    repetition_temp = min(c_trial_num(:)); % 除了coherence=0外的conherence最少重复次数
    zero_cond_num = round(zero_trial_num/repetition_temp / length(unique_stimcondition)); % 因为TEMPO中condition list的设置，coherence=0往往有两个，coherence+0和coherence-0 （其实是一样的）
    one_repetition = (zero_cond_num + length(unique_coherence_exc0) * length(unique_switch_index)) * length(unique_stim_type)*length(unique_stimcondition); % include ctrl and stim
end

trials_num = one_repetition * repetition_temp; %完整的repetition含有的trial数


%% Select start/included/excluede tiral
% If length(BegTrial) > 1 and all elements are positive, they are trials to be included.
% Else, if all elements are negative, they are trials to be excluded.
% This enable us to exclude certain trials ** DURING ** the recording more easily. HH20150410
raw_trials = 1:length(temp_coherence);	% a vector of trial indices
select_trials = false(size(raw_trials));
if length(BegTrial) == 1 && BegTrial > 0 % Backward compatibility
    if BegTrial==1 && EndTrial==trials_num
        select_trials(BegTrial:EndTrial) = true;
        repetition = repetition_temp;
        a=0;b=0;
        
    elseif BegTrial==1 && EndTrial~=trials_num
        b = ceil((trials_num -EndTrial)/one_repetition);
        select_trials(BegTrial:trials_num) = true;
        select_trials( (trials_num-one_repetition*b + 1 ):trials_num) = false;
        EndTrial = trials_num-one_repetition*b ;
        repetition = repetition_temp-b;
        a=0;
        
    elseif BegTrial~=1 && EndTrial==trials_num
        a = ceil(BegTrial/one_repetition);
        select_trials( (one_repetition*a + 1):trials_num) = true;
        repetition = repetition_temp-a;
        b=0;
    else
        a = ceil(BegTrial/one_repetition);
        select_trials( (one_repetition*a + 1):trials_num) = true;
        b = ceil((trials_num-EndTrial)/one_repetition);
        select_trials( (trials_num-one_repetition*b + 1 ):trials_num) = false;
        EndTrial = trials_num-one_repetition*b ;
        repetition = repetition_temp-a-b;
    end
elseif all(BegTrial > 0) % To be included
    select_trials(BegTrial) = true;
    a=0;b=0;
elseif all(BegTrial < 0) % To be excluded
    b = ceil((trials_num-EndTrial)/one_repetition);
    select_trials(1:trials_num) = true;
    select_trials( (trials_num-one_repetition*b + 1 ):trials_num) = false;
    select_trials(-BegTrial) = false;
    EndTrial = trials_num-one_repetition*b ;
    repetition = repetition_temp-b;
    a=0;
    %     temp_trials(outline)= false;
else
    disp('Trial selection error...');
    keyboard;
end


%% %%%%%%%%%%%%%%%%%%%%   traditional Psy curve   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 方法1：sliding window = 6rep，step = 1 rep，不同阶段电刺激效果
% window_width = 6;
% step_width = 1;
% step_num = repetition+1-window_width;
% step_num = step_num+1;  % 最后一次算全部
%
% if step_num < 1 % 重复次数不够window_width = 6次，全部计算，不分windows了
%     step_num = 1;
% end

% 方法2： 每次多增加一个rep，累积rep时候电刺激效果 （指导实验大概做到那个rep）

start_rep = 5; % 从第5个rep开始(第一个循环是1-5rep)
step_num = repetition - start_rep + 1;

if step_num < 1 % 不满5个rep
    step_num = 1;
    not_enough_rep = 1;
else
    not_enough_rep = 0;
end

if Is_step % go step analyse
    step_count = [1:step_num]; % for loop
    last_step = step_num;
else % not step
    step_count = step_num;
    last_step = step_num;
end

for step = step_count
    slide_trials = select_trials;
    if ~not_enough_rep
        slide_trials((start_rep+step-1)*one_repetition+1:end) = false;
    end

    %% choice related
    % determine for each trial whether monkey chooses leftward(target1) or
    % rightward(tarket2) or up(target3,CW) or down(target4,CCW)
    LEFT = 1;
    RIGHT = 2;
    CW = 3;
    CCW = 4;
    
    for i = 1:sum(slide_trials)
        temp = data.event_data(1,:,i+one_repetition*a);
        events = temp(temp>0); % all non-zero entries
        if (sum(events == IN_T1_WIN_CD) > 0) % 右target，stim: from right to left, self-motion: right, 对应heading > 0
            choice(i) = RIGHT;
        elseif (sum(events == IN_T2_WIN_CD) > 0)  % 左target，stim: from left to right, self-motion: Left , 对应heading < 0
            choice(i) = LEFT;
        elseif(sum(events == IN_T3_WIN_CD) > 0)  % 下target, stim: CCW,  self-motion:CW , 对应MOOG中rotation > 0，现校正为<0
            choice(i) = CW;
        elseif(sum(events == IN_T4_WIN_CD) > 0)  % 上target, stim: CW, self-motion:CCW, 对应MOOG中rotation < 0，现校正为>0
            choice(i) = CCW;
        else
            disp('Neither T1 or T2 or T3 or T4 chosen.  This should not happen!.  File must be bogus.');
        end
    end
    
    % group choice for L,CW,Heading  or R, CCW, Rotation and separated by motion type
    choice_LEFT_2T{1} = logical(choice==LEFT);
    choice_RIGHT_2T{1} = logical(choice==RIGHT);
    choice_LEFT_2T{2} = logical(choice==CW);
    choice_RIGHT_2T{2} = logical(choice==CCW);
    choice_LEFT_2T{3} = logical((choice==CCW)|(choice==CW));
    choice_RIGHT_2T{3} = logical((choice==RIGHT)|(choice==LEFT));
    
    % group choice by motion type
    choice_motiontype = [];
    choice_motiontype{1} = logical((choice==RIGHT)|(choice==LEFT)); % Heading
    choice_motiontype{2} = logical((choice==CCW)|(choice==CW)); % Rotation
    choice_motiontype{3} = logical((choice==RIGHT)|(choice==LEFT)|(choice==CCW)|(choice==CW)); % MT
    
%     % group choice by Heading or Rotation
%     choiceHR_temp = [];
%     choiceHR_temp(choice==LEFT) = 1;
%     choiceHR_temp(choice==RIGHT) = 1;
%     choiceHR_temp(choice==CW) = 2;
%     choiceHR_temp(choice==CCW) = 2;
%     
%     valuest = 1:2;
%     catnames = {'H' 'R'};
%     choiceHR = categorical(choiceHR_temp,valuest,catnames);
    
    %% condition related
    % reselect some parameters from slide tiral
    heading = temp_heading(slide_trials);
    rotation = temp_rotation(slide_trials);
    stim_type = temp_stim_type(slide_trials);
    switch_index = temp_switch_index(slide_trials);
    stimcondition = temp_stimcondition(slide_trials);
    coherence = temp_coherence(slide_trials);
    
    unique_heading = munique(heading');
    unique_rotation = munique(rotation'); % 实际上rotation没有0度
    unique_stim_type = munique(stim_type');
    unique_switch_index = munique(switch_index');
    unique_stimcondition = munique(stimcondition');
    unique_coherence = munique(coherence');
    
    unique_heading_exc0 = nonzeros(unique_heading);
    unique_rotation_exc0 = nonzeros(unique_rotation); 
    
    stars_type = munique(temp_stars_type(slide_trials)');

    % group heading/rotation/coherence degree by motion type
    unique_degree = []; degree = [];
    if length(unique_coherence) == 1 % fine task
        unique_degree{1} = unique_heading;
        unique_degree{2} = unique_rotation;
        degree{1} = heading;
        degree{2} = rotation;
        degree{3} = zeros(1,length(heading));
        
        % mt unique_cond: normalize heading and rotation; from rotation to transaltion
        if length(unique_switch_index)==2
            nor_heading_temp = unique_heading(ceil(length(unique_heading)/2)+1:end);
            nor_rotation_temp = unique_rotation(ceil(length(unique_rotation)/2)+1:end);
            %             nor_heading = normc(nor_heading_temp);
            %             nor_rotation = normc(nor_rotation_temp);
            
            % normalize to -1 to 1
            nor_heading_temp = [0; nor_heading_temp];
            nor_heading = mapminmax(nor_heading_temp',0,1);
            nor_rotation_temp = [0; nor_rotation_temp];
            nor_rotation = mapminmax(nor_rotation_temp',0,1);
            nor_heading = nor_heading(2:end)'; % remove zero
            nor_rotation = nor_rotation(2:end)'; % remove zero
            
            if zero_cond_num > 0 % havev zero degree
                unique_degree{3} = [flipud(-nor_rotation);0;nor_heading];
            else
                unique_degree{3} = [flipud(-nor_rotation);nor_heading];
            end
            
            for i = 1:length(unique_heading_exc0)/2
                degree{3}(heading == unique_heading(i) | heading == unique_heading(length(unique_heading)+1-i)) = unique_degree{3}(end+1-i);
            end
            for j = 1:length(unique_rotation_exc0)/2
                degree{3}(rotation == unique_rotation(j) | rotation == unique_rotation(length(unique_rotation)+1-j)) = unique_degree{3}(j);
            end
        else
            unique_degree{3} = 0;
        end
    else % coarse task
        if length(unique_switch_index)==1
            if unique_switch_index == 1 % coarse translation task
                unique_degree{1} = unique_coherence;
                unique_degree{2} = 0;
            elseif unique_switch_index == 2 % coarse rotation task
                unique_degree{1} = 0;
                unique_degree{2} = unique_coherence;
            end
            unique_degree{3} = zeros(1,length(unique_coherence));
        elseif length(unique_switch_index)==2 % coarse 4-T
            unique_degree{1} = unique_coherence;
            unique_degree{2} = unique_coherence;
            unique_degree{3} = mapminmax(unique_coherence',-1,1)';
        end
        degree{1} = coherence;
        degree{2} = coherence;
        degree{3} = mapminmax(coherence,-1,1); % normalize to -1 to 1
    end
    
    % save raw data for condition, choice
    if step == last_step
        raw_data = table((1:trials_num)', heading',rotation',degree{3}',coherence',stimcondition',switch_index',choice',...
            'VariableName',{'Trial' 'Translation' 'Rotation' 'MotionType' 'Coherence' 'Stimcondition' 'SwitchIndex' 'Choice'});
    end

    %%
    % fake switch_index for separating heading and rotation anlysze
    % 1 for heading , 2 for rotation
    if length(unique_switch_index)==1
        if unique_switch_index==1
            fake_switch_index=[1];
        elseif unique_switch_index==2
            fake_switch_index=[2];
        end
    else
        fake_switch_index = [1;2;3]; % 3 for motion type
    end
    
    %% Psychometric function data base: base on monkey choice and task type
    for n = 1:length(unique_stimcondition) % microstior normal n=1:normal   n=2:microsti
        for k=1:length(unique_stim_type)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 1. 选出0度task和非0度task
            if length(unique_coherence) == 1 % fine task
                % 当仅有heading task或者仅有rotation task时候，各自有0度；当h+r task时候，一般情况0度仅存在heading task中
                if length(unique_switch_index)==1 %仅heading或者仅rotation
                    if unique_switch_index==1 %仅heading
                        select_trial_0{n,k} = logical(heading==0 & switch_index==1 & stim_type == unique_stim_type(k) & stimcondition == unique_stimcondition(n)); %所有0度task
                        rep_0{n,k}= sum(select_trial_0{n,k}); %0度task个数（区分ctrl和stim，注意与zero_trial_num有区别）
                    elseif unique_switch_index==2 %仅rotationa
                        select_trial_0{n,k} = logical(rotation==0 & switch_index==2 & stim_type == unique_stim_type(k) & stimcondition == unique_stimcondition(n)); %所有0度task
                        rep_0{n,k} = sum(select_trial_0{n,k}); %0度task个数
                    end
                elseif length(unique_switch_index)==2 %h+r
                    select_trial_0{n,k} = logical(heading==0 & rotation==0 & stim_type == unique_stim_type(k) & stimcondition == unique_stimcondition(n)); %所有0度task
                    rep_0{n,k} = sum(select_trial_0{n,k}); %0度task个数
                end
            else % coarse task
                select_trial_0{n,k} = logical(coherence==0 & stim_type == unique_stim_type(k) & stimcondition == unique_stimcondition(n)); %所有0度task
                rep_0{n,k} = sum(select_trial_0{n,k}); %0度task个数
            end
            select_trial_non_0{n,k} = ~(select_trial_0{n,k});
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 2. select trial for each heading/rotation/MT degree
            for j = 1:length(fake_switch_index)
                ind = fake_switch_index(j);
                select = logical(stim_type == unique_stim_type(k) & stimcondition == unique_stimcondition(n));
                for i = 1:length(unique_degree{ind})
                    if length(unique_coherence) == 1 % fine task
                        if unique_degree{ind}(i)~=0 % 非0度
                            select_degree{n,k,ind,i} = logical(degree{ind} == unique_degree{ind}(i) & select_trial_non_0{n,k} & select); % motion type内
                        else % 0度
                            select_degree{n,k,ind,i} = select_trial_0{n,k};  % 不再区分系统认为的0度正确task type （task_type）
                        end
                        
                    else % coarse task
                        if ind < 3 % heading or rotation task
                            if unique_degree{ind}(i)~=0 % 非0度
                                select_degree{n,k,ind,i} = logical(coherence == unique_coherence(i) & switch_index == unique_switch_index(j) & select);
                            else % 0度
                                select_degree{n,k,ind,i} = logical( coherence == 0 & select); 
                            end
                        else  % mt: from more rotation to more heading, need to seperate same coherence from heading or rotation
                            if i <= floor(length(unique_degree{3})/2) % rotation part
                                select_degree{n,k,ind,i} = logical((degree{ind} == unique_degree{ind}(i) | degree{ind} == unique_degree{ind}(length(unique_degree{ind})+1-i)) & switch_index == 2 & rotation ~= 0 & select);
                            elseif i > ceil(length(unique_degree{3})/2) % heading part
                                select_degree{n,k,ind,i} = logical((degree{ind} == unique_degree{ind}(i) | degree{ind} == unique_degree{ind}(length(unique_degree{ind})+1-i)) & switch_index == 1 & heading ~= 0 & select);
                            else % 0 degree
                                select_degree{n,k,ind,i} = logical( coherence == 0 & select);
%                                 select_degree_4T{n,k,ind,i} = logical( coherence == 0 & select);
                            end
                        end
                    end
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 3. Psychometric function data base: base on monkey choice and degree
            % choose L or R and Heading task → Heading psy function (ignore choose CW or CCW or rotation task trial)
            % choose CW or CCW and Rotation task → Rotation psy function
            %
            %                ↗  choose L or R    →  plot psy function
            % Heading task →
            %                ↘ choose CW or CCW
            
            Psy_correct{n,k,1,step} = nan(1,length(unique_degree{ind})); % for heading
            Psy_correct{n,k,2,step} = nan(1,length(unique_degree{ind})); % for rotation
            Psy_correct{n,k,3,step} = nan(1,length(unique_degree{ind})); % for MT
            
            % for true error   choose R / (L+R+CW+CCW)
            Psy_correct_true{n,k,1,step} = nan(1,length(unique_degree{ind})); % for heading
            Psy_correct_true{n,k,2,step} = nan(1,length(unique_degree{ind})); % for rotation
            Psy_correct_true{n,k,3,step} = nan(1,length(unique_degree{ind})); % for MT
            
            % for positive value asign to preferred direction
            Psy_correct_pd{n,k,1,step} = nan(1,length(unique_degree{ind})); % for heading
            Psy_correct_pd{n,k,2,step} = nan(1,length(unique_degree{ind})); % for rotation
            Psy_correct_pd{n,k,3,step} = nan(1,length(unique_degree{ind})); % for MT

            % 当前motion type内的choice
            for j = 1:length(fake_switch_index)
                ind = fake_switch_index(j);
                for i = 1:length(unique_degree{ind}) % x axis from left to right, from CW to CCW, from rotation to heading
                    % 当前motion type内的choice
                    trials_select_2T{n,k,ind,i} = logical(select_degree{n,k,ind,i} & choice_motiontype{ind}); % Heading task中选heading target， R task中选R target的总数
                    trials_num_2T{n,k,ind,step}(i) = sum(trials_select_2T{n,k,ind,i}); % Save it
                    
                    % make 'S' curve by using the rightward choice / CCW choice / Heading choice for y-axis
                    Psy_correct{n,k,ind,step}(i) = sum(trials_select_2T{n,k,ind,i} & choice_RIGHT_2T{ind}) / trials_num_2T{n,k,ind,step}(i);
                    
                    %                     % save raw choice and trial num
                    %                     Psy_correct_raw{n,k,ind}(1,i) = sum(trials_select_2T{n,k,ind,i} & choice_RIGHT_2T{ind}); % select Right
                    %                     Psy_correct_raw{n,k,ind}(2,i) = sum(trials_select_2T{n,k,ind,i}); % select all this motion type
                    
                    % true error for rightward choice / CCW choice,  R/(L+R+CW+CCW), CCW/(L+R+CW+CCW)
                    % Lwh 20220323
                    trials_num_4T{n,k,ind,step}(i) = sum(select_degree{n,k,ind,i}); % Save it
                    Psy_correct_true{n,k,ind,step}(i) = sum(trials_select_2T{n,k,ind,i} & choice_RIGHT_2T{ind}) / trials_num_4T{n,k,ind,step}(i); % Heading task choose R / heading task all this condition
                    
                    % if have input preferred direction: global_pref_code: 1: 左，CW，rotation；  2: 右，CCW，heading
                    if global_pref_code(ind) == 1 % Left/ALL
                        Psy_correct_pd{n,k,ind,step}(i) = sum(trials_select_2T{n,k,ind,i} & choice_LEFT_2T{ind}) / trials_num_2T{n,k,ind,step}(i);
                    elseif global_pref_code(ind) == 2 % same x axis, R/ALL
                        Psy_correct_pd{n,k,ind,step}(i) = Psy_correct{n,k,ind,step}(i);
                    else
                        Psy_correct_pd{n,k,ind,step}(i) = nan;
                    end
                end
                
                % calculate from left to right, but if global_pref_code(ind) == 1 (preferred left), reverse x axis
                if global_pref_code(ind) == 1 % from non-preferred to preferred direction
                    Psy_correct_pd{n,k,ind,step} = fliplr(Psy_correct_pd{n,k,ind,step});
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 4. fitting
            method = 0; % 0=Maximum likelihood (default), 1=Square error
            tolerance = 0; % 0=No tolerance (default), other = error tolerance of the max abs (heading)
            for j = 1:length(fake_switch_index)
                ind = fake_switch_index(j);
                fit_data_psy_cum{n,k,ind}(:,1) = unique_degree{ind};
                fit_data_psy_cum{n,k,ind}(:,2) = Psy_correct{n,k,ind,step}';
                fit_data_psy_cum{n,k,ind}(:,3) = trials_num_2T{n,k,ind,step}';
                [bb,tt] = cum_gaussfit_max1(fit_data_psy_cum{n,k,ind},method,tolerance);
                Bias_psy{n,k,ind,step} = bb;
                Thresh_psy{n,k,ind,step} = tt;
                psy_perf{n,k,ind,step} = [bb, tt];
               
                % fit the ture error: R/(L+R+CW+CCW), CCW/(L+R+CW+CCW)
                % 得到的PSE和threshold是基于选择为50%，实际上应该基于25%？？？
                fit_data_psy_cum2{n,k,ind}(:,1) = unique_degree{ind};
                fit_data_psy_cum2{n,k,ind}(:,2) = Psy_correct_true{n,k,ind,step}';
                fit_data_psy_cum2{n,k,ind}(:,3) = trials_num_4T{n,k,ind,step}';
                [bb,tt] = cum_gaussfit_max1(fit_data_psy_cum2{n,k,ind},method,tolerance);
                Bias_psy_true50{n,k,ind,step} = bb;
                Thresh_psy_true50{n,k,ind,step} = tt;
                psy_perf_true50{n,k,ind,step} = [bb, tt];
                % 根据累计高斯函数反推25%时的PSE
                Bias_psy_true25{n,k,ind,step} = norminv(0.25,bb,tt);
                
                fit_data_psy_cum{n,k,ind}(:,1) = unique_degree{ind}; % actual we should reverse x axis when global_pref_code(ind) == 1, but it is same
                fit_data_psy_cum{n,k,ind}(:,2) = Psy_correct_pd{n,k,ind,step}';
                fit_data_psy_cum{n,k,ind}(:,3) = trials_num_2T{n,k,ind,step}';
                [bb,tt] = cum_gaussfit_max1(fit_data_psy_cum{n,k,ind},method,tolerance);
                Bias_psy_pd{n,k,ind,step} = bb;
                Thresh_psy_pd{n,k,ind,step} = tt;
                psy_perf_pd{n,k,ind,step} = [bb, tt];
            end
            
            if length(fake_switch_index) == 1
                if fake_switch_index == 1
                    Bias_psy{n,k,2,step} = nan;
                    Thresh_psy{n,k,2,step} = nan;
                    psy_perf{n,k,2,step} = [nan, nan];
                    Bias_psy{n,k,3,step} = nan;
                    Thresh_psy{n,k,3,step} = nan;
                    psy_perf{n,k,3,step} = [nan, nan];
                    
                    Bias_psy_true50{n,k,2,step} = nan;
                    Thresh_psy_true50{n,k,2,step} = nan;
                    psy_perf_true50{n,k,2,step} = [nan, nan];
                    Bias_psy_true50{n,k,3,step} = nan;
                    Thresh_psy_true50{n,k,3,step} = nan;
                    psy_perf_true50{n,k,3,step} = [nan, nan];
                    
                    Bias_psy_true25{n,k,2,step} = nan;
                    Bias_psy_true25{n,k,3,step} = nan;
                    
                    Bias_psy_pd{n,k,2,step} = nan;
                    Thresh_psy_pd{n,k,2,step} = nan;
                    psy_perf_pd{n,k,2,step} = [nan, nan];
                    Bias_psy_pd{n,k,3,step} = nan;
                    Thresh_psy_pd{n,k,3,step} = nan;
                    psy_perf_pd{n,k,3,step} = [nan, nan];
                else
                    Bias_psy{n,k,1,step} = nan;
                    Thresh_psy{n,k,1,step} = nan;
                    psy_perf{n,k,1,step} = [nan, nan];
                    Bias_psy{n,k,3,step} = nan;
                    Thresh_psy{n,k,3,step} = nan;
                    psy_perf{n,k,3,step} = [nan, nan];
                    
                    Bias_psy_true50{n,k,1,step} = nan;
                    Thresh_psy_true50{n,k,1,step} = nan;
                    psy_perf_true50{n,k,1,step} = [nan, nan];
                    Bias_psy_true50{n,k,3,step} = nan;
                    Thresh_psy_true50{n,k,3,step} = nan;
                    psy_perf_true50{n,k,3,step} = [nan, nan];
                    
                    Bias_psy_true25{n,k,1,step} = nan;
                    Bias_psy_true25{n,k,3,step} = nan;
                    
                    ias_psy_pd{n,k,1,step} = nan;
                    Thresh_psy_pd{n,k,1,step} = nan;
                    psy_perf_pd{n,k,1,step} = [nan, nan];
                    Bias_psy_pd{n,k,3,step} = nan;
                    Thresh_psy_pd{n,k,3,step} = nan;
                    psy_perf_pd{n,k,3,step} = [nan, nan];
                end
            end  
            
            % count the correct and wrong trial number
            for j = 1:length(fake_switch_index)
                ind = fake_switch_index(j);
                for i = 1:length(unique_degree{ind})
                    
                    % 1. related to perceptual PSE ponit (instead of stimulu 0 degree)
                    % "find dead-ahead" if the bias is hugh
                    if abs(Bias_psy{1,k,ind,last_step}) > min(nonzeros(abs(unique_degree{ind})))
                        [~,h_0] = min(abs(Bias_psy{1,k,ind,last_step}-unique_degree{ind})); % HH20130905
                    else
                        h_0 = find(unique_degree{ind} >= 0,1);
                    end
                    dead_head = unique_degree{ind}(h_0);
                    
                    % do not correct the dead head
                    dead_head = 0;
                    
                    if unique_degree{ind}(i) < dead_head  % if unique_degree{ind}(i)<0
                        correct_trial_number_related2PSE{n,k,ind}(i) = sum(trials_select_2T{n,k,ind,i} & choice_LEFT_2T{ind}); % choose L
                        wrong_trial_number_related2PSE{n,k,ind}(i) = sum(trials_select_2T{n,k,ind,i} & choice_RIGHT_2T{ind}); % choose R
                    else % 等于dead_head的时候因为是随机决定正确答案是左还是右，暂时默认正确答案是右
                        correct_trial_number_related2PSE{n,k,ind}(i) = sum(trials_select_2T{n,k,ind,i} & choice_RIGHT_2T{ind}); % choose R
                        wrong_trial_number_related2PSE{n,k,ind}(i) = sum(trials_select_2T{n,k,ind,i} & choice_LEFT_2T{ind}); % choose L
                    end            
                    % correct rate between degree
                    correct_rate_relate2PSE{n,k,ind}(i) = correct_trial_number_related2PSE{n,k,ind}(i) / (correct_trial_number_related2PSE{n,k,ind}(i) +  wrong_trial_number_related2PSE{n,k,ind}(i));
                    
                    % 2. related to stimulu 0 degree
                    if unique_degree{ind}(i) < 0
                        correct_trial_number_related2Stim{n,k,ind}(i) = sum(trials_select_2T{n,k,ind,i} & choice_LEFT_2T{ind}); % choose L
                        wrong_trial_number_related2Stim{n,k,ind}(i) = sum(trials_select_2T{n,k,ind,i} & choice_RIGHT_2T{ind}); % choose R
                    else % 等于dead_head的时候因为是随机决定正确答案是左还是右，暂时默认正确答案是右
                        correct_trial_number_related2Stim{n,k,ind}(i) = sum(trials_select_2T{n,k,ind,i} & choice_RIGHT_2T{ind}); % choose R
                        wrong_trial_number_related2Stim{n,k,ind}(i) = sum(trials_select_2T{n,k,ind,i} & choice_LEFT_2T{ind}); % choose L
                    end
                    % correct rate between degree
                    correct_rate_relate2Stim{n,k,ind}(i) = correct_trial_number_related2Stim{n,k,ind}(i) / (correct_trial_number_related2Stim{n,k,ind}(i) +  wrong_trial_number_related2Stim{n,k,ind}(i));
                    
                    % 3. relate to monkey choice: count the "Left" and "Right" choice trial number in one motion type
                    left_trial_number_2T{n,k,ind}(i) = sum(trials_select_2T{n,k,ind,i} & choice_LEFT_2T{ind}); % choose L in heading, choose CW in rotation, choose rotation in mt
                    right_trial_number_2T{n,k,ind}(i) = sum(trials_select_2T{n,k,ind,i} & choice_RIGHT_2T{ind}); % choose R in heading, choose CCW in rotation, choose heading in mt
                end
                
                % correct rate in specific motion type, ignore zero
                correct_rate_2T{n,k,ind} = sum(correct_trial_number_related2PSE{n,k,ind}(unique_degree{ind}~=0)) / sum(trials_num_2T{n,k,ind,step}(unique_degree{ind}~=0));
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 5. P value for microstimulation (GLE method)
    if length(unique_stimcondition) == 2
        for j=1:length(fake_switch_index)
            ind = fake_switch_index(j);
            
            % 1. for traditional psychometric function
            % 2. for ture error: R/(L+R+CW+CCW), CCW/(L+R+CW+CCW)
            for k=1:length(unique_stim_type)
                % solve the problem of extreme behavior in animals to obtain an acurate P value
                % [X,Y] = AdjustXY_for_probit(x,yctrl,ystim, repctrl, repstim)
                [X,Y] = AdjustXY_for_probit(unique_degree{ind},Psy_correct{1,k,ind,step},Psy_correct{2,k,ind,step}, trials_num_2T{1,k,ind,step}, trials_num_2T{2,k,ind,step});
                [X_true,Y_true] = AdjustXY_for_probit(unique_degree{ind},Psy_correct_true{1,k,ind,step},Psy_correct_true{2,k,ind,step}, trials_num_4T{1,k,ind,step}, trials_num_4T{2,k,ind,step});
                
                % no stim case: 用X与(B1:B2)拟合Y
                % bias: 50 pct PD threshoold: bias = (norminv(0.5)-B(1)/B(2))
                % threshold: 84% gaussian threshold: threshold = abs(norminv(0.84)/B(2))
                %
                % stim case：用X*(B1+B3;B2+B4)拟合Y
                % bias = (norminv(0.5)-(B(1)+B(3))) / (B(2)+B(4))
                % threshold = abs(norminv(0.84)/(B(2)+B(4)))
                opts = statset('glmfit');
                opts.MaxIter = 1000; % default value for glmfit is 100.
                [B,DEV,STATS,iterWarn] = glmfit_Lwh(X,Y,'binomial','link','probit','options',opts); % Lwh 20191214
                iterLimit_warn{k,ind,step} = iterWarn;  % 1: Iteration limit reached. ;  0: no warning   % Lwh 20191214
                
                % p_vallues
                p_mu{k,ind,step} = STATS.p(3); %B3的p值，代表stim case与control组在bias上有没有显著差别?
                p_theta{k,ind,step} = STATS.p(4);  %B4的p值，代表stim case与control组在threshold上有没有显著差别?
                
                % comparison for ture error
                % ????????????????????????????????
                [~,~,STATS_true,iterWarn] = glmfit_Lwh(X_true,Y_true,'binomial','link','probit','options',opts); % Lwh 20191214
                iterLimit_warn{k,ind,step} = iterWarn;  % 1: Iteration limit reached. ;  0: no warning   % Lwh 20191214
                p_mu_true{k,ind,step} = STATS_true.p(3); %B3的p值，代表stim case与control组在bias上有没有显著差别?
                p_theta_true{k,ind,step} = STATS_true.p(4);  %B4的p值，代表stim case与control组在threshold上有没有显著差别?
            end
        end
        
        if length(fake_switch_index) == 1
            if fake_switch_index == 1
                p_mu{k,2,step} = nan;
                p_theta{k,2,step} = nan;
                p_mu{k,3,step} = nan;
                p_theta{k,3,step} = nan;
                
                p_mu_true{k,2,step} = nan;
                p_theta_true{k,2,step} = nan;
                p_mu_true{k,3,step} = nan;
                p_theta_true{k,3,step} = nan;
            else
                p_mu{k,1,step} = nan;
                p_theta{k,1,step} = nan;
                p_mu{k,3,step} = nan;
                p_theta{k,3,step} = nan;
                
                p_mu_true{k,1,step} = nan;
                p_theta_true{k,1,step} = nan;
                p_mu_true{k,3,step} = nan;
                p_theta_true{k,3,step} = nan;
            end
        end
        
    else
        p_mu = nan;
        p_theta = nan;
        iterLimit_warn = nan;
        
        p_mu_true = nan;
        p_theta_true = nan;
    end 
end % end step

%% %%%%%%%%%%%%%%%%%%%%    plot traditional Psy curve   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% modified by HH20130901
% modified by Lwh 20170916
% modified by Lwh 20180702
% modified by Lwh 202103
% {n,m,k,j} → {n,k,j}
% n:condition (1 for normal, 2 for micro_stim)
% m:coherence (for coarse task)  % remove this by Lwh 202103
% k:stim_type(1 for visual, 2 for vestibular)
% ind:switch_index (1 for translation, 2 for rotation)
    

% color
% 2AFC,control
c{1,1,1} = [248 28 83]/255; % 红
c{1,1,2} = [41 89 204]/255; % 蓝
c{1,1,3} = [14 153 46]/255; % 绿

% 2AFC,stim
c{2,1,1} = [247 147 30]/255; % 橙色
c{2,1,2} = [247 147 30]/255; % 橙色
c{2,1,3} = [247 147 30]/255; % 橙色

% 4AFC
c_true{1} = [0.3 0.3 0.3]; % control, dark black, true error
c_true{2} = [0.7 0.7 0.7]; % stim, light black, true error

% symbol
sy{1} = 'o'; % control
sy{2} = '^'; % stim

% line
lin{1} = '-'; % control
% lin{2} = '--'; % stim
lin{2} = '-'; % stim

if IF_FIGURE
    for step = step_num %step_num
        figure_name = 10+step;
        if ishandle(figure_name); close(figure_name); end; figure(figure_name);
        set(figure_name,'Position', [10,25, 1500,900], 'Name', 'Psychmetic function', 'color','w');
        orient landscape;
        
       %% plot psychometric function
        nr=2;
        nc=3;
        % first row first column (1)
        % 2nd row first column (nc+1)
        % 3rd row first column (2*nc+1)
        
        ha = tight_subplot(nr,nc,[.08 .04],[.1 .04],[.04 .04]);
        for n = 1:length(unique_stimcondition) % microstior normal n=1:normal   n=2:microsti
            for k=1:length(unique_stim_type)
                for j=1:length(fake_switch_index)
                    ind = fake_switch_index(j);
                    
                    % Psychometric
                    axes(ha(ind));
                    if fake_switch_index(j)==1
                        gap = 0.1;
                        xlabel('Heading Angle(deg)');
                        ylabel('% Rightward choice');
                        title('Psychometric function');
                        if length(unique_coherence)>1 % coarse task
                            textx = max(unique_coherence)-4;
                        else % fine task
                            textx = max(unique_heading)-2;
                        end
                    elseif fake_switch_index(j)==2
                        gap = 0.01;
                        xlabel('Rotation Angle(%)');
                        ylabel('% CCW choice');
                        title('Psychometric function');
                        if length(unique_coherence)>1 % coarse task
                            textx = max(unique_coherence)-4;
                        else % fine task
                            textx = max(unique_rotation)-14;
                        end
                    else
                        gap = 0.1;
                        xlabel('Motion Type');
                        ylabel('Heading choice');
                        title('Psychometric function');
                        if length(unique_coherence)>1 % coarse task
                            textx = 0.5;
                        else % fine task
                            textx = 0.5;
                        end
                    end
                    
                    if length(unique_coherence) > 1 % coarse task
                        xlabel('Coherence');
                    end
                    
                    xh = min(unique_degree{ind}) : gap : max(unique_degree{ind});
                    hold on;
                    %                     h1=plot(unique_degree{ind}, Psy_correct{n,k,ind,step}, cdot{n,k,ind},'markerfacecolor',c{n,k,ind}(1),'markersize',10); % plot dot
                    %                     h2=plot(xh, cum_gaussfit(psy_perf{n,k,ind,step},xh), cline{n,k,ind},'linewidth',2.5); % plot fitting line
                    
                    h1=plot(unique_degree{ind}, Psy_correct{n,k,ind,step}, sy{n},'color',c{n,k,ind},'markerfacecolor',c{n,k,ind},'markersize',10); % plot dot
                    h2=plot(xh, cum_gaussfit(psy_perf{n,k,ind,step},xh), lin{n},'color',c{n,k,ind},'linewidth',2.5); % plot fitting line
                    
                    % plot true error
                    if ind~=3
                        h1=plot(unique_degree{ind}, Psy_correct_true{n,k,ind,step}, '+','color',c_true{n},'markerfacecolor',c_true{n},'markersize',10); % plot dot
                        h2=plot(xh, cum_gaussfit(psy_perf_true50{n,k,ind,step},xh), lin{n},'color',c_true{n},'linewidth',2.5); % plot fitting line
                    end


                    h3=plot([0 0],[0 1],'--k','linewidth',1.0);
                    h4=plot([min(xh) max(xh)],[0.5 0.5],'--k','linewidth',1.0);
                    xlim([min(unique_degree{ind})*1.1 max(unique_degree{ind})*1.1]);
                    ylim([0,1]);
                    set(gca,'xtickmode','auto');
                    set(gca,'xticklabelmode','auto');
                    set(gca,'yticklabelmode','auto');
                    text(min(unique_degree{ind}),1.1-n*0.1,sprintf('\\it\\mu_{psy} = \\rm%6.3g\\circ',psy_perf{n,k,ind,step}(1)),'color',c{n,k,ind},'fontsize',10);
                    text(min(unique_degree{ind}),0.9-n*0.1,sprintf('\\it\\sigma_{psy} = \\rm%6.3g\\circ',psy_perf{n,k,ind,step}(2)),'color',c{n,k,ind},'fontsize',10);
                    text(min(unique_degree{ind}),0.7-n*0.1,sprintf('Correct Rate = %0.3g',correct_rate_2T{n,k,ind}),'color',c{n,k,ind},'fontsize',10);
                    text(min(unique_degree{ind}),0.5-n*0.1,sprintf('\\mu = %0.3g', Bias_psy_true25{n,k,ind,step}),'color',c_true{n},'fontsize',10);
                    
                    try
                        if (length(unique_stimcondition) == 2 && n==1)
                            text(textx,0.3,sprintf('p_u = %0.5g',p_mu{k,ind,step}),'color','k');
                            text(textx,0.2,sprintf('p_{sig} = %0.5g', p_theta{k,ind,step}),'color','k');
                        end
                    end
                    
                    axis square
                    
                    axes(ha(ind+3));
                    hold on
                    %  plot preferred direction
                    if ind~=3
                        h1=plot(unique_degree{ind}, Psy_correct_pd{n,k,ind,step}, sy{n},'color',c{n,k,ind},'markerfacecolor',c{n,k,ind},'markersize',10); % plot dot
                        h2=plot(xh, cum_gaussfit(psy_perf_pd{n,k,ind,step},xh), lin{n},'color',c{n,k,ind},'linewidth',2.5); % plot fitting line
                    end
                    h3=plot([0 0],[0 1],'--k','linewidth',1.0);
                    h4=plot([min(xh) max(xh)],[0.5 0.5],'--k','linewidth',1.0);
                    xlim([min(unique_degree{ind})*1.1 max(unique_degree{ind})*1.1]);
                    ylim([0,1]);
                    set(gca,'xtick',[min(unique_degree{ind}),0,max(unique_degree{ind})]);
                    set(gca,'xticklabel',{'non-prefer','0','prefer'})
                    set(gca,'xticklabelmode','auto');
                    set(gca,'yticklabelmode','auto');
                    xlabel('non-preferred - preferred');
                    ylabel('% PD choice');
                    axis square
                end
            end
        end
        suptitle(sprintf('%s, %s %s %s, %g %s\n',FILE,'Co=',num2str(unique_coherence'),'%',repetition,'repetitions'));
        SetFigure(10)
        
%         % Figure Saving
%         figure(figure_name);saveas(gcf,['Z:\Data\MOOG\FigResults\CueS_Psy&Tuning\' FILE(1:end-4) '.fig']);
%         figure(figure_name);saveas(gcf,['Z:\Data\MOOG\FigResults\CueS_Psy&Tuning\' FILE(1:end-4) '.png'],'png');
    end % end step
end

%% %%%%%%%%%%%%%%%%%%%%   choice selection and motion type selection   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. every condition with all four targets (traditional psy curve only with related 2 targets)
% 2. every condition with correct/wrong MT choice
for k = 1:length(unique_stim_type)
    for n = 1:length(unique_stimcondition) % microstior normal n=1:normal   n=2:microsti
        for ind = 1:2 % 1 for H and 2 for R
            %             if length(fake_switch_index)==3 % only for 4 target task
            if length(unique_degree{ind})>1 % have done this task
                for i = 1:length(unique_degree{ind})
                    % select_cond默认为step最后的一个
                    % 每个condition中选4个选项的数目
                    choose_L{n,k,ind}(i) = sum(select_degree{n,k,ind,i} & choice == LEFT);
                    choose_R{n,k,ind}(i) = sum(select_degree{n,k,ind,i} & choice == RIGHT);
                    choose_CW{n,k,ind}(i) = sum(select_degree{n,k,ind,i} & choice == CW);
                    choose_CCW{n,k,ind}(i) = sum(select_degree{n,k,ind,i} & choice == CCW);
                    choose_all{n,k,ind}(i) =  sum(select_degree{n,k,ind,i});
                    
                    % 每个degree中选对、错motion type的数目
                    trials_type_wrong{n,k,ind,i} = logical(select_degree{n,k,ind,i} & choice_motiontype{3-ind}); % 选错motion type：heading task中选rotation，或者，rotation task中选heading   choice_type{3-ind}
                    % trials_select_2T{n,k,ind,i} = logical(select_degree_2T{n,k,ind,i} & choice_motiontype{ind}); % Heading task中选heading target， R task中选R target的总数
                    mt_right_num{n,k,ind}(i) = sum(trials_select_2T{n,k,ind,i});  % Heading task中选Heading数目， R task中选R数目
                    mt_wrong_num{n,k,ind}(i) = sum(trials_type_wrong{n,k,ind,i}); % Heading task中选rotation 数目， rotation task中选H数目
                    
                    mt_right_pro{n,k,ind}(i) = mt_right_num{n,k,ind}(i) / choose_all{n,k,ind}(i);  % Heading task中选Heading数目， R task中选R数目
                    mt_wrong_pro{n,k,ind}(i) = mt_wrong_num{n,k,ind}(i) / choose_all{n,k,ind}(i); % Heading task中选rotation 数目， rotation task中选H数目
                end
                
                % 生成一个choice矩阵， target selection
                % Heading task中R choice, L choice, CCW, CW个数（横坐标是heading）
                % Rotation task中R choice, L choice, CCW, CW个数 （横坐标是rotation angle,从CW-0-CCW）
                choice_array{n,k,ind}(1,:) = choose_R{n,k,ind};
                choice_array{n,k,ind}(2,:) = choose_L{n,k,ind};
                choice_array{n,k,ind}(3,:) = choose_CCW{n,k,ind};
                choice_array{n,k,ind}(4,:) = choose_CW{n,k,ind};
                
                % 按照两条轴逆时针顺序R-CW-L-CCW
                choice_array2{n,k,ind}(1,:) = choose_R{n,k,ind};
                choice_array2{n,k,ind}(2,:) = choose_CW{n,k,ind};
                choice_array2{n,k,ind}(3,:) = choose_L{n,k,ind};
                choice_array2{n,k,ind}(4,:) = choose_CCW{n,k,ind};
                
                % note! in one task type task, e.g. fine rotation task,
                % monkey may choose the translation target, we ignore them
                % in analyses, but we can see the choice in choice_array 
                
            else % does not do this task, e.g. fine rotation only task, fill with nan, keep same size with the other task
                choice_array{n,k,ind} = nan(4,length(unique_degree{3-ind}));
                choice_array2{n,k,ind} = nan(4,length(unique_degree{3-ind}));
                mt_right_pro{n,k,ind} = nan;
                mt_wrong_pro{n,k,ind} = nan;
                choose_all{n,k,ind} = nan;
            end
        end
    end
end

%% Correct rate, (not consider 0 degree)
% 1. correct rate for 2T
% 2. correct rate for 4T
for k = 1:length(unique_stim_type)
    for n = 1:length(unique_stimcondition) % microstior normal n=1:normal   n=2:microsti
        % T or R
        for ind = 1:2
            if length(unique_degree{ind})>1 % have done this task
                TRdegree0 = unique_degree{ind}' == 0;
                degree_num = length(unique_degree{ind});
                mid_point = (degree_num+1)/2; % due with whether with zero degree
                temp = 1:degree_num;
                left_part = temp(temp<mid_point);
                right_part = temp(temp>mid_point);
                
                % first row, "right" choice + second row, "left" choice
                correct_trial_num =  sum( choice_array{n,k,ind}(2*(ind-1)+1,right_part) ) + sum( choice_array{n,k,ind}(2*(ind-1)+2,left_part) );
                all_trial_num = sum(sum(choice_array{n,k,ind}([2*(ind-1)+1 2*(ind-1)+2],~TRdegree0))); % two row for this type
                all_trial_true_num = sum(sum(choice_array{n,k,ind}(:,~TRdegree0))); % four row for this type, true error
                correct_rate(n,ind) = correct_trial_num / all_trial_num * 100; % 1. correct rate for 2T
                correct_rate_true(n,ind) = correct_trial_num / all_trial_true_num * 100; % 2. correct rate for 4T
                
            else % does not do this task, e.g. fine rotation only task, fill with nan, keep same size with the other task
                correct_rate(n,ind) = nan;
                correct_rate_true(n,ind) = nan;
            end
        end
        
        % MT
        ind = 3;
        if length(fake_switch_index) == 3 % 4T task
            MTdegree0 = unique_degree{ind}' == 0;
            trans_choose_trans = choice_array{n,k,1}([1 2],~MTdegree0);  % correct choice
            trans_choose_all = choice_array{n,k,1}(:,~MTdegree0); % all choice
            rot_choose_rot = choice_array{n,k,2}([3 4],~MTdegree0);  % correct choice
            rot_choose_all = choice_array{n,k,2}(:,~MTdegree0); % all choice
            correct_rate(n,ind) = sum([trans_choose_trans(:);rot_choose_rot(:)]) / sum([trans_choose_all(:);rot_choose_all(:)]) * 100;
            correct_rate_true(n,ind) = correct_rate(n,ind); % same!
        else % 2T task
            correct_rate(n,ind) = nan;
            correct_rate_true(n,ind) = nan;
        end
    end
end

    
%% choice changed: (stim-control)/control
% 1. Heading/Rotation/MT task中，电刺激前后，每一个degree与非电刺激0度的correct rate做chi-square test to fine the similar difficulty trial for monkey (Salzman, et al, 1992)
% 2. 比较电刺激前后choice的变化
% 3. p by chi-square test between control and microstimulation trial
choice_change_pro_diff = arrayfun(@(x) nan(1,2),[1:3],'UniformOutput',0);
choice_change_pro_in0 = arrayfun(@(x) nan(1,2),[1:3],'UniformOutput',0);
choice_change_pro_around0 = arrayfun(@(x) nan(1,2),[1:3],'UniformOutput',0);

choice_change_index_diff = arrayfun(@(x) nan(1,3),[1:length(unique_stim_type)]','UniformOutput',0);
choice_change_index_in0 = arrayfun(@(x) nan(1,3),[1:length(unique_stim_type)]','UniformOutput',0);
choice_change_index_around0 = arrayfun(@(x) nan(1,3),[1:length(unique_stim_type)]','UniformOutput',0);

chi_square_p_diff = arrayfun(@(x) nan(1,3),[1:length(unique_stim_type)]','UniformOutput',0);
chi_square_p_in0 = arrayfun(@(x) nan(1,3),[1:length(unique_stim_type)]','UniformOutput',0);
chi_square_p_around0 = arrayfun(@(x) nan(1,3),[1:length(unique_stim_type)]','UniformOutput',0);

for k = 1:length(unique_stim_type)
    for j=1:length(fake_switch_index)
        ind = fake_switch_index(j);
        
        % 1. find difficult degree in control and stim trial
        % 1.1 find 0 degree in no-stim, change to "find dead-ahead" if the bias is hugh
        if abs(Bias_psy{1,k,ind,last_step}) > min(nonzeros(abs(unique_degree{ind})))
            [~,h_0] = min(abs(Bias_psy{1,k,ind,last_step}-unique_degree{ind})); % HH20130905
        else
            h_0 = find(unique_degree{ind} >= 0,1);
        end
        dead_head = unique_degree{ind}(h_0);
        
        if zero_cond_num > 0 % 如果有0度
            correct_control = correct_trial_number_related2PSE{1,k,ind}(unique_degree{ind}' == dead_head); % 非电刺激trial 0度
            wrong_control = wrong_trial_number_related2PSE{1,k,ind}(unique_degree{ind}' == dead_head); % 非电刺激trial 0度
        else % 没有0度选取中间两个角度的平均
            correct_control = mean(correct_trial_number_related2PSE{1,k,ind}([length(unique_degree{ind})/2 length(unique_degree{ind})/2+1]));
            wrong_control = mean(wrong_trial_number_related2PSE{1,k,ind}([length(unique_degree{ind})/2 length(unique_degree{ind})/2+1]));
        end
        
        % 1.2 chi-squre with other degree
        for n = 1:length(unique_stimcondition) % microstior normal n=1:normal   n=2:microsti
            for i = 1:length(unique_degree{ind})
                correct_test = correct_trial_number_related2PSE{n,k,ind}(i);
                wrong_test = wrong_trial_number_related2PSE{n,k,ind}(i);
%                 chi_square_Lwh([wrong_control correct_control; wrong_test correct_test],'right')
                chi_square_p{n,k,ind}(i) = chi_square_Lwh([wrong_control correct_control; wrong_test correct_test],'right'); % 联列表 p<0.05，left：wrong_text>>wrong control: 更难了
            end
        end
        % 1.3 p > 0.1: not significantly different from random performance on nonstimulated trials (find the difficult degree)
        % p > 0, h = 0: 不拒绝0假设，右侧单尾：wrong control并没有比wrong test多（也就是wrong contol<=wrong test，也即是test条件下wrong 更多了，更难了）
        difficult_degree{k,ind} = logical(chi_square_p{1,k,ind}>0.05 | chi_square_p{2,k,ind}>0.05);
            
        % 2. count left/CW/rotation choice and right/CCW/heading choice with or without microstimulation
        % choice_change{k,ind} =
        %          L/CW/rot    R/CCW/trans
        % nonstim    a            b
        %  stim      c            d
        % 注意ind=3时，第一列时Rot choice变化，第二列时Trans choice的变化
        degree_num = length(unique_degree{ind});
        for n = 1:length(unique_stimcondition) % microstior normal n=1:normal   n=2:microstim
            % 2.1 choice_change in difficult degree
            choice_change_in_diff{k,ind}(n,:) = [sum(left_trial_number_2T{n,k,ind}(difficult_degree{k,ind})) sum(right_trial_number_2T{n,k,ind}(difficult_degree{k,ind}))];
            
            % 2.2 choice_change in zero degree
            if zero_cond_num > 0 % 如果有0度
                choice_change_in_0{k,ind}(n,:) = [left_trial_number_2T{n,k,ind}(unique_degree{ind}' == 0) right_trial_number_2T{n,k,ind}(unique_degree{ind}' == 0)]; % in zero degree
            else % 没有0度选取中间两个角度的平均
                choice_change_in_0{k,ind}(n,:) = [mean(left_trial_number_2T{n,k,ind}([floor(degree_num/2):(floor((degree_num+1)/2)+1)])),...
                    mean(right_trial_number_2T{n,k,ind}([floor(degree_num/2):(floor((degree_num+1)/2)+1)]))]; % 此时与choice_change_around_0类似，只是取平均
            end
            
            % 2.3 choice_change around zero degree (0度左右各多加一个角度)
            choice_change_around_0{k,ind}(n,:) = [sum(left_trial_number_2T{n,k,ind}([floor(degree_num/2):(floor((degree_num+1)/2)+1)])),...
                sum(right_trial_number_2T{n,k,ind}([floor(degree_num/2):(floor((degree_num+1)/2)+1)]))]; % e.g. [3 4 5]in 7 degrees or [3 4] in 6 degrees
        end
        
        
        % choice_change_around_0{k,3}只有T和R的变化，choice_change_around_0{k,1/2}只有motiontype内的变化，
        % 现在总结ind=3时4个choice的变化
        if ind == 3
            % use choice_array2 to count
            % because in choice_array2{n,k,2} (rotation degress axis), the zero
            % degree is same in choice_array2{n,k,1} (translation degress axis),
            % adjust the choice_array2{n,k,2} zero degree to all ZERO for sum up
            choice_array_adj = choice_array2;
            for n = 1:length(unique_stimcondition) % microstior normal n=1:normal   n=2:microstim
                if zero_cond_num > 0 % 如果有0度
                    choice_array_adj{n,k,2}(:,unique_degree{2}' == 0) = zeros(4,1);
                end
            end
            % *************** only around zero ****************
            % 每一行：R-CW-L-CCW; 第一列ctrl，第二列stim
            for n = 1:length(unique_stimcondition) % microstior normal n=1:normal   n=2:microstim
                Four_choice_num(:,n) = sum(choice_array_adj{n,k,1}(:,[floor(degree_num/2):(floor((degree_num+1)/2)+1)]) + ...
                    choice_array_adj{n,k,2}(:,[floor(degree_num/2):(floor((degree_num+1)/2)+1)]),2);  
            end
            % use case/N
            Four_choice_num_pro = Four_choice_num./sum(Four_choice_num(:,1));

            % plot
            % ****************************************************************
             IF_FIGURE_this = true;
            if IF_FIGURE_this % IF_FIGURE
            figure;
            subplot(1,2,1) % four choice change
            hold on
            for n = 1:length(unique_stimcondition)
                vec1 = [Four_choice_num_pro(1,n),0]; % 右，right
                vec2 = [0 -Four_choice_num_pro(2,n)]; % 下，CW
                vec3 = [-Four_choice_num_pro(3,n),0]; % 左，left
                vec4 = [0 Four_choice_num_pro(4,n)]; % 上，CCW
                if n == 1
                    plot([vec1(1) vec2(1) vec3(1) vec4(1) vec1(1)],[vec1(2) vec2(2) vec3(2) vec4(2) vec1(2)],'k.-');
                else
                    plot([vec1(1) vec2(1) vec3(1) vec4(1) vec1(1)],[vec1(2) vec2(2) vec3(2) vec4(2) vec1(2)],'g.-');
                end
            end
            xlim([-0.5 0.5]);
            ylim([-0.5 0.5]);
            axis square
            
            set(gca,'xtick',[-0.25 0 0.25],'xticklabel',{'Left' '' 'Right'})
            set(gca,'ytick',[-0.25 0 0.25],'yticklabel',{'CW' '' 'CCW'})
            plot([-0.5 0.5],[0 0],'k-')
            plot([0 0],[-0.5 0.5],'k-')
            SetFigure(15)
            legend('ctrl','stim');
            
            subplot(1,2,2) % t/r choice change, 顺序:ctrlR,stimR,ctrlT,stimT
            temp = choice_change_around_0{k,3}(:) ./ sum(Four_choice_num(:,1));
            bar(temp);
            set(gca, 'xticklabel',{'ctrlR','stimR','ctrlT','stimT'});
            ylabel('Choice frequency (%)');
            hold on
            plot([0 5],[0.5 0.5],'k-');
            
            % chi-square test, 下面code也求p value (chi_square_p_around0)了，这里画图需要暂时先算一次
            chis_p = chi_square_Lwh(choice_change_around_0{k,3}); % around zero degree
            text(1, 0.8, sprintf('chi-square test, p = %.2g', chis_p));
            ylim([0 1]);
            SetFigure(15)
            end
        end
        
        % 2.4 choice change proportion (chocie-change index, CCI) = (stim-control)/control
        % heading task: (stimL-ctrlL)/ctrlL, (stimR-ctrlR)/ctrlR             (within motion type)
        % rotation task: (stimCW-ctrlCW)/ctrlCW, (stimCCW-ctrlCCW)/ctrlCCW   (within motion type)
        % mt task: (stimR-ctrlR)/ctrlR, (stimT-ctrlT)/ctrlT                  (between motion type)
        choice_change_pro_diff{k,ind} = (choice_change_in_diff{k,ind}(2,:) - choice_change_in_diff{k,ind}(1,:)) ./ choice_change_in_diff{k,ind}(1,:);
        choice_change_pro_in0{k,ind} = (choice_change_in_0{k,ind}(2,:) - choice_change_in_0{k,ind}(1,:)) ./ choice_change_in_0{k,ind}(1,:);
        choice_change_pro_around0{k,ind} = (choice_change_around_0{k,ind}(2,:) - choice_change_around_0{k,ind}(1,:)) ./ choice_change_around_0{k,ind}(1,:);
        
        % 2.5 turn chi2 to scale from "left" to 'right' and used as choice change index
        % (ad-bc)>0: choose more "Right" or less "Left" or both(but much more Right) in stim  / more CCW / more trans
        % (ad-bc)<0: choose more "Left" or less "Right" or both in stim
        % 卡方值 chi2 = n*(ad-bc)^2 / [(a+b)(c+d)(a+c)(b+d)] 
        % choice_change_index = sign * sqrt(chi2)  (直接用上式子计算，不考虑频数是否适合不校正,函数chi_square_Lwh得出的chi2是修正过的)
        % 对于MT(ind=3)时候，choice_change_index越负代表电刺激后选择Rotation target变多了，越正代表电刺激后选择Translation target变多了，=0代表电刺激没有作用
        aa1 = choice_change_in_diff{k,ind}(1,1); aa2 = choice_change_in_0{k,ind}(1,1); aa3 = choice_change_around_0{k,ind}(1,1);
        bb1 = choice_change_in_diff{k,ind}(1,2); bb2 = choice_change_in_0{k,ind}(1,2); bb3 = choice_change_around_0{k,ind}(1,2);
        cc1 = choice_change_in_diff{k,ind}(2,1); cc2 = choice_change_in_0{k,ind}(2,1); cc3 = choice_change_around_0{k,ind}(2,1);
        dd1 = choice_change_in_diff{k,ind}(2,2); dd2 = choice_change_in_0{k,ind}(2,2); dd3 = choice_change_around_0{k,ind}(2,2);
        
        sign1 = aa1*dd1-bb1*cc1;
        sign2 = aa2*dd2-bb2*cc2;
        sign3 = aa3*dd3-bb3*cc3;
        
        choice_change_index_diff{k}(ind) = sign(sign1) * sqrt((aa1+bb1+cc1+dd1)*sign1^2 / ((aa1+bb1)*(cc1+dd1)*(aa1+cc1)*(bb1+dd1)) );
        choice_change_index_in0{k}(ind) = sign(sign2) * sqrt((aa2+bb2+cc2+dd2)*sign2^2 / ((aa2+bb2)*(cc2+dd2)*(aa2+cc2)*(bb2+dd2)) );
        choice_change_index_around0{k}(ind) = sign(sign3) * sqrt((aa3+bb3+cc3+dd3)*sign3^2 / ((aa3+bb3)*(cc3+dd3)*(aa3+cc3)*(bb3+dd3)) );

        % 3. chi-squre between stim and nostim
        [chi_square_p_diff{k}(ind)] = chi_square_Lwh(choice_change_in_diff{k,ind}); % in difficult degree
        [chi_square_p_in0{k}(ind)] = chi_square_Lwh(choice_change_in_0{k,ind}); % in zero degree
        [chi_square_p_around0{k}(ind)] = chi_square_Lwh(choice_change_around_0{k,ind}); % around zero degree
    end
end

%% 随时间(trial)电刺激的变化 Salzman,et al. 1992
% 1. 选取电刺激PSE shift方向为正
% 2. slide_window
% 3. only plot positive proportion

% 注意，在4-Target task中slide-window是在实际的trial序列中进行，而不是在当前motion
% type内的trial中进行，也就是figure中的每个子图的横坐标是实际的trial数目，而不是当前motion
% type内trial的数目

% slide_window = one_repetition; % actual trial num, not this specific motion type trial num
% slide_step = round(one_repetition/3);
slide_window = 30; % actual trial num, not this specific motion type trial num
slide_step = 10;

% prepare 2*1*3 cell with nan
slide_mid_trial = arrayfun(@(x) nan(1,1),nan(2,1,3),'un',0);
slide_positive_pro = arrayfun(@(x) nan(1,1),nan(2,1,3),'un',0);
slide_negative_pro = arrayfun(@(x) nan(1,1),nan(2,1,3),'un',0);
slide_CR = arrayfun(@(x) nan(1,1),nan(2,1,3),'un',0);

for k = 1:length(unique_stim_type)
    for j=1:length(fake_switch_index)
        ind = fake_switch_index(j);
        for n = 1:length(unique_stimcondition) % microstior normal n=1:normal   n=2:microsti
            if ind < 3 % heading or rotation
                select = logical(stim_type == unique_stim_type(k) & stimcondition == unique_stimcondition(n) & switch_index == unique_switch_index(j)); %所有0度task
            else % MT
                select = logical(stim_type == unique_stim_type(k) & stimcondition == unique_stimcondition(n));
            end
            
            % only in slide_window
            slide_start = 1;
            slide_end = slide_start + slide_window - 1;
            slide_num = 1;
            while slide_end <= length(coherence)
                select_slide_trials = false(size(coherence));
                select_slide_trials(slide_start:slide_end) = true;
                select_this_slide = logical(select & select_slide_trials); % only select choice in this motion tpye

                slide_mid_trial{n,k,ind}(slide_num) = (slide_start + slide_end) / 2 ; % for marker
                
                % 电刺激trial选择最终PSE偏移方向的目标点的比例
                if Bias_psy{2,k,ind,last_step} > Bias_psy{1,k,ind,last_step} % bias to "Left", set "left" to positive
                    slide_positive_pro{n,k,ind}(slide_num) = sum(select_this_slide & choice_LEFT_2T{ind}) / sum(select_this_slide & (choice_LEFT_2T{ind}| choice_RIGHT_2T{ind}));
                    slide_negative_pro{n,k,ind}(slide_num) = sum(select_this_slide & choice_RIGHT_2T{ind}) / sum(select_this_slide & (choice_LEFT_2T{ind}| choice_RIGHT_2T{ind}));
                else % bias to "right", set "right" to positive
                    slide_positive_pro{n,k,ind}(slide_num) = sum(select_this_slide & choice_RIGHT_2T{ind}) / sum(select_this_slide & (choice_LEFT_2T{ind}| choice_RIGHT_2T{ind}));
                    slide_negative_pro{n,k,ind}(slide_num) = sum(select_this_slide & choice_LEFT_2T{ind}) / sum(select_this_slide & (choice_LEFT_2T{ind}| choice_RIGHT_2T{ind}));
                end
                
                % 电刺激trial正确率变化
                choose_left = sum(choice_LEFT_2T{ind}(select_this_slide) & degree{ind}(select_this_slide)<0); % 左刺激选左
                set_left = sum(degree{ind}(select_this_slide)<0); % 左刺激
                choose_right = sum(choice_RIGHT_2T{ind}(select_this_slide) & degree{ind}(select_this_slide)>0); % 右刺激选右
                set_right = sum(degree{ind}(select_this_slide)>0); % 右刺激
                slide_CR{n,k,ind}(slide_num) = (choose_left + choose_right) / (set_left + set_right);
                
                slide_start = slide_start + slide_step;
                slide_end = slide_start + slide_window - 1;
                slide_num = slide_num + 1;
            end
        end
    end
end

% % plot
% figure
% for k = 1:length(unique_stim_type)
%     for j=1:length(fake_switch_index)
%         ind = fake_switch_index(j);
%         subplot(1,3,ind)
%         for n = 1:length(unique_stimcondition) % microstior normal n=1:normal   n=2:microsti
%             hold on
%             plot(slide_mid_trial{n,k,ind},slide_positive_pro{n,k,ind},c{n,k,ind}); % only plot positive pproportion
%         end
%     end
% end

%% get eyetrace data
% time information
h = data.htb_header{EYE_DB};	
eye_timeWin = 1000 * (h.skip + 1) / (h.speed_units / h.speed); % in ms % added by YXF from HH's memory guided saccade protocol;20140905
h2 = data.htb_header{EVENT_DB};	
event_timeWin = 1000 * (h2.skip + 1) / (h2.speed_units / h2.speed);  % in ms % added by YXF from HH's memory guided saccade protocol;20140905

% event data
event_in_bin = squeeze(data.event_data(:,:,select_trials));   % TrialNum * 5000
TrialBegT = find(event_in_bin==TRIAL_START_CD);  %01
FPOnT = find(event_in_bin == FP_ON_CD);  %02
fixInT = find(event_in_bin == IN_FIX_WIN_CD);  %03
stimOnT = find(event_in_bin == VSTIM_ON_CD);  %04
stimOffT = find(event_in_bin == VSTIM_OFF_CD);  %05
targetOnT =find(event_in_bin==TARGS_ON_CD);  %06
sacOnT = find(event_in_bin == SACCADE_BEGIN_CD); %07
sacOffT = find((event_in_bin == IN_T1_WIN_CD)|(event_in_bin == IN_T2_WIN_CD)|(event_in_bin == IN_T3_WIN_CD)|(event_in_bin == IN_T4_WIN_CD));
TrialEndT = find(event_in_bin == TRIAL_END_CD); %0F


% 按照eye trace取样率，bin(5ms = 1bin)  （eye trace取样率是200，我们的data是1000）
FPOnBin = ceil(FPOnT*event_timeWin/eye_timeWin);
stimOnBin = ceil(stimOnT*event_timeWin/eye_timeWin);
stimOffBin = ceil(stimOffT*event_timeWin/eye_timeWin);
fixInBin = ceil(fixInT*event_timeWin/eye_timeWin);
% sacOnBin = ceil(sacOnT*event_timeWin/eye_timeWin);
% sacOffBin = ceil(sacOffT*event_timeWin/eye_timeWin);
TrialEndBin = ceil(TrialEndT*event_timeWin/eye_timeWin);

% eye data, in bin
if data.one_time_params(LEFT_EYE_X_CHANNEL) > 0
    LEYE_H = data.one_time_params(LEFT_EYE_X_CHANNEL);
end

if data.one_time_params(LEFT_EYE_Y_CHANNEL) > 0
    LEYE_V = data.one_time_params(LEFT_EYE_Y_CHANNEL);
end

if data.one_time_params(RIGHT_EYE_X_CHANNEL) > 0
    REYE_H = data.one_time_params(RIGHT_EYE_X_CHANNEL);
end

if data.one_time_params(RIGHT_EYE_Y_CHANNEL) > 0
    REYE_V = data.one_time_params(RIGHT_EYE_Y_CHANNEL);
end

% eye_data = data.eye_data(:,:,select_trials); % bin
% use calibration eye data
eye_data = data.eye_positions_calibrated(:,:,select_trials); % bin
 
eye_data_LH = squeeze(data.eye_data(LEYE_H,:,select_trials));
eye_data_LV = squeeze(data.eye_data(LEYE_V,:,select_trials));
eye_data_RH = squeeze(data.eye_data(REYE_H,:,select_trials));
eye_data_RV = squeeze(data.eye_data(REYE_V,:,select_trials));

stimOnBin_trial = mod(stimOnBin,1000);
stimOffBin_trial = mod(stimOffBin,1000);
plot_LR = 2; % 1 = left; 2 = right
IF_eyetrace_FIGURE = 0;
for n = 1:sum(select_trials)
    %%% Calibrate eye offset (average 200 ms after VSTIM_ON_CD(04)
    %     eye_offsetLX = mean(eye_data_LH(stimOnBin_trial(n):stimOnBin_trial(n)+200/eye_timeWin,n));
    %     eye_offsetLY = mean(eye_data_LV(stimOnBin_trial(n):stimOnBin_trial(n)+200/eye_timeWin,n));
    eye_offsetRX = mean(eye_data_RH(stimOnBin_trial(n):stimOnBin_trial(n)+200/eye_timeWin,n));
    eye_offsetRY = mean(eye_data_RV(stimOnBin_trial(n):stimOnBin_trial(n)+200/eye_timeWin,n));
    eye_offsetRX = 0;
    eye_offsetRY = 0;
    
    %     eye_data_LH(:,n) = eye_data_LH(:,n) - eye_offsetLX;
    %     eye_data_LV(:,n) = eye_data_LV(:,n) - eye_offsetLY;
    eye_data_RH(:,n) = eye_data_RH(:,n) - eye_offsetRX;
    eye_data_RV(:,n) = eye_data_RV(:,n) - eye_offsetRY;
    
    if IF_eyetrace_FIGURE %%% plot 'visual on' period
        figure(542);
        hold on
        color_this = linspace(0,1,length(stimOnBin_trial:stimOffBin_trial));
        if plot_LR == 1
            scatter(eye_data_LH(stimOnBin_trial:stimOffBin_trial,n),eye_data_LV(stimOnBin_trial:stimOffBin_trial,n),[],color_this,'o','filled'); % from yellow to blue
            plot(eye_data_LH(stimOnBin_trial:stimOffBin_trial,n),eye_data_LV(stimOnBin_trial:stimOffBin_trial,n),'k-');
        elseif plot_LR == 2
            % plot all trial
            plot(eye_data_RH(stimOnBin_trial-500/eye_timeWin:stimOffBin_trial+200/eye_timeWin,n),eye_data_RV(stimOnBin_trial-500/eye_timeWin:stimOffBin_trial+200/eye_timeWin,n),'-'); % 刺激前500ms到刺激后200ms
%             scatter(eye_data_RH(stimOnBin_trial:stimOffBin_trial,n),eye_data_RV(stimOnBin_trial:stimOffBin_trial,n),[],color_this,'o','filled'); % from yellow to blue
        end
    end
    
    %%% save eye data for each trial during 'visual on' period
    eye_data_VisualOn_X(n,:) = eye_data_RH(stimOnBin_trial(n):stimOffBin_trial(n),n)'; % 每一行一个trial，每一列一个bin
    eye_data_VisualOn_Y(n,:) = eye_data_RV(stimOnBin_trial(n):stimOffBin_trial(n),n)';
end


%% Save and Output
% Reorganized. HH20141124
config.batch_flag = batch_flag;

% % Output information for test. HH20160415
% if isempty(batch_flag)
%     config.batch_flag = 'test.m';
%     disp('Saving results to \batch\test\ ');
% end

unique_coherence = unique_coherence';
unique_heading = unique_heading';
unique_rotation = unique_rotation';

% save last step PSE shift to excel: stim - ctrl
for ind = 1:3 % h r mt
    PSE_shift{ind} = Bias_psy{2,:,ind,last_step} - Bias_psy{1,:,ind,last_step};
    Threshold_shift{ind} = Thresh_psy{2,:,ind,last_step} - Thresh_psy{1,:,ind,last_step};
end


%%%%%%%%%%%%%%%%%%%%% Change here %%%%%%%%%%%%%%%%%%%%%%%%%%%%
result = PackResult(FILE, SpikeChan, repetition, rep_0, unique_stimcondition, unique_coherence, unique_stim_type, ... % Obligatory!!
    eyepos_x, eyepos_y, zero_cond_num,last_step,stars_type,...
    unique_heading, unique_rotation, unique_degree,unique_switch_index,...
    Psy_correct, Bias_psy, Thresh_psy, psy_perf, ...
    p_mu, p_theta,PSE_shift, Threshold_shift,...
    Psy_correct_true, Bias_psy_true50, Thresh_psy_true50, psy_perf_true50, Bias_psy_true25,...
    p_mu_true, p_theta_true,...
    choice_array,choice_array2,...
    correct_rate, correct_rate_true, ...
    mt_right_pro,mt_wrong_pro,choose_all,...
    iterLimit_warn,Is_step,step_num,...
    choice_change_pro_diff,choice_change_pro_in0,choice_change_pro_around0,...
    choice_change_index_diff,choice_change_index_in0,choice_change_index_around0,...
    chi_square_p_diff, chi_square_p_in0, chi_square_p_around0,...
    slide_mid_trial,slide_positive_pro,slide_negative_pro,slide_CR,...
    plot_LR,eye_data_VisualOn_X, eye_data_VisualOn_Y,stimOnBin,stimOffBin,stimOnBin_trial,stimOffBin_trial,...
    raw_data);

    
config.batch_flag = batch_flag;
config.suffix = 'CueS_performance';
config.xls_column_begin = 'T_PSE_shift';
config.xls_column_end = 'MT_threshold_p';

% Figures to save
config.save_figures = gcf;

% Only once   %不需要循环输出到excel的数据
config.sprint_once_marker = '';
config.sprint_once_contents = '';

% Loop across each stim_type  %按照几列为循环输出到excel中
config.sprint_loop_marker = {'gg';
    'gg';};

config.sprint_loop_contents = {
    'result.PSE_shift{ind},result.p_mu{1,ind,last_step}',;
    'result.Threshold_shift{ind},result.p_theta{1,ind,last_step}';};

config.append = 1; % Overwrite or append

IF_EXCEL = 0;
SaveResult(config, result,IF_FIGURE, IF_EXCEL);

if isempty(config.batch_flag)
    disp('This is Done~ ');
end

return;
