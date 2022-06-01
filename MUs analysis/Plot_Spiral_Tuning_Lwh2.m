% % Plot_Spiral_Tuning_Lwh.m -- Plots response as a function of azimuth and rotation for
% % spiral tuning task
% %--	YG, 07/12/04
% %     HH, 2013
% %  Lwh 201707
% %  Lwh 201710 Changed Spial tuning potocal by added spiral condition(CW+Exp/Con, CCW+Exp/Con)
% %  Lwh 201806 combine heading tuning and spiral tuning analyses
% %  Lwh 201812, pref is mean preference self-motion, not object motion
%        Heading pref L, heading = -90, azmuth = 180, dot motion:RIGHT!!!!
%        Heading pref R, heading = 90, azmuth = 0, dot motion:LEFT!!!!
%        Rotation pref CW, rotation = 90, dot motion:CCW!!!!
%        Rotation pref CCW, rotation = 270, dot motion:CW!!!!
%    原程序dot motion CW设定为pref CW是错误的！！！改正
%
% % Lwh  1, 2021  session 319, m7c1274(含)以后的task改用圆点视觉刺激，并且调整旋转中心为fixation point，用stars_type定义
%
%
%% Lwh 20190410 重新整理 spiral tuning code
% azimuth   rotation_axis     stimulus         self-motion     polar angle(azi)
%    0         0         from Right to Left       Right                0
%   180        0         from Left to Right       Left                180
%   90         0                Exp              Forward               90
%   270        0                Con             Backward              270
%   0          90               CCW                CW                 180
%   0         270               CW                 CCW                 0
%   90        90              Exp+CCW           Forward+CW            135
%   90        270             Exp+CW            Forward+CCW            45
%   270       90              Con+CCW           Backward+CW           225
%   270       270             Con+CW            Backward+CCW          315

% 最终定义prefer tuning为prefer self-motion 方向

%% Lwh 20201102 modify spiral stimuli, 统一translation和rotation，使其能连续变化
% 原来程序中得rotation part中，rotation polar
% angle不管怎么样（除了纯旋转），向前分量一直是distance=0.2m，旋转分量一直时最大旋转角度RotAmplitude 45°/s （参考下表，全是固定值）

% Translation part, Distance = 2m
% Azimuth(polar angle)   RotAzimuth       stimulus         self-motion(camera)
%       0                   888       from Right to Left         Right
%      45                   888
%      90                   888            Expasion             Forward
%     135                   888                                 Forward+L
%     180                   888       from Left to Right          Left
%     225                   888                                 Backward+L
%     270                   888           Contraction            Backward
%     315                   888                                 Backward+R
%
% 向前分量 = Distance * sin(Azi)
% 左右分量 = Distance * cos(Azi)
%
% Rotation part, RAmplitude = 45°
% Azimuth(polar angle)   RotAzimuth       stimulus         self-motion(camera)
%     888                     0              CW                   CCW
%     888                    45            Exp+CW             Forward+CCW
%     888                    90              Exp                Forward       → 与Translation part重复，Tempo中删掉次条件
%     888                   135            Exp+CCW            Forward+CW
%     888                   180              CCW                  CW
%     888                   225            Con+CCW            Backward+CW
%     888                   270              Con               Backward       → 与Translation part重复，Tempo中删掉次条件
%     888                   315            Con+CW             Backward+CCW
%
% 向前分量 = Distance * sin(RotAzi)
% 旋转分量 = RAmplitude * cos(RotAzi)
%
% 共8个translation条件+6 rotatin条件+1 null = 15 condition


% %-----------------------------------------------------------------------------------------------------------------------
function Plot_Spiral_Tuning_Lwh2(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, batch_flag);
TEMPO_Defs;
Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

IF_FIGURE = 1;  % 1：画图； 0：不画图
IF_Fit = 1;
IF_DrawRecLoc = 0;
IF_EXCEL = 1;

%% get data
%get the column of values for azimuth and elevation and stim_type
temp_azimuth = data.moog_params(AZIMUTH,:,MOOG);
temp_rotation_axis = data.moog_params(ROT_AZIMUTH ,:,MOOG);
temp_elevation = data.moog_params(ELEVATION,:,MOOG);
temp_stim_type = data.moog_params(STIM_TYPE,:,MOOG);
temp_amplitude = data.moog_params(AMPLITUDE,:,CAMERAS);
temp_motion_coherence = data.moog_params(COHERENCE,:,MOOG);
% temp_num_sigmas = data.moog_params(NUM_SIGMAS,:,MOOG);
temp_spike_rates = data.spike_rates(SpikeChan, :);% get the firing rates for all the trials
temp_stars_type = data.moog_params(STARS_TYPE,:,MOOG); % 更换刺激从三角形变成圆点，STARS_TYPE从1变成0
temp_eyepos_x = data.targ_params(1,:,MOOG);
temp_eyepos_y = data.targ_params(3,:,MOOG);

% some file miss ROT_POLARAZI and ROT_AMPLITUDE, but why?
try
    temp_rotazimuth = data.moog_params(ROT_POLARAZI ,:,MOOG);
    temp_rotamplitude = data.moog_params(ROT_AMPLITUDE,:,CAMERAS);
catch
    temp_rotazimuth = zeros(size(temp_azimuth));
    temp_rotamplitude = zeros(size(temp_azimuth));
end

null_value =  data.one_time_params(NULL_VALUE);

% transform to heading and rotation (0 = Expansion)
% 在极坐标中
%  90   0   -90  -180
% CCW  Exp  CW   Con (self-motion)
% Left Exp  Right
% include null trial (-9999)
temp_heading = aziToHeading(temp_azimuth,null_value);
if unique(temp_rotazimuth)==0                              % 1. old version，use rot axis to transform to rotation (20201102 修改前)
    temp_rotation = axisToRotation(temp_azimuth,temp_rotation_axis,null_value);
else                                                       % 2. new version, use rot azimuth to transform to rotation (20201102 修改后)
    temp_rotation = aziToHeading(temp_rotazimuth,null_value);
    % 此时radial motion只在heading part中，此时rotation=888，在rotation condition list 中补上radial motion
    temp_rotation(temp_heading==0)=0;
    temp_rotation(temp_heading==-180)=-180;
end

%get indices of any NULL conditions (for measuring spontaneous activity
trials = 1:length(temp_azimuth);

% If length(BegTrial) > 1 and all elements are positive, they are trials to be included.
% Else, if all elements are negative, they are trials to be excluded.
% This enable us to exclude certain trials ** DURING ** the recording more easily. HH20150410
select_trials = false(size(trials));
if length(BegTrial) == 1 && BegTrial > 0 % Backward compatibility
    select_trials(BegTrial:EndTrial) = true;
elseif all(BegTrial > 0) % To be included
    select_trials(BegTrial) = true;
elseif all(BegTrial < 0) % To be excluded
    temp_spike_rates(-BegTrial) = NaN; % changed by Lwh 201907 ： the trial we excluded only need to nan the spike rate
    %     select_trials(-BegTrial) = true;
    select_trials = ~ select_trials;
else
    disp('Trial selection error...');
    keyboard;
end

temp_heading = temp_heading(select_trials);
temp_rotation = temp_rotation(select_trials);
temp_elevation = temp_elevation(select_trials);
temp_stim_type = temp_stim_type(select_trials);
temp_amplitude = temp_amplitude(select_trials);
temp_rotamplitude = temp_rotamplitude(select_trials);
temp_motion_coherence = temp_motion_coherence(select_trials);
temp_stars_type = temp_stars_type(select_trials);

% for some mismatch smr and htp file (more htp trial/condition than smr)
if strcmp(FILE,'m7c1074r1.htb') || strcmp(FILE,'m7c1153r1.htb')  || strcmp(FILE,'m7c1022r1.htb') || strcmp(FILE,'m16c1009r1.htb')|| strcmp(FILE,'m16c1039r3.htb')
    temp_spike_rates = temp_spike_rates(1:end);
else
    temp_spike_rates = temp_spike_rates(select_trials);
end

%暂时求出unique判断condition数目和重复数
unique_heading = munique(temp_heading'); %included null
unique_heading(unique_heading==-9999)=[];% not included zero and null(-9999)
unique_heading(unique_heading==888)=[];% not included 888

unique_rotation = munique(temp_rotation'); %included null
unique_rotation(unique_rotation==-9999)=[];% not included zero and null(-9999)
unique_rotation(unique_rotation==888)=[];% not included 888

unique_elevation = munique(temp_elevation');
unique_stim_type = munique(temp_stim_type');
unique_stim_type(unique_stim_type==-9999)=[];% not included null(-9999)
unique_amplitude = munique(nonzeros(temp_amplitude')); % m/s
unique_rotamplitude = munique(nonzeros(temp_rotamplitude')); % degree/s
unique_coherence = munique(temp_motion_coherence');
unique_coherence(unique_coherence==-9999)=[];% not included null(-9999)

% for rescue file data
if sum(isnan(unique_elevation))>0
    unique_elevation = null_value;
end

% Now I change the code to cope with situations where different angles have different repetitions (unbalanced ANOVA).
% So this "repetitionN" is no longer important. HH20150704
if length(unique_heading) == 8 && length(unique_rotation)==8
    mean_repetitionN = floor( length(temp_spike_rates) / ((length(unique_heading)+length(unique_rotation)-2+1) * length(unique_elevation)* length(unique_stim_type)*length(unique_coherence)) ); % take minimum repetition
    temp_trials = mean_repetitionN * ((length(unique_heading)+length(unique_rotation)-2+1) * length(unique_elevation)* length(unique_stim_type)*length(unique_coherence)); %heading(+8), rotation(+8), repeat Exp+Con(-2), 1 null trial(+1)
else
    keyboard;
end

%根据完整的repetition重新选择trial范围
select_trials = true(size(1:temp_trials));

temp_heading = temp_heading(select_trials);
temp_rotation = temp_rotation(select_trials);
temp_elevation = temp_elevation(select_trials);
temp_stim_type = temp_stim_type(select_trials);
temp_amplitude = temp_amplitude(select_trials);
temp_rotamplitude = temp_rotamplitude(select_trials);
temp_motion_coherence = temp_motion_coherence(select_trials);
temp_spike_rates = temp_spike_rates(select_trials);
temp_stars_type = temp_stars_type(select_trials);
temp_eyepos_x = temp_eyepos_x(select_trials);
temp_eyepos_y = temp_eyepos_y(select_trials);

null_trials = logical((temp_heading == data.one_time_params(NULL_VALUE)) ); %spontaneous
% calculate spontaneous firing rate
for c=1:length(unique_coherence)
    for k=1:length(unique_stim_type)
        spon_found{c,k} = find(null_trials==1);
        resp_spon{c,k} = temp_spike_rates(spon_found{c,k});
        spon_resp{c,k} = nanmean(resp_spon{c,k});
        spon_resp_std{c,k} = nanstd(resp_spon{c,k});
    end
end

%去掉null trial
heading = temp_heading(~null_trials & select_trials);
rotation = temp_rotation(~null_trials & select_trials);
elevation = temp_elevation(~null_trials & select_trials);
stim_type = temp_stim_type(~null_trials & select_trials);
amplitude = temp_amplitude(~null_trials & select_trials);
rotamplitude = temp_rotamplitude(~null_trials & select_trials);
spike_rates = temp_spike_rates(~null_trials & select_trials);
motion_coherence = temp_motion_coherence(~null_trials & select_trials);

% **********************************************************************************************************
% remove time trend from spike_rate
% spike_rates = spike_rate_adjust(spike_rates_raw,methodCode)
% (空值或0：默认不作处理；1:detrend处理, 2多项式拟合，3 local regression)
spike_rates = spike_rate_adjust(spike_rates,0);
% **********************************************************************************************************

%unique parameters
unique_heading = munique(heading');
unique_heading(unique_heading==888)=[];% not included 888
unique_rotation = munique(rotation');
unique_rotation(unique_rotation==888)=[];% not included 888

unique_elevation = munique(elevation');
unique_stim_type = munique(stim_type');
unique_amplitude = munique(amplitude');
unique_rotamplitude = munique(rotamplitude');
unique_coherence = munique(motion_coherence');
stars_type = munique(temp_stars_type');
eyepos_x = munique(temp_eyepos_x');
eyepos_y = munique(temp_eyepos_x');

% separate heading trial and rotation trial
rad_trials = logical((heading==0 & rotation==0) | (heading==-180 & rotation==-180) );
ro_trials = logical(rotation~=888 & ~rad_trials);
ho_trials = logical(~ro_trials & ~rad_trials);

% overlap with exp and con
h_trials = logical(ho_trials | rad_trials);
r_trials = logical(ro_trials | rad_trials);

%% combine for loop analyze
condition{1} = heading;
condition{2} = rotation;

unique_cond{1} = unique_heading;
unique_cond{2} = unique_rotation;

trial_type{1,1} = h_trials;
trial_type{1,2} = r_trials;

trial_type{2,1} = ho_trials;
trial_type{2,2} = ro_trials;

all_trial = true(1,length(heading));

%% creat basic matrix represents for each response vector
resp = [];
for c=1:length(unique_coherence)
    for k=1:length(unique_stim_type)
        for s=1:2 % heading or rotation
            resp{c,s} = nan(length(unique_cond{s}),length(unique_stim_type));
            resp_std{c,s} = resp{c,s};
            resp_err{c,s} = resp{c,s};
            
            resp_trial{c,k,s} = nan(100,length(unique_cond{s})); % A large enough matrix (if you're not TOO CRAZY). HH20150704
            
            resp_L_exc{c,s} = [];
            resp_R_exc{c,s} = [];
            
            for i=1:length(unique_cond{s})
                select = logical(trial_type{1,s} & condition{s}==unique_cond{s}(i) & stim_type==unique_stim_type(k) & motion_coherence==unique_coherence(c));
                
                % Cope with situations where different angles have different numbers of repetition (unbalanced ANOVA).  HH20150704
                rep_this = sum(select);
                repetitions(c,k,s,i) = rep_this; % Save it
                
                if (rep_this > 0)  % there are situations where visual/combined has >2 coherence and vestibular only has one coherence
                    resp{c,s}(i,k) = nanmean(spike_rates(select));
                    resp_std{c,s}(i,k) = nanstd(spike_rates(select));
                    resp_err{c,s}(i,k) = nanstd(spike_rates(select)) / sqrt(rep_this);
                    resp_sse{c,s}(i,k) = sum((spike_rates(select)-resp{c,s}(i,k)).^2);
                    resp_trialnum{c,s}(i,k) = length(spike_rates(select));
                    
                    spike_temp = spike_rates(select);
                    try
                        resp_trial{c,k,s}(1:rep_this,i) = spike_temp;
                    catch
                        keyboard
                    end
                end
                
                % group by L and R; CW and CCW
                if unique_cond{s}(i)<0 && unique_cond{s}(i)~=-180 % L or CW
                    resp_L_exc{c,s} = [resp_L_exc{c,s} spike_temp];
                elseif unique_cond{s}(i)>0 % R or CCW
                    resp_R_exc{c,s} = [resp_R_exc{c,s} spike_temp];
                end
            end
            resp_trial{c,k,s}(all(isnan(resp_trial{c,k,s}),2),:) = []; % Trim the large matrix by removing unnecessary NaNs. HH20150704
        end
        
        % normalize FR: max to min - 0 to 1 (需要减去spon？)
        nor_temp = mapminmax([resp{c,1};resp{c,2}]',0,1);
        resp_norm{c,1} = nor_temp(1:size(resp{c,1},1))';
        resp_norm{c,2} = nor_temp(size(resp{c,1},1)+1:end)';
    end
end

non_nan{1} = ~isnan(resp{c,s}(:,k));
non_nan{2} = non_nan{1};
non_nan{2}([1 5]) = 0; % remove exp/con

% *** 此时resp第一个为-180，第二个为-135，。。。。以此类推，从小到大排列
unique_rep = unique(repetitions);
unique_rep = unique_rep(unique_rep>0);
if length(unique_rep)>1  % m7c816r1.htb
    unique_rep = unique_rep(1);
    disp('**************  repetition err ******************');
end

%% 判断是否相对spon有反应（高于spon）   分2种情况：包括和不包括exp/con
for c=1:length(unique_coherence)
    for k = 1: length(unique_stim_type)
        non_nan{1} = ~isnan(resp{c,s}(:,k));
        non_nan{2} = non_nan{1};
        non_nan{2}([1 5]) = 0;
        % r=1 :没有排除Exp，Con
        % r=2: 排除Exp,Con
        for r = 1:2
            for s=1:2
                non_nan_resp_trial =resp_trial{c,k,s}(:,non_nan{r});
                n = sum(non_nan{r});
                
                p_spon_stim_trial = [];
                for i = 1:n
                    p_spon_stim_trial(i) = ranksum(non_nan_resp_trial(:,i),resp_spon{c,k},'tail','right');
                end
                p_sti(c,k,r,s) = min(p_spon_stim_trial(:));
            end
        end
    end
end
p_sti_min = min([p_sti(1,1,1,1), p_sti(1,1,1,2)]); % h or r 中最小值，只要有一个显著，则判断为该细胞显著反应

%% Preference direction of self-motion
unique_elevation_temp = [];
for c=1:length(unique_coherence)
    for k = 1: length(unique_stim_type)
        for s=1:3  %1 for heading tuning, 2 for rotation tuning   3 for spiral space
            
            % VectorSum
            if s~=3
                unique_elevation_temp{s} = zeros(length(unique_cond{s}),1);
                
                if length(resp{c,s}) >= length(unique_cond{s})
                    % [az(c,k,s), el(c,k,s), amp(c,k,s)] = vectorsumAngle(resp{c,s}(:,k), unique_cond{s}, unique_elevation{s});
                    % vectorsumAngle 只能对极坐标生效，现在heading和rotation的0点为极坐标的90°，故调整以适用vectorsumAngle
                    azi{s} = HeadingToazi(unique_cond{s});
                    [az, el, amp,xx,yy,zz] = vectorsumSpiral(resp{c,s}(:,k), azi{s}, unique_elevation_temp{s});
                    [az_norm, el_norm, amp_norm] = vectorsumSpiral(resp_norm{c,s}(:,k), azi{s}, unique_elevation_temp{s}); % normalize FR
                    
                    % 返回为heading/rotation
                    pref_direc(c,k,s) = aziToHeading(az);
                    pref_el(c,k,s) = el;
                    pref_amp(c,k,s) = amp;
                    single_pref_x(s) = xx;
                    single_pref_y(s) = yy;
                    single_pref_z(s) = zz;
                    
                    pref_direc_norm(c,k,s) = aziToHeading(az_norm);
                    pref_el_norm(c,k,s) = el_norm;
                    pref_amp_norm(c,k,s) = amp_norm;
                else
                    pref_direc(c,k,s) = NaN;   % this hard-coded method cannot handle < 8 azims  -CRF
                end
                
                % Tuning Index  分2种情况：包括和不包括exp/con
                
                % D-prime / discriminability / sensitivity index / discrimination index
                % d' = (Rright - Rleft) / sqrt(1/2*(STDright^2+STDleft^2))
                
                % fine task之前使用的是左右22.5°计算，我的实验中没有这两个角度，故使用左右45°
                % coarse task使用左右90°
                d_cond(1) = 45;
                d_cond(2) = 90;
                
                Rspon = spon_resp{c,k};
                STDspon = spon_resp_std{c,k};
                
                for i = 1:2 % fine or coarse task
                    Rright = resp{c,s}(unique_cond{s} == d_cond(i));
                    Rleft = resp{c,s}(unique_cond{s} == -d_cond(i));
                    
                    STDright = resp_std{c,s}(unique_cond{s} == d_cond(i));
                    STDleft = resp_std{c,s}(unique_cond{s} == -d_cond(i));
                    
                    d_prime(c,k,s,i) = (Rright - Rleft) / (sqrt(1/2 * (STDright^2 + STDleft^2)));
                    
                    % 每个特定角度与平均FR求d’
                    if s == 1
                        d_prime_single{i}(1) = (Rright - Rspon) / (sqrt(1/2 * (STDright^2 + STDspon^2))); % R
                        d_prime_single{i}(2) = (Rleft - Rspon) / (sqrt(1/2 * (STDleft^2 + STDspon^2))); % L
                    else
                        d_prime_single{i}(3) = (Rright - Rspon) / (sqrt(1/2 * (STDright^2 + STDspon^2))); % CCW
                        d_prime_single{i}(4) = (Rleft - Rspon) / (sqrt(1/2 * (STDleft^2 + STDspon^2))); % CW
                    end
                    
                    % neural sensitivity = lg(|d'|)
                    neu_sensitivity(c,k,s,i) = log10(abs(d_prime(c,k,s,i)));
                end
                d_prime_fine(c,k,s) = d_prime(c,k,s,1);
                d_prime_coarse(c,k,s) = d_prime(c,k,s,2);
                
                neu_sensitivity_fine(c,k,s) = neu_sensitivity(c,k,s,1);  % translation中的左右d'，rotation中的CW，CCW d'
                neu_sensitivity_coarse(c,k,s) = neu_sensitivity(c,k,s,2);
                
                % 最大 d prime
                % 1. FR差值最大的两个点
                [maxFR(c,k,s),I_max] = max(resp{c,s});
                [~,I_min] = min(resp{c,s});
                d_prime_maxDiff(c,k,s) = ( resp{c,s}(I_max) - resp{c,s}(I_min) ) / sqrt(1/2 * (resp_std{c,s}(I_max)^2 + resp_std{c,s}(I_min)^2));
                
                % 2. FR差值最大的一条轴（相隔180度）
                % 计算所有(4条)对侧轴的差值（180~0前后，-135~45，-90~90，-45~-135）
                for ii = 1:4
                    all_axis(ii) = abs(resp{c,s}(ii) - resp{c,s}(ii+4));
                end
                [~,Index] = max(all_axis);
                I_max = [];I_min = [];
                if resp{c,s}(Index) > resp{c,s}(Index+4)
                    I_max = Index;
                    I_min = Index+4;
                else
                    I_max = Index+4;
                    I_min = Index;
                end
                d_prime_maxAxis(c,k,s) = ( resp{c,s}(I_max) - resp{c,s}(I_min) ) / sqrt(1/2 * (resp_std{c,s}(I_max)^2 + resp_std{c,s}(I_min)^2));
                
                % 3. 最大FR及其180度的对侧
                % d' = (Rmax - Rc-lat) / sqrt(1/2*(STDmax^2+STDc-lat^2))    c-lat:contralateral
                I_max = [];I_min = [];
                [~,I_max] = max(resp{c,s});
                I_min = I_max + 4;
                if I_min > 8
                    I_min = I_min-8;
                end
                d_prime_maxFR(c,k,s) = ( resp{c,s}(I_max) - resp{c,s}(I_min) ) / sqrt(1/2 * (resp_std{c,s}(I_max)^2 + resp_std{c,s}(I_min)^2));
                
                % DDI
                for r = 1:2 % 1.包括exp/con; 2.不包括
                    non_nan_resp = resp{c,s}(non_nan{r},k);
                    non_nan_resp_std = resp_std{c,s}(non_nan{r},k);
                    
                    % DDI: delta / (delta + 2 * sqrt(SSE/(N-M)))
                    %                     DDI(c,k,r,s) = ( max(non_nan_resp)-min(non_nan_resp) ) / ( max(non_nan_resp)-min(non_nan_resp)+ 2 * sqrt(sum(resp_sse{c,s}) / (sum(resp_trialnum{c,s}) - sum(non_nan{r}))) ); % 和下面一样的
                    DDI(c,k,r,s) = ( max(non_nan_resp)-min(non_nan_resp) ) / ( max(non_nan_resp)-min(non_nan_resp)+ 2 * sqrt( sum(non_nan_resp_std.^2 / sum(non_nan{r})))) ;   % Lwh 202011
                end
                
                % HTI = |ΣRi| / Σ|Ri|   Ri=mean-spon
                %     =  AMPvectorsum / sum(|mean-spon|)
                HTI(c,k,s) = pref_amp(c,k,s)/sum(abs(resp{c,s}(:,k)-spon_resp{c,k}));
                
                % modulation index, DeAnglis and UKA, 2002, add by Lwh202101
                % direction/disparity tuning index, DeAngelis1, Newsome, 2004. tuning strength
                modulation_index(c,k,s) = (max(resp{c,s}) - min(resp{c,s})) / (max(resp{c,s}) - Rspon); % modulation index = (Rmax-Rmin)/(Rmax-Spon)
                
                % 0度（向前）斜率
                all_slope(s,:) = gradient(resp{c,s}(:,k)) ./ gradient(unique_cond{s});
                ahead_slope(s) = all_slope(s,unique_cond{s}==0); % 正号表示pref right，负号表示pref left
                
                % anova test
                for r = 1:2
                    p_value(c,k,r,s) = anova1(resp_trial{c,k,s}(:,non_nan{r}),'','off');
                end
                
                % pair discrimination index ratio:比较与0°对称角度下的反应，正负45°，正负90°，正负135°( not consider forward and backward)
                % change from YunYang, Sheng Liu, Dora, 2011 :DSDI = 1/n * Σ( (Rfar(i)-Rnear(i)) / |Rfar(i)-Rnear(i)| + std)
                % Lwh 20210730
                pair_unique_cond = unique_cond{s}(unique_cond{s}~=0 & unique_cond{s}~=-180 & unique_cond{s}~=180);
                for i = 1:length(pair_unique_cond)/2
                    pair_1 = spike_rates(condition{s} == pair_unique_cond(i));
                    pair_2 = spike_rates(condition{s} == -pair_unique_cond(i));
                    DI_temp(i) = abs(mean(pair_1) - mean(pair_2)) / (abs(mean(pair_1) - mean(pair_2)) + std([pair_1 pair_2])); % from 0 to 1 分子取绝对值，不考虑左还是右
                end
                pair_DI(s) = sum(DI_temp) / (length(pair_unique_cond)/2); % from 0 to 1 :1: strong preference; 0: no preference
                
                
                % non-pair discrimination index, 参考DDI，不考虑exp、con
                % DDI: delta / (delta + 2 * sqrt(SSE/(N-M)))
                % group all left heading and right heading / CW and CCW rotation
                meanA = nanmean(resp_L_exc{c,s});
                meanB = nanmean(resp_R_exc{c,s});
                stdA = nanstd(resp_L_exc{c,s});
                stdB = nanstd(resp_R_exc{c,s});
                sseA = sum((resp_L_exc{c,s} - meanA).^2);
                sseB = sum((resp_R_exc{c,s} - meanB).^2);
                nonpair_DI(s) = abs(meanA-meanB) / (abs(meanA-meanB)+2*sqrt((sseA+sseB)/(length(resp_L_exc{c,s})+length(resp_R_exc{c,s}) - 2)) ); % DDI: delta / (delta + 2 * sqrt(SSE/(N-M)))
                
                % variability调整不一样
                nonpair_DI2(s) = abs(meanA-meanB) / (abs(meanA-meanB)+ std([resp_L_exc{c,s} resp_R_exc{c,s}])); % DSDI = 1/n Σ(delta/(|delta|+σavg)), just like pair discrimination index
                
            elseif s==3 % spiral space: 将rotation plane 立起来，rotation90在最上面/北极点（dot:CW; self-motion:CCW）
                
                % rotation plane
                % 立起来后的顺序：-180(后) -90（下） 0（前） 90（上） 180（后）
                %                 spi_elevation = unique_rotation;
                %                 spi_elevation(spi_elevation>90) = 180 - spi_elevation(spi_elevation>90);
                %                 spi_elevation(spi_elevation<-90) = -180 -
                %                 spi_elevation(spi_elevation<-90);、
                spi_elevation = [0 -45 -90 -45 0 45 90 45]';
                spi_azimuth = [270 270 888 90 90 90 888 270]'; % 888: 在z轴上，azi没有意义，随便给个值
                
                % 去掉Exp和Con两个值，heading plane上做vectorsum的时候已经用了
                spi_azimuth_temp = spi_azimuth;
                spi_azimuth_temp([1 5]) = [];
                
                spi_elevation_temp = spi_elevation;
                spi_elevation_temp([1 5]) = [];
                
                spi_resp_temp = resp{c,2}(:,k);
                spi_resp_temp([1 5]) = [];
                
                % 加上heading plane
                spi_azimuth_all = [spi_azimuth_temp;azi{1}]; % azi{1}:con-left-exp-right
                spi_elevation_all = [spi_elevation_temp;unique_elevation_temp{1}]; % 与heading plane的夹角
                spi_resp_all = [spi_resp_temp;resp{c,1}(:,k)];
                maxFR(c,k,s) = max(spi_resp_all);
                
                [az, el, amp,x_pref,y_pref,z_pref] = vectorsumSpiral(spi_resp_all, spi_azimuth_all, spi_elevation_all);
                pref_direc(c,k,s) = aziToHeading(az);
                pref_el(c,k,s) = el;
                pref_amp(c,k,s) = amp;
                
                % normalize FR
                spi_resp_all_norm = [resp_norm{c,2}([2:4,6:8],k);resp_norm{c,1}(:,k)];
                [az_norm, el_norm, amp_norm] = vectorsumSpiral(spi_resp_all_norm, spi_azimuth_all, spi_elevation_all);
                pref_direc_norm(c,k,s) = aziToHeading(az_norm);
                pref_el_norm(c,k,s) = el_norm;
                group_result_norm(c,k,s) = amp_norm;
                
                % 此时el为总vectro-sum与translation平面夹角，与rotation平面夹角未知！！！二者相加不是90°!!!  Lwh 20201116
                vector_sum = [x_pref,y_pref,z_pref];
                trans_normal_vector = [0,0,1]; % translation plane 法向量
                spiral2Tran = cosd(el); % el:-90~90, cosd(el):0~1  总vectro-sum与translation平面夹角的cos值
                
                rot_normal_vector = [-1,0,0]; % rotation plane 法向量
                spiral2Rot_norm = dot(vector_sum,rot_normal_vector)/(norm(vector_sum)*norm(rot_normal_vector)); % cosA = dot(a,b)/(norm(a)*norm(b)) % 总vectro-sum与rotation平面法向量夹角的cos
                if spiral2Rot_norm<0
                    rot_normal_vector = [1,0,0]; % rotation plane 法向量
                    spiral2Rot_norm = dot(vector_sum,rot_normal_vector)/(norm(vector_sum)*norm(rot_normal_vector)); % cosA = dot(a,b)/(norm(a)*norm(b)) % 总vectro-sum与rotation平面法向量夹角的cos = 与rotation平面夹角的sin
                end
                spiral2Rot = sqrt(1-spiral2Rot_norm^2); % sinA^2+cosA^2 = 1，总vectro-sum与rotation平面夹角的cos值
                
                % spiral_index does not consider neural response variability, try to use model fitting to find preferred spiral direction
                % 1: bias to translation; -1: bias to rotation
                spiral_index = (spiral2Tran-spiral2Rot)/(spiral2Rot+spiral2Tran); % spiral2Rot、spiral2Tran (cosA)越接近1说明A越小，合向量越靠近（偏好）该平面；越接近0说明A越大，合向量越原理该平面
                spiral_index_angle = (acosd(spiral2Rot) - acosd(spiral2Tran)) / (acosd(spiral2Tran) + acosd(spiral2Rot));
                
                % commented by Lwh 20201016
                % Lwh 20191217, 压缩Exp-Con的轴 (投影到 Left-Right-CW-CCW plane, x-z plane)
                % 当pref
                % Exp/Con时，y_pref(纵深)值较大，x_pref(横轴)和z_pref(纵轴)值较小，此时虽然x_pref与z_pref绝对值较小但是比值就能影响theta,弃用
                
                % ANOVA for spiral space
                spi_rot_trial = [];
                condtion_num = length(unique_heading) + length(unique_rotation) - 2;
                
                spi_rot_trial = resp_trial{c,k,2}(:,[2:4,6:8]);
                resp_trial_all = [resp_trial{c,k,1}(:);spi_rot_trial(:)];
                resp_trial_all = reshape(resp_trial_all,[],condtion_num);
                p_value(c,k,1,s) = anova1(resp_trial_all,'','off');  % include exp/con (所有角度，包括heading与rotation)
                
                resp_trial_all_noexpcon = resp_trial_all(:,[2:4,6:end]);
                p_value(c,k,2,s) = anova1(resp_trial_all_noexpcon,'','off');  % exclude exp/con
                
                %  Motion type的DDI：计算区分最大最小的反应能力貌似没有意义？ 应该是区分两种motion
                %  type的能力，故直接用DDI_index/d prime index 去衡量，所以注释掉有关mt DDI / d'的计算
            end
        end
        
        % motion type tuning
        % heading与rotation plane有显著区别
        p_value_mt(c,k) = anova1([resp_trial{c,k,1}(:),resp_trial{c,k,2}(:)],'','off');
        
    end
end


% %% circular analysis
% for c=1:length(unique_coherence)
%     for k=1:length(unique_stim_type)
%         for s=1:2 % heading or rotation
%             rep_this = unique(repetitions(c,k,s,:));
%             non_nan_resp_trial = resp_trial{c,k,s}(:,non_nan{1});
%             alpha = deg2rad(repmat(unique_cond{s}(non_nan{1})',rep_this,1));
%
%             mu = circ_mean(alpha(:), non_nan_resp_trial(:));  % angular mean
%
%             pref_direc_cir(c,k,s) = rad2deg(mu);  % 结果和vector-sum一样的
%         end
%     end
% end


%% Modulation Index for motion type stimulation: >0: pref Heading;  <0: pref Rotation
% 1. 0度时候的斜率，取绝对值，只比较大小
ahead_slope_index = (abs(ahead_slope(1)) - abs(ahead_slope(2))) / (abs(ahead_slope(1)) + abs(ahead_slope(2)));

% 2. d'
d_prime_index_fine = (abs(d_prime_fine(1,1,1)) - abs(d_prime_fine(1,1,2)) ) ./ (abs(d_prime_fine(1,1,1)) + abs(d_prime_fine(1,1,2)) ) ;
d_prime_index_coarse = (abs(d_prime_coarse(1,1,1)) - abs(d_prime_coarse(1,1,2)) ) ./ (abs(d_prime_coarse(1,1,1)) + abs(d_prime_coarse(1,1,2)) ) ;
d_prime_index(1) = d_prime_index_fine;
d_prime_index(2) = d_prime_index_coarse;

d_prime_maxDiff_index = (abs(d_prime_maxDiff(1,1,1)) - abs(d_prime_maxDiff(1,1,2))) ./ (abs(d_prime_maxDiff(1,1,1)) + abs(d_prime_maxDiff(1,1,2)));
d_prime_maxAxis_index = (abs(d_prime_maxAxis(1,1,1)) - abs(d_prime_maxAxis(1,1,2))) ./ (abs(d_prime_maxAxis(1,1,1)) + abs(d_prime_maxAxis(1,1,2)));
d_prime_maxFR_index = (abs(d_prime_maxFR(1,1,1)) - abs(d_prime_maxFR(1,1,2))) ./ (abs(d_prime_maxFR(1,1,1)) + abs(d_prime_maxFR(1,1,2)));

d_prime_single_fine = d_prime_single{1};
d_prime_single_coarse = d_prime_single{2};

% 3. DDI
DDI_inc_index = (DDI(1,1,1,1) - DDI(1,1,1,2)) ./ (DDI(1,1,1,1) + DDI(1,1,1,2)) ;
DDI_exc_index = (DDI(1,1,2,1) - DDI(1,1,2,2)) ./ (DDI(1,1,2,1) + DDI(1,1,2,2)) ;

% 4. HTI
HTI_index = (HTI(1,1,1) - HTI(1,1,2)) ./ (HTI(1,1,1) + HTI(1,1,2)) ;

% 5. Neural sensitivity = log10(abs(d'))
neu_sen_fine = squeeze(neu_sensitivity_fine)';
neu_sen_coarse = squeeze(neu_sensitivity_coarse)';

neu_sen_index_fine = (neu_sen_fine(1) - neu_sen_fine(2))/(neu_sen_fine(1) + neu_sen_fine(2));
neu_sen_index_coarse = (neu_sen_coarse(1) - neu_sen_coarse(2))/(neu_sen_coarse(1) + neu_sen_coarse(2));
neu_sen_index(1) = neu_sen_index_fine; % >0: pref H; <0: pref R
neu_sen_index(2) = neu_sen_index_coarse;

% 6. 
pair_TRDI = (pair_DI(1)-pair_DI(2))/(pair_DI(1)+pair_DI(2)); % from -1 to 1, -1:Rot, 1:Translation
nonpair_TRDI = (nonpair_DI(1)-nonpair_DI(2))/(nonpair_DI(1)+nonpair_DI(2)); % from -1 to 1, -1:Rot, 1:Translation
nonpair_TRDI2 = (nonpair_DI2(1)-nonpair_DI2(2))/(nonpair_DI2(1)+nonpair_DI2(2)); % from -1 to 1, -1:Rot, 1:Translation

% 7. new index for T-R index, 20210826
% 7.1. LR firing rate different / CWCCW different
LR_diff = abs(resp{1,1}(unique_heading==90) - resp{1,1}(unique_heading==-90));
CWCCW_diff = abs(resp{1,2}(unique_rotation'==90) - resp{1,2}(unique_rotation'==-90));
fr_axis_ratio = LR_diff / CWCCW_diff; % >1: more T; <1: more R

% 7.2. vector-sum projection on translation axis and rotatino axis
axis_proj = squeeze(pref_amp(c,k,:)) .* cosd(HeadingToazi(squeeze(pref_direc(c,k,:))));
vec_axis_ratio = abs(axis_proj(1) / axis_proj(2)); % >1: more T; <1: more R
                

%% separable predictions
% y = spike_rates;
x{1} = nan(size(heading));
x{2} = nan(size(rotation));
% remove 0 and -180 （forward and backward）
x{1}(heading<0 & ho_trials) = -1; % left
x{1}(heading>0 & ho_trials) = 1; % right
x{2}(rotation<0 & ro_trials) = -1; % CW
x{2}(rotation>0 & ro_trials) = 1; % CCW

angle{1} = nan(size(heading));
angle{2} = nan(size(rotation));
angle{1}(ho_trials) = heading(ho_trials);
angle{2}(ro_trials) = rotation(ro_trials);

% % linear regression = L/R + abs(heading) + L/R*asb(heading)
% ,比较第一项的系数，不同次regress之间不饿能比较！！！
% for hr = 1:2
% %     [b{hr},bint{hr},r,rint,stats] = regress(y', [x{hr}' abs(angle{hr})' (x{hr}.*abs(angle{hr}))']);
%     [b,dev,stats] = glmfit([x{hr}' abs(angle{hr})' (x{hr}.*abs(angle{hr}))'],y');
% end


% SVM to heading(L/R) or rotation(CW/CCW), 一维。。
for hr = 1:2
    X = [resp_L_exc{1,hr},resp_R_exc{1,hr}]';
    if hr == 1
        y = [repmat("Left",length(resp_L_exc{1,hr}),1); repmat("Right",length(resp_R_exc{1,hr}),1)];
    else
        y = [repmat("CW",length(resp_L_exc{1,hr}),1); repmat("CCW",length(resp_R_exc{1,hr}),1)];
    end
    
    % 训练SVM model
    SVMModel{hr} = fitcsvm(X,y,'standardize',true); % 标准化原始数据
    
    %     %  训练情况
    %     sv = SVMModel{hr}.SupportVectors;
    %     figure
    %     gscatter(zeros(size(X)),X,y)
    %     hold on
    %     plot(zeros(size(sv)),sv,'ko','MarkerSize',10)
    %     legend('versicolor','virginica','Support Vector')
    %     hold off
    
    % 验证SVM model
    CVMdl{hr} = crossval(SVMModel{hr},'kfold',10); % 对SVM model进行交叉验证:将原始数据分割成k组，每个单独的数据集作为验证集，其余的k-1个数据集用来训练，交叉验证重复k次
    classLoss(hr)  = kfoldLoss(CVMdl{hr}); % 估计out-of-sample misclassification rate
    SVM_correct_rate(hr) = 1 - classLoss(hr); % preformance
end
% SVM对heaidng/rotation的分类能力作为衡量标准
SVM_index = (SVM_correct_rate(1)-SVM_correct_rate(2))/(SVM_correct_rate(1)+SVM_correct_rate(2)); % 1: only heaidng correct, -1: only rotation correct; 0: correct rate equal


% ROC
for hr = 1:2
    AuROC(hr) =  abs(rocN(resp_L_exc{1,hr},resp_R_exc{1,hr})-0.5)+0.5;
end
AUC_index = (AuROC(1)-AuROC(2))/(AuROC(1)+AuROC(2)); % 1:heaidng seperate, -1: rotation seperate; 0: separability equal


% spiral space in 2-D plane 二维插值yyy
for s = 1:2
    unique_cond_plot{s} = [unique_cond{s};180];
    % make it circular
    circular_add{s} = [1:length(unique_cond{s}) 1];
end
xx = -180:15:180;
yy = -180:15:180;
[xq,yq] = meshgrid(xx,yy);

xx_ori = [-180:45:180,zeros(size(-180:45:180))];
yy_ori = [zeros(size(-180:45:180)),-180:45:180];
zz_ori = [resp{1,1}(circular_add{1})', resp{1,2}(circular_add{2})'];
% remove the repeated one zero for translation or rotation plane
xx_ori(5) = [];
yy_ori(5) = [];
zz_ori(5) = [];

vq1 = griddata(xx_ori,yy_ori,zz_ori,xq,yq,'v4');
vq2 = griddata(xx_ori,yy_ori,zz_ori,xq,yq,'natural');

% save the result for 2-D v4 interpolate
x_v4interp = xq;
y_v4interp = yq;
z_v4interp = vq1;

if IF_FIGURE
    figure(12);clf; % 2D
    subplot(1,2,1)
    plot3(xx_ori,yy_ori,zz_ori,'o');
    hold on
    mesh(xq,yq,vq1);
    title('v4');
    xlabel('Translation');
    ylabel('Rotation');
    
    subplot(1,2,2)
    plot3(xx_ori,yy_ori,zz_ori,'o');
    hold on
    mesh(xq,yq,vq2);
    title('natural');
    xlabel('Translation');
    ylabel('Rotation');
    
%     figure(13);clf; % back to 3d in ball
%     az1 = deg2rad(HeadingToazi(xq));
%     el1 = yq;
%     for i = 1:size(el1,1)
%         for j = 1:size(el1,2)
%             if el1(i,j) > 90
%                 el1(i,j) = 90 - (el1(i,j)-90);
%             elseif el1(i,j) < -90
%                 el1(i,j) = -90-(el1(i,j)+90);
%             end
%         end
%     end
%     el1 = deg2rad(el1);
%     amp = vq2;
%     [xint,yint,zint] = sph2cart(az1,el1,amp);
%     plot3(xint,yint,zint,'ko');
%     hold on
%     surf(xint,yint,zint)
    
end

% SVD, singular value decomposition: separable or inseparable
% 因为矩阵vq1是插值的来的，以下奇异值分解结果不太准确。。。。
[U,S,V] = svd(vq1); % S:奇异值singular values矩阵（对角线），U：AA^T的所有特征向量组成的矩阵-左奇异向量singular vector，V：A^TA的所有特征向量组成的矩阵-右奇异向量
singular_values = diag(S);
first_singular_value = singular_values(1);
first_left_singular_vector = U(:,1);
first_right_singular_vector = V(:,1);

% M = k + σU(i)V(j)^T


% separability index, DSI(direction separability index in Ryo, Dora, Gregory, 2017)
nor_singular_value = normc(singular_values);
SI = 1 - (nor_singular_value(1))^2 / sum(nor_singular_value.^2); % 0: separable; greater: more inseparable
% population 比较不同cell type（h,r,spiral cell）的SI是否有区别
% 猜想：h,r cell中SI较小，spiral cell中SI较大。 显著性分析: ANCOVA


%% -------------------------------------------------------------------
%check significance of DDI and calculate p value, do bootstrap at
%the same time to test value varience (permutation test)
perm_num = 1000;
bin = 0.005;
spike_rates_perm = [];
for n=1: perm_num
    % this is permuted based on trials
    
    % 计算heading和rotation DDI p value：
    %******* 不能将heading tuning task 和 rotation tuning task的spike一起混合随机抽，应该各自混合随机抽！！！！！！*********
    %             spike_rates_perm = spike_rates;
    %             spike_rates_perm = spike_rates_perm( randperm(length(spike_rates_perm)) );
    for c=1:length(unique_coherence)
        for k = 1: length(unique_stim_type)
            for r = 1:2
                for s = 1:3
                    if s~=3
                        temp = spike_rates(trial_type{r,s});
                    else
                        temp = spike_rates(all_trial);
                    end
                    spike_rates_perm{r,s} = temp(randperm(length(temp)));
                end
            end
        end
    end
    
    % re-creat a matrix similar as resp_mat  (由于spike打乱是根据两个tuning task来的，与之前的构造resp不一样，之前是分3种情况：heading，rotation，spiral，再合并，现在是直接两种情况heading_plane，rotation_plane)
    resp_perm_temp = [];
    resp_std_perm_temp = [];
    
    %1 for heading tuning, 2 for rotation tuning
    for c=1:length(unique_coherence)
        for k = 1: length(unique_stim_type)
            for r = 1:2
                for s = 1:2
                    perm_temp = [];
                    trial_perm = [];
                    
                    perm_temp = logical(trial_type{r,s}==1);
                    trial_perm = trial_type{r,s}(perm_temp);
                    stim_type_perm = stim_type(perm_temp);
                    motion_coherence_perm = motion_coherence(perm_temp);
                    cond_perm = condition{s}(perm_temp);
                    unique_cond_perm = unique(cond_perm);
                    
                    for i=1:length(unique_cond_perm)
                        select = logical(trial_perm & cond_perm==unique_cond_perm(i) & stim_type_perm==unique_stim_type(k) & motion_coherence_perm==unique_coherence(c));
                        if (sum(select) > 0)
                            resp_perm{c,r,s}(i,k) = nanmean(spike_rates_perm{r,s}(select));
                            resp_std_perm{c,r,s}(i,k) = nanstd(spike_rates_perm{r,s}(select));
                        else
                            resp_perm{c,r,s}(i,k) = 0;
                            resp_std_perm{c,r,s}(i,k) = 0;
                        end
                    end
                end
            end
            
            % for spiral tuning s==3
            resp_perm{c,1,3} = [resp_perm{c,1,1}(:,k);resp_perm{c,1,2}([2:4,6:8],k)];
            resp_std_perm{c,1,3} = [resp_std_perm{c,1,1}(:,k);resp_std_perm{c,1,2}([2:4,6:8],k)];
            resp_perm{c,2,3} = NaN;   % spiral space tuning 不需要区分r1,r2
            resp_std_perm{c,2,3} = NaN;
        end
    end
    
    % re-calculate DDI now
    for c=1:length(unique_coherence)
        for k = 1: length(unique_stim_type)
            for r = 1:2
                for s = 1:3 %1 for heading tuning, 2 for rotation tuning
                    if r==2 && s==3
                        DDI_perm(c,k,r,s,n) = NaN;
                    else
                        %Modulation Index
                        non_nan_perm = ~isnan(resp_perm{c,r,s}(:,k));
                        non_nan_resp_perm = resp_perm{c,r,s}(non_nan_perm,k);
                        non_nan_resp_std_perm = resp_std_perm{c,r,s}(non_nan_perm,k);
                        %                         DDI_perm(c,k,r,s,n) = ( max(non_nan_resp_perm)-min(non_nan_resp_perm) ) / ( max(non_nan_resp_perm)-min(non_nan_resp_perm)+ ...
                        %                             2 * sqrt( sum(non_nan_resp_std_perm.^2 / (unique_rep-1) / sum(non_nan_perm)))) ;    % delta / (delta + 2 * sqrt(SSE/(N-M))) = delta / (delta + 2 * sqrt(SSE/(rep-1)/M)))
                        DDI_perm(c,k,r,s,n) = ( max(non_nan_resp_perm)-min(non_nan_resp_perm) ) / ( max(non_nan_resp_perm)-min(non_nan_resp_perm)+ ...
                            2 * sqrt( sum(non_nan_resp_std_perm.^2 / sum(non_nan_perm)))) ;    % delta / (delta + 2 * sqrt(SSE/(N-M))) = delta / (delta + 2 * sqrt(SSE/M))) Lwh 202011
                    end
                end
            end
        end
    end
end

% now calculate p value or significant test
for c=1:length(unique_coherence)
    for k = 1: length(unique_stim_type)
        for r = 1:2
            for s=1:2
                p_DDI(c,k,r,s) = length(find(DDI_perm(c,k,r,s,:)>=DDI(c,k,r,s)) )/perm_num;
            end
        end
    end
end

%% Xu Hong preference direction analyses    added in 20190422
% vector analysis, only use the rotation plane

% Tuning index (TI) = FRmax - FRmin / FRmax + FR min
% Decompose the neural respons onto the Vradial and Vrotation
% Vradial = FR x cos(alpha)   Vrotation = FR x cos(beta)
% alpha: angle between the spiral motion and the expansion motion
% beta:  angle between the spiral motion and the Vrotation

Vrad_temp = resp{1,2} .* cos(abs(unique_rotation));
Vrot_temp = resp{1,2} .* cos(abs(abs(unique_rotation)-90));
Vrad = sum(Vrad_temp);
Vrot = sum(Vrot_temp);


%% Mutual information between neural responses and heading stimulus or rotation stimulus
% 看特定的unit对两种刺激的编码程度，MI越大，说明神经元的反应与刺激条件关系越密切 （FR与角度关系越密切）
% Strong et al. 1998
for c = 1:length(unique_coherence)
    for k = 1: length(unique_stim_type)
        for s = 1:2
            muinf_fr{c,k,s} = resp_trial{c,k,s}(:);
            muinf_cond{c,k,s} = repmat([-180:45:135],mean_repetitionN,1);
            muinf_cond{c,k,s} = muinf_cond{c,k,s}(:);
            
            [MI(c,k,s),normalMI(c,k,s)] = nmi(muinf_cond{c,k,s},muinf_fr{c,k,s});
        end
    end
end


%% Fitting
% for s = 1:2
%     unique_cond_plot{s} = [unique_cond{s};180];
%     % make it circular
%     circular_add{s} = [1:length(unique_cond{s}) 1];
% end

if IF_Fit
    for k=1: length(unique_stim_type)
        for c=1:length(unique_coherence)
            for s = 1:2
                if sum(abs(unique_cond_plot{s}) == 180) > 0
                    fo_ = fitoptions('method','NonlinearLeastSquares','Robust','On','Algorithm','T',...
                        'StartPoint',[10 10 0 1.4284],'Lower',[0 0 -pi 0],'Upper',[Inf 500 pi 500],'MaxFunEvals',5000);
                    ft_ = fittype('a*exp(-2*(1-cos(x-pref))/sigma^2)+b',...
                        'dependent',{'y'},'independent',{'x'},...
                        'coefficients',{'a', 'b', 'pref', 'sigma'});
                    
                    not_nan = ~isnan(resp{c,s}(:,k));
                    [cf_{s}, goodness, output]= fit(unique_cond_plot{s}(not_nan)*pi/180,resp{c,s}(not_nan,k),ft_,fo_);
                    
                    GaussfitRsquare(s) = goodness.rsquare;
                    
                    % if output.exitflag > 0 && cf_.sigma >= 0.3 changed byLwh 20181107
                    % 增加约束goodness.rsquare>0 Lwh 202108
                    if output.exitflag > 0 && cf_{s}.sigma >= 0.3 && cf_{s}.sigma < 500 && goodness.rsquare>0% output.exitflag > 0 收敛
                        fitOK(s) = 1;
                    else
                        fitOK(s) = 0;
                    end
                end
                
                if fitOK(s)
                    prefGauss_temp = [];
                    prefGauss_temp = cf_{s}.pref/pi*180;
                    if prefGauss_temp<-180
                        prefGauss(c,k,s) = 360 + prefGauss_temp;
                    elseif prefGauss_temp>180
                        prefGauss(c,k,s) = prefGauss_temp-360;
                    else
                        prefGauss(c,k,s) = prefGauss_temp;
                    end
                    
                    maxFRGauss(c,k,s) = feval(cf_{s},cf_{s}.pref);
                    
                    %         widGauss(c,k) = 2*acos(1-cf_.sigma^2*log(2)/2)/pi*180; % HH20140425
                    widGauss(c,k,s) = abs(2*acos(abs(1-cf_{s}.sigma^2*log(2)/2))/pi*180);
                else
                    prefGauss(c,k,s) = nan;
                    widGauss(c,k,s) = nan; % HH20140425
                    maxFRGauss(c,k,s) = nan;
                end
            end
        end
    end
    % else
    %     fitOK = nan;
    %     widGauss = nan;
    %     prefGauss = nan;
    %     maxFRGauss = nan;
    %     GaussfitRsquare = nan;
end

%% cell type 20200723
% 1: translation cell; 2: rotation cell; 3: spiral cell;
% 4: no tuning but above spon; 5: no active response
p_cri = 0.001;
select = true;
data.p_stim = min([p_sti(1,1,1,1), p_sti(1,1,1,2)]); % h or r 中最小值，只要有一个显著，则判断为该细胞显著反应
data.p_stim_h = p_sti(1,1,2,1); % 排除exp/con
data.p_stim_r = p_sti(1,1,2,2); % 排除exp/con
data.p_turn_h = p_value(1,1,2,1); % 排除exp/con
data.p_turn_r = p_value(1,1,2,2); % 排除exp/con
data.p_turn_spi = p_value(1,1,1,3); % 包括exp/con
data.spiral_index = spiral_index;

data.p_turn_mt = p_value_mt;

[cell_type,mt_type] = cell_type_classification(data,select,p_cri);


%% plot figures
if IF_FIGURE
    % temporarily hard coded, will be probematic if there are more than 3*3 conditions
    % repetitions -GY
    f{1,1}='ro-'; f{1,2}='bo-'; f{1,3}='mo-';
    f{2,1}='ro-'; f{2,2}='bo-'; f{2,3}='mo-';
    f{3,1}='ro-'; f{3,2}='bo-'; f{3,3}='mo-';
    figure(10);
    clf;
    set(gcf,'color','white');
    set(10,'Position', [30,40 1300,950], 'Name', 'Spiral Tuning');
    
    title(sprintf('%s, %s %g , %s %g ,%s %g\n',FILE,'Unit',SpikeChan,'Rep',mean_repetitionN,'Stim type',unique_stim_type(1)),'fontsize',15);
    axis off;
    
    %% errobar plot
    ylimt = max([resp{1,1};resp{1,2}])+max([resp_err{1,1};resp_err{1,2}]);
    
    for k=1: length(unique_stim_type)
        for c=1:length(unique_coherence)
            for s = 1:2
                if s==1
                    axes('position',[0.05 0.64 0.26 0.27]);
                    title('Heading Plane','fontsize',12);
                elseif s==2
                    axes('position',[0.05 0.34 0.26 0.27]);
                    title('Rotation Plane','fontsize',12);
                    %             elseif s==4
                    %                 axes('position',[0.05 0.04 0.26 0.27]);
                end
                
                non_nan = [];
                non_nan = ~isnan(resp{c,s}(circular_add{s},k));
                set(errorbar(unique_cond_plot{s}(non_nan), resp{c,s}(circular_add{s}(non_nan),k), resp_err{c,s}(circular_add{s}(non_nan),k), f{c,s} ),'linewidth',2);
                
                xlim( [min(unique_cond_plot{s}), max(unique_cond_plot{s})] );
                ylim([0 ylimt]);
                set(gca,'xtick',[-180.0000 -135.0000  -90.0000  -45.0000  0  45.0000  90.0000  135.0000  180.0000]);
                if s==2
                    set(gca,'xticklabel',{'Con','','CW','','Exp','','CCW','','Con'});
                end
                
                % Fitting
                hold on;
                y_limits = ylim;
                
                plot([ min(unique_cond_plot{s})  max(unique_cond_plot{s})] ,[spon_resp{c,k} spon_resp{c,k}], 'k--');
                
                % Forward and backward
                
                axis manual;
                plot([0 0], y_limits,'k--');
                
                %Preferred
                line([pref_direc(c,k,s) pref_direc(c,k,s)],y_limits,'Color',f{1,s}(1),'linestyle','-','LineWidth',2); % VectorSum %Lwh20170106 change the color to red
                
                if IF_Fit
                    if sum(abs(unique_cond_plot{s}) == 180) > 0
                        if fitOK(s)
                            xx = -180:180;
                            yy = feval(cf_{s},[-180:180]/180*pi);
                            plot(xx,yy,[f{1,3}(1) '--'],'linewidth',2); %Lwh20170207
                        end
                    end
                    text(-170,y_limits(2)-y_limits(2)/20,sprintf('Rsquare = %0.3g', GaussfitRsquare(s)),'color',[f{1,3}(1)]);
                    
                    % Transferred from azimuth to heading. HH20140415
                    line([prefGauss(c,k,s) prefGauss(c,k,s)],y_limits,'color',f{c,3}(1),'linestyle','--','linewidth',2);   % Gaussian fitting %Lwh20170106 change the color by add "unique_stim_type(k)" with "+1"
                end
                
                if k == 1 ;ylabel('spikes/s'); end
            end
        end
    end
    
    %% polar coordinates plot
    for k=1: length(unique_stim_type)
        for c=1:length(unique_coherence)
            for s = 1:2
                max_scale_temp{k,c,s} = max([resp{c,s}(circular_add{s},k)] + [resp_err{c,s}(circular_add{s},k)]);
            end
        end
    end
    
    max_scale = max(cell2mat(squeeze(max_scale_temp)));
    
    for k=1: length(unique_stim_type)
        for c=1:length(unique_coherence)
            for s = 1:2
                if s==1
                    axes('position',[0.33 0.64 0.27 0.27]);
                elseif s==2
                    axes('position',[0.33 0.34 0.27 0.27]);
                elseif s==4
                    axes('position',[0.33 0.04 0.27 0.27]);
                end
                
                %             non_nan = [];
                %             non_nan = ~isnan(resp{c,s}(:,k));
                
                non_nan = [];
                non_nan = ~isnan(resp{c,s}(circular_add{s},k));
                
                %             polarwitherrorbar函数输入为极坐标，不能直接用heading or rotation，转换成azimuth
                azi{s} = HeadingToazi(unique_cond_plot{s});
                
                % for spon firing rate higher than stim.
                if max(resp{c,s}(circular_add{s}(non_nan),k))>spon_resp{c,k}
                    polarwitherrorbar([azi{s}(non_nan)/180*pi], [resp{c,s}(circular_add{s}(non_nan),k)],[resp_err{c,s}(circular_add{s}(non_nan),k)],f{c,s}(1),max_scale);
                    hold on;
                    h_p = polar((0:10:360)/180*pi, spon_resp{c,k} * ones(1,length(0:10:360)),'k');
                    set(h_p,'LineWidth',2.5);
                    h_p = polar([HeadingToazi(pref_direc(c,k,s)) HeadingToazi(pref_direc(c,k,s))]/180*pi,[0 pref_amp(c,k,s)],f{c,s}(1));
                    set(h_p,'LineWidth',2.5);
                else
                    h_p = polar((0:10:360)/180*pi, spon_resp{c,k} * ones(1,length(0:10:360)),'k');
                    set(h_p,'LineWidth',2.5);
                    hold on;
                    polarwitherrorbar([azi{s}(non_nan)/180*pi], [resp{c,s}(circular_add{s}(non_nan),k)],[resp_err{c,s}(circular_add{s}(non_nan),k)],f{c,s}(1));
                    hold on;
                    h_p = polar([HeadingToazi(pref_direc(c,k,s)) HeadingToazi(pref_direc(c,k,s))]/180*pi,[0 pref_amp(c,k,s)],f{c,s}(1));
                    set(h_p,'LineWidth',2.5);
                end
                
                % annotation for rotation polar
                axes('position',[0.33 0.34 0.27 0.27]);
                text(0.47,0.99,'Exp','fontsize',10,'background','w');
                text(0.1,0.5,'CW','fontsize',10,'background','w');
                text(0.84,0.5,'CCW','fontsize',10,'background','w');
                text(0.47,0.01,'Con','fontsize',10,'background','w');
                axis off;
                
            end
        end
    end
    
    %% annotation
    for s = 1:2
        if s==1
            axes('position',[0.6,0.6, 0.3,0.3] );
            xlim( [0,100] );
            ylim( [0,length(unique_stim_type)*length(unique_coherence)+2] );
        elseif s==2
            axes('position',[0.6,0.3, 0.3,0.3] );
            xlim( [0,100] );
            ylim( [0,length(unique_stim_type)*length(unique_coherence)+2] );
        elseif s==4
            axes('position',[0.6,0.0, 0.3,0.3] );
            xlim( [0,100] );
            ylim( [0,length(unique_stim_type)*length(unique_coherence)+2] );
        end
        
        %     text(2,2.5,['prefer', repmat(' ',1,12), 'DDI',repmat(' ',1,12), 'HTI']);
        text(2,2.5,['prefer', repmat(' ',1,12), 'DDI']);
        text(5,1.7,['p', repmat(' ',1,12),'half-width']);
        for k=1:length(unique_stim_type)
            for c=1:length(unique_coherence)
                text(-1,2.2, num2str(pref_direc(c,k,s)) ); %vectorsum prefer
                text(20,2.2, num2str(DDI(c,k,1,s)) );
                %             text(37,2.2, num2str(HTI(c,k,s)) );
                text(-1,1.4, num2str(p_value(c,k,2,s)) ); % anova p without exp/con
                if IF_Fit
                    text(22,1.4, num2str(widGauss(c,k,s)) );
                end
            end
        end
        axis off;
    end
    
    % show DDI without exp and con
    axes('position',[0.40 0.07 0.1 0.2]);
    y = [DDI(1,1,2,1) DDI(1,1,2,2)];
    y_temp = diag(y);
    h_b = bar(y_temp,'stacked');
    h_b(1).FaceColor = 'r';
    h_b(2).FaceColor = 'b';
    title('DDI comparision (excluded Exp/Con)');
    ylim([0 1]);
    box off
    
    % show Spiral Index
    text(6,1,sprintf('Spiral Index = %0.3g',spiral_index));
    text(6,0.8,sprintf('D prime Index (±45) = %0.3g',d_prime_index_fine));
    text(6,0.6,sprintf('D prime Index (±90) = %0.3g',d_prime_index_coarse));
    text(6,0.4,sprintf('DDI Index (exclude Exp/Con) = %0.3g',DDI_exc_index));
    
    % show cell type
    text(6,0.2,sprintf('cell type = %g (1:T cell; 2: R cell; 3: S cell; 4: No tuning; 5: not actived)',cell_type));
    
    %% Plot RF
    xlRange{1} = 'A1:DG3000';
    xlRange{2} = 'A1:DG3000';
    
    if exist('Z:\Data\MOOG\Results\Result_LJY.xlsm')
        XlsData = ReadXls('Z:\Data\MOOG\Results\Result_LJY.xlsm',[2],3,xlRange); % 第一个sheet默认ringbell，第二个sheet默认Arthas
    elseif exist('Z:\Data\MOOG\Results\Result_Lwh.xlsm')
        XlsData = ReadXls('Z:\Data\MOOG\Results\Result_Lwh.xlsm',[1 2],3,xlRange);
    end
    num = XlsData.num;
    txt = XlsData.txt;
    raw = XlsData.raw;
    header = XlsData.header;
    FileNo = txt(:,header.FileNo);
    RF = num(:,header.X : header.H); % X	Y	W	H

    id = strmatch(FILE,FileNo);
    if ~isnan(RF(id,1))
        rf_x = RF(id,1);
        rf_y = RF(id,2);
        rf_w = RF(id,3);
        rf_h = RF(id,4);
        axes('position',[0.73,0.55, 0.3,0.3] )
        axis equal; box on; axis([-60 60 -60 60]);
        line([-60 60],[0 0],'color','k','LineStyle',':'); hold on;
        line([0 0],[-60 60],'color','k','LineStyle',':');
        set(gca,{'xtick','ytick'},{-60:20:60,-60:20:60});
        xlabel('degree');
        rectangle('position',[rf_x-rf_w/2 rf_y-rf_h/2, rf_w rf_h],...
            'Curvature',[0.3 0.3],'EdgeColor','b','LineWidth',1.5);
        title(sprintf('(%g,%g)',eyepos_x,eyepos_y),'fontsize',15);
    end

    
    %% Plot example neuron preferred direction in 3-d spiral space

    % make a matric in ploar axis
    amp_s = [];el_s=[];az_s=[];
    x=[];y=[];z=[];
    
    amp_s{1} = resp{1,1};
    el_s{1} = zeros(size(resp{1,1}));
    az_s{1} = HeadingToazi(unique_cond{1});
    
    amp_s{2} = resp{1,2};
    az_s{2} = spi_azimuth;
    % el_s{2} = spi_elevation;
    el_s{2} = [-180:45:135]';
    
    % plot the circle
    % amp_s{1} = repmat(100,8,1);
    % amp_s{2} = repmat(100,8,1);
    
    x{1} = amp_s{1}.*cosd(az_s{1});
    y{1} = amp_s{1}.*sind(az_s{1});
    z{1} = zeros(size(x{1}));
    
    y{2} = amp_s{2}.*cosd(el_s{2});
    z{2} = amp_s{2}.*sind(el_s{2});
    x{2} = zeros(size(y{2}));
    
    %加上preferred direction的点
    x{3} = x_pref;
    y{3} = y_pref;
    z{3} = z_pref;
    
    figure(11);clf;
    set(gcf,'color','white');
    set(11,'Position', [990,400 500,500], 'Name', '3-D Spiral Preference');
    % figure
    % % for spon
    % [x_spon,y_spon,z_spon]  = ellipsoid(0,0,0,spon_resp{1,1},spon_resp{1,1},spon_resp{1,1});
    % surf(x_spon,y_spon,z_spon,'facecolor',[.7 .7 .7]);
    % alpha(0.1);
    hold on;
    
    for s = 1:2
        % 画点
        scatter3(x{s},y{s},z{s},'MarkerEdgeColor',f{c,s}(1),'MarkerFaceColor',f{c,s}(1));
        
        % 连线
        plot3([x{s};x{s}(1)],[y{s};y{s}(1)],[z{s};z{s}(1)],'-','color',f{c,s}(1));
    end
    
    % vector-sum preference
    x_pref_plot = [0;x_pref];
    y_pref_plot = [0;y_pref];
    z_pref_plot = [0;z_pref];
    % plot3(x_pref_plot, y_pref_plot, z_pref_plot,'k.-','MarkerSize',10, 'LineWidth',2);
    h4 = quiver3(0,0,0, x_pref,y_pref,z_pref);
    set(h4,'maxheadsize',0.5,'LineWidth',1.5,'color','k');  %set the size
    
    % vector-sum preference在两个面上的投影
    h7 = quiver3(0,0,0, x_pref,y_pref,0); % heading plane
    set(h7,'maxheadsize',0,'LineWidth',1.5,'color','k','linestyle','--');
    h8 = quiver3(0,0,0, 0,y_pref,z_pref); % rotation plane
    set(h8,'maxheadsize',0,'LineWidth',1.5,'color','k','linestyle','--');
    
    % vector-sum preference在x-z plane 投影
    proj_vec = [x_pref 0 z_pref]; % 直接将vector_sum的y值取0即可
    h9 = quiver3(0,0,0, proj_vec(1),proj_vec(2),proj_vec(3)); % heading plane
    set(h9,'maxheadsize',0,'LineWidth',1.5,'color','g','linestyle','--');
    
    
    % single plane preference
    % heading plane:
    x_pref_h = single_pref_x(1);
    y_pref_h = single_pref_y(1);
    z_pref_h = single_pref_z(1);
    hold on
    h5 = quiver3(0,0,0, x_pref_h,y_pref_h,z_pref_h);
    set(h5,'maxheadsize',0.5,'LineWidth',1.5,'color','r');  %set the size
    
    % rotation plane
    x_pref_r = single_pref_z(2);
    y_pref_r = single_pref_y(2);
    z_pref_r = single_pref_x(2);
    hold on
    h6 = quiver3(0,0,0, x_pref_r,y_pref_r,z_pref_r);
    set(h6,'maxheadsize',0.5,'LineWidth',1.5,'color','b');  %set the size
    
    % annotation
    x_all = [x{1};x{2};x{3}];
    y_all = [y{1};y{2};y{3}];
    z_all = [z{1};z{2};z{3}];
    
    lims = max(abs([x_all;y_all;z_all]));
    lims = [-lims lims];
    
    plot3(lims,[0 0],[0 0],'LineWidth',1.5,'color','k') % for x-axis
    hold on
    plot3([0 0], lims,[0 0],'LineWidth',1.5,'color','k') % for y-axis
    plot3([0 0],[0 0], lims,'LineWidth',1.5,'color','k') % for z-axis
    set(gca,'ticklength',[0;0]);
    
    set(gca,'xlim',lims);
    set(gca,'xtick',lims,'xticklabel',{'Left','Right'});
    set(gca,'ylim',lims);
    set(gca,'ytick',lims,'yticklabel',{'Con','Exp'});
    set(gca,'zlim',lims);
    set(gca,'ztick',lims,'zticklabel',{'CW','CCW'}); % self-motion direction
    
    title(sprintf('Spiral Index = %g',spiral_index));
    
    view(20,20)
    axis equal
    
    % % 3-D 可视化 use raw data
    % x = cat(1,x{:});
    % y = cat(1,y{:});
    % z = cat(1,z{:});
    %
    % y = y(~isnan(y));
    % z = z(~isnan(z));
    % x = x(~isnan(z));
    %
    % X=[x,y,z];
    % C=convhulln(X);
    % for i = 1:size(C,1)
    %    j = C(i,[1 2 3 1]);
    %    patch(X(j,1),X(j,2),X(j,3),rand,'FaceAlpha',0.6);
    % end
    %
    % % Modify the view.
    % view(3), axis equal off tight vis3d; camzoom(1.5)
    % colormap(spring)
    
%     % unit sphere visulization
%     [x,y,z] = sphere(8);
%     figure;
%     surf(x,y,z)
%     axis equal
%     hold on
%     lims = [-1.5 1.5];
%     plot3(lims,[0 0],[0 0],'LineWidth',1.5,'color','k') % for x-axis
%     plot3([0 0], lims,[0 0],'LineWidth',1.5,'color','k') % for y-axis
%     plot3([0 0],[0 0], lims,'LineWidth',1.5,'color','k') % for z-axis
%     set(gca,'xlim',lims);
%     set(gca,'xtick',[-1.5 -1 0 1 1.5],'xticklabel',{'Left','-1','0','1','Right'});
%     set(gca,'ylim',lims);
%     set(gca,'ytick',[-1.5 -1 0 1 1.5],'yticklabel',{'Con','-1','0','1','Exp'});
%     set(gca,'zlim',lims);
%     set(gca,'ztick',[-1.5 -1 0 1 1.5],'zticklabel',{'CW','-1','0','1','CCW'});
%     
%     %     % heading marker
%     %     text(x(5,1),y(5,1),z(5,1),'-90');
%     %     text(x(5,2),y(5,2),z(5,2),'-135');
%     %     text(x(5,3),y(5,3),z(5,3),'180');
%     %     text(x(5,4),y(5,4),z(5,4),'135');
%     %     text(x(5,5),y(5,5),z(5,5),'90');
%     %     text(x(5,6),y(5,6),z(5,6),'45');
%     %     text(x(5,7),y(5,7),z(5,7),'0');
%     %     text(x(5,8),y(5,8),z(5,8),'-45');
%     %     %     text(x(5,9),y(5,9),z(5,9),'-90'); % = (5,1)
%     
%     %     % rotation marker: 间隔为22.5°，22.5通过插值获取
%     %     text(x(5,7),y(5,7),z(5,7),'Exp');
%     %     text(x(5,3),y(5,3),z(5,3),'Con');
%     %     text(x(7,7),y(7,7),z(7,7),'45');
%     %     text(x(7,3),y(7,3),z(7,3),'135');
%     %     text(x(3,7),y(3,7),z(3,7),'-45');
%     %     text(x(3,3),y(3,3),z(3,3),'-135');
%     %     text(x(1,:),y(1,:),z(1,:),'-90');
%     %     text(x(1,:),y(1,:),z(9,:),'90');
%     
%     
%     col = spon_resp{1,1}.*ones(size(x)); % spon response
%     % replace heading response
%     col(5,:) = resp{1,1}([3 2 1 8 7 6 5 4 3]);
%     % replace rotation response
%     col(7,7) = resp{1,2}(unique_cond{2}==45);
%     col(7,3) = resp{1,2}(unique_cond{2}==135);
%     col(3,7) = resp{1,2}(unique_cond{2}==-45);
%     col(3,3) = resp{1,2}(unique_cond{2}==-135);
%     col(1,:) = resp{1,2}(unique_cond{2}==-90);
%     col(9,:) = resp{1,2}(unique_cond{2}==90);
%     surf(x,y,z,col)
    
end % end IF_FIGURE


%% 2D gaussian fitting
% xg = [];
% yg = [];
% zg = [];
% 
% % heading plane
% for i = 1:length(unique_cond{1})
%     for j = 1:length(resp_trial{1,1,1}(:,i))
%         xg = [xg; unique_cond{1}(i).*pi/180];
%         yg = [yg; 0;];
%         zg =[zg; resp_trial{1,1,1}(j,i)];
%     end
% end
% 
% % rotation plane
% % remove the duplicate of exp and con
% unique_cond_r = unique_cond{2}([2:4,6:8]).*pi/180;
% resp_trial_r = resp_trial{1,1,2};
% for i = 1:length(unique_cond_r)
%     for j = 1:length(resp_trial_r(:,i))
%         xg = [xg; 0;];
%         yg = [yg; unique_cond_r(i)];
%         zg =[zg; resp_trial_r(j,i)];
%     end
% end
% 
% c = [xg,yg,zg];
% 
% 
% 
% % 1. fit 'NonlinearLeastSquares': % Algorithm: 'levenberg-marquardt' or 'Trust-region'
% fo_ = fitoptions('method','NonlinearLeastSquares','Robust','On','Algorithm','T',...
%     'StartPoint',[10 10 0 1.4284 0 1.4284],'Lower',[0 0 -pi 0 -pi 0],'Upper',[Inf 1000 pi 500 pi 500],'MaxFunEvals',5000);
% ft_ = fittype('a*exp(-(x-prefx).^2/2/sigmax^2 -(y-prefy).^2/2/sigmay^2) + b',...
%     'dependent',{'z'},'independent',{'x','y'},...
%     'coefficients',{'a', 'b', 'prefx', 'sigmax', 'prefy', 'sigmay'});
% 
% [cf2D_, goodness, output]= fit([xg,yg],zg,ft_,fo_);
% [xx,yy] = meshgrid(-pi:pi/6:pi,-pi:pi/6:pi);
% 
% zz = feval(cf2D_,xx,yy);
% figure
% plot3(xx,yy,zz,'ko'); %Lwh20170207
% 
% % 2. lsqcurvefit
% fun1 = @(a,x) a(1)*exp(-(x(1,:)-a(3)).^2/2/a(4)^2 -(x(2,:)-a(5)).^2/2/a(6)^2) + a(2);
% % 拟合初始值
% a0 = [10 10 0 1.4 0 1.4];
% lb = [-Inf -1000 -pi 0 -pi 0];
% ub = [Inf 1000 pi 500 pi 500];
% 
% options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt',... % 'trust-region-reflective'
%     'FiniteDifferenceType','central',... % 'forward'
%     'MaxIterations',1000); 
% % [a,resnorm]=lsqcurvefit(fun1,ones(1, 6),[xg';yg'],zg',lb,ub); %resnorm残差平方和  X是参数
% [a,resnorm]=lsqcurvefit(fun1,ones(1, 6),[x_v4interp(:)';y_v4interp(:)'],z_v4interp(:)',lb,ub); %resnorm残差平方和  X是参数
% 
% figure
% Zfit = a(1)*exp(-(xx-a(3)).^2/2/a(4)^2 -(yy-a(5)).^2/2/a(6)^2) + a(2);
% surf(xx, yy, Zfit);
% hold on
% % scatter3(xg, yg, zg, 50);
% scatter3(deg2rad(x_v4interp(:)), deg2rad(y_v4interp(:)), deg2rad(z_v4interp(:)), 50);


%% 3D fitting
% von Mises–Fisher distribution
% Spherical Harmonics


%% plot recorsing location on MRI
if IF_DrawRecLoc
    [session, hemi, xloc, yloc, guidetube, offset, depth_temp, file_name, rf_temp, fixation_x, fixation_y, ab_stim] = GetBatchfileInfo(1);
    if ~isnan(session{1})
        depth_temp = str2double(depth_temp);
        ind = strmatch(FILE,file_name);
        if length(ind) == 1
            
            % plot recorsing location on MRI
            if ~isnan(depth_temp)
                xSelect = str2double(xloc(ind));
                ySelect = str2double(yloc(ind));
                %                 try
                %                     temp = get(803);
                %                 catch
                DrawMapping_Lwh_2020(1, [xSelect ySelect]);
                %                 end
                
                jframe=getJFrame(803);
                jframe.setAlwaysOnTop(1); %将figure置顶
                
                offSet = round((str2double(guidetube(ind))  + str2double(offset(ind)) - 2.0) * 100);  % Related to Guide Tube 2.0 cm!!
                depth = depth_temp(ind);
                xx = ySelect;
                yy = offSet + round(depth/100);
                hold on
                if ~isempty(batch_flag)
                    try
                        lh = findall(803,'type','line');
                        set(lh,'XData',[],'YData',[]); % 在Batch时去掉之前画的location marker
                        th = findall(803,'tag','file_marker');
                        set(th,'string',''); % 在Batch时去掉之前画的location marker
                    end
                end
                figure(803);
                plot(xx,yy,'go','markersize',10,'markerfacecolor','g');
                file_marker = text(12,250,sprintf('%s',FILE),'fontsize',13,'color','w');
                file_marker.Tag = 'file_marker';
            end
        end
    end
end






%% 只plot forward motion 平面（对应 CueS task）
az_aroundF = ([0:90:270]);

% 45 degree around forward motion R-CW-L-CCW
resp_aroundF{1} = [resp{1}(6) resp{2}(4) resp{1}(4) resp{2}(6)];
resp_trial_aroundF{1} = [resp_trial{1}(:,6) resp_trial{2}(:,4) resp_trial{1}(:,4) resp_trial{2}(:,6)];

% 90 degree around forward motion
resp_aroundF{2} = [resp{1}(7) resp{2}(3) resp{1}(3) resp{2}(7)];
resp_trial_aroundF{2} = [resp_trial{1}(:,7) resp_trial{2}(:,3) resp_trial{1}(:,3) resp_trial{2}(:,7)];

for ii = 1:2 % fine and coarse
    % vector sum
    [az, ~, amp,~,~,~] = vectorsumSpiral(resp_aroundF{ii},az_aroundF);
    pref_direc_aroundF(ii) = deg2rad(az);
    pref_amp_aroundF(ii) = amp;
    
    % anova
    [p_aroundF(ii),~,s] = anova1(resp_trial_aroundF{ii},'','off');
    
    % 两两比较
    [comp,~] = multcompare(s,'Display','off'); % 第6列时p value
    
    % rearrangement com
    [~,Indx] = sort(resp_aroundF{ii}); % R-CW-L-CCW
    Indx = fliplr(Indx); % 从大到小
    nn = 1;
    for i = 1:3
        for j = i+1:4
            comp_indx(nn) = find((comp(:,1)==Indx(i)& comp(:,2)==Indx(j)) | (comp(:,2)==Indx(i)& comp(:,1)==Indx(j)));
            nn = nn+1;
        end
    end
    comp_new = comp(comp_indx',:);
    
    % 四个轴的mean FR从大到小命名为A，B，C，D
    % 暂时分四类
    if comp_new(1,6)<0.05 && comp_new(2,6)<0.05 && comp_new(3,6)<0.05  % 单轴preffer：A>>B，C，D
        pref_axis_num(ii) = 1;
        pref_axis{ii} = Indx(1);
    elseif comp_new(2,6)<0.05 && comp_new(3,6)<0.05 && comp_new(5,6)<0.05 && comp_new(1,6)>0.05 % 双轴pref：A>>CD,B>>CD,A~B
        pref_axis_num(ii) = 2;
        pref_axis{ii} = [Indx(1) Indx(2)];
    elseif comp_new(3,6)<0.05 && comp_new(5,6)<0.05 && comp_new(6,6)<0.05 % A>>D,B>>D,C>>D
        pref_axis_num(ii) = 3;
        pref_axis{ii} = [Indx(1) Indx(2) Indx(3)];
    else
        pref_axis_num(ii) = 0;
        pref_axis{ii} = nan;
    end
end

%% %% 两个plane围成的面积
% s_h = polyarea(x{1},y{1});
% s_r = polyarea(y{2},z{2});

%% Data Saving
% Reorganized. HH20141124
config.batch_flag = batch_flag;

% % Output information for test. HH20160415
% if isempty(batch_flag)
%     config.batch_flag = 'test.m';
%     disp('Saving results to \batch\test\ ');
% end

%%%%%%%%%%%%%%%%%%%% Change here %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if IF_Fit
    result = PackResult(FILE, SpikeChan, mean_repetitionN, unique_stim_type, ... % Obligatory!!
        eyepos_x, eyepos_y, unique_coherence,stars_type,unique_heading,unique_rotation,unique_cond,...
        p_sti, p_sti_min, pref_direc, pref_el, pref_amp, p_value,p_value_mt, ...
        pref_direc_norm, pref_el_norm, pref_amp_norm,...
        all_slope,ahead_slope,d_prime, d_prime_fine, d_prime_coarse, d_prime_maxDiff, d_prime_maxAxis,d_prime_maxFR,...
        DDI, p_DDI, HTI, modulation_index, neu_sensitivity, neu_sen_fine, neu_sen_coarse,pair_DI,nonpair_DI,nonpair_DI2,...
        ahead_slope_index, d_prime_index,d_prime_index_coarse,d_prime_index_fine, d_prime_maxDiff_index, d_prime_maxAxis_index, d_prime_maxFR_index,...
        d_prime_single, d_prime_single_fine, d_prime_single_coarse,...
        DDI_inc_index, DDI_exc_index, HTI_index, neu_sen_index,pair_TRDI,nonpair_TRDI,nonpair_TRDI2,...
        fitOK, widGauss, prefGauss,maxFR,maxFRGauss,GaussfitRsquare, cf_, ...
        resp,resp_std,spon_resp,spon_resp_std,resp_trial,resp_norm,...
        MI, normalMI,spiral2Tran, spiral2Rot, spiral_index,spiral_index_angle,SVM_index,SVM_correct_rate, AuROC,AUC_index,SI,...
        fr_axis_ratio, vec_axis_ratio,...
        pref_direc_aroundF,pref_amp_aroundF,p_aroundF,pref_axis_num,pref_axis,cell_type,mt_type,...
        xx_ori, yy_ori, zz_ori, x_v4interp, y_v4interp, z_v4interp);
else
    result = PackResult(FILE, SpikeChan, mean_repetitionN, unique_stim_type, ... % Obligatory!!
        eyepos_x, eyepos_y, unique_coherence,stars_type,unique_heading,unique_rotation,unique_cond,...
        p_sti, p_sti_min, pref_direc, pref_el, pref_amp, p_value,p_value_mt, ...
        pref_direc_norm, pref_el_norm, pref_amp_norm,...
        all_slope,ahead_slope,d_prime, d_prime_fine, d_prime_coarse, d_prime_maxDiff, d_prime_maxAxis,d_prime_maxFR,...
        DDI, p_DDI, HTI, modulation_index, neu_sensitivity, neu_sen_fine, neu_sen_coarse,pair_DI,nonpair_DI,nonpair_DI2,...
        ahead_slope_index, d_prime_index,d_prime_index_coarse,d_prime_index_fine, d_prime_maxDiff_index, d_prime_maxAxis_index, d_prime_maxFR_index,...
        d_prime_single, d_prime_single_fine, d_prime_single_coarse,...
        DDI_inc_index, DDI_exc_index, HTI_index, neu_sen_index,pair_TRDI,nonpair_TRDI,nonpair_TRDI2,...
        maxFR, ...
        resp,resp_std,spon_resp,spon_resp_std,resp_trial,resp_norm,...
        MI, normalMI,spiral2Tran, spiral2Rot, spiral_index,spiral_index_angle,SVM_index,SVM_correct_rate, AuROC,AUC_index,SI,...
        fr_axis_ratio, vec_axis_ratio,...
        pref_direc_aroundF,pref_amp_aroundF,p_aroundF,pref_axis_num,pref_axis,cell_type,mt_type,...
        xx_ori, yy_ori, zz_ori, x_v4interp, y_v4interp, z_v4interp);
end

config.suffix = 'SpiT';
config.xls_column_begin = 'p_stim';
config.xls_column_end = 'p_value_all_cond';

% Figures to save
% config.save_figures = gcf;
if IF_FIGURE
    config.save_figures(1) = figure(10);  % changed by Lwh
    config.save_figures(2) = figure(11);  % changed by Lwh
end

% % % Only once
config.sprint_once_marker = 'gggggg';
config.sprint_once_contents = 'result.p_sti_min, result.cell_type, result.mean_repetitionN, result.d_prime_index_coarse, result.spiral_index, result.p_value_mt';

% % % Loop across each stim_type
config.sprint_loop_marker = {
    'gg';};
config.sprint_loop_contents = {
    'result.pref_direc(1,1,ind), result.p_value(1,1,r,ind)';
    };


% % % Only once   %不需要循环输出到excel的数据

% % % Loop across each stim_type  %按照几列为循环输出到excel中
% config.sprint_loop_marker = {'ggggg'};
% config.sprint_loop_contents = {'result.p_sti(c,k,s), result.az(c,k,s),result.DDI(c,k,s), result.HTI(c,k,s), result.p_value(c,k,s)'};

config.append = 1; % Overwrite or append

SaveResult(config, result, IF_FIGURE, IF_EXCEL);

if isempty(config.batch_flag)
    disp('This is Done~ ');
end


