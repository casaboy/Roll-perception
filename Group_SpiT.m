function function_handles_spi = Group_SpiT(loadM_data,handles,p_cri)

global group_result;

% move Group_Classification.m code here:
%% Extract some data to global value for easier access
% 将常用的参数提取出来
% global group_result;  % 只选取main depth位置的spiral tuning数据 (Depth_code == 0)
for n = 1:length(loadM_data)
    main_depth_loc = [];
    try
        main_depth_loc = [(loadM_data(n).Data.Depth_code)]==0;
    catch
        keyboard
    end
    
    % if length(loadM_data(n).unique_SpikeChan)==1 % 只有一个MU
    % ******* 只考虑第一个unit ********
    select_tuning_data = loadM_data(n).Data(main_depth_loc).SpikeChan(1).centerEye_Protocol{1}; % 在main depth中盯原点(0,0)做了spiral tuning task
    
    if ~(isfield(select_tuning_data,'FILE') && size(select_tuning_data,2)) % main depth 没有做spiral tuning
        % 找出是否在别的深度也做了tuning task
        depth_ori = [loadM_data(n).Data.Depth_code];
        % 非main depth
        depth_up = fliplr(depth_ori(depth_ori<0));
        depth_down = depth_ori(depth_ori>0);
        
        depth_temp = nan(2,max(length(depth_up),length(depth_down)));
        depth_temp(1,1:length(depth_up)) = depth_up;
        depth_temp(2,1:length(depth_down)) = depth_down;
        depth_temp = depth_temp(:);  % 查找顺序-1, 1, -2, 2, -3, 3, .....
        depth_temp(isnan(depth_temp)) = [];
        
        if ~isempty(depth_temp)
            for dep = 1:length(depth_temp)
                select_tuning_data = loadM_data(n).Data(depth_ori == depth_temp(dep)).SpikeChan(1).centerEye_Protocol{1};
                
                if ~isempty(select_tuning_data)
                    break % 找到有
                end
            end
        end
    end
    
    % basic properties
    group_result(n).Date = loadM_data(n).Date;
    group_result(n).monkey = loadM_data(n).monkey;
    group_result(n).cellID = loadM_data(n).cellID;
    group_result(n).Area = loadM_data(n).Area;
    group_result(n).RF = loadM_data(n).RF;
    group_result(n).unique_SpikeChan = loadM_data(n).unique_SpikeChan;
    group_result(n).GridX = loadM_data(n).GridX;
    group_result(n).GridY = loadM_data(n).GridY;
    group_result(n).main_depth = loadM_data(n).main_depth;
    group_result(n).guidetube = loadM_data(n).guidetube;
    group_result(n).offset = loadM_data(n).offset;
    

    if isfield(select_tuning_data,'FILE') && size(select_tuning_data,2) % 找到该unit做了spiral tuning
        group_result(n).SpiT_FILE = select_tuning_data.FILE;
        group_result(n).SpiT_repN = select_tuning_data.mean_repetitionN;
        group_result(n).stars_type = select_tuning_data.stars_type; % 0:圆点，1：三角形
        group_result(n).SpiT_coherence = select_tuning_data.unique_coherence;
        group_result(n).unique_cond = select_tuning_data.unique_cond;
        group_result(n).unique_heading = select_tuning_data.unique_heading;
        group_result(n).unique_rotation = select_tuning_data.unique_rotation;
        group_result(n).MI = squeeze(select_tuning_data.MI)';
        group_result(n).normalMI = squeeze(select_tuning_data.normalMI)';
        group_result(n).pref_direc = squeeze(select_tuning_data.pref_direc)';
        group_result(n).pref_el = squeeze(select_tuning_data.pref_el)';
        group_result(n).pref_amp = squeeze(select_tuning_data.pref_amp)';
        group_result(n).pref_direc_norm = squeeze(select_tuning_data.pref_direc_norm)';
        group_result(n).pref_el_norm = squeeze(select_tuning_data.pref_el_norm)';
        group_result(n).pref_amp_norm = squeeze(select_tuning_data.pref_amp_norm)';
        group_result(n).fitOK = select_tuning_data.fitOK;
        group_result(n).widGauss = squeeze(select_tuning_data.widGauss)';
        group_result(n).prefGauss = squeeze(select_tuning_data.prefGauss)';
        group_result(n).maxFRGauss = squeeze(select_tuning_data.maxFRGauss)';
        group_result(n).GaussfitRsquare = select_tuning_data.GaussfitRsquare;
        group_result(n).maxFR = squeeze(select_tuning_data.maxFR)';
        group_result(n).sponFR = cell2mat(select_tuning_data.spon_resp);
        
        % compare to spon: ranksum. one tail (whether response above spon)
        group_result(n).p_stim = select_tuning_data.p_sti_min;  % min
        group_result(n).p_stim_raw = select_tuning_data.p_sti; % p_sti(c,k,r,s)  r=1:没有排除exp/Con,  r=2:排除exp/con
        group_result(n).p_stim_h = group_result(n).p_stim_raw(1,1,2,1); % exclude exp/con, translation response above spon
        group_result(n).p_stim_r = group_result(n).p_stim_raw(1,1,2,2); % exclude exp/con, rotation response above spon
        
        % ANOVA test, tuning or not
        p_turn_temp = select_tuning_data.p_value; % p_turn(c,k,r,s)
        group_result(n).p_turn_h = p_turn_temp(1,1,2,1); % 排除exp/con
        group_result(n).p_turn_r = p_turn_temp(1,1,2,2); % 排除exp/con
        group_result(n).p_turn_spi = p_turn_temp(1,1,1,3); % 包括exp/con
        group_result(n).p_turn_mt = squeeze(select_tuning_data.p_value_mt); % % heading与rotation plane有显著区别(n).p_turn_mt = squeeze(select_tuning_data.p_value_mt); % % heading与rotation plane有显著区别
        group_result(n).p_value_includeExpCon = squeeze(select_tuning_data.p_value(:,:,1,:))'; % 第一列T，第二列R，第三列Spiral Space
        group_result(n).p_value_excludeExpCon = squeeze(select_tuning_data.p_value(:,:,2,:))'; % 第一列T，第二列R
        group_result(n).spiral_index = select_tuning_data.spiral_index;
        group_result(n).spiral_index_angle = select_tuning_data.spiral_index_angle;
        
        % global_pref_code 1: 左，CW，rotation；  2: 右，CCW，heading
        % 1. h and r
        for tt = 1:2
            if group_result(n).pref_direc(tt) < 0
                group_result(n).global_pref_code(tt) = 1; % 1: 左，CW，rotation；  2: 右，CCW，heading
            elseif group_result(n).pref_direc(tt) > 0
                group_result(n).global_pref_code(tt) = 2;
            else
                group_result(n).global_pref_code(tt) = nan;
            end
        end
        group_result(n).global_pref_p(1) = group_result(n).p_turn_h; % ANOVA test
        group_result(n).global_pref_p(2) = group_result(n).p_turn_r;
        group_result(n).global_pref_p(3) = group_result(n).p_turn_spi;
        
        % 2. mt: how to define ?
        if group_result(n).spiral_index < 0
            group_result(n).global_pref_code(3) = 1; % 1: 左，CW，rotation；  2: 右，CCW，heading
        elseif group_result(n).spiral_index > 0
            group_result(n).global_pref_code(3) = 2; % 1: 左，CW，rotation；  2: 右，CCW，heading
        else
            group_result(n).global_pref_code(3) = nan;
        end
        
        % modulation index
        group_result(n).ahead_slope = select_tuning_data.ahead_slope;
        group_result(n).d_prime_fine = squeeze(select_tuning_data.d_prime_fine)'; % ±45,第一列T，第二列R
        group_result(n).d_prime_coarse = squeeze(select_tuning_data.d_prime_coarse)'; % ±90,第一列T，第二列R
        group_result(n).d_prime_maxDiff = squeeze(select_tuning_data.d_prime_maxDiff)'; % FR差值最大两点
        group_result(n).d_prime_maxAxis = squeeze(select_tuning_data.d_prime_maxAxis)'; % 差值最大的一条轴
        group_result(n).d_prime_maxFR = squeeze(select_tuning_data.d_prime_maxFR)'; % FR最大点及其对侧
        group_result(n).DDI_includeExpCon = squeeze(select_tuning_data.DDI(:,:,1,:))'; % 第一列T，第二列R
        group_result(n).DDI_excludeExpCon = squeeze(select_tuning_data.DDI(:,:,2,:))'; % 第一列T，第二列R
        group_result(n).p_DDI_includeExpCon = squeeze(select_tuning_data.p_DDI(:,:,1,:))'; % 第一列T，第二列R
        group_result(n).p_DDI_excludeExpCon = squeeze(select_tuning_data.p_DDI(:,:,2,:))'; % 第一列T，第二列R
        group_result(n).HTI = squeeze(select_tuning_data.HTI)';
        group_result(n).MDI = squeeze(select_tuning_data.modulation_index)';
        group_result(n).neu_sen_fine = select_tuning_data.neu_sen_fine;
        group_result(n).neu_sen_coarse = select_tuning_data.neu_sen_coarse;
        group_result(n).pair_DI = select_tuning_data.pair_DI;
        group_result(n).nonpair_DI = select_tuning_data.nonpair_DI;
        group_result(n).nonpair_DI2 = select_tuning_data.nonpair_DI2;
        group_result(n).SI = select_tuning_data.SI;
        group_result(n).d_prime_single_fine = select_tuning_data.d_prime_single_fine; % 每个特定角度与平均FR求d’, fine中是L,R,CW,CCW 45°
        group_result(n).d_prime_single_coarse = select_tuning_data.d_prime_single_coarse; % 每个特定角度与平均FR求d’，coarse中是L,R,CW,CCW 90°
        
        % T-R index
        group_result(n).ahead_slope_index = select_tuning_data.ahead_slope_index;
        group_result(n).d_prime_index_fine = select_tuning_data.d_prime_index_fine;
        group_result(n).d_prime_index_coarse = select_tuning_data.d_prime_index_coarse;
        group_result(n).d_prime_maxDiff_index = select_tuning_data.d_prime_maxDiff_index;
        group_result(n).d_prime_maxAxis_index = select_tuning_data.d_prime_maxAxis_index;
        group_result(n).d_prime_maxFR_index = select_tuning_data.d_prime_maxFR_index;
        group_result(n).DDI_inc_index = select_tuning_data.DDI_inc_index;
        group_result(n).DDI_exc_index = select_tuning_data.DDI_exc_index;
        group_result(n).HTI_index = select_tuning_data.HTI_index;
        group_result(n).neu_sen_index_fine = select_tuning_data.neu_sen_index(1);
        group_result(n).neu_sen_index_coarse = select_tuning_data.neu_sen_index(2);
        group_result(n).pair_TRDI = select_tuning_data.pair_TRDI;
        group_result(n).nonpair_TRDI = select_tuning_data.nonpair_TRDI;
        group_result(n).nonpair_TRDI2 = select_tuning_data.nonpair_TRDI2;

        % separable predictions, spiral index
        group_result(n).SVM_correct_rate = select_tuning_data.SVM_correct_rate;
        group_result(n).SVM_index = select_tuning_data.SVM_index;
        group_result(n).AuROC = select_tuning_data.AuROC;
        group_result(n).AUC_index = select_tuning_data.AUC_index;
            
        % raw data
        % first row: translation tuning
        % second row: rotation tuning
        %      -180,  -135,  -90,  -45,  0,  45,  90,  135
        % H  backward       left      forward    right
        % R  backward        CW       forward    CCW
        group_result(n).spit_resp = cell2mat(select_tuning_data.resp)';
        group_result(n).resp_norm = cell2mat(select_tuning_data.resp_norm)';
        group_result(n).spit_resp_trial = squeeze(select_tuning_data.resp_trial);
        
        %         % spike count mean and variance across cell
        %         group_result(n).mean_fr(1) = nanmean(group_result(n).spit_resp_trial{1}(:)); % translation
        %         group_result(n).mean_fr(2) = nanmean(group_result(n).spit_resp_trial{2}(:)); % rotation
        %         group_result(n).std_fr(1) = nanstd(group_result(n).spit_resp_trial{1}(:)); % translation
        %         group_result(n).std_fr(2) = nanstd(group_result(n).spit_resp_trial{2}(:)); % rotation
        
        % spike count mean and variance across trial
        % 该细胞在每个condition下所有trial的平均反应，排除Exp，Con
        % tranlsation trial
        select_t_trial = logical(group_result(n).unique_heading ~=0 & group_result(n).unique_heading~=180 & group_result(n).unique_heading~=-180)';
        group_result(n).mean_fr_trial(:,1) = nanmean(group_result(n).spit_resp_trial{1}(:,select_t_trial));
        group_result(n).std_fr_trial(:,1) = nanstd(group_result(n).spit_resp_trial{1}(:,select_t_trial));
        
        % rotation trial
        select_r_trial = logical(group_result(n).unique_rotation ~=0 & group_result(n).unique_rotation~=180 & group_result(n).unique_rotation~=-180)';
        group_result(n).mean_fr_trial(:,2) = nanmean(group_result(n).spit_resp_trial{2}(:,select_r_trial));
        group_result(n).std_fr_trial(:,2) = nanstd(group_result(n).spit_resp_trial{2}(:,select_r_trial));
        
        % other
        % around forward motion, % first is fine, second is coarse
        group_result(n).pref_direc_aroundF = squeeze(select_tuning_data.pref_direc_aroundF);
        group_result(n).pref_amp_aroundF = squeeze(select_tuning_data.pref_amp_aroundF);
        group_result(n).p_aroundF = squeeze(select_tuning_data.p_aroundF);
        group_result(n).pref_axis_num = squeeze(select_tuning_data.pref_axis_num);
        group_result(n).pref_axis = squeeze(select_tuning_data.pref_axis);
        
        % 提取forward motion附近的点的RF raw data，R-CW-L-CCW
        unique_heading = group_result(n).unique_heading;
        unique_rotation = group_result(n).unique_rotation;
        h_0 = find(unique_heading==0);
        r_0 = find(unique_rotation==0);
        h_p90 = unique_heading==90; % R
        h_n90 = unique_heading==-90; % L
        r_p90 = unique_rotation==90; % CCW
        r_n90 = unique_rotation==-90; % CW
        
        % fine, 45 degree(或者靠近0度的角度) around forward motion
        temp_L_fine = group_result(n).spit_resp_trial{1}(:,h_0-1);
        temp_R_fine = group_result(n).spit_resp_trial{1}(:,h_0+1);
        temp_CW_fine = group_result(n).spit_resp_trial{2}(:,r_0-1);
        temp_CCW_fine = group_result(n).spit_resp_trial{2}(:,r_0+1);
        
        % coarse, 90 degree around forward motion
        temp_R_coarse = group_result(n).spit_resp_trial{1}(:,h_p90);
        temp_L_coarse = group_result(n).spit_resp_trial{1}(:,h_n90);
        temp_CCW_coarse = group_result(n).spit_resp_trial{2}(:,r_p90);
        temp_CW_coarse = group_result(n).spit_resp_trial{2}(:,r_n90);
        
        group_result(n).spit_resp_0_fine = [temp_R_fine,temp_CW_fine,temp_L_fine,temp_CCW_fine];
        group_result(n).spit_resp_0_coarse = [temp_R_coarse,temp_CW_coarse,temp_L_coarse,temp_CCW_coarse];
        
        % new index for T-R index, 20210826
        % remove the max FR lower than 10 hz
        if max([mean(temp_R_coarse), mean(temp_L_coarse), mean(temp_CCW_coarse), mean(temp_CW_coarse)]) < 5
            keyboard
            group_result(n).d_prime_ratio = nan;
            group_result(n).fr_axis_ratio = nan;
            group_result(n).vec_axis_ratio = nan;
        else
            % 1. d prime ratio
            group_result(n).d_prime_ratio = abs(group_result(n).d_prime_coarse(1) / group_result(n).d_prime_coarse(2)); 
            
            % 2. LR firing rate different / CWCCW different
            LR_diff = abs(mean(temp_R_coarse) - mean(temp_L_coarse));
            CWCCW_diff = abs(mean(temp_CCW_coarse) - mean(temp_CW_coarse));
            group_result(n).fr_axis_ratio = LR_diff / CWCCW_diff; % >1: more T; <1: more R
            
            % 3. vector-sum projection on translation axis and rotatino axis
            axis_proj = group_result(n).pref_amp .* cosd(HeadingToazi(group_result(n).pref_direc));
            group_result(n).vec_axis_ratio = abs(axis_proj(1) / axis_proj(2)); % >1: more T; <1: more R
        end
        
        % gaussian fitting result
        try
            group_result(n).fittingresult = select_tuning_data.cf_;
        catch
            group_result(n).fittingresult = [];
        end
    else
        group_result(n).SpiT_FILE = nan;
        group_result(n).SpiT_repN = nan;
        group_result(n).stars_type = nan;
        group_result(n).SpiT_coherence = nan;
        group_result(n).unique_cond = {[nan],[nan]};
        group_result(n).unique_heading = nan;
        group_result(n).unique_rotation = nan;
        group_result(n).MI = nan(1,2);
        group_result(n).normalMI = nan(1,2);
        
        group_result(n).pref_direc = nan(1,3);
        group_result(n).pref_el = nan(1,3);
        group_result(n).pref_amp = nan(1,3);
        group_result(n).pref_direc_norm = nan(1,3);
        group_result(n).pref_el_norm = nan(1,3);
        group_result(n).pref_amp_norm = nan(1,3);
        group_result(n).fitOK = nan(1,2);
        group_result(n).widGauss = nan(1,2);
        group_result(n).prefGauss = nan(1,2);
        group_result(n).maxFRGauss = nan(1,2);
        group_result(n).GaussfitRsquare = nan(1,2);
        group_result(n).maxFR = nan(1,3);
        group_result(n).sponFR = nan;
        
        % compare to spon: ranksum. one tail (whether response above spon)
        group_result(n).p_stim = nan; % min
        group_result(n).p_stim_raw = nan; % p_stim_raw(c,k,r,s)  r=1:没有排除exp/Con,  r=2:排除exp/con
        group_result(n).p_stim_h = nan; % exclude exp/con, translation response above spon
        group_result(n).p_stim_r = nan; % exclude exp/con, rotation response above spon
        
        % ANOVA test, tuning or not
        group_result(n).p_turn_h = nan;
        group_result(n).p_turn_r = nan;
        group_result(n).p_turn_spi = nan;
        group_result(n).p_turn_mt = nan;
        group_result(n).p_value_includeExpCon = nan(1,3);
        group_result(n).p_value_excludeExpCon = nan(1,3);
        
        % modulation index
        group_result(n).ahead_slope = nan(1,2);
        group_result(n).d_prime_fine = nan(1,2);
        group_result(n).d_prime_coarse = nan(1,2);
        group_result(n).d_prime_maxDiff = nan(1,2);
        group_result(n).d_prime_maxAxis = nan(1,2);
        group_result(n).d_prime_maxFR = nan(1,2);
        group_result(n).DDI_includeExpCon = nan(1,2);
        group_result(n).DDI_excludeExpCon = nan(1,2);
        group_result(n).p_DDI_includeExpCon = nan(1,2);
        group_result(n).p_DDI_excludeExpCon = nan(1,2);
        group_result(n).HTI = nan(1,2);
        group_result(n).MDI = nan(1,2);
        group_result(n).neu_sen_fine = nan(1,2);
        group_result(n).neu_sen_coarse =nan(1,2);
        group_result(n).pair_DI = nan(1,2);
        group_result(n).nonpair_DI = nan(1,2);
        group_result(n).nonpair_DI2 = nan(1,2);
        group_result(n).SI = nan;
        group_result(n).d_prime_single_fine = nan(1,4);
        group_result(n).d_prime_single_coarse = nan(1,4);
        
        % T-R index
        group_result(n).ahead_slope_index = nan;
        group_result(n).d_prime_index_fine = nan;
        group_result(n).d_prime_index_coarse = nan;
        group_result(n).d_prime_maxDiff_index = nan;
        group_result(n).d_prime_maxAxis_index = nan;
        group_result(n).d_prime_maxFR_index = nan;
        group_result(n).DDI_inc_index = nan;
        group_result(n).DDI_exc_index = nan;
        group_result(n).HTI_index = nan;
        group_result(n).neu_sen_index_fine = nan;
        group_result(n).neu_sen_index_coarse = nan;
        group_result(n).pair_TRDI = nan;
        group_result(n).nonpair_TRDI = nan;
        group_result(n).nonpair_TRDI2 = nan;
        
        % separable predictions, spiral index
        group_result(n).spiral_index = nan;
        group_result(n).spiral_index_angle = nan;
        group_result(n).SVM_correct_rate = nan(1,2);
        group_result(n).SVM_index = nan;
        group_result(n).AuROC = nan(1,2);
        group_result(n).AUC_index = nan;
        
        % global_pref_code 1: 左，CW，rotation；  2: 右，CCW，heading
        if loadM_data(n).cellID == 1124
            group_result(n).global_pref_code = [2 1 nan]; % 该文件ao数据丢失，但是仍有tuning
        else
            group_result(n).global_pref_code = nan(1,3);
        end
        
        group_result(n).global_pref_p = nan(1,3);
        
        % raw data
        group_result(n).spit_resp = nan(2,8);
        group_result(n).resp_norm = nan(2,8);
        group_result(n).spit_resp_trial = {nan(1,8);nan(1,8)};
        
%         % spike count mean and variance across cell
%         group_result(n).mean_fr = nan(1,2);
%         group_result(n).std_fr = nan(1,2);
        
        % spike count mean and variance across trial
        group_result(n).mean_fr_trial = nan(6,2);
        group_result(n).std_fr_trial = nan(6,2);
        
        % other
        % around forward motion, % first is fine, second is coarse
        group_result(n).pref_direc_aroundF = nan(1,2);
        group_result(n).pref_amp_aroundF = nan(1,2);
        group_result(n).p_aroundF = nan(1,2);
        group_result(n).pref_axis_num = nan(1,2);
        group_result(n).pref_axis = [{nan} {nan}];
        
        % 提取forward motion附近的点的RF raw data，R-CW-L-CCW
        group_result(n).spit_resp_0_fine = nan(8,4);
        group_result(n).spit_resp_0_coarse = nan(8,4);
        
        group_result(n).d_prime_ratio = nan;
        group_result(n).fr_axis_ratio = nan;
        group_result(n).vec_axis_ratio = nan;
        
        group_result(n).fittingresult = [];
    end
    
    
    % find if do AfterStim tuning test Lwh 20211231
    select_AS_tuning_data = loadM_data(n).Data(main_depth_loc).SpikeChan(1).AfterStim_Tuning;
    if size(select_AS_tuning_data,2) && isfield(select_AS_tuning_data,'FILE') % 电刺激前做在centerEye_Protocol下做了spiral tuning, 有FILE字段且不为空值
        group_result(n).AS_SpiT_FILE = select_AS_tuning_data.FILE;
        group_result(n).AS_pref_direc = squeeze(select_AS_tuning_data.pref_direc)';
        group_result(n).AS_p_stim = select_AS_tuning_data.p_sti_min;  % min
        p_turn_temp = select_AS_tuning_data.p_value; % p_turn(c,k,r,s)
        group_result(n).AS_p_turn_h = p_turn_temp(1,1,2,1); % 排除exp/con
        group_result(n).AS_p_turn_r = p_turn_temp(1,1,2,2); % 排除exp/con
        group_result(n).AS_p_turn_spi = p_turn_temp(1,1,1,3); % 包括exp/con
        group_result(n).AS_spiral_index = select_AS_tuning_data.spiral_index;
        group_result(n).AS_spit_resp = cell2mat(select_AS_tuning_data.resp)';
        group_result(n).AS_resp_norm = cell2mat(select_AS_tuning_data.resp_norm)';
        group_result(n).AS_spit_resp_trial = squeeze(select_AS_tuning_data.resp_trial);
    else
        group_result(n).AS_SpiT_FILE = nan;
        group_result(n).AS_pref_direc = nan(1,3);
        group_result(n).AS_p_stim = nan;  % min
        group_result(n).AS_p_turn_h = nan;
        group_result(n).AS_p_turn_r =  nan;
        group_result(n).AS_p_turn_spi =  nan;
        group_result(n).AS_spiral_index =  nan;
        group_result(n).AS_spit_resp = nan(2,8);
        group_result(n).AS_resp_norm = nan(2,8);
        group_result(n).AS_spit_resp_trial = {nan(1,8);nan(1,8)};
    end

end


% Add flexible monkey mask here (but I've still decided to choose monkey for analysis below). HH20150723
monkey_included_for_loading = [7 16];



%% cluster related
kk = 1;
for n = 1:length(loadM_data)
    for m = 1:length(loadM_data(n).Data) % 不同深度
        % 1. 提取vary eye position数据,同一深度下同时做了vary eye position和 center eye
        if isfield(loadM_data(n).Data(m).SpikeChan(1).varyEye_Protocol{1},'FILE') && size(loadM_data(n).Data(m).SpikeChan(1).varyEye_Protocol{1},2) && ...
                isfield(loadM_data(n).Data(m).SpikeChan(1).centerEye_Protocol{1},'FILE') && size(loadM_data(n).Data(m).SpikeChan(1).centerEye_Protocol{1},2) % 相同深度下有没有同时做eye center和vary eye 的
            group_result(n).pref_direc_centerEye(m,:) = squeeze([loadM_data(n).Data(m).SpikeChan(1).centerEye_Protocol{1}.pref_direc])';
            group_result(n).pref_direc_varyEye(m,:) = squeeze([loadM_data(n).Data(m).SpikeChan(1).varyEye_Protocol{1}.pref_direc])';
            
            % save for analysis
            pref_direc_centerEye(kk,:) = squeeze([loadM_data(n).Data(m).SpikeChan(1).centerEye_Protocol{1}.pref_direc])';
            pref_direc_varyEye(kk,:) = squeeze([loadM_data(n).Data(m).SpikeChan(1).varyEye_Protocol{1}.pref_direc])';
            RF_varyEye(kk,:) = loadM_data(n).RF;
            kk = kk + 1;
        end
    end
    
    % 2. cluster信息，提取不同深度下tuning resp, preferred direction
    Depth_code_num = length([loadM_data(n).Data.Depth_code]);
    if Depth_code_num > 1 % 在两个深度做了task
        Depth_code_SpiT = ones(Depth_code_num,1);
        % 判断是否测了tuning
        for m = 1:Depth_code_num
            file_name = [];
            try
                file_name = loadM_data(n).Data(m).SpikeChan(1).centerEye_Protocol{1}.FILE;
            catch
                Depth_code_SpiT(m) = 0;
            end
        end
        
        if sum(Depth_code_SpiT)>1 % 至少在两个不同的深度做了SpiT
            select_depth = find(Depth_code_SpiT);
            cluster_filename = [];
            for mm = 1:sum(Depth_code_SpiT)
                group_result(n).cluster(mm).depth = loadM_data(n).Data(select_depth(mm)).Depth;
                group_result(n).cluster(mm).depth_code = loadM_data(n).Data(select_depth(mm)).Depth_code;
                
                temp_resp = loadM_data(n).Data(select_depth(mm)).SpikeChan(1).centerEye_Protocol{1}.resp;
                group_result(n).cluster(mm).resp = cell2mat(temp_resp);
                
                temp_spiral_index = loadM_data(n).Data(select_depth(mm)).SpikeChan(1).centerEye_Protocol{1}.spiral_index;
                group_result(n).cluster(mm).spiral_index = temp_spiral_index;
                
                temp_pref_direc = loadM_data(n).Data(select_depth(mm)).SpikeChan(1).centerEye_Protocol{1}.pref_direc;
                group_result(n).cluster(mm).pref_direc = squeeze(temp_pref_direc)';
                
                cluster_filename = [cluster_filename,{loadM_data(n).Data(select_depth(mm)).SpikeChan(1).centerEye_Protocol{1}.FILE}]; 
            end
            group_result(n).cluster_FILE = cluster_filename;   
        else
            group_result(n).cluster(1).depth = nan;
            group_result(n).cluster(1).depth_code = nan;
            group_result(n).cluster(1).resp = nan;
            group_result(n).cluster(1).spiral_index = nan;
            group_result(n).cluster(1).pref_direc = [nan nan nan];
            group_result(n).cluster_FILE = nan;
        end
    else
        group_result(n).cluster(1).depth = nan;
        group_result(n).cluster(1).depth_code = nan;
        group_result(n).cluster(1).resp = nan;
        group_result(n).cluster(1).spiral_index = nan;
        group_result(n).cluster(1).pref_direc = [nan nan nan];
        group_result(n).cluster_FILE = nan;
    end
end

% 1. prefered direction在同一个cluster中不同深度时的standard deviation，可用作比较电刺激效果
% 2. 直接用该位点相邻100um的cluster index平均值作为指标
                
% 深度相隔100，200，300，400um， 允许误差20um
depth100 = [100:100:400];
tol = 25;
for n = 1:length(loadM_data)
    if ~isnan(group_result(n).cluster(1).depth)
        group_result(n).depth_all = [group_result(n).cluster.depth];
        
        if length(group_result(n).depth_all)>1 
            depth_two = nchoosek(group_result(n).depth_all,2); % [1 2 3]两两组合,例如[1 2],[1 3],[2 3]
            depth_diff = depth_two(:,2) - depth_two(:,1);
            
            % prepare
            depth_diff_plot_temp = arrayfun(@(x) nan(1,2), 1:4,'UniformOutput',0);
            
            % 去掉不合要求的深度, 例如相隔100时只留下[1 2],[2 3]，相隔200时只留下[1 3]，相隔300时为[nan nan];
            for dep = 1:4 % 100 200 300 400
                for kk = 1:length(depth_diff)
                    if (depth_diff(kk) < depth100(dep) + tol && depth_diff(kk) > depth100(dep) - tol)
                        if isnan(depth_diff_plot_temp{dep}(1))
                            depth_diff_plot_temp{dep} = []; % 去掉默认的nan
                        end
                        depth_diff_plot_temp{dep} = [depth_diff_plot_temp{dep};depth_two(kk,:)];
                    end
                end
            end
            group_result(n).depth_diff_plot = depth_diff_plot_temp; % use for calculate cluster correlation in @f2p9
            
            if ~isnan(group_result(n).depth_diff_plot{1}(1)) % 如果连相隔100时都为[]，那就不用继续运行了, 可能做了两个深度，但是相距不足100±25um
                % 计算刺激位点上下100um深度（tol=25）下的prefered direction的stardand deviation，DeAngelis and Newsome, 2004, 文章计算的时400um内的PD的SD
                % group_result(n).cluster_FILE
                
                % fine the stim depth and ±100um
                depth_code = [group_result(n).cluster.depth_code];
                depth =  [group_result(n).cluster.depth];
                stim_depth = depth(depth_code == 0);
                select_stim = logical(depth_code == 0);
                % select_depth = logical(depth_code == -1 | depth_code == 0 | depth_code == 1);
                select_depth = logical(abs(depth - stim_depth) < 100 + tol) ; % 直接计算深度比depth_code准确点
                
                % 1. PD and PD std
                % PD
                depth_pd = cell2mat({group_result(n).cluster.pref_direc}');
                select_pd_temp = depth_pd(select_depth',:);
                
                % PD std
                select_pd = HeadingToazi(select_pd_temp); % in degree, range from 0 to 2pi
                [~, pd_std_temp] = circ_std(deg2rad(select_pd)); % circular standard deviation, in rad
                group_result(n).pd_std = rad2deg(pd_std_temp); % in degree
                group_result(n).pd_std2 = std(select_pd); % just simply std, in degree
                
                % 计算刺激位点上下100um深度（tol=25）下的平均pearson correlation
                % Gu: cluster index可定义为每个位点的tuning和临近两个位点的pearson correlation 的平均值（上下各一个位点）
                stim_resp = group_result(n).cluster(select_stim).resp; % stim depth resp
                
                % up 100 and down 100 um
                select_up = logical( (depth-stim_depth)<-tol & (depth-stim_depth)>-100-tol); % -125~-25 um
                select_down = logical( (depth-stim_depth)< 100+tol & (depth-stim_depth)>tol); % 25~125 um
                if sum(select_up)>0
                    up_resp = group_result(n).cluster(select_up).resp; % stim depth resp
                else
                    up_resp = nan(size(stim_resp));
                end
                
                if sum(select_down)>0
                    down_resp = group_result(n).cluster(select_down).resp; % stim depth resp
                else
                    down_resp = nan(size(stim_resp));
                end
                

                % 2. average of corr coe for ±100um
                cluster_index_temp = [];
                for TR = 1:2
                    [cluster_index_temp(1,TR), ~] = corr(stim_resp(:,TR),up_resp(:,TR),'type','Pearson');
                    [cluster_index_temp(2,TR), ~] = corr(stim_resp(:,TR),down_resp(:,TR),'type','Pearson');
                end
                group_result(n).cluster_index = nanmean(cluster_index_temp,1); % T and R
                
            else
                group_result(n).depth_all = nan;
                group_result(n).depth_diff_plot = arrayfun(@(x) nan(1,2), 1:4,'UniformOutput',0);
                group_result(n).pd_std = nan(1,3);
                group_result(n).pd_std2 = nan(1,3);
                group_result(n).cluster_index = nan(1,2);
            end
        end
    else
        group_result(n).depth_all = nan;
        group_result(n).depth_diff_plot = arrayfun(@(x) nan(1,2), 1:4,'UniformOutput',0);
        group_result(n).pd_std = nan(1,3);
        group_result(n).pd_std2 = nan(1,3);
        group_result(n).cluster_index = nan(1,2);
    end
end

%% cell classification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 此步骤在Plot_Spiral_Tuning_Lwh2.m 中也有，存于select_tuning_data中，不错默认了p criteria=0.001，
% 此时直接从Group_GUI设定 p criteria
[cell_type_all,mt_type_all,~] = cell_type_classification(group_result,[],p_cri);

% save back to group_result
for n = 1:length(loadM_data)
    group_result(n).cell_type = cell_type_all(n);
    group_result(n).mt_type = mt_type_all(n);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Cell Selection and Cell Counter
select_all = [];
select_typical = [];
select_no_typical = [];
loop_num = [];
selection_num = [];
monkey_included_for_analysis = [];
cell_selection();

    function cell_selection(typical_cell_selection_num)  % Cell Selection and Cell Counter
        if nargin < 1
            typical_cell_selection_num = 2; % Default,select all cell type
        end
        
        if typical_cell_selection_num == 2 % select all cell type (T+R+S)
            loop_num = [1:3];
        elseif typical_cell_selection_num == 3 || typical_cell_selection_num == 4 || typical_cell_selection_num == 5 % T,R,Spi cell
            loop_num = typical_cell_selection_num - 2; % 后面用作 = cell_type
        elseif typical_cell_selection_num == 1 % select all cell
            loop_num = 0; 
        end
        
        selection_num = typical_cell_selection_num;
        
        area = {loadM_data.Area}';
        %         select_all_all_monkey =  cellfun(@(x) strcmp(x,'MST'),area) & [group_result.SpiT_repN]' >= 3; % include area "MST"
        select_all_all_monkey =  cellfun(@(x) strncmp(x,'MST',3),area) & [group_result.SpiT_repN]' >= 3; % include area "MST" and "MST?"
        
        typical_cell_selection_criteria = ...
            {%  Logic                                   Notes
            % Bottom-line
            'Bottom-line (all)', select_all_all_monkey; % Just all bottom-line cells
            % Cell type base
            'T+R+Spi cell', (([group_result.cell_type] == 1)' | ([group_result.cell_type] == 2)' | ([group_result.cell_type] == 3)' ) & select_all_all_monkey;
            'Translation cell', ([group_result.cell_type] == 1)' & select_all_all_monkey;
            'Rotation cell', ([group_result.cell_type] == 2)' & select_all_all_monkey;
            'Spiral cell', ([group_result.cell_type] == 3)' & select_all_all_monkey;
            'Motion type cell', ([group_result.mt_type] == 1)' & select_all_all_monkey;
            };
        
        select_typical_cells_all_monkey = select_all_all_monkey & typical_cell_selection_criteria{typical_cell_selection_num,2};
        select_no_typical_cells_all_monkey = select_all_all_monkey & ~ typical_cell_selection_criteria{typical_cell_selection_num,2};
        t_criterion_txt = typical_cell_selection_criteria{typical_cell_selection_num,1};
        
        % -------- Count cell numbers for each monkey. HH20150723 --------
        n_monkey = length(monkey_included_for_loading);
        cell_nums = zeros(n_monkey + 1,2); % All units / Typical Cells
        for mm = 1:n_monkey
            select_monkey{mm} = [loadM_data.monkey]' == monkey_included_for_loading(mm);
            cell_nums(mm,1) = sum(select_all_all_monkey & select_monkey{mm});
            cell_nums(mm,2) = sum(select_typical_cells_all_monkey & select_monkey{mm});
        end
        
        % -------- Update actual dataset for analysis. HH20150723 --------
        % 根据sheet选择强制限定monkey（只选择一个sheet/一个monkey的时候） Lwh 20210601
        monkey_choose = str2num(get(handles.sheetN,'string'));
        if length(monkey_choose)==1
            monkey_included_for_analysis = monkey_included_for_loading(monkey_choose);
            
            % update check box, 默认都是勾选=1
            if monkey_choose == 1
                set(handles.Ringbell_data,'value',1);
                set(handles.Arthas_data,'value',0);
            else
                set(handles.Ringbell_data,'value',0);
                set(handles.Arthas_data,'value',1);
            end
        else
            monkey_included_for_analysis = monkey_included_for_loading(logical([get(findall(gcbf,'tag','Ringbell_data'),'value') get(findall(gcbf,'tag','Arthas_data'),'value')]));
        end
        
        monkey_mask_for_analysis = false(length(group_result),1);
        for mm = 1:length(monkey_included_for_analysis)
            monkey_mask_for_analysis = monkey_mask_for_analysis | ([loadM_data.monkey]' == monkey_included_for_analysis(mm));
        end
        
        % -------- Affect all analysis below --------
        select_all = select_all_all_monkey & monkey_mask_for_analysis;
        select_typical = select_typical_cells_all_monkey & monkey_mask_for_analysis;
        select_no_typical = select_no_typical_cells_all_monkey & monkey_mask_for_analysis;
        
        cell_nums(end,:) = [sum(select_all) sum(select_typical)];
        
        % -------- Update cell counter ---------
        h_all = findall(gcbf,'tag','num_all_units');
        set(h_all,'string',sprintf('%10d%13d\n',cell_nums'),'fontsize',13);
        h_t_criterion = findall(gcbf,'tag','t_criterion');
        set(h_t_criterion,'string',{typical_cell_selection_criteria{:,1}});
        set(h_t_criterion,'value',typical_cell_selection_num);
        
        if sum(select_all)~=0
            [~,~,cell_type_num] = cell_type_classification(group_result,select_all,p_cri);
        else
            cell_type_num = [];
        end
        
    end

% -------- Update/refresh some related datasets that are influenced by cell_selection ---------

%% Final Preparation

% =========== Data for common use ============
function_handles_spi = [];

% ================ Miscellaneous ===================
set(0,'defaultAxesColorOrder',[41 89 204; 248 28 83; 14 153 46]/255);
c{1} = [248 28 83]/255; % 红
c{2} = [41 89 204]/255; % 蓝
c{3} = [14 153 46]/255; % 绿

c{4}='c'; % 青色
c{5}='w';
c{11} = [0.8 0.5 0.5]; % 浅红色
c{12} = [0.5 0.5 0.8]; % 浅蓝色
c{13} = [0.5 0.8 0.5]; % 浅绿色

c{98} = [100 70 140]/255; % 紫色
c{99} = [0.6 0.6 0.6]; % 灰色，数字越小颜色越深

transparent = get(findall(gcbf,'tag','transparent'),'value');
p_critical = 0.05;
figN = 100;

%% ====================================== Function Handles =============================================%
% Today I reconstruct the group analysis codes in a more logical and effcient way:
% 1. I divided each figure into small nested functions in Group_HD
% 2. Group_HD returns a cell called "function_handles" back to GROUP_GUI (with all the data available for each nested function)
% 3. Now I can choose which figure(s) to plot or which function to debug in GROUP_GUI
% 4. I feel so happy.
% @HH20150425

function_handles_spi = {
    'Basic properties',{
    'Classification',@f1p1;
    'Response Strength (Mutual Information)',@f1p2;
    'RF distribution of MST',@f1p3
    'RF across site location (Grid-xy)',@f1p4
    'Compare Max FiringRate with Gaussian Fit maxFR',@f1p5
    'Response reliability in T and R',@f1p6
    'Max Response in T and R',@f1p7
    };
    'Spiral Tuning',{
    'DDI comparison among cell type',@f2p1;
    'D prime comparison among cell type',@f2p2;
    'D prime comparison (±45 and ±90)',@f2p7;
    'Preferred Direction',@f2p3;
    'PD vs DDI/d prime',@f2p23;
    'Spiral space preference',@f2p4;
    'Translation, rotation prefererence correlation with cell type',@f2p5;
    'Translation, rotation prefererence correlation with RF location',@f2p12;
    'Center Eye position and Vary Eye position',@f2p6;
    'RF effect in preferred direction',@f2p8
    'Calculate Prefer OpticFlow (local planar motion)',@f2p15
    'Save Prefer OpticFlow similarity (local planar motion)',@f2p19
    'Prefer OpticFlow similarity vs RF',@f2p20
    'Clustering (correlation and PD SD)',@f2p9 % cluster
    'Clustering for different similarity sites', @f2p24 % cluster
    'Preferred Direction around forward motion(compare to CueS task)',@f2p10
    'Tuning Index',@f2p11
    'Gaussian Fit Result in T and R',@f2p16
    'Gaussian Bandwidth in T and R',@f2p17
    'Bandwidth and DDI',@f2p18
    'AfterStim Spiral Tuning',@f2p21
    'ROC analysis and compare to d prime',@f2p22
    };
    'Spiral Space Tuning (Motion type tuning)',{
    'Spiral Index',@f3p1
    'Some ratio for T-R',@f3p2
    '',@f3p3
    };
    'PCA for spiral tuning',{
    'PCA',@f4p1
    };
    
    
    'NoShow',{@cell_selection};
    };

%% ====================================== Function Definitions =============================================%
    function f1p1(debug)      % Cell Classification
        %  饼图比较Translation cell和Rotationg cell和spiral cell and no response cell (4类)
        % 【2】根据p_value, spiral_index 进行分类
        % 此步在Plot_Spiral_Tuning_Lwh2.m 中已经进行(调用cell_type_classification.m)，可以直接从select_tuning_data中获得，无需重新获取
        
        if debug  ; dbstack;   keyboard;      end
        methods_of_select = {select_typical, 'Typical cells'};
        
        set(figure(figN),'name','Cell Classifcation','pos',[200,300 800,500]); clf
        
        heading_cell_num = sum([group_result(select_all).cell_type] == 1);
        rotation_cell_num = sum([group_result(select_all).cell_type] == 2);
        spiral_cell_num = sum([group_result(select_all).cell_type] == 3);
        
        % combine cell type 4 and 5: no tuning cell
        no_tuning_cell_num = sum([group_result(select_all).cell_type] == 4) + sum([group_result(select_all).cell_type] == 5);
        
        cell_all_num = heading_cell_num+rotation_cell_num+spiral_cell_num+no_tuning_cell_num;
        
        mt_tuning_num = sum([group_result(select_all).mt_type] == 1);
        all_active_num = heading_cell_num + rotation_cell_num + spiral_cell_num + no_tuning_cell_num;
        
%         cm = [c{1};c{2};c{3};c{99}/10;c{99}];
%         colormap(cm)
        pp = pie([heading_cell_num rotation_cell_num spiral_cell_num no_tuning_cell_num],[0 0 0 0]);
        legend('Translation cell','Rotation cell','Spiral cell','No tuning cell','location','bestoutside');
        
        axes('position',[0.67 0.05 0.3 0.5]);
        
        text(0,1,sprintf('Translation Cell = %g (%2.4g%%)',heading_cell_num, heading_cell_num/cell_all_num*100));
        text(0,0.8,sprintf('Rotation Cell = %g (%2.4g%%)',rotation_cell_num, rotation_cell_num/cell_all_num*100));
        text(0,0.6,sprintf('Spiral Cell = %g (%2.4g%%)',spiral_cell_num, spiral_cell_num/cell_all_num*100));
        text(0,0.4,sprintf('No tuning = %g (%2.4g%%)',no_tuning_cell_num, no_tuning_cell_num/cell_all_num*100));
        text(0,0.2,sprintf('N = %g',cell_all_num));
        axis off;
        
        SetFigure(15);figN = figN+1;
        
%         % motion type tuning: Excitory, 且heading与rotation plane有显著区别
%         axes('position',[0.1 0.03 0.3 0.05]);
%         x1 = 1;
%         x2 = all_active_num / mt_tuning_num;
%         rectangle('Position',[0 0 x1 1],'facecolor','k');
%         hold on
%         rectangle('Position',[x1 0 x2 1],'facecolor','w');
%         text((x1+x2)/2-2.7,1.7,sprintf('Motion type tuning cell / excitory cell = %.2f%%',mt_tuning_num / all_active_num*100));
% %         text((x1+x2)-6.2,0.6,'Excitory cell');
%         axis off

        SetFigure(10);figN = figN+1;
        
        %         % PCA
        %
        %         %%% group raw data
        %         T_response_all = [];
        %         response_all = [];
        %         response_trial = [];
        %         for n = 1:length(group_result)
        %             if methods_of_select{1}(n)
        %                 % mean fr
        %                 TR_component = repmat([-180:45:135]',2,1);
        %                 spit_resp_temp = group_result(n).spit_resp; % mean firing rate: -180,  -135,  -90,  -45,  0,  45,  90,  135
        %                 T_response = spit_resp_temp(1,:);
        %                 R_response = spit_resp_temp(2,:);
        %                 response_all_this = [T_response';R_response'];
        %                 %             response_all_this([3,13]) = []; % remove the repeat 0 and 180 in rotation
        %                 remove_this = logical(TR_component==0 | TR_component==-180);
        %                 %             response_all_this(remove_this) = []; % remove the repeat 0 and 180 in both translation and rotation
        %
        %                 response_all = [response_all,response_all_this];
        %             end
        %         end
        %
        %         % remove nan
        %         remove_this1 = any(isnan(response_all));
        %         response_all(:,remove_this1) = [];
        %
        %         % Normalize
        %         response_all_zscore = zscore(response_all,[],2); % normalize
        %         [row,col]=size(response_all_zscore);
        %
        %         % cell type
        %         cell_type = [group_result.cell_type]';
        %         cell_type = cell_type(methods_of_select{1});
        %
        %         % n * feature
        %         response_all_zscore = response_all_zscore';
        %         C = cov(response_all_zscore); % 直接调用cov直接计算协方差矩阵即可
        %         [W, Lambda] = eig(C); % W是特征向量组成的矩阵（4×4），Lambda是特征值组成的对角矩阵
        %         ev = (diag(Lambda))'; % 提取特征值
        %         ev = ev(:, end:-1:1); % eig计算出的特征值是升序的，这里手动倒序（W同理）
        %         W = W(:, end:-1:1);
        %         sum(W.*W, 1) % 可以验证每个特征向量各元素的平方和均为1
        %         Wr = W(:, 1:3); % 提取前3个主成分的特征向量
        %         Tr = response_all_zscore * Wr; % 新坐标空间的数据点
        %
        %         % 作图
        %         figure;
        %         stairs(cumsum(ev)/sum(ev), 'LineWidth',1.5);
        %         axis([1 4 0 1]);
        %         xlabel('$ k $', 'Interpreter', 'latex');
        %         ylabel('$ f(k)=\frac{\sum _{i=1}^i \lambda_k}{\sum_{i=1}^m \lambda_i} $','Interpreter', 'latex');
        %         hold on;
        %         plot([1 4], [0.95 0.95], '--'); % 从图中可以看出，取r = 3即可
        %
        %         figure;
        %         scatter3(Tr(:,1), Tr(:,2), Tr(:,3), 130, categorical(cell_type), '.');
        %         colormap(winter);
        %         xlabel('Principal Component 1');
        %         ylabel('Principal Component 2');
    end

    function f1p2(debug)  % response strength, mutual information
        % MI：mutual information
        % MI越大说明FR与当前刺激(heading or rotation)的depedence越高
        if debug  ; dbstack;   keyboard;      end
        methods_of_select = {
            % select_all, 'All cells';
            select_typical, 'Typical cells';
            % select_no_typical, 'Non-typical cells'
            };
        
        % Mutual Information
        MI = cell2mat({group_result(methods_of_select{1,1}).MI}');
        MI_max = max(MI(:));
        MI_min = min(MI(:));
        xbins = linspace(MI_min,MI_max,10);
        MI_tN = hist(MI(:,1),xbins);
        MI_rN = hist(MI(:,2),xbins);
        
        % normal MI
        normalMI = cell2mat({group_result(methods_of_select{1,1}).normalMI}');
        normMI_max = max(normalMI(:));
        normMI_min = min(normalMI(:));
        xbins2 = linspace(normMI_min,normMI_max,10);
        normMI_tN = hist(normalMI(:,1),xbins2);
        normMI_rN = hist(normalMI(:,2),xbins2);
        
        figure(figN);
        set(figN,'Position', [300,200 750,700], 'Name', 'response of MST for Translation and Rotation stimulus (MI)');
        subplot(2,2,1)
        h = bar(xbins,[MI_tN',MI_rN'],1,'grouped','facecolor','r');
        set(h(2),'facecolor','b');
        xlim([MI_min-MI_max/10 MI_max+MI_max/10]);
        title('Mutual Information distribution');
        
        % unit-by-unit MI
        subplot(2,2,2)
        plot(MI(:,1),MI(:,2),'ko');
        hold on
        plot([0 MI_max+MI_max/10],[0 MI_max+MI_max/10],'k:');
        xlim([0 MI_max+MI_max/10]);
        ylim([0 MI_max+MI_max/10]);
        title('unit-by-unit MI');
        xlabel('Translation MI');
        ylabel('Rotation MI');
        
        nan_this = any(isnan([MI(:,1),MI(:,2)]),2);
        
        [r,p] = corr(MI(~nan_this,1),MI(~nan_this,2));
        
        
        % normalize MI
        subplot(2,2,3)
        h = bar(xbins2,[normMI_tN',normMI_rN'],1,'grouped','facecolor','r');
        set(h(2),'facecolor','b');
        xlim([normMI_min-normMI_max/10 normMI_max+normMI_max/10]);
        title('normal MI distribution');
        
        % unit-by-unit-normalize MI
        subplot(2,2,4)
        plot(normalMI(:,1),normalMI(:,2),'ko');
        hold on
        plot([0 normMI_max+normMI_max/10],[0 normMI_max+normMI_max/10],'k:');
        xlim([0 1]);
        ylim([0 1]);
        %         xlim([normMI_min-normMI_max/10 normMI_max+normMI_max/10]);
        %         ylim([normMI_min-normMI_max/10 normMI_max+normMI_max/10]);
        %
        title('unit-by-unit normal MI');
        xlabel('Translation normal MI');
        ylabel('Rotation normal MI');
        
        SetFigure(15);figN = figN+1;
    end

    function f2p7(debug)  % d prime comparison  fine and coarse
        if debug  ; dbstack;   keyboard;      end
        figure(figN);
        set(figN,'Position', [300,0 600,1000], 'Name', 'd prime in ±45 and ±90');
        
        methods_of_select = {
            % select_all, 'All cells';
            select_typical, 'Typical cells';
            % select_no_typical, 'Non-typical cells'
            };
        
        d_prime_fine = cell2mat({group_result(methods_of_select{1,1}).d_prime_fine}');
        d_prime_coarse = cell2mat({group_result(methods_of_select{1,1}).d_prime_coarse}');
        neu_sen_fine = cell2mat({group_result(methods_of_select{1,1}).neu_sen_fine}');
        neu_sen_coarse = cell2mat({group_result(methods_of_select{1,1}).neu_sen_coarse}');
        d_prime_index_fine = cell2mat({group_result(methods_of_select{1,1}).d_prime_index_fine}');
        d_prime_index_coarse = cell2mat({group_result(methods_of_select{1,1}).d_prime_index_coarse}');
        
        % m7c1004r1 rotation tuning, CW CCW Firing Rate 完全一样，导致d'=0，
        % neu_sen_coarse为一个很小的值，去掉
        % 还发现别的file也存在此问题，neu_sen_coarse=-inf
        remove_file = [];
        for d = 1:2
            remove_file1 = logical(neu_sen_fine(:,d) < -15);
            remove_file2 = logical(neu_sen_coarse(:,d) < -15);
            remove_file{d} = logical(remove_file1 | remove_file2);
            neu_sen_fine(remove_file{d},d) = nan;
            neu_sen_coarse(remove_file{d},d) = nan;
        end
        
        % Plot d prime
        for d = 1:2 % h or r
            subplot(3,2,d)
            plot(d_prime_fine(:,d),d_prime_coarse(:,d),'o','color',c{d});
            xlimit = xlim;
            ylimit = ylim;
            xylimit = max(abs([xlimit ylimit]));
            xlim([-xylimit xylimit]);
            ylim([-xylimit xylimit]);
            xlabel('d prime(fine)');
            ylabel('d prime(coarse)');
            
            hold on
            plot([-xylimit xylimit],[-xylimit xylimit],'k:');
            [R,p_cor] = corrcoef(d_prime_fine(:,d),d_prime_coarse(:,d));
            R = R(1,2);
            p_cor = p_cor(1,2);
            text(-xylimit+xylimit/10,xylimit-xylimit/10,sprintf('R = %2g',R));
            text(-xylimit+xylimit/10,xylimit-xylimit/5,sprintf('p = %2g',p_cor));
            
            if d == 1
                title('Translation tuning');
            else
                title('Rotation tuning');
            end
        end
        
        % plot neural sensitivity
        for d = 1:2 % h o r
            subplot(3,2,d+2)
            plot(neu_sen_fine(:,d),neu_sen_coarse(:,d),'o','color',c{d});
            xlimit = xlim;
            ylimit = ylim;
            xylimit = max(abs([xlimit ylimit]));
            xlim([-xylimit xylimit]);
            ylim([-xylimit xylimit]);
            xlabel('neuron sensitivity(fine)');
            ylabel('neuron sensitivity(coarse)');
            
            hold on
            plot([-xylimit xylimit],[-xylimit xylimit],'k:');
            [R,p_cor] = corrcoef(neu_sen_fine(:,d),neu_sen_coarse(:,d),'Rows','complete'); % 排除nan
            R = R(1,2);
            p_cor = p_cor(1,2);
            text(-xylimit+xylimit/10,xylimit-xylimit/10,sprintf('R = %2g',R));
            text(-xylimit+xylimit/10,xylimit-xylimit/5,sprintf('p = %2g',p_cor));
        end
        
        % 比较fine和coarse index之间的关系
        subplot(3,2,5)
        plot(d_prime_index_fine,d_prime_index_coarse,'ko');
        xlim([-1 1]);
        ylim([-1 1]);
        xlabel('d prime index (Fine)');
        ylabel('d prime index (Coarse)');
        hold on
        plot([-1 1],[-1 1],'k:');
        plot([0 0],[-1 1],'k:');
        plot([-1 1],[0 0],'k:');
        [R,p_cor] = corr(d_prime_index_fine,d_prime_index_coarse,'type','Spearman');
        text(-0.9,1.2,sprintf('R = %2g',R));
        text(-0.9,1.1,sprintf('p = %2g',p_cor));
        %         text(0.2,0.5,'pref Trans','fontsize',15,'color','r','fontweight','bold');
        %         text(-0.8,-0.5,'pref Rot','fontsize',15,'color','b','fontweight','bold');
        
        % 单独看看fine and coarse task中d prime index分布
        d_prime_index(:,1) = d_prime_index_fine;
        d_prime_index(:,2) = d_prime_index_coarse;
        xbins = linspace(-1,1,10);
        for d = 1:2 % fine or coarse
            subplot(6,2,(d-1)*2+10)
            nbins = hist(d_prime_index(:,d),xbins);
            bar(xbins,nbins,1);
            top = max(nbins);
            xlim([-1.2 1.2]);
            ylim([0 top+top/4]);
            hold on
            mean_d = nanmean(d_prime_index(:,d));
            [h, p_t] = ttest(d_prime_index(:,d));
            plot(mean_d,top+top/5,'kv');
            text(0.4,top,sprintf('p = %2g',p_t));
            
            if d==1
                title('Fine task');
            else
                title('Coarse task');
            end
        end
        
        SetFigure(10);figN = figN+1;
        
        subplot(3,2,5)
        text(0.2,0.5,'pref Trans','fontsize',15,'color','r','fontweight','bold');
        text(-0.8,-0.5,'pref Rot','fontsize',15,'color','b','fontweight','bold');
    end

    function f2p1(debug)  % DDI compasison
        if debug  ; dbstack;   keyboard;      end
        methods_of_select = {select_typical, 'Typical cells' };
        
        DDI{1} = cell2mat({group_result.DDI_includeExpCon}');
        DDI{2} = cell2mat({group_result.DDI_excludeExpCon}');
        
        DDI{1} = DDI{1}(methods_of_select{1},:);
        DDI{2} = DDI{2}(methods_of_select{1},:);
        
        % 带Histogram的散点图
        for r = 1:2 % 包括or排除exp/con的
            if r == 1
                set(figure(figN),'Position', [200,100 1200,500], 'Name', 'DDI (included Exp/Con)');
            else
                set(figure(figN),'Position', [200,100 1200,500], 'Name', 'DDI (excluded Exp/Con)');
            end
            
            h1 = LinearCorrelation({
                DDI{r}([group_result(select_typical).cell_type]' == 1,1);
                DDI{r}([group_result(select_typical).cell_type]' == 2,1);
                DDI{r}([group_result(select_typical).cell_type]' == 3,1);
                },...
                {
                DDI{r}([group_result(select_typical).cell_type]' == 1,2);
                DDI{r}([group_result(select_typical).cell_type]' == 2,2);
                DDI{r}([group_result(select_typical).cell_type]' == 3,2);
                },...
                'CombinedIndex',[7],'PlotCombinedOnly',1,...
                'Xlabel','DDI for Translation tuning','Ylabel','DDI for Rotation tuning',...
                'FaceColors',{c{1},c{2},c{3}},'Markers',{'o','o','o'},...
                'LineStyles',{'k:','k:','k:','k-'},'MarkerSize',10,...
                'figN',figN,'XHist',11,'YHist',11,...
                'XHistStyle','grouped','YHistStyle','grouped','SameScale',1,...
                'Method','Pearson','FittingMethod',2,...
                'dot_leg',{'T cell', 'R cell',' Spi cell'}); figN = figN + 1;
        end

        % 不带Histogram的散点图 + DDI_index    
        for r = 1:2 % 包括or排除exp/con的
            if r == 1
                set(figure(figN),'Position', [200,100 1200,500], 'Name', 'DDI (included Exp/Con)');
            else
                set(figure(figN),'Position', [200,100 1200,500], 'Name', 'DDI (excluded Exp/Con)');
            end
            
            % DDI scatter plot
            subplot(1,2,1)
            for d = 1:3 % h r spi cell type
                hold on
                h(d) = plot(DDI{r}([group_result(select_typical).cell_type]' == d,1), DDI{r}([group_result(select_typical).cell_type]' == d,2),'o','color',c{d});
            end
            xlim([0 1]);
            ylim([0 1]);
            plot([0 1],[0 1],'k:');
            plot([0 1],[0.5 0.5],'k:');
            plot([0.5 0.5],[0 1],'k:');
            
            [r,pp] = plot_corr_line(DDI{r}(:,1), DDI{r}(:,2),'MethodOfCorr','Pearson','FittingMethod',2,'LineStyle','k-');
            text(0.05,0.6,sprintf('r = %3.3g',r));
            text(0.05,0.55,sprintf('p = %3.3g',pp));
            
            xlabel('DDI for translation tuning');
            ylabel('DDI for rotation tuning');
            hl = legend(h([1 2 3]),'Heading cell','Rotation cell','Spiral cell','location','northwest');
            box off;
            if r==1
                title('DDI included Exp/Con');
            else
                title('DDI excluded Exp/Con');
            end
            
            % mark the largest DDI cell
            [~,largerone] = max(DDI{r}(:,1).*DDI{r}(:,2));
            plot(DDI{r}(largerone,1),DDI{r}(largerone,2),'*','color','k');
            
            % DDI index
            DDI_index = (DDI{r}(:,1) - DDI{r}(:,2))./ (DDI{r}(:,1) + DDI{r}(:,2));
            xbins = -0.5:0.1:0.5;
            subplot(1,2,2)
            xbins_c = [];
            xbins_c{1} = linspace(-0.25,0.25,12);
            xbins_c{2} = linspace(-0.25,0.25,12);
            xbins_c{3} = linspace(-0.25,0.25,12);
            
            hold on
            for d = 1:3 % 先画spiral
                i_temp = 4-d;
                nbins_type = hist(DDI_index([group_result(select_typical).cell_type]' == i_temp),xbins_c{i_temp});
                nbins_type = log10(nbins_type);
                bar(xbins_c{i_temp},nbins_type,1,'facecolor',c{i_temp});
                alpha(0.5);
                
                mean_DDI_index = nanmean(DDI_index([group_result(select_typical).cell_type]' == i_temp));
                % 与0比较
                if d == 1 % spiral
                    [~,p] = ttest(DDI_index([group_result(select_typical).cell_type]' == i_temp),0);
                elseif d == 2 % rotation
                    [~,p] = ttest(DDI_index([group_result(select_typical).cell_type]' == i_temp),0,'tail','left'); % 单边检验，小于0
                elseif  d==3 % heading cell
                    [~,p] = ttest(DDI_index([group_result(select_typical).cell_type]' == i_temp),0,'tail','right'); % 单边检验，大于0
                end
                xlim([-0.4 0.4]);
                plot(mean_DDI_index,1.6,'v','color',c{i_temp});
                text(-0.3,0.8+d*0.15,sprintf('p = %2.2g',p),'color',c{i_temp});
            end
            title('DDI Index');
            SetFigure(15);
            figN = figN+1;
        end
    end

    function f2p2(debug)  % D prime comparison among cell type, d prime index
        if debug  ; dbstack;   keyboard;      end
        set(figure(figN),'Position', [200,100 800,800], 'Name', 'D prime');clf;
        
        d_prime_fine = cell2mat({group_result.d_prime_fine}');
        d_prime_coarse = cell2mat({group_result.d_prime_coarse}');
        neu_sen_fine = cell2mat({group_result.neu_sen_fine}');
        neu_sen_coarse = cell2mat({group_result.neu_sen_coarse}');
        d_prime_index_fine = cell2mat({group_result.d_prime_index_fine}');
        d_prime_index_coarse = cell2mat({group_result.d_prime_index_coarse}');
        
        % d prime越大或越小，说明细胞越能区分左右（CW,CCW），所以只需要考虑绝对值的大小
        d_prime_abs{1} = (abs(d_prime_fine));
        d_prime_abs{2} = (abs(d_prime_coarse));
        
%         % ROC
%         AuROC = cell2mat({group_result.AuROC}');
%         d_prime_abs{2} = AuROC;

        methods_of_select = {
            % select_all, 'All cells';
            select_typical, 'Typical cells';
            % select_no_typical, 'Non-typical cells'
            };
        
        islogscale = 1;
        for j = 1:2 % fine or coarse
            
            % 1. d prime
            for d = loop_num % h r spi cell type
                subplot(2,2,j)
                hold on
                x_temp = []; y_temp = [];norm = [];norm_temp = [];
                
                if length(loop_num)==3
                    x_temp = d_prime_abs{j}(([group_result.cell_type]' == d & methods_of_select{1,1}),1);
                    y_temp = d_prime_abs{j}(([group_result.cell_type]' == d & methods_of_select{1,1}),2);
                else
                    x_temp = d_prime_abs{j}(methods_of_select{1,1},1)';
                    y_temp = d_prime_abs{j}(methods_of_select{1,1},2)';
                end
                
                x = x_temp;
                y = y_temp;
                
                h(d) = plot(x, y,'o','color',c{d},'markersize',8);
                if islogscale
                    set(gca,'xscale','log','yscale','log'); % log scale
                end
                
                %                 % fitting separately
                %                 [r,ppp] = plot_corr_line(x,y,'MethodOfCorr','pearson','FittingMethod',2,'LineStyles','k-','LineWidth',1);
                
                xlabel('Translation tuning');
                ylabel('Rotation tuning');
                
                if j==1
                    title('Fine |d prime| (±45°)');
                else
                    title('Coarse |d prime| (±90°)');
                end
            end

            % fitting and correlation all cell type
            select_all_celltype = logical([group_result.cell_type]' == 1 | [group_result.cell_type]' == 2 | [group_result.cell_type]' == 3);
            xx = d_prime_abs{j}((select_all_celltype & methods_of_select{1,1}),1);
            yy = d_prime_abs{j}((select_all_celltype & methods_of_select{1,1}),2);
            
            if islogscale
                % correlation
                % xx,yy中任意一项有NAN则去除
                noxx = isnan(xx);
                noyy = isnan(yy);
                remove_xy = logical(noxx | noyy);
                xx = xx(~remove_xy);
                yy = yy(~remove_xy);
                [rrr,ppp] = corr(xx,yy,'type','pearson');
                r = rrr^2;
                
                % PCA fitting in log scale
                % remove Inf/-Inf if have
                xxlog = log10(xx);
                yylog = log10(yy);
                noxx2 = logical(xxlog == inf | xxlog == -inf);
                noyy2 = logical(yylog == inf | yylog == -inf);
                remove_xy2 = logical(noxx2 | noyy2);
                xxlog = xxlog(~remove_xy2);
                yylog = yylog(~remove_xy2);

                fitType1 = regress_perp(xxlog,yylog,1,0.05);
                linPara(1) = fitType1.k;
                linPara(2) = fitType1.b;
                linParaSE = [std(fitType1.bootK) std(fitType1.bootB)];
                xxx = linspace(min(xxlog),max(xxlog),150);
                Y = linPara(1) * xxx + linPara(2);
                xxx = xxx(min(yylog) <= Y & Y <= max(yylog));
                Y = Y(min(yylog) <= Y & Y <= max(yylog));
                if ppp<0.05
                    h = plot(10.^xxx,10.^Y,'k-','linewidth',1);
                end
                xlim([0.001 100]);
                ylim([0.001 100]);
                
                % ttest
                [~,p_ttest] = ttest(xxlog,yylog);
                text(0.01,0.01,sprintf('t test, p = %.2g',p_ttest));
            else
                [r,ppp] = plot_corr_line(xx,yy,'MethodOfCorr','pearson','FittingMethod',2,'LineStyles','k-','LineWidth',1);
            end
            axis square
            xlimit = xlim;
            ylimit = ylim;
            str = [sprintf('r = %.3g', r) newline sprintf('p = %.3g', ppp)];
            text(xlimit(1),ylimit(2),str);
            
            
            
%             % Show individual cell selected from the figure. HH20150424
%             if length(loop_num)==3
%                 select_actual_plot = logical(([group_result.cell_type]'==1 | [group_result.cell_type]'==2 | [group_result.cell_type]' == 3) & select_all);
%             else
%                 select_actual_plot = logical([group_result.cell_type]'==loop_num & select_all);
%             end
%             
%             x_temp = d_prime_abs{j}(select_actual_plot,1);
%             y_temp = d_prime_abs{j}(select_actual_plot,2);
%             
%             h_all = plot(x_temp,y_temp,'visible','off');hold on;
%             set([gca h_all],'ButtonDownFcn',{@Show_individual_cell, h_all, select_actual_plot});
            
            % 2. d prime Index
            d_prime_index = (d_prime_abs{j}(:,1) - d_prime_abs{j}(:,2)) ./ (d_prime_abs{j}(:,1) + d_prime_abs{j}(:,2));
            
            % bar plot
            xbins_c = [];
            xbins_c{1} = -1:0.2:1;
            xbins_c{2} = -1:0.2:1;
            xbins_c{3} = -1:0.2:1;
            
            for d = loop_num % h cell, r cell, spiral cell
                hold on
                temp = [];
                temp = d_prime_index(([group_result.cell_type]' == d & methods_of_select{1,1}));
                nbins_type = hist((temp),xbins_c{d});
                subplot(2,2,2+j)
                bar(xbins_c{d},nbins_type,1,'facecolor',c{d});
                alpha(0.5);
                
                % 与0比较
                if d == 1 % spiral
                    [~,p] = ttest(temp,0);
                elseif d == 2 % rotation
                    [~,p] = ttest(temp,0,'tail','left'); % 单边检验，小于0
                elseif  d==3 % heading cell
                    [~,p] = ttest(temp,0,'tail','right'); % 单边检验，大于0
                end
                
                text(-1.3,15+d*2,sprintf('p = %2.2g',p),'color',c{d});
                xlim([-1.5 1.5]);
                
                
                % 拟合正态分布
                x = xbins_c{d}';
                y = nbins_type';
                
                % Set up fittype and options.
                ft_ = fittype('a/(sqrt(2*pi)*sigma)*exp( -((x-u)^2) / (2*sigma^2) )', 'independent', 'x', 'dependent', 'y','coefficients',{'a','u','sigma'});
                fo_ = fitoptions('method','NonlinearleastSquares','Robust','On','Algorithm','T',...
                    'StartPoint',[2 0 1.4284],'Lower',[0 -1.5 0],'upper',[Inf 1.5 Inf],'MaxFunEvals',5000);
                
                [cf_, goodness, output]= fit(x,y,ft_,fo_);
                if output.exitflag > 0 && cf_.sigma >= 0.3 && cf_.sigma < 1000
                    xx = -1.5:0.001:1.5;
                    yy = feval(cf_,xx);
                    hold on
                    plot(xx,yy,'k-');
                end
                
                xlabel('d prime index');
                ylabel('Cases');
            end
        end
        
        SetFigure(13);
        figN = figN+1;
    end

    function f2p3(debug) % Preferred Direction
        if debug  ; dbstack;   keyboard;      end
        set(figure(figN),'Position',[90,80 1000,700], 'Name', 'Preferred Direction');
        
        methods_of_select = {
            % select_all, 'All cells';
            select_typical, 'Typical cells';
            % select_no_typical, 'Non-typical cells'
            };
        
        % VectorSum result
        pref_direc = cell2mat({group_result.pref_direc}');
        %         motion_type_theta = cell2mat({group_result.motion_type_theta}');
        
        xbins{1} = linspace(-180,180,11);
        xbins{2} = linspace(-180,180,11);
        xbins{3} = linspace(-90,90,11);
        
        xtl{1} = {[-180:90:180]};
        xtl{2} = {'Con' 'CW' 'Exp' 'CCW' 'Con'};
        
        yt(1) = {'Translation preference'};
        yt(2) = {'Rotation preference'};
        yt(3) = {'Elevation preference in Spiral Space'};
        
        t(1) = {'Bottom-line(all)'};
        t(2) = {'All cell type'};
        t(3) = {'Translation cell'};
        t(4) = {'Rotation cell'};
        t(5) = {'Spiral cell'};
        t(6) = {'Motion type cell'};
        
        % 预处理决定distribution ylim 范围
        for d = 1:2 % heading pref and rotation pref
            nbins(d,:) = hist(pref_direc(methods_of_select{1,1},d),xbins{d});
        end
        ylimit = ceil(max(nbins(:))/10)*10;
        
        % 预处理决定polarplot rlim 范围
        h_temp = figure(9999);
        for d = 1:2 % heading pref and rotation pref
            h = polarhistogram(HeadingToazi(pref_direc(methods_of_select{1,1},d))*pi/180,16);
            rbins(d,:) = h.Values;
        end
        rlimit = max(rbins(:));
        close(h_temp);
        
        for d = 1:2 % heading pref and rotation pref
            subplot(2,3,d)
            
            % bar plot
            nbins = hist(pref_direc(methods_of_select{1,1},d),xbins{d});
            hold on
            h = bar(xbins{d},nbins,'facecolor',c{d});
            xlim([-205 205]);
            ylim([0 ylimit+ylimit/3]);
            set(gca,'XTick',[-180:90:180]);
            set(gca,'xticklabel',xtl{d});
            
            xlabel(yt(d));
            title(t(selection_num));
            
            % 检验tuning的分布情况
            % 1. uniformity test,p>0.05为符合该mode
            range  = 0 : 360;
            num_perm =1000;
            bin_num=16;
            
            az = HeadingToazi(pref_direc(methods_of_select{1,1},d));
            az = az';
            
            dim = size(az);
            V1 = az;
            V1 = V1';
            dim = size(V1);
            hist_raw = hist(V1,bin_num);
            R = sum( (hist_raw - mean(hist_raw)).^2 );
            
            for b = 1:num_perm
                if min(range) == 0
                    diff_angle_perm = rand(1,dim(1))*max(range);
                elseif min(range) == -1
                    diff_angle_perm = rand(1,dim(1))*2-1;
                else
                    warning = 'reset data range properly';
                end
                hist_raw_perm = hist(diff_angle_perm,bin_num);
                R_perm(b) = sum( (hist_raw_perm - mean(hist_raw)).^2 );
                %                 plot(hist_raw_perm,'b-');
                %                 hold on;
            end
            p_uniform(d) =length(find(R_perm>R)) / num_perm;
            
            % KS test
            [h,p,~,~] = kstest(az',[az',unifcdf(az')]); % h=1,p<0.05 不是均匀分布
            
            % circular data uniformity test: Hodges-Ajne test
            p_HA(d) = circ_otest(deg2rad(az));  % h=1,p<0.05 不是均匀分布
            
            % 2. modality test   p>0.05为符合该mode
            [mode(d), p_mod{d}] = modalityTestForYong(az,[1 2],1000);
            
            text(-170,ylimit+ylimit/3.2,sprintf('Hodges-Ajne test,p=%0.3g',p_HA(d)),'fontsize',10); % circular data uniformity test
            text(-170,ylimit+ylimit/5,sprintf('modality test, p_u_n_i=%0.2g',p_mod{d}(1)),'fontsize',10);
            text(-170,ylimit+ylimit/10,sprintf('modality test, p_b_i=%0.2g',p_mod{d}(2)),'fontsize',10);
            
            % uniformly distribution?
            % fit Bingham function, ANOVA, Bonferroni-Holm corrected
            % Ting-Yu Chang, 2020. Functional links between sensory
            % representations, .....
            
            % rose plot
            subplot(2,3,3+d)
            h = polarhistogram(HeadingToazi(pref_direc(methods_of_select{1,1},d))*pi/180,16,'facecolor',c{d},'facealpha',0.5);
            if d == 1
                set(gca,'ThetaTick',[0:45:315],'ThetaTickLabel',{'90','45','0','-45','-90','-135','180','135'});
            else
                set(gca,'ThetaTick',[0:45:315],'ThetaTickLabel',{'CCW','45','0','-45','CW','-135','180','135'});
            end
            rlim([0 rlimit]);
        end
        
        % spiral space preferrence: show the preferred elevation represent the heading or rotation pref
        pref_el = cell2mat({group_result.pref_el}');
        pref_el_3D = pref_el(:,3);
        
        % bar plot
        subplot(2,3,3)
        nbins = hist(pref_el_3D(methods_of_select{1,1}),xbins{3});
        hold on
        bar(xbins{3},nbins,'facecolor',c{4});
        set(gca,'XTick',[-90:30:90]);
        xlabel(yt(3));
        title(t(selection_num));
        
        % rose plot: -90 南极点，cw self motion， ccw stim
        subplot(2,3,6)
        polarhistogram(pref_el_3D(methods_of_select{1,1})*pi/180,12,'facecolor',c{4},'facealpha',0.5);
        set(gca,'thetalim',[-90 90],'ThetaTick',[-90:30:90]);
        
        SetFigure(13);
        figN = figN+1;
        
        % preferred spiral directions as Graziano, Andersen, Snowden, 1994
        % compass plot: x axis: rotation; y axis: expasion-contraction;
        % intermediate: spiral motion
        % same as rotation prefer above
    end

    function f2p23(debug) % 'PD vs DDI/d prime'
        if debug  ; dbstack;   keyboard;      end
        methods_of_select = {select_typical, 'Typical cells';};
        
        % preferred direction
        pref_direc = cell2mat({group_result.pref_direc}');
        DDI = cell2mat({group_result.DDI_excludeExpCon}');
        d_prime_coarse = cell2mat({group_result.d_prime_coarse}');
        
        pref_direc = pref_direc(methods_of_select{1,1},:);
        DDI = DDI(methods_of_select{1,1},:);
        d_prime = log10(abs(d_prime_coarse(methods_of_select{1,1},:)));
        
        % PD vs DDI
        set(figure(figN),'Position',[36 81 1500 650], 'Name', 'PD vs DDI');clf
        for hr = 1:2
            subplot(1,2,hr)
            select{1} = logical(abs(pref_direc(:,hr))<=45 & abs(pref_direc(:,hr))>=0); % 0-45
            select{2} = logical(abs(pref_direc(:,hr))<=90 & abs(pref_direc(:,hr))>45); % 45-90
            select{3} = logical(abs(pref_direc(:,hr))<=135 & abs(pref_direc(:,hr))>90); % 90-135
            select{4} = logical(abs(pref_direc(:,hr))<=180 & abs(pref_direc(:,hr))>135); % 135-180
            
            for i = 1:4
                DDI_group{hr,i} = DDI(select{i},hr); % 0-45
                mean_DDI_group(hr,i) = nanmean(DDI_group{hr,i});
                sem_DDI_group(hr,i) = nanstd(DDI_group{hr,i}) / sqrt(sum(select{i}));
                
                dprime_group{hr,i} = d_prime(select{i},hr); % 0-45
                mean_dprime_group(hr,i) = nanmean(dprime_group{hr,i});
                sem_dprime_group(hr,i) = nanstd(dprime_group{hr,i}) / sqrt(sum(select{i}));
                
                hold on
                % bar
                bar(i,mean_DDI_group(hr,i),0.6,'facecolor','none','edgecolor','k','linewidth',1.5);
                
                % scatter
                plot(rand(size(DDI_group{hr,i}))*0.4-0.2+i,DDI_group{hr,i},'.','color',c{hr});
                
                % mean±std
                errorbar(i,mean_DDI_group(hr,i),sem_DDI_group(hr,i),'color','k','LineStyle','none','linewidth',1.5);
            end
            ylim([0 1]);
            xticks([1:4]);
            xticklabels({'0-45','45-90','90-135','135-180'});
            xlabel('|Preferred direction from forward| (°)');
            ylabel('DDI');
            
            % anova
%             [p,~,stats] = anova1([DDI_group{hr,1};DDI_group{hr,2};DDI_group{hr,3};DDI_group{hr,4}],[ones(size(DDI_group{hr,1}));2*ones(size(DDI_group{hr,2}));3*ones(size(DDI_group{hr,3}));4*ones(size(DDI_group{hr,4}))]);
%             [c1, m1, h, mnms] = multcompare(stats) ;
            
            % 两两ttest
            [~,p_t(1)] = ttest2(DDI_group{hr,1},DDI_group{hr,2});
            [~,p_t(2)] = ttest2(DDI_group{hr,1},DDI_group{hr,3});
            [~,p_t(3)] = ttest2(DDI_group{hr,1},DDI_group{hr,4});
            [~,p_t(4)] = ttest2(DDI_group{hr,2},DDI_group{hr,3});
            [~,p_t(5)] = ttest2(DDI_group{hr,2},DDI_group{hr,4});
            [~,p_t(6)] = ttest2(DDI_group{hr,3},DDI_group{hr,4});
            str = [sprintf('12,%2.3g; 13,%2.3g; 14,%2.3g',p_t(1),p_t(2),p_t(3)) newline sprintf('23,%2.3g; 24,%2.3g; 34,%2.3g',p_t(4),p_t(5),p_t(6))];
            text(0.5,1,str);
        end
        suptitle('PD vs DDI');
        SetFigure(13);
        figN = figN + 1;
        
        % PD vs d'
        set(figure(figN),'Position',[100 81 1500 650], 'Name', 'PD vs d prime');clf
        for hr = 1:2
            subplot(1,2,hr)
            select{1} = logical(abs(pref_direc(:,hr))<=45 & abs(pref_direc(:,hr))>=0); % 0-45
            select{2} = logical(abs(pref_direc(:,hr))<=90 & abs(pref_direc(:,hr))>45); % 45-90
            select{3} = logical(abs(pref_direc(:,hr))<=135 & abs(pref_direc(:,hr))>90); % 90-135
            select{4} = logical(abs(pref_direc(:,hr))<=180 & abs(pref_direc(:,hr))>135); % 135-180
            
            for i = 1:4
                d_prime_group{hr,i} = d_prime(select{i},hr); % 0-45
                mean_d_prime_group(hr,i) = nanmean(d_prime_group{hr,i});
                sem_d_prime_group(hr,i) = nanstd(d_prime_group{hr,i}) / sqrt(sum(select{i}));
                
                dprime_group{hr,i} = d_prime(select{i},hr); % 0-45
                mean_dprime_group(hr,i) = nanmean(dprime_group{hr,i});
                sem_dprime_group(hr,i) = nanstd(dprime_group{hr,i}) / sqrt(sum(select{i}));
                
                hold on
                % bar
%                 bar(i,mean_d_prime_group(hr,i),0.6,'facecolor','none','edgecolor','k','linewidth',1.5);
                
                % scatter
                plot(rand(size(d_prime_group{hr,i}))*0.4-0.2+i,d_prime_group{hr,i},'.','color',c{hr});
                
                % mean±std
                errorbar(i,mean_d_prime_group(hr,i),sem_d_prime_group(hr,i),'color','k','LineStyle','none','linewidth',1.5);
            end
%             ylim([0 1]);
            xticks([1:4]);
            xticklabels({'0-45','45-90','90-135','135-180'});
            xlabel('|Preferred direction from forward| (°)');
            ylabel('log(|d prime|)');

            % 两两ttest
            [~,p_t(1)] = ttest2(d_prime_group{hr,1},d_prime_group{hr,2});
            [~,p_t(2)] = ttest2(d_prime_group{hr,1},d_prime_group{hr,3});
            [~,p_t(3)] = ttest2(d_prime_group{hr,1},d_prime_group{hr,4});
            [~,p_t(4)] = ttest2(d_prime_group{hr,2},d_prime_group{hr,3});
            [~,p_t(5)] = ttest2(d_prime_group{hr,2},d_prime_group{hr,4});
            [~,p_t(6)] = ttest2(d_prime_group{hr,3},d_prime_group{hr,4});
            str = [sprintf('12,%2.3g; 13,%2.3g; 14,%2.3g',p_t(1),p_t(2),p_t(3)) newline sprintf('23,%2.3g; 24,%2.3g; 34,%2.3g',p_t(4),p_t(5),p_t(6))];
            text(0.5,1,str);
        end
        suptitle('PD vs d prime');
        SetFigure(13);
        figN = figN + 1;
        
    end

    function f2p4(debug)  % Spiral space preference
        if debug  ; dbstack;   keyboard;      end
        set(figure(figN),'Position',[36 81 1098 729], 'Name', '2-D VectorSum Spiral Space Preference');
        
        methods_of_select = {
            % select_all, 'All cells';
            select_typical, 'Typical cells';
            % select_no_typical, 'Non-typical cells'
            };
        
        % VectorSum result
        pref_direc = cell2mat({group_result.pref_direc}'); % 0 = expasion = az90 (in deg)
        pref_el = cell2mat({group_result.pref_el}');
        pref_amp = cell2mat({group_result.pref_amp}');
        pref_direc_3D = pref_direc(:,3); % 0 = expasion = az90 (in deg)
        pref_az_3D = deg2rad(HeadingToazi(pref_direc_3D)); % az0 = right (in rad)
        pref_el_3D = pref_el(:,3);
        pref_amp_3D = pref_amp(:,3);
        
        %         % normalize FR data
        %         pref_direc_norm = cell2mat({group_result.pref_direc_norm}');
        %         pref_el_norm = cell2mat({group_result.pref_el_norm}');
        %         pref_amp_norm = cell2mat({group_result.pref_amp_norm}');
        %         pref_direc_3D_norm = pref_direc_norm(:,3);
        %         pref_el_3D_norm = pref_el_norm(:,3);
        %         pref_amp_3D_norm = pref_amp_norm(:,3);
        
        cell_type = [group_result.cell_type]';
        cell_type = cell_type(methods_of_select{1,1});
        
        for d = loop_num
            if length(loop_num)==3
                az_3D{d} = pref_az_3D([group_result.cell_type]' == d & methods_of_select{1,1});
                direc_3D{d} = pref_direc_3D([group_result.cell_type]' == d & methods_of_select{1,1});
                el_3D{d} = pref_el_3D([group_result.cell_type]' == d & methods_of_select{1,1});
                amp_3D{d} = pref_amp_3D([group_result.cell_type]' == d & methods_of_select{1,1});
            else
                az_3D{d} = pref_az_3D(methods_of_select{1,1});
                direc_3D{d} = pref_direc_3D(methods_of_select{1,1});
                el_3D{d} = pref_el_3D(methods_of_select{1,1});
                amp_3D{d} = pref_amp_3D(methods_of_select{1,1});
            end
        end
        
        %% 2D visualization
        ymult = 1098/ 729;
        h.ax_raw = axes('Position',[0.07 0.07 0.4 0.4*ymult]);
        h.ax_xhist = axes('position',[0.07 0.75 0.4 0.206]);
%         set(h.ax_xhist,'xtick',[]);
        ylabel('Cases');  hold on;
        h.ax_yhist = axes('position',[0.57 0.07 0.206/ymult 0.4*ymult]);
%         set(h.ax_yhist,'xtick',[]);
        ylabel('Cases');  hold on;
        view([90 -90])
                
        % Lambert cylindrical equal-area projection (Snyder, 1987)  将球形的数据投影到直角坐标中
        ele_sin = [sind(-90),sind(-45),sind(0),sind(45),sind(90)];
        % dot plot
        for d = loop_num
            axes(h.ax_raw);
            plot(h.ax_raw, direc_3D{d},sind(el_3D{d}),'ko','markersize',11,'markerfacecolor',c{d});
            hold on
        end
        set(gca,'xlim',[-180 180]);
        set(gca,'ylim',[-1 1]);
        set(gca,'xtick',[-180:45:180],'xticklabel',{'-180','-135','Left' '-45' '0' '45' 'Right','135','180'}); % 实际上只有ele=0，direc=±90时才是Left和right！！
        set(gca,'ytick',ele_sin,'yticklabel',{'-90 (CW)' '-45' '0' '45' '90 (CCW)'});
        xlabel('VectorSum preferred direction (°)');
        ylabel('VectorSum preferred elevation (°)');
        
        % Histogram
        binnum = 8; group_num = 3;
        % 1. x hist
        axes(h.ax_xhist);
        plotGroupedBar(direc_3D,binnum,[-180 180],c,'probability');
        
        % 2. y hist
        axes(h.ax_yhist);
        y_this = cellfun(@(x) sind(x),el_3D,'uniformoutput',0);
        plotGroupedBar(y_this,binnum,[ele_sin(1) ele_sin(end)],c,'probability');
        set(gca,'xticklabel',{'-90' '-45' '0' '45' '90'});
        SetFigure(10)
        figN = figN+1;
        
        
        %% 3D visualization 
        set(figure(figN),'Position',[70,70 1500,1000], 'Name', '3D visualization');
        
        % 3维球面上表示spiral preference, draw cell type spiral space preference on the sphere
        % ignore origin amplitude, use r = 1
        amp_3D_one = cellfun(@(x) ones(size(x)),amp_3D,'uniformoutput',0);
        subplot(2,3,1)
        [X,Y,Z] = sphere(16); % 单位球
        surf(X,Y,Z,'FaceAlpha',0,'facecolor','w');
        hold on
        for d = loop_num
            % 转为笛卡尔坐标系
            [XX,YY,ZZ] = sph2cart(az_3D{d},deg2rad(el_3D{d}),amp_3D_one{d});
            plot3(XX,YY,ZZ,'.','markersize',15,'Color',c{d});
        end
        hold on
        
        % axis
        lims = [-1.5 1.5];
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
        view(20,20)
        axis equal
        
        
        % %         下面的表示方式太复杂了，直接去elevation打得绝对值画hist就好。。。
        %                 % 3维球面上表示spiral preference, 俯视压缩，侧视压缩
        %                 % 1. 俯视图, 主要看点的径向分布，越靠近圆心说明ele越大，越rotation，越靠近圆边说明ele越小，越translation
        %                 subplot(2,3,2)
        %                 hold on % 单位圆
        %                 rectangle('position',[-1 -1 2 2],'curvature',[1 1]); % ele = 0
        %                 rectangle('position',[-cosd(30) -cosd(30) 2*cosd(30) 2*cosd(30)],'curvature',[1 1],'linestyle','--'); % ele = 30
        %                 rectangle('position',[-cosd(60) -cosd(60) 2*cosd(60) 2*cosd(60)],'curvature',[1 1],'linestyle','--'); % ele = 60
        %                 % dot marker
        %                 XX = []; YY = [];
        %                 for d = loop_num
        %                     [XX{d},YY{d},~] = sph2cart(az_3D{d},deg2rad(el_3D{d}),amp_3D_one{d});
        %                     plot(XX{d},YY{d},'.','markersize',15,'Color',c{d});
        %                 end
        %                 plot([0 1],[0 0],'k-');
        %                 axis equal
        %                 ylim([-1 1]); xlim([-1 1]);
        %                 SameScale
        %                 set(gca,'xtick',[-1 1],'xticklabel',{'Left','Right'});
        %                 set(gca,'ytick',[-1 1],'yticklabel',{'Con','Exp'});
        %                 title('俯视图 Elevation distribution (2D)');
        %
        %                 % 1.1 只取圆的半径表示：从圆心到圆边
        %                 subplot(2,3,5)
        %                 % dot
        %                 for d = loop_num
        %                     hold on
        %                     r_dis{d} = sqrt(XX{d}.^2 + YY{d}.^2); % from 0 to 1
        %                     % normalize the gab between ele0(=1), ele30, ele60, ele90(=0)
        %                     r_new{d} = pi/2 - acos(r_dis{d}); % from 0 to pi/2
        %                     plot(r_new{d},-(-d*10+40),'.','markersize',15,'Color',c{d});
        %                 end
        %                 xlim([0 pi/2]);
        %
        %                 % hist, first plot spiral cell
        %                 xbins = [];
        %                 xbins = linspace(0,pi/2,10);
        %                 y_all = [];
        %                 for d = loop_num
        %                     h = histogram(r_new{4-d},xbins,'facecolor',c{4-d});
        %                     % mean, median
        %                     mean_r(4-d) = nanmean(r_new{4-d});
        %                     meadian_r(4-d) = nanmedian(r_new{4-d});
        %                     y_all = [y_all h.Values];
        %                 end
        %                 ymax = max(y_all) + max(y_all)/8;
        %                 h = plot(mean_r,ymax,'v','color',c{1});
        %                 set(h(2),'color',c{2}); set(h(3),'color',c{3});
        %                 plot([pi/6,pi/6],[0 ymax],'k--');
        %                 plot([pi/3,pi/3],[0 ymax],'k--');
        %                 plot([pi/2,pi/2],[0 ymax],'k-');
        %                 set(gca,'xtick',[0 pi/6 pi/3 pi/2],'xticklabel',{'Ele90','Ele60','Ele30','Ele0'});
        %
        %                 p(1) = ranksum(r_new{1},r_new{2}); % median different?
        %                 p(2) = ranksum(r_new{1},r_new{3});
        %                 p(3) = ranksum(r_new{2},r_new{3});
        %                 text(0.1,ymax,sprintf('T vs R, p = %0.3g',p(1)));
        %                 text(0.1,ymax-5,sprintf('T vs S, p = %0.3g',p(2)));
        %                 text(0.1,ymax-10,sprintf('R vs S, p = %0.3g',p(3)));
        %                 title('Elevation distribution (1D)');
        %
        %
        %                 % 2. 侧视图
        %                 subplot(2,3,3)
        %                 binedge = linspace(-pi/2,pi/2,11);
        %                 for d = loop_num
        %                     % circle marker
        %                     polarplot(deg2rad(el_3D{d}),0.7+d*0.1,'o','color',c{d},'markerfacecolor',c{d});
        %                     hold on
        %
        %                     % hist, first plot spiral cell
        %                     nbins = histcounts(deg2rad(el_3D{4-d}),binedge);
        %                     nbins = [nbins 0];
        %                     nbins = mapminmax(nbins,0,0.6);
        %                     nbins = nbins(1:end-1);
        %                     polarhistogram('BinEdges',binedge,'BinCounts',nbins,'facecolor',c{4-d},'edgecolor','k');
        %                 end
        %                 set(gca,'thetalim',[-90 90],'ThetaTick',[-90:45:90],'Rlim',[0 1],'RTick',[]);
        %                 title('侧视图 Elevation distribution');
        %                 figN = figN+1;
        
        % elevation distribution
        subplot(2,3,4)
        ele_abs = cellfun(@(x) abs(x),el_3D,'UniformOutput',0);
        
        xbins = [];
        xbins = linspace(0,90,10);
        y_all = [];
        hold on
        %         for d = loop_num
        %             h = histogram(ele_abs{4-d},xbins,'facecolor',c{4-d});
        %             y_all = [y_all h.Values];
        %
        %             mean_ele(4-d) = nanmean(ele_abs{4-d});
        %             meadian_ele(4-d) = nanmedian(ele_abs{4-d});
        %         end
        h = histogram([el_3D{1};el_3D{2};el_3D{3}],xbins,'facecolor','k');
        y_all = [y_all h.Values];
        
        mean_ele(4-d) = nanmean(ele_abs{4-d});
        meadian_ele(4-d) = nanmedian(ele_abs{4-d});
        
        
        xlim([0 90]);
        set(gca,'xtick',[0 30 60 90],'xticklabel',{'0 (Trans.)','30','60','90 (Rot.)'});
        xlabel('Elevation'); ylabel('Cases');
        ymax = max(y_all)+max(y_all)/6;
        
%         h = plot(meadian_ele,ymax,'v','color',c{1});
%         set(h(2),'color',c{2}); set(h(3),'color',c{3});
%         p(1) = ranksum(ele_abs{1},ele_abs{2}); % median different?
%         p(2) = ranksum(ele_abs{1},ele_abs{3});
%         p(3) = ranksum(ele_abs{2},ele_abs{3});
%         text(65,ymax,sprintf('T vs R, p = %0.3g',p(1)));
%         text(65,ymax-3,sprintf('T vs S, p = %0.3g',p(2)));
%         text(65,ymax-6,sprintf('R vs S, p = %0.3g',p(3)));
        title('Vector-sum Elevation distribution');
        
        % azimuth distribution
        subplot(2,3,5)
        xbins = [];
        xbins = linspace(-180,180,13);
        y_all = [];
        hold on
%         for d = loop_num
%             h = histogram(direc_3D{4-d},xbins,'facecolor',c{4-d});
%             y_all = [y_all h.Values];
%         end
                h = histogram([direc_3D{1};direc_3D{2};direc_3D{3}],xbins,'facecolor','k');
        y_all = [y_all h.Values];
        
        
        xlim([-180 180]);
        %         set(gca,'xtick',[-180:45:180],'xticklabel',{'-180 (Con)','-135','-90 (Left)','-45','0 (Exp)','45','90 (Right)','135','180 (Con)'});
        set(gca,'xtick',[-180:45:180]);
        text(-180,-4,'Con','horizontalAlignment','center');
        text(-90,-4,'Left','horizontalAlignment','center');
        text(0,-4,'Exp','horizontalAlignment','center');
        text(90,-4,'Right','horizontalAlignment','center');
        text(180,-4,'Con','horizontalAlignment','center');
        xlabel('Direction'); ylabel('Cases');
        title('Vector-sum Direction distribution');
        
        % amplitude distribution
        % 因为时vector-sum的结果，spiral细胞因为有两个相互垂直的分量，所以做vector-sum的时候amplitude变大
        subplot(2,3,6)
        amp_all = [];
        for d = loop_num
            amp_all = [amp_all;amp_3D{d}];
        end
        xbins = [];
        xbins = linspace(min(amp_all),max(amp_all),11);
        hold on
        for d = loop_num
            histogram(amp_3D{4-d},xbins,'facecolor',c{4-d});
        end
        xlabel('Amplitude'); ylabel('Cases');
        title('Vector-sum Amplitude distribution');
        
        SetFigure(10);
        figN = figN + 1;
        
        
    end

    function f2p5(debug) % Translation, rotation prefererence correlation
        if debug  ; dbstack;   keyboard;      end
        set(figure(figN),'Position',[30,50 800,800], 'Name', 'Translation, rotation prefererence correlation');
        
        methods_of_select = {
            % select_all, 'All cells';
            select_typical, 'Typical cells';
            % select_no_typical, 'Non-typical cells'
            };
        
        % VectorSum result
        pref_direc_temp = cell2mat({group_result.pref_direc}');
        
        for d = loop_num
            if length(loop_num)==3
                pref_direc = pref_direc_temp([group_result.cell_type]' == d & methods_of_select{1,1},:);
            else
                pref_direc = pref_direc_temp(methods_of_select{1,1},:);
            end
            
            plot(pref_direc(:,1),pref_direc(:,2),'ko','markersize',15,'markerfacecolor',c{d});
            hold on
            
            xlim([-180 180]);ylim([-180 180]);
            set(gca,'xtick',-180:90:180);
            set(gca,'ytick',-180:90:180);
            set(gca,'yticklabel',{'Con (-180)' 'CW (-90)' 'Exp (0)' 'CCW (90)' 'Con (180)'});
            
            xlabel('Translation preference');
            ylabel('Rotation preference');
        end
        if loop_num==1
            title('Translation cell');
        elseif loop_num==2
            title('Rotation cell');
        elseif loop_num==3
            title('Spiral cell');
        else
            title('All cell');
        end
        % 4 * 4等分
        dir_range = [-180:90:180];
        hold on
        for dis = 2:4
            h = plot([dir_range(dis) dir_range(dis)],[-180 180],'k:');
            h = plot([-180 180],[dir_range(dis) dir_range(dis)],'k:');
        end
        SetFigure(13);
        figN = figN+1;
        
        % Show individual cell selected from the figure. HH20150424
        if length(loop_num)==3
            select_actual_plot = logical(([group_result.cell_type]'==1 | [group_result.cell_type]'==2 | [group_result.cell_type]' == 3) & select_all);
        else
            select_actual_plot = logical([group_result.cell_type]'==loop_num & select_all);
        end
        pref_direc = pref_direc_temp(select_actual_plot,:);
        h_all = plot(pref_direc(:,1),pref_direc(:,2),'visible','off');hold on;
        set([gca h_all],'ButtonDownFcn',{@Show_individual_cell, h_all, select_actual_plot});
        
        
        % DDI （include exp,con）表示每个unit的长轴和短轴，画椭圆/十字
        %         % 改用d prime coarse
        % 用DDI exc
        DDI_temp = cell2mat({group_result.DDI_excludeExpCon}');
        DDI = DDI_temp(methods_of_select{1,1},:);
        
        d_prime_coarse_temp = cell2mat({group_result.d_prime_coarse}');
        d_prime_coarse = abs(d_prime_coarse_temp(methods_of_select{1,1},:)); % 绝对值 |d'|
        
        pref_direc_temp = cell2mat({group_result.pref_direc}');
        pref_direc = pref_direc_temp(methods_of_select{1,1},:);
        
        cell_type = [group_result.cell_type]';
        cell_type = cell_type(methods_of_select{1,1});
        
%         DDI = d_prime_coarse; % 改用d prime coarse
        
        % 画椭圆表示
        set(figure(figN),'Position',[400,150 800,800], 'Name', '2-D Spiral Space Preference with d prime (Ellipse)');
        for n = 1:length(DDI)
            DDI_this = DDI(n,:);
            if ~isnan(DDI_this(1))
                if DDI_this(1) > DDI_this(2) % heading tuning 较强,焦点在x轴上
                    semimajor = DDI_this(1)*10; % 半长轴,a, 长轴=2a
                    semiminor = DDI_this(2)*10; % 半短轴,b，短轴=2b
                    azi = 0; % 方向角
                else % 焦点在y轴上
                    semimajor = DDI_this(2)*10;
                    semiminor = DDI_this(1)*10;
                    azi = pi/2; % 方向角
                end
                hold on
                if ~isnan(cell_type(n))
                    if cell_type(n)==1 || cell_type(n)==2 || (cell_type(n)==3 && randi(3)==1) % spiral cell太多了，随机只画1/3的细胞
                        PlotEllipse(pref_direc(n,1),pref_direc(n,2),2*semimajor,2*semiminor,azi,{'color',c{cell_type(n)}});
                    end
                end
            end
        end
%         % 原始数据点
%         hold on
%         plot(pref_direc(:,1),pref_direc(:,2),'k.');
        
        xlim([-185 185]);ylim([-185 185]);
        set(gca,'xtick',-180:90:180);
        set(gca,'ytick',-180:90:180);
        set(gca,'yticklabel',{'Con' 'CW' 'Exp' 'CCW' 'Con'});
        xlabel('Translation preference');
        ylabel('Rotation preference');
        
        % marker the largest DDI cell
        hold on
        [~,largeone] = max(DDI(:,1).*DDI(:,2));
        if DDI(largeone,1)>DDI(largeone,2)
            semimajor = DDI(largeone,1)*10;
            semiminor = DDI(largeone,2)*10;
            azi = 0; % 方向角
        else
            semimajor = DDI(largeone,2)*10;
            semiminor = DDI(largeone,1)*10;
            azi = pi/2; % 方向角
        end
        PlotEllipse(pref_direc(largeone,1),pref_direc(largeone,2),2*semimajor,2*semiminor,azi,{'color','c'});
        
        % 辅助线4 * 4等分
        dir_range = [-180:90:180];
        hold on
        for dis = 2:4
            h = plot([dir_range(dis) dir_range(dis)],[-180 180],'k:');
            h = plot([-180 180],[dir_range(dis) dir_range(dis)],'k:');
        end
%         title('Semimaj/Semiminor = d prime coarse');
        title('Semimaj/Semiminor = DDI');
        SetFigure(20);figN = figN+1;
        
        % 画十字
        set(figure(figN),'Position',[800,150 800,800], 'Name', '2-D Spiral Space Preference with d prime (Cross)');
        for n = 1:length(DDI)
            DDI_this = DDI(n,:); % from 0.5-1
            if ~isnan(DDI_this(1))
                xlength = DDI_this(1)*20;
                ylength = DDI_this(2)*20;
                %                 if xlength > ylength
                %                     xlength = 2 * xlength; % 较大的值乘以2倍显示，区分两者差别
                %                 else
                %                     ylength = 2 * ylength;
                %                 end
                xx1 = pref_direc(n,1)-xlength/2;
                xx2 = pref_direc(n,1)+xlength/2;
                yy1 = pref_direc(n,2)-ylength/2;
                yy2 = pref_direc(n,2)+ylength/2;
                hold on
                
                if ~isnan(cell_type(n))
                    plot([xx1 xx2],[pref_direc(n,2) pref_direc(n,2)],'color',c{cell_type(n)});
                    plot([pref_direc(n,1) pref_direc(n,1)],[yy1 yy2],'color',c{cell_type(n)});
                end
            end
        end
        xlim([-185 185]);ylim([-185 185]);
        set(gca,'xtick',-180:90:180);
        set(gca,'ytick',-180:90:180);
        set(gca,'yticklabel',{'Con' 'CW' 'Exp' 'CCW' 'Con'});
        xlabel('Translation preference');
        ylabel('Rotation preference');
        
        % 辅助线4 * 4等分
        dir_range = [-180:90:180];
        for d = 2:4
            h = plot([dir_range(d) dir_range(d)],[-180 180],'k:');
            h = plot([-180 180],[dir_range(d) dir_range(d)],'k:');
        end
%         title('Cross length = d prime coarse');
        title('Cross length = DDI');
        SetFigure(20);figN = figN+1;
        
    end

    function f2p12(debug) % 'Preferred direction correlation with RF location',
        if debug  ; dbstack;   keyboard;      end
        
        methods_of_select = {select_typical, 'Typical cells';};
        
        % show the cell with RF location
        % RF: X Y H W (X,Y是RF中心点位置)
        rf_temp =  {group_result.RF}';
        have_rf = ~cellfun(@isempty,rf_temp);
        rf_temp(~have_rf) = {[nan nan nan nan]};
        rf = cell2mat(rf_temp); % RF: X Y H W (X,Y是RF中心点位置)
        rf_area = rf(:,3) .* rf(:,4);
        rf_dia = 2.* sqrt(rf_area / pi); % 直径
        rf_ecc = sqrt(rf(:,1).^2 + rf(:,2).^2); % RF中心点距离原点的距离
        rf_ang = atan(rf(:,2) ./ rf(:,1)) * 180/pi; % rf_ang>0 中心点在上半视野，rf_ang<0 中心点在下半视野
        
        select_withRF = logical(~isnan(rf(:,1)) & methods_of_select{1,1});
        
        % 挑选RF较小的cell。 理论上ecc范围是0-45°（左右半屏幕45度视角，实际比这个要大一点）
        select_RFsmall = logical((rf_ecc - rf_dia/2) > -min(rf_dia)/2 & rf_dia<45); % RF不过fovea或者只过fovea很小距离的cell
        
        % 挑选RF location
        select_RFupper = logical(rf_ang > 0); % 上半视野
        select_RFlower = logical(rf_ang < 0); % 下半视野
        
        % VectorSum result
        pref_direc = [];
        pref_direc = cell2mat({group_result.pref_direc}');
        
        set(figure(figN),'Position',[50,50 700,700], 'Name', 'Prefer Direction with RF location');
        h1 = plot(pref_direc(:,1),pref_direc(:,2),'ko','markersize',11); % all unit
        hold on
        h2 = plot(pref_direc(select_withRF,1),pref_direc(select_withRF,2),'ko','markersize',11,'markerfacecolor',[0.6 0.6 0.6]); % unit with RF
        
        % unit with RF location, 暂时用颜色表示rf位于上/下视野，完善：在圆中填充上半/下半部分颜色表示RF位置
        h3 = plot(pref_direc(select_withRF & select_RFsmall & select_RFupper,1),pref_direc(select_withRF & select_RFsmall & select_RFupper,2),'ko','markersize',11,'markerfacecolor','r');
        h4 = plot(pref_direc(select_withRF & select_RFsmall & select_RFlower,1),pref_direc(select_withRF & select_RFsmall & select_RFlower,2),'ko','markersize',11,'markerfacecolor','b');
        
        xlim([-180 180]);ylim([-180 180]);
        set(gca,'xtick',-180:90:180);
        set(gca,'ytick',-180:90:180);
        set(gca,'yticklabel',{'Con' 'CW' 'Exp' 'CCW' 'Con'});
        
        xlabel('Translation preference');
        ylabel('Rotation preference');
        
        % 4 * 4等分
        dir_range = [-180:90:180];
        for d = 2:4
            h = plot([dir_range(d) dir_range(d)],[-180 180],'k:');
            h = plot([-180 180],[dir_range(d) dir_range(d)],'k:');
        end
        set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        
        legend([h1 h2 h3 h4],'no RF','large RF','upper small RF','lower small RF','location','east');
        SetFigure(13);
        title('Diagonal match in Lower RF, Ant-idiagonal match in Upper RF');
        
        % Show individual cell selected from the figure. HH20150424
        select_actual_plot = select_withRF;
        h_all = plot(pref_direc(select_actual_plot,1),pref_direc(select_actual_plot,2),'visible','off');hold on;
        set(get(get(h_all,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        set([gca h_all],'ButtonDownFcn',{@Show_individual_cell, h_all, select_actual_plot});
        figN = figN+1;
        
        % 根据RF判断Heading preference与rotation preference的 local optic flow是否相似
        
        % !!!!!!!!!!!!!!!!!!!!!!!
        % !!!  注意 !!!!!!!!!!!!!
        % 仅仅靠RF的位置判断local optic flow不太准确，最好的办法是获取RF内的optic flow的向量，求和，并比较！
        
        
        % 对pref_direction编号1-8（等分为4*4，主要集中在对角线，所以只取对角线的8块区域）
        % 1. 正对角线的4个block
        for d = 1:4
            heading_pref_select = logical(pref_direc(:,1) > dir_range(d) & pref_direc(:,1) < dir_range(d+1) ); % 从-180到180：1 2 3 4
            rotation_pref_select = logical(pref_direc(:,2) > dir_range(d) & pref_direc(:,2) < dir_range(d+1) ); % 1 2 3 4
            select_prefDir{d} = logical(heading_pref_select & rotation_pref_select);
        end
        
        % 2. 反对角线的4个block
        for d = 1:4
            heading_pref_select = logical(pref_direc(:,1) > dir_range(d) & pref_direc(:,1) < dir_range(d+1) ); % 1 2 3 4
            rotation_pref_select = logical(pref_direc(:,2) > dir_range(5-d) & pref_direc(:,2) < dir_range(6-d) ); % 4 3 2 1
            select_prefDir{d+4} = logical(heading_pref_select & rotation_pref_select);
        end
        
        %         % 这8个block的RF位置预期(match case)为：（假设细胞prefer相似的 local optic  flow）
        %         % ！！ 此方法不准确！！！！！！！！！！！！！！
        %         % ***********************
        %         % ***  1-4. Lower RF  ***
        %         % ***  5-8: Upper RF  ***
        %         % ***********************
        %         set(figure(figN),'Position',[200,200 700,700], 'Name', 'Prefer Direction match or not');
        %         h1 = plot(pref_direc(select_withRF,1),pref_direc(select_withRF,2),'ko','markersize',11,'markerfacecolor',[0.6 0.6 0.6]); % unit with RF
        %         hold on
        %         for i = 1:8
        %             if i < 5
        %                 try
        %                     h2 = plot(pref_direc(select_prefDir{i} & select_RFsmall & select_RFlower,1),...
        %                         pref_direc(select_prefDir{i} & select_RFsmall & select_RFlower,2),'ko','markersize',11,'markerfacecolor','r'); % match
        %                 end
        %                 try
        %                     h3 = plot(pref_direc(select_prefDir{i} & select_RFsmall & select_RFupper,1),...
        %                         pref_direc(select_prefDir{i} & select_RFsmall & select_RFupper,2),'ko','markersize',11,'markerfacecolor','b'); % not match
        %                 end
        %             else
        %                 try
        %                     h3 = plot(pref_direc(select_prefDir{i} & select_RFsmall & select_RFlower,1),...
        %                         pref_direc(select_prefDir{i} & select_RFsmall & select_RFlower,2),'ko','markersize',11,'markerfacecolor','b'); % not match
        %                 end
        %                 try
        %                     h2 = plot(pref_direc(select_prefDir{i} & select_RFsmall & select_RFupper,1),...
        %                         pref_direc(select_prefDir{i} & select_RFsmall & select_RFupper,2),'ko','markersize',11,'markerfacecolor','r'); % match
        %                 end
        %             end
        %         end
        %         xlim([-180 180]);ylim([-180 180]);
        %         set(gca,'xtick',-180:90:180);
        %         set(gca,'ytick',-180:90:180);
        %         set(gca,'yticklabel',{'Con' 'CW' 'Exp' 'CCW' 'Con'});
        %
        %         % 4 * 4等分
        %         dir_range = [-180:90:180];
        %         for i = 2:4
        %             h = plot([dir_range(i) dir_range(i)],[-180 180],'k:');
        %             h = plot([-180 180],[dir_range(i) dir_range(i)],'k:');
        %         end
        %
        %         xlabel('Translation preference');
        %         ylabel('Rotation preference');
        %         title('Diagonal match in Lower RF, Ant-idiagonal match in Upper RF');
        %         SetFigure(13);
        %         figN = figN+1;
        
        %         set(figure(figN),'Position',[1000,300 800,400], 'Name', 'Prefer Direction match or not number summary');
        %         for i = 1:8
        %             subplot(2,4,i)
        %             if i < 5
        %                 match = sum(select_prefDir{i} & select_RFlower & select_RFsmall);
        %                 nomatch = sum(select_prefDir{i} & select_RFupper & select_RFsmall);
        %             else
        %                 match = sum(select_prefDir{i} & select_RFupper & select_RFsmall);
        %                 nomatch = sum(select_prefDir{i} & select_RFlower & select_RFsmall);
        %             end
        %             bar([match nomatch]);
        %         end
        %         SetFigure(13);
        %         figN = figN+1;
    end


    function f2p6(debug)  % Center Eye position and Vary Eye position
        if debug  ; dbstack;   keyboard;      end
        set(figure(figN),'Position',[100,100 800,700], 'Name', 'Vary Eye position');
        
        % 比较两个eye position preference的改变
        diffpref = abs(pref_direc_centerEye - pref_direc_varyEye);
        diffpref(diffpref>180) = 360 - diffpref(diffpref>180);
        
        xbins = [];
        xbins = linspace(0,180,15);
        
        for d = 1:2 % heading or rotation pref
            subplot(2,2,d)
            nbins = hist(diffpref(:,d),xbins);
            hold on
            bar(xbins,nbins,1,'facecolor',c{d});
            xlim([-10 190]);
            set(gca,'XTick',[0:45:180]);
            xlabel('|Preference change (deg)|');
            ylabel('Number of cells');
            
            top = max(nbins);
            
            if d==1
                text(130,top,['N = ',num2str(length(diffpref))],'fontsize',15);
            end
            
            % mean
            mean_diff(d) = mean(diffpref(:,d));
            plot(mean_diff(d),top+top/10,'kv','markerfacecolor',c{d});
            text(mean_diff(d)-25,top+top/5,sprintf('Mean = %2.2g',mean_diff(d)),'fontsize',10);
        end
        
        % pref change 与 RF离心度 关系
        ecc = sqrt(RF_varyEye(:,1).^2 +  RF_varyEye(:,2).^2 );
        
        subplot(2,2,3)
        for d = 1:2
            hold on
            plot(ecc,diffpref(:,d),'o','color',c{d});
            
            %             % Type II regression analysis  % 运行很慢很慢
            %             option = 1; % free intercept (default)
            %             corr_option = 1; % spearman correlation
            %             xh = 13:0.01:21;
            %             [slope, intercept, bint, aint, r, p] = regress_perp(ecc, diffpref(:,i), 0.05,option,corr_option);
            %             plot(xh,intercept+slope*xh,'--','color',c{i},'linewidth',2.5);
            % %             text(25,ylimit(2)-2*i,sprintf('R^2 = %5.3g, p = %5.3g',r, p),'color',c{4-i});
        end
        xlabel('RF eccentricity (°)');
        ylabel('Prefference difference');
        xlim([4 31]);
        
        % pref change 与 RF大小关系
        rf_size = RF_varyEye(:,3) .* RF_varyEye(:,4);
        rf_dia = 2*sqrt(rf_size/pi);
        
        subplot(2,2,d+2)
        for d = 1:2
            hold on;
            plot(rf_dia,diffpref(:,d),'o','color',c{d});
            
            %             % Type II regression analysis  % 运行很慢很慢
            %             option = 1; % free intercept (default)
            %             corr_option = 1; % spearman correlation
            %             xh = 0:0.01:100;
            %             [slope, intercept, bint, aint, r, p] = regress_perp(rf_dia, diffpref(:,i), 0.05,option,corr_option);
            %             plot(xh,intercept+slope*xh,'--','color',c{i},'linewidth',2.5);
            %             %             text(25,ylimit(2)-2*i,sprintf('R^2 = %5.3g, p = %5.3g',r, p),'color',c{i}
            
        end
        xlabel('RF diameter (°)');
        ylabel('Prefference difference');
        xlim([20 110]); 
        
        SetFigure(13);
        figN = figN+1;
    end

    function f1p3(debug) % RF, all MST cell, 不一定每个cell都测了RF
        if debug  ; dbstack;   keyboard;      end
        
        % ------------------ default parameter -----------------------
        %         methods_of_select = {select_typical, 'Typical cells'}; % GUI 选 Bottom-line, 所有recording cell， 但是不是每个cell都测了RF
        methods_of_select = {select_all}; % 默认all cell
        
        % RF size: 转化为圆的半径比较
        % RF eccentricity: rf中心距离（0，0）距离
        rf_temp =  {group_result.RF}'; % X,Y,W,H  (X,Y为RF中心点坐标)
        rf_temp = cell2mat(rf_temp(methods_of_select{1,1}));
        have_rf = ~isnan(rf_temp(:,1));
        rf = rf_temp(have_rf,:); % 不是每个cell都测了RF
        rf_size = rf(:,3) .* rf(:,4);
        rf_dia = 2.* sqrt(rf_size / pi);
        rf_sqrt = sqrt(rf_size / pi);
        rf_ecc = sqrt(rf(:,1).^2 + rf(:,2).^2); % Eccentricity is the distance from the fixation point to the center of the receptive field.
        
        % plot MST RF
        set(figure(figN),'Position',[100,250, 1500,500], 'Name', 'MST RF');
        subplot(1,3,1)
        axis equal; box on; axis([-60 60 -60 60]);
        line([-60 60],[0 0],'color','k','LineStyle',':'); hold on;
        line([0 0],[-60 60],'color','k','LineStyle',':');
        set(gca,{'xtick','ytick'},{-60:20:60,-60:20:60});
        xlabel('degree');
        ylabel('degree');
        title('MST RF');
        
        % 用灰色画，避免太多了看不清
        for d = 1:length(rf)
            if ~isnan(rf(d,1))
                rectangle('position',[rf(d,1)-rf(d,3)/2 rf(d,2)-rf(d,4)/2, rf(d,3) rf(d,4)],...
                    'Curvature',[0.3 0.3],'EdgeColor',[0.6 0.6 0.6],'LineWidth',1);
                hold on
            end
        end
        %         % 每隔几个RF大小细胞用黑色画RF
        %         % 重新排列RF：从小到大
        %         [~,Index] = sort(rf_dia);
        %         rf = rf(Index,:);
        %         rf_num = length(rf);
        %         specific_bin = floor(rf_num/10);
        %         specific_this = [1:specific_bin:length(rf)];
        %
        %         for d = 1:length(rf)
        %             if ~isnan(rf(d,1)) && sum(d == specific_this) > 0
        %                 rectangle('position',[rf(d,1)-rf(d,3)/2 rf(d,2)-rf(d,4)/2, rf(d,3) rf(d,4)],...
        %                     'Curvature',[0.3 0.3],'EdgeColor','k','LineWidth',2);
        %                 hold on
        %             end
        %         end
        
        % included fovea
        criti_value = 1; % in visual degree
        includefovea = false(length(rf),1);
        criti_value1 = 3; % in visual degree
        closefovea = false(length(rf),1);
        
        for d = 1:length(rf)
            % include fovea
            if (rf(d,1)-rf(d,3)/2) < -criti_value && (rf(d,1)+rf(d,3)/2) > criti_value && (rf(d,2)-rf(d,4)/2) < -criti_value && (rf(d,2)+rf(d,4)/2) > criti_value
                includefovea(d) = true;
                %                 rectangle('position',[rf(d,1)-rf(d,3)/2 rf(d,2)-rf(d,4)/2, rf(d,3) rf(d,4)],...
                %                     'Curvature',[0.3 0.3],'EdgeColor','r','LineWidth',1.5);
            end
            
            % came close to the fovea but not include(Proximal edge ≤ 3° from the fixation point), Roy, et al., 1992
            if (abs(rf(d,1)-rf(d,3)/2) <= criti_value1) || (abs(rf(d,1)+rf(d,3)/2) <= criti_value1) || (abs(rf(d,2)-rf(d,4)/2) <= criti_value1) || (abs(rf(d,2)+rf(d,4)/2) <= criti_value1) || ... % (Proximal edge ≤ 3° from the fixation point)
                    ((rf(d,1)-rf(d,3)/2) < 0 && (rf(d,1)+rf(d,3)/2) > 0 && (rf(d,2)-rf(d,4)/2) < 0 && (rf(d,2)+rf(d,4)/2) > 0 && ~includefovea(d)) % includefovea中包含fovea但是被criti_value排除的细胞
                
                closefovea(d) = true;
                %                 rectangle('position',[rf(d,1)-rf(d,3)/2 rf(d,2)-rf(d,4)/2, rf(d,3) rf(d,4)],...
                %                     'Curvature',[0.3 0.3],'EdgeColor','b','LineWidth',1.5);
            end
            
            
        end

        text(-55,55,sprintf('All MST neuron N=%g',sum(have_rf))); % all recording neurons, not selected
        text(-55,50,sprintf('Included Fovea N=%g',sum(includefovea)),'color','r');
        text(-55,45,sprintf('Close to Fovea N=%g',sum(closefovea)),'color','b');
        
        % RF distribution
        subplot(1,3,2)
        xbins = [0:10:110];
        h = hist(rf_dia,xbins);
        h2 = hist(rf_dia(includefovea),xbins);
        h3 = hist(rf_dia(closefovea),xbins);
        
        hold on;
        h_h = bar(xbins,h,1,'facecolor','w');
        h_h2 = bar(xbins,[h2;h3]','stacked','barwidth',1);
        set(h_h2(1),'facecolor','r');
        set(h_h2(2),'facecolor','b');
        xlabel('Equivalent diameter of receptive field (degree)');
        ylabel('Number of cases');
        set(gca,'xlim',[0 120]);
        % mean
        mean_d = nanmean(rf_dia);
        top = max(h);
        plot([mean_d mean_d],[0 top+1],'k--');
        plot(mean_d,top+2,'kv','markerfacecolor','k');
        text(mean_d+4,top+2,num2str(mean_d),'fontsize',15);
        hl = legend('All','Included fovea','Close to fovea','location','best');
        set(hl,'Box','off');
        
        % squart root of the area - exxentricity
        subplot(1,3,3)
        plot(rf_ecc,rf_sqrt,'ko');
        xlabel('eccentricity');
        ylabel('sqrt(area)');
        [r,pp] = plot_corr_line(rf_ecc,rf_sqrt,'MethodOfCorr','Pearson','FittingMethod',2,'LineStyle','k-');
        text(30,50,sprintf('r = %g',r));
        text(30,45,sprintf('p = %g',pp));
        
        SetFigure(13);
        figN = figN+1;
        
        
        % plot RF等效圆
        set(figure(figN),'Position',[100,250, 1500,500], 'Name', 'RF equivalent cicle');clf
        for j = 1:2 % contain fovea or not contain
            subplot(1,2,j) % containing fovea
            hold on
            if j == 1
                select = includefovea;
            else
                select = ~includefovea;
            end
            rf_this = rf(select,:);
            rf_dia_this = rf_dia(select);
            for i = 1:length(rf_dia_this)
                rectangle('Position',[rf_this(i,1)-rf_dia_this(i)/2,rf_this(i,2)-rf_dia_this(i)/2,rf_dia_this(i),rf_dia_this(i)],'curvature',[1,1],'EdgeColor',[0.6 0.6 0.6],'LineWidth',1);
            end
            axis square;
            line([-60 60],[0 0],'color','k','LineStyle','-');
            line([0 0],[-60 60],'color','k','LineStyle','-');
            if j == 1
                title('include fovea');
            else
                title('exclude fovea');
            end
        end
        SetFigure(13);
        figN = figN+1;
    end

    function f1p4(debug) % 'RF across site location (Grid-xy)'
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'}; % GUI 选 Bottom-line, 所有recording cell， 但是不是每个cell都测了RF
        
        % RF size: 转化为圆的半径比较
        % RF eccentricity: rf中心距离（0，0）距离
        rf =  {group_result.RF}'; % X,Y,W,H  (X,Y为RF中心点坐标)
        rf = cell2mat(rf);
        rf(~methods_of_select{1,1},:) = nan(sum(~methods_of_select{1,1}),4); % set nan not select site
        rf_size = rf(:,3) .* rf(:,4);
        rf_dia = 2.* sqrt(rf_size / pi);
        rf_sqrt = sqrt(rf_size / pi);
        rf_ecc = sqrt(rf(:,1).^2 + rf(:,2).^2); % Eccentricity is the distance from the fixation point to the center of the receptive field.
        
        % Grid x-y
        GridX = cell2mat({group_result.GridX}');
        GridY = cell2mat({group_result.GridY}');
        GridX(~methods_of_select{1,1}) = nan(sum(~methods_of_select{1,1}),1);
        GridY(~methods_of_select{1,1}) = nan(sum(~methods_of_select{1,1}),1);
        
        uniqueX = unique(GridX);
        uniqueY = unique(GridY);
        uniqueX(isnan(uniqueX)) = [];
        uniqueY(isnan(uniqueY)) = [];
        
        
        % ------------------------ Plot ------------------------------
        set(figure(figN),'name','RF across x-y position','pos',[150,0,800,1000]); clf
        ha = tight_subplot(length(uniqueY), length(uniqueX),[.03 .03], [.05 .05], [.05 .05]);
        mean_d_this = [];
        for xx = 1:length(uniqueX)
            for yy = 1:length(uniqueY)
                axes(ha(xx+(yy-1)*length(uniqueX)))
                % do not show the axis and label
                %                 ha(xx+(yy-1)*length(uniqueX)).YAxis.Visible = 'off';
                %                 ha(xx+(yy-1)*length(uniqueX)).XAxis.Visible = 'off';
                
                % show the x and y number
                if xx == 1 && yy~=1
                    text(-120,0.5,num2str(uniqueY(length(uniqueY)+1-yy)));
                end
                if yy == length(uniqueY)
                    text(-20,-80,num2str(uniqueX(xx)));
                end
                
                axis equal; box on; axis([-60 60 -60 60]);
                line([-60 60],[0 0],'color','k','LineStyle',':'); hold on;
                line([0 0],[-60 60],'color','k','LineStyle',':');
                
                % x-y note
                if xx == 1 && yy == 1
                    text(-110,0,'Middle');
                    text(20,00,'Lateral');
                    text(-50,-65,'Anterior');
                    text(-50,60,'Posterior');
                    ha(1).YAxis.Visible = 'off';
                    ha(1).XAxis.Visible = 'off';
                end
                
                select = logical(GridX == uniqueX(xx) & GridY == uniqueY(length(uniqueY)+1-yy)); % Y 轴从大到小排列，原点(0,0)位于左下角
                select_this = find(GridX == uniqueX(xx) & GridY == uniqueY(length(uniqueY)+1-yy)); % Y 轴从大到小排列，原点(0,0)位于左下角
                
                if sum(select)~=0
                    for d = 1:sum(select)
                        select_temp = select_this(d);
                        if ~isnan(rf(select_temp,1))
                            rectangle('position',[rf(select_temp,1)-rf(select_temp,3)/2 rf(select_temp,2)-rf(select_temp,4)/2, rf(select_temp,3) rf(select_temp,4)],...
                                'Curvature',[0.3 0.3],'EdgeColor',[0.6 0.6 0.6],'LineWidth',1);
                            hold on
                            
                            r_dia_save{xx,yy}(d) = rf_dia(select_temp);
                        end
                    end
                    mean_d_this(xx,yy) = nanmean(rf_dia(select)); % Y轴从小到大，后面需要翻转符合实际Grid分布 （原点(0,0)位于左下角）
                    text(-50,80,num2str(mean_d_this(xx,yy)));
                end
            end
        end
        figN = figN+1;
        
        set(figure(figN),'name','RF across x-y position (headt map)','pos',[1000,150,500,500]); clf
        mean_d = rot90(mean_d_this,3);
        mean_d = fliplr(mean_d);
        
        cmap = colormap('hot');
        heatmap(mean_d);
        colorbar
        figN = figN+1;
        
        set(figure(figN),'name','RF across AP position','pos',[1500,150,500,400]); clf % 前后变化
        r_dia_save = rot90(r_dia_save,3);
        r_dia_save = fliplr(r_dia_save);
        
        % 每一行统计,从大到小（后到前）
        r_dia_row_all = nan(size(r_dia_save,1),100);
        for d = 1:size(r_dia_save,1)
            r_dia_row = cell2mat(r_dia_save(d,:));
            r_dia_row(r_dia_row==0) = [];
            mead_d_row(d,1) = nanmean(r_dia_row);
            sem_d_row(d,1) = nanstd(r_dia_row)/sqrt(sum(~isnan(r_dia_row)));
            
            % save data across row
            r_dia_row_all(d,1:length(r_dia_row)) = r_dia_row;
        end
        
        % 去掉nan
        select_nan = isnan(mead_d_row);
        mead_d_row(select_nan) = [];
        sem_d_row(select_nan) = [];
        r_dia_row_all(:,~any(~isnan(r_dia_row_all))) = [];
        
        uniqueY(flipud(select_nan)) = [];
        
        % 从前到后
        mead_d_row = flipud(mead_d_row);
        sem_d_row = flipud(sem_d_row);
        r_dia_row_all = flipud(r_dia_row_all);
        
        % plot
        hold on
        for i = 1:length(uniqueY)
            if sum(~isnan(r_dia_row_all(i,:)))>0
                for j = 1: sum(~isnan(r_dia_row_all(i,:)))
                    %                 plot(uniqueY(i)+ rand*0.4-0.2,r_dia_row_all(i,j),'ko');
                    plot(i+ rand*0.4-0.2,r_dia_row_all(i,j),'ko');
                end
            end
        end
        barwitherr(sem_d_row,mead_d_row,'facecolor','none');
        set(gca,'xtick',[1:length(uniqueY)],'xticklabel',uniqueY);
        xlabel('Row (A -- P location)');xticks([13:22]);
        ylabel('RF diameter');
        SetFigure(13);
        figN = figN+1;
    end

    function f1p5(debug) % 'Comparison Max FiringRate with Gaussian Fit maxFR'
        if debug  ; dbstack;   keyboard;      end
        
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'}; % 默认T cell+R cell+ Spi cell
        
        % ---------------------- get data ----------------------------
        % 只对fitOK的细胞
        fitOK = cell2mat({group_result.fitOK}'); %  if output.exitflag > 0 && cf_{s}.sigma >= 0.3 && cf_{s}.sigma < 500; fitOK = 1;
        fitOK = fitOK(methods_of_select{1,1},:);
        
        cell_type = [group_result.cell_type]';
        cell_type =  cell_type(methods_of_select{1,1},:);
        
        % %%%%%%%%%%%%%%%%%%%%
        % maxFR and gaussian-maxFR
        maxFR = cell2mat({group_result.maxFR}');
        maxFRGauss = cell2mat({group_result.maxFRGauss}');
        maxFR = maxFR(methods_of_select{1,1},:);
        maxFRGauss = maxFRGauss(methods_of_select{1,1},:);
        
        set(figure(figN),'Position',[100,70,800,800], 'Name', ' max firing rate vs Gaussian Fitting max FR');
        tl{1} = 'Translation tuning';
        tl{2} = 'Rotation tuning';
        for i = 1:2 % h or r
            select = fitOK(:,i) == 1;
            
            % correlation
            subplot(2,2,i)
            plot(maxFR(select,i),maxFRGauss(select,i),'ko');
            SameScale(1);
            [r,pp] = plot_corr_line(maxFR(select,i),maxFRGauss(select,i),'MethodOfCorr','pearson','FittingMethod',2,'LineStyles','r-','LineWidth',2);
            xlimit = xlim;
            ylimit = ylim;
            text((xlimit(2)-xlimit(1))/20+xlimit(1),ylimit(2)-(ylimit(2)-ylimit(1))/20,sprintf('r=%0.3g, p=%0.3g',r,pp));
            xlabel('max FR');
            ylabel('max GaussianFitting FR');
            title(tl{i});
            
            % difference
            subplot(2,2,i+2)
            diffFR = maxFR(select,i)-maxFRGauss(select,i);
            histogram(diffFR,30);
            meandiff = mean(diffFR);
            [~,p] = ttest(diffFR);
            xlimit = xlim;
            ylimit = ylim;
            text((xlimit(2)-xlimit(1))/20+xlimit(1),ylimit(2)-(ylimit(2)-ylimit(1))/20,sprintf('mean = %0.3g',meandiff));
            text((xlimit(2)-xlimit(1))/20+xlimit(1),ylimit(2)-(ylimit(2)-ylimit(1))/7,sprintf('t test, p = %0.3g',p));
            xlabel('maxFR - Gaussian-maxFR');
            ylabel('Cases');
        end
        SetFigure(13)
        figN = figN + 1; 
    end

    function f1p6(debug) %'Response reliability in T and R' / Fanor Factor
        % Armstrong and Moore 2007
        if debug  ; dbstack;   keyboard;      end
        
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'}; % 默认T cell+R cell+ Spi cell
        
        % ---------------------- get data ----------------------------
        % mean response magnitude (spike count)
        % variance across trials
        mean_fr_temp = {group_result.mean_fr_trial}';
        mean_fr_temp = mean_fr_temp(methods_of_select{1,1});
        mean_fr = cell2mat(mean_fr_temp);
        
        std_fr_temp = {group_result.std_fr_trial}';
        std_fr_temp = std_fr_temp(methods_of_select{1,1});
        std_fr = cell2mat(std_fr_temp);
        
        % calculate for each neuron (Bichot,  Thompson, Rao, Schall, 2001)
        progressbar('Fitting for each neuron...');
        for i = 1:length(mean_fr_temp)
            mean_fr_this = mean_fr_temp{i};
            std_fr_this = std_fr_temp{i};
            
            if ~isnan(mean_fr_this(1))
                % fitting： power function 幂函数
                for hr = 1:2
                    f = [];
                    f = fit(mean_fr_this(:,hr),std_fr_this(:,hr),'b*x^m','StartPoint',[1.5,0.5]);
                    power(i,hr) = f.m;
                    coefficient(i,hr) = f.b;
                end
            else
                power(i,hr) = nan;
                coefficient(i,hr) = nan;
            end
            
            progressbar(i/length(mean_fr_temp));
        end
        
        % ttest for power and cofficient
        [~,ppow] = ttest(power(:,1),power(:,2));
        [~,pcoe] = ttest(coefficient(:,1),coefficient(:,2));
        
        
        % 去掉nan and zero
        removethis1 = logical(sum(isnan(mean_fr),2)>0 | sum(isnan(std_fr),2)>0);
        removethis2 = logical(sum(mean_fr==0,2)>0 | sum(std_fr==0,2)>0);
        removethis = logical(removethis1 | removethis2);
        mean_fr(removethis,:) = [];
        std_fr(removethis,:) = [];
        
        xlimit = [min(min(mean_fr)) max(max(mean_fr))];
        ylimit = [min(min(std_fr)) max(max(std_fr))];
        
        % plot
        set(figure(figN),'Position',[610,70,800,800], 'Name', 'Response reliability in T and R');clf
        hold on
        for hr = 1:2
            plot(mean_fr(:,hr),std_fr(:,hr),'.','color',c{hr});
            
            % fitting： power function 幂函数
            f = [];
            [f,gof(hr),output] = fit(mean_fr(:,hr),std_fr(:,hr),'b*x^m','StartPoint',[1.5,0.5]);
            h = plot(f,xlimit,ylimit);
            set(h,'color',c{hr});
            legend off
            
            set(gca,'XScale','log') ;
            set(gca,'YScale','log') ;
            
            power(hr) = f.m;
            coefficient(hr) = f.b;
            text(xlimit(1),ylimit(2)-hr*(ylimit(2)-ylimit(1))/5,sprintf('y = %0.3g x^{%0.3g}',coefficient(hr),power(hr)),'FontSize',20,'color',c{hr});
        end
        xlabel('Mean Spike Count');
        ylabel('Spike Count Variance');
        str = ['Response reliability' newline sprintf('p_{power}=%.2g, p_{coe}=%.2g',ppow,pcoe)];
        title(str);
        
        
        %         % ANCOVA 协方差分析两个group是否具有显著差异
        %         xx = [mean_fr(:,1);mean_fr(:,2)];
        %         yy = [std_fr(:,1);std_fr(:,2)];
        %         group = [ones(length(mean_fr(:,1)),1); ones(length(mean_fr(:,1)),1)*2];
        %         [~,a] = aoctool((xx),(yy),group,0.05,'','','','off');
        %         group_p_value = a{2,6};
        %         text(xlimit(1),ylimit(2)-3*(ylimit(2)-ylimit(1))/5,sprintf('ANCOVA p = %0.3g',group_p_value),'FontSize',20);
        
        % how to pair t-test to power and coefficient ???
        
        SetFigure(15);
        figN = figN + 1;
    end

    function f1p7(debug) % Max Response in T and R
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        %          methods_of_select = {select_typical, 'Typical cells'}; % 默认T cell+R cell+ Spi cell
        % select all cell
        
        % ---------------------- get data ----------------------------
        maxFR = cell2mat({group_result.maxFR}');
        maxFR = maxFR(select_all,:); % select all cell
        sponFR = [group_result.sponFR]';
        sponFR = sponFR(select_all);
        
        maxFR_group{1} = maxFR;
        maxFR_group{2} = maxFR-sponFR;
        
         % plot
         set(figure(figN),'Position',[100,70,800,800], 'Name', ' maxFR in T and R');clf
         for i = 1:2 % max and max-spon
             subplot(2,2,i)
             plot(maxFR_group{i}(:,1),maxFR_group{i}(:,2),'ko');
             [r,pp] = plot_corr_line(maxFR_group{i}(:,1),maxFR_group{i}(:,2),'MethodOfCorr','pearson','FittingMethod',2,'LineStyles','k-','LineWidth',2);
             SameScale(1);
             xlimit = xlim;
             ylimit = ylim;
             xlim([-20 xlimit(2)]);
             ylim([-20 ylimit(2)]);
             text(0,ylimit(2)-(ylimit(2)-ylimit(1))/20,sprintf('r=%0.3g, p=%0.3g',r,pp));
             xlabel('Translation');
             ylabel('Rotation');
             cellnum = sum(select_all);
             text(0,ylimit(2)-(ylimit(2)-ylimit(1))/5,sprintf('N=%0.3g',cellnum));
             if i == 1
                 title('Max Response');
             elseif i == 2
                 title('Max Response - Spontaneous Response');
             end
         end
         
         % difference
         subplot(2,2,3)
         diffFR = maxFR_group{i}(:,1)-maxFR_group{i}(:,2);
         histogram(diffFR,30);
         meandiff = nanmean(diffFR);
         [~,p] = ttest(diffFR);
         xlimit = xlim;
         ylimit = ylim;
         text((xlimit(2)-xlimit(1))/20+xlimit(1),ylimit(2)-(ylimit(2)-ylimit(1))/20,sprintf('mean = %0.3g',meandiff));
         text((xlimit(2)-xlimit(1))/20+xlimit(1),ylimit(2)-(ylimit(2)-ylimit(1))/7,sprintf('t test, p = %0.3g',p));
         xlabel('T-R');
         ylabel('Cases');
         title('Max Response difference');
         
         SetFigure(13)
         figN = figN + 1;
         
         
         % maxFR seperate cell type, % select all cell
         set(figure(figN),'Position',[950,70,800,800], 'Name', ' maxFR in T and R (cell type)');
         cell_type1 = [group_result.cell_type]';
         cell_type1 =  cell_type1(select_all,:);
         this_celltype{1} = logical(cell_type1 == 1); % translation cell
         this_celltype{2} = logical(cell_type1 == 2); % rotation cell
         this_celltype{3} = logical(cell_type1 == 3); % spiral cell
         this_celltype{4} = logical(cell_type1~=1 & cell_type1~=2 & cell_type1~=3); % non-excitory or no tuning
         
         for ct = 1:4
             subplot(2,2,ct)
             plot(maxFR_group{1}(this_celltype{ct},1),maxFR_group{1}(this_celltype{ct},2),'o','color',c{ct});
             if sum(this_celltype{ct})>5
                 [r,pp] = plot_corr_line(maxFR_group{1}(this_celltype{ct},1),maxFR_group{1}(this_celltype{ct},2),'MethodOfCorr','pearson','FittingMethod',2,'LineStyles',{'-','color',c{ct}},'LineWidth',2);
             end
             SameScale(1);
             xlimit = xlim;
             ylimit = ylim;
             xlim([-20 xlimit(2)]);
             ylim([-20 ylimit(2)]);
             text(0,ylimit(2)-(ylimit(2)-ylimit(1))/20,sprintf('r=%0.3g, p=%0.3g',r,pp));
             cellnum = sum(this_celltype{ct});
             text(0,ylimit(2)-(ylimit(2)-ylimit(1))/5,sprintf('n=%g',cellnum));
             xlabel('Translation');
             ylabel('Rotation');
             if ct == 1
                 title('Translation cell');
             elseif ct == 2
                 title('Rotation cell');
             elseif ct == 3
                 title('Spiral cell');
             else
                 title('Other cell');
             end
         end
         suptitle('max FR');
         SetFigure(13)
         figN = figN + 1;   
    end

    function f2p8(debug) % RF effect in preferred direction
        if debug  ; dbstack;   keyboard;      end
        
        % RF size: 转化为圆的半径比较
        % RF eccentricity: rf中心距离（0，0）距离
        rf_temp =  {group_result.RF}';
        have_rf = ~cellfun(@isempty,rf_temp);
        rf_temp(~have_rf) = {[nan nan nan nan]};
        rf = cell2mat(rf_temp);
        rf_area = rf(:,3) .* rf(:,4);
        rf_dia = 2.* sqrt(rf_area / pi); % 直径
        rf_ecc = sqrt(rf(:,1).^2 + rf(:,2).^2);
        rf_ang = atan(rf(:,2) ./ rf(:,1)) * 180/pi;
        pref_direc = cell2mat({group_result.pref_direc}');
        
        methods_of_select = {
            % select_all, 'All cells';
            select_typical, 'Typical cells';
            % select_no_typical, 'Non-typical cells'
            };
        
        % 极坐标区分，同amp→同rf_ecc，同az→同rf_ang
        thetabinNum = 9; % rf_ang
        RbinNum = 7; % rf_ecc
        
        theta_edge = linspace(-90,90,thetabinNum); % 范围暂时确定 (rf_ang,都在右半视野，-90~90 degree)
        R_edge = linspace(0,50,RbinNum);  % 范围暂时确定 (eccentricity)
        Rrange = max(R_edge) - min(R_edge);
        
        % 标准化R为0-1
        rNorm = (R_edge)/max(R_edge)*max(R_edge/Rrange);
        
        % meshgrid, matrix (RbinNum × thetabinNum)  （网格为RbinNum × thetabinNum，data实际为RbinNum-1 × thetabinNum-1   （两条线之间才有一个data！！））
        YY = (rNorm)'*cosd(theta_edge);
        XX = (rNorm)'*sind(theta_edge);
        
        % 区分RF大小 暂时按照mean分成两组, 只包括有tuning的细胞
        mean_dia = nanmean(rf_dia);
        rf_size{1} = logical(rf_dia <= mean_dia & methods_of_select{1,1});
        rf_size{2} = logical(rf_dia > mean_dia & methods_of_select{1,1});
        rf_size{3} = logical(~isnan(rf_dia) & methods_of_select{1,1});
        
        % heading preference 与 rotation preference 作差
        for r = 1:3 % small rf, big rf, all
            for x = 1:length(theta_edge)-1 % （两条线之间才有一个data！！）
                for y = 1:length(R_edge)-1 % （两条线之间才有一个data！！）
                    loc{r,x,y} = logical(rf_size{r} & rf_ang >=theta_edge(x) & rf_ang<theta_edge(x+1) & rf_ecc>=R_edge(y) & rf_ecc<R_edge(y+1));  % RbinNum-1 × thetabinNum-1 （两条线之间才有一个data！！）
                    diff_pref{r,x,y} =  abs(pref_direc(loc{r,x,y},2) - pref_direc(loc{r,x,y},1));
                    diff_pref{r,x,y}(diff_pref{r,x,y}>180) = 360 -  diff_pref{r,x,y}(diff_pref{r,x,y}>180);
                    mean_diff_pref{r}(x,y) = mean(diff_pref{r,x,y});
                end
            end
        end
        
        set(figure(figN),'Position',[350,100, 1300,550], 'Name', 'RF effect in preference');
        for r = 1:3
            Z = [];
            subplot(1,3,r)
            Z = mean_diff_pref{r}';
            % 伪造最后一行和最后一列
            [last_x,last_y] = size(Z);
            Z(last_x+1,:) = zeros(1,last_y);
            Z(:,last_y+1) = zeros(last_x+1,1);
            h = pcolor(XX,YY,Z);
            axis equal tight
            caxis([0,180]); % 统一colorbar
            colorbar
            % pref : 每一行：不同rf_angle  →  极坐标从-90到+90
            % 每一列：不同rf_ecc  → 极坐标从内往外
            
            view([90 90]);
            set(gca, 'XDir', 'reverse');
            axis off;
            text(-0.97,0.28,'-90');
            text(-0.7,0.8,'-60');
            text(-0.3,1,'-30');
            text(0,1.1,'0');
            text(0.97,0.28,'90');
            text(0.7,0.8,'60');
            text(0.3,1,'30');
            
            text(0.1,0.2,'10');
            text(0.3,0.35,'20');
            text(0.5,0.55,'30');
            text(0.65,0.7,'40');
            
            if r == 1
                title('small RF');
            elseif r==2
                title('big RF');
            else
                title('All');
            end
        end
        suptitle('|T pref - R pref|');
        
        SetFigure(13);
        figN = figN+1;
    end

    function f2p15(debug) %RF effect in preferred direction (local planar motion)
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        
        % ---------------------- get data ----------------------------
        % select above spon and significant tuning cell, translation and
        % rotation tuning
        cell_type = [group_result.cell_type]';
        select_temp = logical(methods_of_select{1,1} & (cell_type==1 | cell_type==2 | cell_type==3));
        
        pref_direc = cell2mat({group_result.pref_direc}');
        
        % RF size: 转化为圆的半径比较
        % RF eccentricity: rf中心距离（0，0）距离
        rf_temp =  {group_result.RF}';
        have_rf = ~cellfun(@isempty,rf_temp);
        rf_temp(~have_rf) = {[nan nan nan nan]};
        rf = cell2mat(rf_temp);
        rf_area = rf(:,3) .* rf(:,4);
        rf_dia = 2.* sqrt(rf_area / pi); % 直径
        rf_ecc = sqrt(rf(:,1).^2 + rf(:,2).^2);
        rf_ang = atan(rf(:,2) ./ rf(:,1)) * 180/pi;
        pref_direc = cell2mat({group_result.pref_direc}');
        
        % 根据prefer direction，还原出optic flow，并且计算RF内的向量和
        % RF 过大时候，向量和并不准确
        % 1. 将RF细分为4个象限，计算各自的向量和再一一比较
        % 2. 只选取RF较小的unit （RF较大的unit，不同h,r的prefer direction的视觉刺激较难一致）
        
        selectCell = find(select_temp & ~isnan(rf_dia) & ~isnan(pref_direc(:,1)));
        select = logical(select_temp & ~isnan(rf_dia) & ~isnan(pref_direc(:,1)));
        
        TVectorSumALL = nan(length(group_result),2); RVectorSumALL = nan(length(group_result),2);
        TVectorSumCenter = nan(length(group_result),2); RVectorSumCenter = nan(length(group_result),2);
        TVectorSumQuad = []; RVectorSumQuad = [];
        
        TRAngleALL = nan(size(selectCell));TRAngleCenter = nan(size(selectCell));TRAngleQuad = nan(size(selectCell));
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        progressbar('Calculate Prefer OpticFlow...');
        caln = 0;
        for n = 555 % 1:length(group_result)
            if select(n)
                caln = caln+1;
                 %             pause(3)
                rf_this = [];
                rf_this = rf(n,:);
                this_tuning_file = group_result(n).SpiT_FILE;
                
                % randomDots_Lwh202010(pref_direc(translation,rotation,firstframe,secondframe);
                % optic flow field范围：
                
                % translation part
                [TFirFrame, TSecFrame] = randomDots_Lwh202010(pref_direc(n,1),888,30,35);
                % rotation part
                [RFirFrame, RSecFrame] = randomDots_Lwh202010(888,pref_direc(n,2),30,35);
                
                Tus = TSecFrame(1,:)-TFirFrame(1,:);
                Tvs = TSecFrame(2,:)-TFirFrame(2,:);
                Rus = RSecFrame(1,:)-RFirFrame(1,:);
                Rvs = RSecFrame(2,:)-RFirFrame(2,:);
                
                % select local flow dot in firstframe in RF
                TselectDotALL = []; RselectDotALL= [];
                TselectDotALL = logical(TFirFrame(1,:)>(rf_this(1)-rf_this(3)/2) & TFirFrame(1,:)<(rf_this(1)+rf_this(3)/2) & TFirFrame(2,:)>(rf_this(2)-rf_this(4)/2) & TFirFrame(2,:)<(rf_this(2)+rf_this(4)/2));
                RselectDotALL = logical(RFirFrame(1,:)>(rf_this(1)-rf_this(3)/2) & RFirFrame(1,:)<(rf_this(1)+rf_this(3)/2) & RFirFrame(2,:)>(rf_this(2)-rf_this(4)/2) & RFirFrame(2,:)<(rf_this(2)+rf_this(4)/2));
                
                % select local flow dot in firstframe in small field of RF:  
                % 1. one middle small field (长宽变小一半)
                TselectDotCenter = [];RselectDotCenter = [];
                TselectDotCenter = logical(TFirFrame(1,:)>(rf_this(1)-rf_this(3)/4) & TFirFrame(1,:)<(rf_this(1)+rf_this(3)/4) & TFirFrame(2,:)>(rf_this(2)-rf_this(4)/4) & TFirFrame(2,:)<(rf_this(2)+rf_this(4)/4));
                RselectDotCenter = logical(RFirFrame(1,:)>(rf_this(1)-rf_this(3)/4) & RFirFrame(1,:)<(rf_this(1)+rf_this(3)/4) & RFirFrame(2,:)>(rf_this(2)-rf_this(4)/4) & RFirFrame(2,:)<(rf_this(2)+rf_this(4)/4));
                
                % 4 quadrant RF
                TselectDotQuad = [];RselectDotQuad=[];
                TselectDotQuad{2} = logical(TFirFrame(1,:)>rf_this(1) & TFirFrame(1,:)<(rf_this(1)+rf_this(3)/2) & TFirFrame(2,:)>rf_this(2) & TFirFrame(2,:)<(rf_this(2)+rf_this(4)/2)); % 第一象限
                RselectDotQuad{2} = logical(RFirFrame(1,:)>rf_this(1) & RFirFrame(1,:)<(rf_this(1)+rf_this(3)/2) & RFirFrame(2,:)>rf_this(2) & RFirFrame(2,:)<(rf_this(2)+rf_this(4)/2));
                
                TselectDotQuad{1} = logical(TFirFrame(1,:)>(rf_this(1)-rf_this(3)/2) & TFirFrame(1,:)<rf_this(1) & TFirFrame(2,:)>rf_this(2) & TFirFrame(2,:)<(rf_this(2)+rf_this(4)/2)); % 第二象限
                RselectDotQuad{1} = logical(RFirFrame(1,:)>(rf_this(1)-rf_this(3)/2) & RFirFrame(1,:)<rf_this(1) & RFirFrame(2,:)>rf_this(2) & RFirFrame(2,:)<(rf_this(2)+rf_this(4)/2));
                
                TselectDotQuad{3} = logical(TFirFrame(1,:)>(rf_this(1)-rf_this(3)/2) & TFirFrame(1,:)<rf_this(1) & TFirFrame(2,:)>(rf_this(2)-rf_this(4)/2) & TFirFrame(2,:)<rf_this(2)); % 第三象限
                RselectDotQuad{3} = logical(RFirFrame(1,:)>(rf_this(1)-rf_this(3)/2) & RFirFrame(1,:)<rf_this(1) & RFirFrame(2,:)>(rf_this(2)-rf_this(4)/2) & RFirFrame(2,:)<rf_this(2));
                
                TselectDotQuad{4} = logical(TFirFrame(1,:)>rf_this(1) & TFirFrame(1,:)<(rf_this(1)+rf_this(3)/2) & TFirFrame(2,:)>(rf_this(2)-rf_this(4)/2) & TFirFrame(2,:)<rf_this(2)); % 第四象限
                RselectDotQuad{4} = logical(RFirFrame(1,:)>rf_this(1) & RFirFrame(1,:)<(rf_this(1)+rf_this(3)/2) & RFirFrame(2,:)>(rf_this(2)-rf_this(4)/2) & RFirFrame(2,:)<rf_this(2));
                
                %% similarity
                % 1. vector sum local optic flow
                % 1.1 overall vector sum
                TVectorSumALL(n,:) = [sum(Tus(TselectDotALL)),sum(Tvs(TselectDotALL))]; % [x,y]
                RVectorSumALL(n,:) = [sum(Rus(RselectDotALL)),sum(Rvs(RselectDotALL))]; % [x,y]
                
                % 1.2 small field vector sum
                TVectorSumCenter(n,:) = [sum(Tus(TselectDotCenter)),sum(Tvs(TselectDotCenter))]; % [x,y]
                RVectorSumCenter(n,:) = [sum(Rus(RselectDotCenter)),sum(Rvs(RselectDotCenter))]; % [x,y]
                
                for j = 1:4 % 四个象限
                    TVectorSumQuad(n,j,:) = [sum(Tus(TselectDotQuad{j})),sum(Tvs(TselectDotQuad{j}))];
                    RVectorSumQuad(n,j,:) = [sum(Rus(RselectDotQuad{j})),sum(Rvs(RselectDotQuad{j}))];
                    TRAngleQuad(n,j) = cossAngle(TVectorSumQuad(n,j,:),RVectorSumQuad(n,j,:));
                end
                
                % 比较Translation与rotation vector sum 夹角
                TRAngleALL(n,1) = cossAngle(TVectorSumALL(n,:),RVectorSumALL(n,:));
                TRAngleCenter(n,1) = cossAngle(TVectorSumCenter(n,:),RVectorSumCenter(n,:));
                
                % Translation, rotation vector sum 与0°(向右)夹角
                T0AngleALL(n,1) = cossAngle(TVectorSumALL(n,:),[1 0]);
                R0AngleALL(n,1) = cossAngle(RVectorSumCenter(n,:),[1 0]);
                
                % 2. 1D orientation distribution
                % 对RF内的点插值，使得T，R的点的起始位置相同
                rf_xlim = [rf_this(1)-rf_this(3)/2 rf_this(1)+rf_this(3)/2];
                rf_ylim = [rf_this(2)-rf_this(4)/2 rf_this(2)+rf_this(4)/2];
                win_xlim = [-45 45];
                win_ylim = [-45 45];
                
                %                 % RF内插值
                %                 figure
                %                 quiver(TFirFrame(1,TselectDotALL),TFirFrame(2,TselectDotALL),Tus(TselectDotALL),Tvs(TselectDotALL));
                %                 [xq,yq] = meshgrid(rf_xlim(1):2:rf_xlim(2), rf_ylim(1):2:rf_ylim(2));
                %                 uq = griddata(TFirFrame(1,TselectDotALL),TFirFrame(2,TselectDotALL),Tus(TselectDotALL),xq,yq);
                %                 vq = griddata(TFirFrame(1,TselectDotALL),TFirFrame(2,TselectDotALL),Tvs(TselectDotALL),xq,yq);
                
                % whole filed内插值
                [xq,yq] = meshgrid(win_xlim(1):2:win_xlim(2), win_ylim(1):2:win_ylim(2)); % same in T and R
                Tuq = griddata(TFirFrame(1,:),TFirFrame(2,:),Tus,xq,yq,'natural');
                Tvq = griddata(TFirFrame(1,:),TFirFrame(2,:),Tvs,xq,yq,'natural');
                
                Ruq = griddata(RFirFrame(1,:),RFirFrame(2,:),Rus,xq,yq,'natural');
                Rvq = griddata(RFirFrame(1,:),RFirFrame(2,:),Rvs,xq,yq,'natural');
                
                % select vector in RF
                select_in_RF = logical(xq>rf_xlim(1) & xq<rf_xlim(2) & yq>rf_ylim(1) & yq<rf_ylim(2));
                RFx = xq(select_in_RF);
                RFy = yq(select_in_RF);
                RFu_T = Tuq(select_in_RF);
                RFv_T = Tvq(select_in_RF);
                RFu_R = Ruq(select_in_RF);
                RFv_R = Rvq(select_in_RF);
                
                % remove Nan
                remove = logical(isnan(RFu_T) | isnan(RFv_T) | isnan(RFu_R) | isnan(RFv_R));
                RFx = RFx(~remove);
                RFy = RFy(~remove);
                RFu_T = RFu_T(~remove);
                RFv_T = RFv_T(~remove);
                RFu_R = RFu_R(~remove);
                RFv_R = RFv_R(~remove);
                
                % group
                T_vec = [RFu_T,RFv_T];
                R_vec = [RFu_R,RFv_R];
                vec{1} = T_vec;
                vec{2} = R_vec;
                
                % zscore
                zscore_all = zscore([RFu_T;RFu_R;RFv_T;RFv_R]);
                RFu_Tz = zscore_all(length(RFu_T)*0+1:length(RFu_T)*1);
                RFu_Rz = zscore_all(length(RFu_T)*1+1:length(RFu_T)*2);
                RFv_Tz = zscore_all(length(RFu_T)*2+1:length(RFu_T)*3);
                RFv_Rz = zscore_all(length(RFu_T)*3+1:length(RFu_T)*4);
                vecZscore{1} = [RFu_Tz RFv_Tz];
                vecZscore{2} = [RFu_Rz RFv_Rz];
                
                
                [Ttheta,Trho] = cart2pol(RFu_T,RFv_T);
                [Rtheta,Rrho] = cart2pol(RFu_R,RFv_R);
                Ttheta = rad2deg(Ttheta);
                Rtheta = rad2deg(Rtheta);
                
                binedges = linspace(-180,180,19);
                
%                 figure
%                 subplot(1,2,1)
%                 histogram(Ttheta,binedges);
%                 xlim([-180 180]);
%                 xlabel('vector orientation');
%                 ylabel('Cases');
%                 subplot(1,2,2)
%                 histogram(Rtheta,binedges);
%                 xlim([-180 180]);

                % 3. cosin similarity
                D_cos = 1 - pdist2(T_vec,R_vec,'cosine');
                meanD_cos = mean(mean(D_cos));
                % Adjusted Cosine Similarity
                
                % 4.correlation distance
                D_cor = pdist2(T_vec,R_vec,'correlation');
                meanD_cor = mean(mean(D_cor));
                
                % 5. 2-D correlation coefficient
                corcoe2D = corr2(T_vec,R_vec);
                
                %                 % 6. SPAEF, The neural basis of predictive pursuit, Seng Bum Michael Yoo, 2020
                %                 % SPAEF=-∞: anticorrelated， 0：uncorrelated  1：perfect matching
                %                 % SPAEF 太过于敏感？相似性都偏低  SPAEF貌似只针对n*n*1数据
                %
                %                 % pearson correlation coefficient
                %                 AA_temp = corrcoef(T_vec,R_vec);
                %                 AA = AA_temp(1,2);
                %                 % fraction of the coefficient of variation
                %                 BB = (std2(R_vec)/mean2(R_vec)) / (std2(T_vec)/mean2(T_vec));
                %                 % activity silimarity measured by histogram profiles: histogram intersection of the T histogram (K) and the R histogram (L), each containing n bins
                %                 
                %
                %                 binedges = linspace(-3,3,21);
                % %                 figure
                % %                 subplot(1,2,1)
                % %                 histogram2(RFu_Tz,RFv_Tz,binedges,binedges);
                % %                 xlabel('u');ylabel('v');
                % %                 subplot(1,2,2)
                % %                 histogram2(RFu_Rz,RFv_Rz,binedges,binedges);
                %
                %                 Thist = histcounts2(RFu_Tz,RFv_Tz,binedges,binedges);
                %                 Rhist = histcounts2(RFu_Rz,RFv_Rz,binedges,binedges);
                %                 Thist = Thist(:);
                %                 Rhist = Rhist(:);
                %                 CC = sum(min(Thist,Rhist))/sum(Thist);
                %                 SPAEF = 1 - sqrt((AA-1)^2 + (BB-1)^2 + (CC-1)^2);
                
                % 7. global distributions, Dinh and Xu, 2008. Measuring the Similarity of Vector Fields Using Global Distributions
                % insensitive to Rotation。。。
                % randomly selecting 10000 pairs of points, computing the euclidean distance and dot product
                pairN = 10000;
                xedges = linspace(0,10,101);
                yedges = linspace(-1,1,101);
                xmid = xedges+(xedges(2)-xedges(1))/2;
                xmid = xmid(1:end-1);
                ymid = yedges+(yedges(2)-yedges(1))/2;
                ymid = ymid(1:end-1);
                
%                 figure(1002);clf
                for tr = 1:2 % T and R
                    firVec = randi(length(T_vec),[10000,1]);
                    secVec = randi(length(T_vec),[10000,1]);
                    
                    for i = 1:pairN
                        D_euc{tr}(i,1) = pdist2(vec{tr}(firVec(i),:), vec{tr}(secVec(i),:));
                        
                        firnorm = vec{tr}(firVec(i),:)/norm(vec{tr}(firVec(i),:));
                        secnorm = vec{tr}(secVec(i),:)/norm(vec{tr}(secVec(i),:));
                        dot_product{tr}(i,1) = dot(firnorm, secnorm); % 只取方向
                    end
                    
                    
                    histN{tr} = histcounts2(D_euc{tr}, dot_product{tr},xedges,yedges,'Normalization','probability'); % normalization, sum=1
%                     %                     subplot(2,2,tr)
%                     %                     histogram2(D_euc{tr}, dot_product{tr},xedges,yedges);
%                     %                     subplot(2,2,tr+2)
%                     subplot(1,2,tr)
%                     imagesc([xmid(1) xmid(end)],[ymid(1) ymid(end)],histN{tr}');
%                     set(gca,'YDir','normal');
%                     xlabel('Euclidean distance');
%                     ylabel('dot product');
                end
                
                
                histN1D = cellfun(@(x) x(:), histN, 'UniformOutput', 0); % a按一列存储
                
                % chi square distance between two histogram
                % https://stats.stackexchange.com/questions/184101/comparing-two-histograms-using-chi-square-distance
                % Σ((xi-yi)^2/(xi+yi))
                D_chisq = 0;
                for i = 1:length(histN1D{1})
                    if histN1D{1}(i)+histN1D{2}(i)~=0 % define 0/0=0
                        D_chisq = D_chisq + (histN1D{1}(i)-histN1D{2}(i))^2 / (histN1D{1}(i)+histN1D{2}(i));
                    end
                end
                           
                %% plot optic flow
                figure(1001)
                set(figure(1001),'Position',[100,50,900,900], 'Name', '');clf
                subplot(2,2,1)
                quiver(TFirFrame(1,:),TFirFrame(2,:),Tus,Tvs,3); % raw, random dot
                hold on
                quiver(xq,yq,Tuq,Tvq); % 插值result
                rectangle('position',[rf_this(1)-rf_this(3)/2 rf_this(2)-rf_this(4)/2, rf_this(3) rf_this(4)],...
                    'Curvature',[0.3 0.3],'EdgeColor',[0.6 0.6 0.6],'LineWidth',1); % show RF
                text(-40,47,sprintf('Preferred Direction = %g',pref_direc(n,1)),'color','r');
                xlim([-45 45]);
                ylim([-45 45]);
                axis square
                % vector sum
                % normalize TVectorSumALL to unit length
                [theta,~] = cart2pol(TVectorSumALL(n,1),TVectorSumALL(n,2));
                [x,y] = pol2cart(theta,25);
                quiver(rf_this(1),rf_this(2),x,y);
                title(this_tuning_file);

                subplot(2,2,2)
                quiver(RFirFrame(1,:),RFirFrame(2,:),Rus,Rvs,3);
                hold on
                quiver(xq,yq,Ruq,Rvq); % 插值
                rectangle('position',[rf_this(1)-rf_this(3)/2 rf_this(2)-rf_this(4)/2, rf_this(3) rf_this(4)],...
                    'Curvature',[0.3 0.3],'EdgeColor',[0.6 0.6 0.6],'LineWidth',1);
                text(-40,47,sprintf('Preferred Direction = %g',pref_direc(n,2)),'color','r');
                xlim([-45 45]);
                ylim([-45 45]);
                axis square
                % vector sum
                [theta,~] = cart2pol(RVectorSumALL(n,1),RVectorSumALL(n,2));
                [x,y] = pol2cart(theta,25);
                quiver(rf_this(1),rf_this(2),x,y);
                title(this_tuning_file);
                
                % In RF
                subplot(2,2,3)
                quiver(RFx,RFy,RFu_T,RFv_T); % 插值result
                xlim(rf_xlim);
                ylim(rf_ylim);
                subplot(2,2,4)
                quiver(RFx,RFy,RFu_R,RFv_R); % 插值result
                xlim(rf_xlim);
                ylim(rf_ylim);
                
                % anotation
                xlimit = xlim;
                ylimit = ylim;
                str = [sprintf('mean cosine distance = %.3g',meanD_cos) newline ...
                    sprintf('mean correlatioin distance = %.3g',meanD_cor) newline ...
                    sprintf('2-D correlation coefficient = %.3g',corcoe2D) newline ...
                    sprintf('chi-square distance = %.3g',D_chisq)];
                text(xlimit(1),ylimit(2)+(ylimit(2)-ylimit(1))/8,str);
                
                % save figure
                savepath = 'Z:\Data\Tempo\Batch\opticflow_in_RF\';
                filename = group_result(n).SpiT_FILE;
                fig_filename = strcat(filename(1:end-4),'.bmp');
                fig_savepath = strcat(savepath,fig_filename);
                saveas(figure(1001),fig_savepath);
                
                % save data to PC
                RFxy = [RFx RFy];
                RFTvec = [RFu_T RFv_T];
                RFRvec = [RFu_R RFv_R];
                data_filename = strcat(filename(1:end-4),'RFxy.mat');
                data_savepath = strcat(savepath,data_filename);
                save(data_savepath,'RFxy','RFTvec','RFRvec');
             
                % save to group_result
                group_result(n).TVectorSumALL = TVectorSumALL(n,:);
                group_result(n).RVectorSumALL = RVectorSumALL(n,:);
                group_result(n).TVectorSumMidd = TVectorSumCenter(n,:);
                group_result(n).RVectorSumMidd = RVectorSumCenter(n,:);
                group_result(n).TRAngleMidd = TRAngleCenter(n,1);
                group_result(n).TRAngleALL = TRAngleALL(n,1);
                
                % raw vector field
                group_result(n).TFirFrame = TFirFrame;
                group_result(n).TSecFrame = TSecFrame;
                group_result(n).RFirFrame = RFirFrame;
                group_result(n).RSecFrame = RSecFrame;
                group_result(n).Tvec = [Tus Tvs];
                group_result(n).Rvec = [Rus Rvs];
                
                % whole filed内插值后RF内vector field
                group_result(n).RFxy = [RFx RFy];
                group_result(n).RFTvec = [RFu_T RFv_T];
                group_result(n).RFRvec = [RFu_R RFv_R];
                
                % similarity
                group_result(n).D_cos = D_cos; % cosine distance
                group_result(n).meanD_cos = meanD_cos; 
                group_result(n).D_cor = D_cor; % correlation distance
                group_result(n).meanD_cor = meanD_cor; 
                group_result(n).corcoe2D = corcoe2D; % 2-D correlation coefficient
                group_result(n).histN = histN; % global distribution
                group_result(n).D_chisq = D_chisq; % global distribution chi-square distance

                progressbar(caln/length(selectCell));
            else
                group_result(n).TVectorSumALL = [nan nan];
                group_result(n).RVectorSumALL = [nan nan];
                group_result(n).TVectorSumMidd = [nan nan];
                group_result(n).RVectorSumMidd = [nan nan];
                group_result(n).TRAngleMidd = nan;
                group_result(n).TRAngleALL = nan;
                group_result(n).TFirFrame = nan;
                group_result(n).TSecFrame = nan;
                group_result(n).RFirFrame = nan;
                group_result(n).RSecFrame = nan;
                group_result(n).Tvec = [nan nan];
                group_result(n).Rvec = [nan nan];
                group_result(n).RFxy = [nan nan];
                group_result(n).RFTvec = [nan nan];
                group_result(n).RFRvec = [nan nan];
                group_result(n).D_cos = nan; % cosine distance
                group_result(n).meanD_cos = nan; 
                group_result(n).D_cor = nan; % correlation distance
                group_result(n).meanD_cor = nan; 
                group_result(n).corcoe2D = nan; % 2-D correlation coefficient
                group_result(n).histN = nan;
                group_result(n).D_chisq = nan; % global distribution chi-square distance
            end
%             progressbar(caln/length(select));
        end
        
        cd('D:\user\desktop\test202010');
        save('RVectorSumALL.mat','RVectorSumALL');
        save('RVectorSumCenter.mat','RVectorSumCenter');
        save('RVectorSumQuad.mat','RVectorSumQuad');
        save('TRAngleALL.mat','TRAngleALL');
        save('TRAngleCenter.mat','TRAngleCenter');
        save('TRAngleQuad.mat','TRAngleQuad');
        save('TVectorSumALL.mat','TVectorSumALL');
        save('TVectorSumCenter.mat','TVectorSumCenter');
        save('TVectorSumQuad.mat','TVectorSumQuad');
        
        save('T0AngleALL.mat','T0AngleALL');
        save('R0AngleALL.mat','R0AngleALL');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
% 
%         cd('Z:\BaiduNetdiskWorkspace\data\Batch\opticflow_in_RF');
%         namelist = dir('Z:\BaiduNetdiskWorkspace\data\Batch\opticflow_in_RF');
%         len = length(namelist);
%         RFxy = [];
%         RFTvec = [];
%         RFRvec = [];
%         
%         
%         for ii = 1:len
%             try
%                 file_name11{ii}=namelist(ii).name;
%                 if file_name11{ii}(1)=='m'
%                     RFxy{ii} = load(file_name11{ii});
%                 end
%             catch
%                 keyboard
%             end
%         end
%         
%         % 插值适得vector个数一样，与RF大小无关。
%         
%         
%         TRAngleALL(TRAngleALL== 0.00000) = nan;
%         
%         % ------------------------ Plot ------------------------------
%         
%         edges = [0:22.5:180];
%         set(figure(figN),'Position',[100,100 700,300], 'Name', '');
%         subplot(1,2,1)
%         histogram(TRAngleALL,edges);
%         set(gca,'xtick',[0:45:180]);
%         title('RF vector-sum');
%         
%         subplot(1,2,2)
%         histogram(TRAngleCenter,edges);
%         set(gca,'xtick',[0:45:180]);
%         title('middle RF vector-sum');
%         figN = figN+1;
%         
%         set(figure(figN),'Position',[800,100 700,700], 'Name', '');
%         for j = 1:4
%             subplot(2,2,j)
%             histogram(TRAngleQuad(:,j),edges);
%             set(gca,'xtick',[0:45:180]);
%         end
%         suptitle('quadrant RF vector-sum');
%         
%         % 与0°(向右）夹角
%         % Translation, rotation vector sum 与0°(向右)夹角
%         T0AngleALL = cossAngle(TVectorSumALL,[1 0]);
%         R0AngleALL = cossAngle(RVectorSumCenter,[1 0]);
%         
%         figure
%         plot(T0AngleALL,R0AngleALL,'k.');
%         xlabel('T0AngleALL');
%         ylabel('R0AngleALL')
%         
%         figure % 与RF位置，大小关系
%         subplot(1,3,1)
%         plot(rf_dia(selectCell),TRAngleALL(selectCell),'k.');
%         title('RF diameter')
%         
%         subplot(1,3,2)
%         plot(rf_ecc(selectCell),TRAngleALL(selectCell),'k.');
%         title('RF eccentricity')
%         
%         subplot(1,3,3)
%         plot(rf_ang(selectCell),TRAngleALL(selectCell),'k.');
%         title('RF angle')
%         
%         % partial correlation
%         Y = [];
%         Y(:,1) = TRAngleALL(selectCell);
%         Y(:,2) = rf_dia(selectCell);
%         Y(:,3) = rf_ecc(selectCell);
%         Y(:,4) = rf_ang(selectCell);
%         
%         [r,p] = partialcorr(Y,'Rows','pairwise');
%         partial_corr_coef = r(1,2:4);
%         partial_corr_p = p(1,2:4);
%         
%         % regress
%         X =[ones(size(Y(:,1))),Y(:,2),Y(:,3),Y(:,4) Y(:,2).*Y(:,3) Y(:,2).*Y(:,4) Y(:,3).*Y(:,4) Y(:,2).^2];
%         [b,bint,r,rint,stats] = regress(Y(:,1),X); % stats: R^2, F, p,误差方差
%         R2 = stats(1);
%         
%         % FR of prefer direction
%         fitOK = cell2mat({group_result.fitOK}');
%         prefFRGauss = cell2mat(cellfun(@squeeze,{group_result.prefFRGauss},'uniformoutput',0))';
%         select_fit = logical(fitOK(:,1) == 1 & fitOK(:,2) == 1 & methods_of_select{1,1} & select);
%         
%         global_pref_p = cell2mat({group_result.global_pref_p}'); % anova for all direction (h and r)
%         
%         %         GaussfitRsquare = cell2mat({group_result.GaussfitRsquare}');
%         %         GaussfitRsquare = GaussfitRsquare(methods_of_select{1,1},:);
%         
%         % 将RF内optic flow vector分为3大类：TR相似<45, TR垂直45~135，TR相反>135
%         RFtype = nan(size(TRAngleALL));
%         RFtype(TRAngleALL<=45) = 1; % similar
%         RFtype(TRAngleALL>45 & TRAngleALL<=135) = 2; % vertical
%         RFtype(TRAngleALL>135) = 3; % opposite
%         
%         % 不同RF type中FR在Translation与rotation间的差别 → 是否只对局部的optic flow响应 
%         figure % 两种刺激在prefer时候产生的FR，于RF内刺激的一致性的关系  （如果FR一致，RF内刺激也一致，说明细胞只对RF内刺激响应）
%         r = [];p = [];
%         for ii = 1:3
%             TprefFR = prefFRGauss(select_fit & RFtype==ii,1);
%             RprefFR = prefFRGauss(select_fit & RFtype==ii,2);
%             
%             subplot(2,3,ii)
%             plot(TprefFR,RprefFR,'k.');
%             xlimit = xlim;
%             ylimit = ylim;
%             limitmax = max(xlimit(2),ylimit(2));
%             hold on
%             plot([0 limitmax],[0 limitmax],'r--');
%             xlim([0 limitmax]);
%             ylim([0 limitmax]);
%             xlabel('FR in T preferred');
%             ylabel('FR in R preferred');
%             
%             [r,p] = corrcoef(TprefFR,RprefFR);
%             
%             % different distribution
%             diffFR = TprefFR - RprefFR;
%             subplot(2,3,ii+3)
%             histogram(diffFR);
%             title('Different FR');
%         end
    end

    function f2p19(debug) % 'Result RF effect in preferred direction (local planar motion)'
        if debug  ; dbstack;   keyboard;      end
        try
            pathname = 'Z:\BaiduNetdiskWorkspace\data\Batch\opticflow_in_RF';
            cd(pathname);
        catch
            pathname = 'Z:\Data\Tempo\Batch\opticflow_in_RF\mat_file';
            cd(pathname);
        end
        
        matfile = dir([pathname,'\*.mat']); %所要读取文件的类型
        dircell = struct2cell(matfile);
        ori_name = dircell(1,:);
        temp_name = sort_nat(ori_name)';  % 文件名按照自然顺序排序
        
        filenum_this = length(temp_name); % this folder file number
        
        If_figure = 0;
        
        RF = cell2mat({group_result.RF}');
        
        progressbar('');
        for n = 1:length(group_result)
            this_tuning_file = group_result(n).SpiT_FILE;
            if ~isnan(this_tuning_file)
                
%                 % example neuron: n = 555, m7c1394r2
%                 if str2double(this_tuning_file(find(this_tuning_file == 'c')+1:find(this_tuning_file == 'r')-1))==1394
%                     If_figure = 1;
%                 else
%                     If_figure = 0;
%                 end
                
                this_RFfile_name = [this_tuning_file(1:end-4) 'RFxy.mat'];
                RF_fileID = strcmp(this_RFfile_name,temp_name);
                if sum(RF_fileID)==1
                    RFxy = []; RFTvec = []; RFRvec = [];
                    vector_data = load(char(this_RFfile_name));
                    RFxy = vector_data.RFxy;
                    RFTvec = vector_data.RFTvec;
                    RFRvec = vector_data.RFRvec;
                    
                    % 为了使得similarity index具有可比性，插值适得vector个数一样，与RF大小无关。是否有必要？
                    xbinnum = 21;
                    ybinnum = 21;
                    xbin = linspace(min(RFxy(:,1)),max(RFxy(:,1)),xbinnum);
                    ybin = linspace(min(RFxy(:,2)),max(RFxy(:,2)),ybinnum);
                    [xq, yq] = meshgrid(xbin,ybin); % same in T and R
                    RFxy_norm(:,1) = xq(:);
                    RFxy_norm(:,2) = yq(:);
                    
                    RFTvec_norm(:,1) = griddata(RFxy(:,1),RFxy(:,2),RFTvec(:,1),xq(:),yq(:),'natural');
                    RFTvec_norm(:,2) = griddata(RFxy(:,1),RFxy(:,2),RFTvec(:,2),xq(:),yq(:),'natural');
                    RFRvec_norm(:,1) = griddata(RFxy(:,1),RFxy(:,2),RFRvec(:,1),xq(:),yq(:),'natural');
                    RFRvec_norm(:,2) = griddata(RFxy(:,1),RFxy(:,2),RFRvec(:,2),xq(:),yq(:),'natural');
                    
                    Tuu = reshape(RFTvec_norm(:,1),size(xq));
                    Tvv = reshape(RFTvec_norm(:,2),size(xq));
                    Ruu = reshape(RFRvec_norm(:,1),size(xq));
                    Rvv = reshape(RFRvec_norm(:,2),size(xq));
                    
                    uu{1} = Tuu;
                    vv{1} = Tvv;
                    uu{2} = Ruu;
                    vv{2} = Rvv;
                    
                    if If_figure
                        rf_this = RF(n,:);
  
                        figure(figN);clf
                        subplot(2,2,1)
                        quiver(RFxy(:,1),RFxy(:,2),RFTvec(:,1),RFTvec(:,2)); % raw
                        hold on
%                         quiver(xq(:),yq(:),Tuu(:),Tvv(:)); % interplot
                        title('Translation');
                        axis tight
                        
                        subplot(2,2,2)
                        quiver(RFxy(:,1),RFxy(:,2),RFRvec(:,1),RFRvec(:,2)); % raw
                        hold on
%                         quiver(xq(:),yq(:),Ruu(:),Rvv(:)); % interplot
                        title('Rotation');
                        axis tight
                        
                        figure(figN);
                        subplot(2,2,3)
                        hold on
                        cla; even_stream_arrow(xq, yq, Tuu, Tvv,0.2,2, 'LineStyle', '-', 'LineWidth', 1, 'Color', 'k', 'ArrowLength', 5, 'ArrowTipAngle', 45, 'ArrowBaseAngle', 20, 'ArrowDensity', 0.3);
                        hold on
                        rectangle('position',[rf_this(1)-rf_this(3)/2 rf_this(2)-rf_this(4)/2, rf_this(3) rf_this(4)],...
                            'Curvature',[0.3 0.3],'EdgeColor',[0.6 0.6 0.6],'LineWidth',1);
                        xlim([-45 45]);
                        ylim([-45 45]);
                        axis square
                        
                        subplot(2,2,4)
                        hold on
                        cla; even_stream_arrow(xq, yq, Ruu, Rvv,0.2,2, 'LineStyle', '-', 'LineWidth', 1, 'Color', 'k', 'ArrowLength', 5, 'ArrowTipAngle', 45, 'ArrowBaseAngle', 20, 'ArrowDensity', 0.3);
                        hold on
                        rectangle('position',[rf_this(1)-rf_this(3)/2 rf_this(2)-rf_this(4)/2, rf_this(3) rf_this(4)],...
                            'Curvature',[0.3 0.3],'EdgeColor',[0.6 0.6 0.6],'LineWidth',1);
                        xlim([-45 45]);
                        ylim([-45 45]);
                        axis square
                    end
                    
                    % slide window to get small field vector-sum
                    slide_wid = 5;
                    slide_step = 2;
                    j=1;
                    while slide_wid+(j-1)*slide_step <= xbinnum
                        k=1;
                        while slide_wid+(k-1)*slide_step <= ybinnum
                            select_x = 1+(j-1)*slide_step : slide_wid+(j-1)*slide_step;
                            select_y = 1+(k-1)*slide_step : slide_wid+(k-1)*slide_step;
                            
                            select = false(size(xq));
                            select(select_x,select_y) = true;
                            
                            for tr = 1:2
                                vector_sum{tr}(1) = sum(uu{tr}(select));
                                vector_sum{tr}(2) = sum(vv{tr}(select));
                                vector_sum_norm{tr} = vector_sum{tr}/norm(vector_sum{tr});
                            end
                            
%                             angle(j,k) = acosd(dot(vector_sum{1},vector_sum{2})/(norm(vector_sum{1})*norm(vector_sum{2}))); % 注意这里的y轴和figure是上下颠倒的
                            angle(j,k) = cossAngle(vector_sum{1},vector_sum{2});
                            
%                             % plot vectr result
%                             select_center = [mean(xq(select)) mean(yq(select))];
%                             figure(figN);
%                             subplot(2,2,1); hold on
%                             quiver(select_center(1), select_center(2),vector_sum{1}(1),vector_sum{1}(2),0.2,'color','k'); 
%                             subplot(2,2,2); hold on
%                             quiver(select_center(1), select_center(2),vector_sum{2}(1),vector_sum{2}(2),0.2,'color','k'); 
                            
                            k = k + 1;
                        end
                        j = j + 1;
                    end
                    
                    % similar vector-sum: angle smaller than 45
                    group_result(n).RFsimilarity = sum(sum(abs(angle)<45)) / (size(angle,1)*size(angle,2)) * 100;
                    if If_figure
                        str = [this_tuning_file newline sprintf('similarity = %.3g', group_result(n).RFsimilarity)];
                        suptitle(str);
                    end
                    
                    % save figure
                    if If_figure
                        try
                            savefold_all_temp = dir('Z:\Data\Tempo\Batch\opticflow_in_RF\streamline_InRF\');
                            savefold_all = {savefold_all_temp.name}';
                            savefold_all(cellfun(@(x) isequal(x,'.'),savefold_all,'UniformOutput',true)) = [];
                            savefold_all(cellfun(@(x) isequal(x,'..'),savefold_all,'UniformOutput',true)) = [];
                            
                            savepath_temp = 'Z:\Data\Tempo\Batch\opticflow_in_RF\streamline_InRF\';
                            save_location = discretize(group_result(n).RFsimilarity,linspace(0,100,11)); % save location
                            
                            savepath = [savepath_temp savefold_all{save_location}];
                            savename = num2str(group_result(n).RFsimilarity);
                            save_bmp = []; save_fig = [];
                            save_bmp = strcat(savepath,'\',savename,'_',this_tuning_file(1:end-4),'.bmp');
                            save_fig = strcat(savepath,'\',savename,'_',this_tuning_file(1:end-4),'.fig');
                            saveas(figure(figN),save_bmp);
                            saveas(figure(figN),save_fig);
                        catch
                            keyboard
                        end
                    end

                    % 计数
                    progressbar(n/length(group_result));
                    
                else
                    group_result(n).RFsimilarity = nan;
                end
            else
                group_result(n).RFsimilarity = nan;
            end
        end
        disp('Finishih');
        figN = figN+1;
    end

    function f2p20(debug) % 'Prefer OpticFlow similarity vs RF'
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        
        try
            RFsimilarity = [group_result.RFsimilarity]';
        catch
            run_this = strcmpi({function_handles_spi{2,2}{:,1}}','Save Prefer OpticFlow similarity (local planar motion)');
            feval(function_handles_spi{2,2}{run_this,2},false); % first run "Spiral Tuning" - "Spiral Tuning" - "Save Prefer OpticFlow (local planar motion)"
            RFsimilarity = [group_result.RFsimilarity]';
        end
        RFsimilarity = RFsimilarity(methods_of_select{1,1});
        
        % RF: X Y H W (X,Y是RF中心点位置)
        rf_temp =  {group_result.RF}';
        rf = cell2mat(rf_temp); % RF: X Y H W (X,Y是RF中心点位置)
        rf = rf(methods_of_select{1,1},:);
        rf_area = rf(:,3) .* rf(:,4);
        rf_dia = 2.* sqrt(rf_area / pi); % 直径
        rf_ecc = sqrt(rf(:,1).^2 + rf(:,2).^2); % RF中心点距离原点的距离
        rf_ang = atand(rf(:,2) ./ rf(:,1)); % rf_ang>0 中心点在上半视野，rf_ang<0 中心点在下半视野
        rf_ang = abs(rf_ang);
        
        % included fovea or not
        criti_value = 0; % in visual degree
        includefovea = false(length(rf),1);
        for d = 1:length(rf)
            if (rf(d,1)-rf(d,3)/2) < -criti_value && (rf(d,1)+rf(d,3)/2) > criti_value && (rf(d,2)-rf(d,4)/2) < -criti_value && (rf(d,2)+rf(d,4)/2) > criti_value
                includefovea(d) = true;
            end
        end
        
        % plot        
        set(figure(figN),'Position',[350,100, 600,800], 'Name', 'Prefer OpticFlow Similarity between T and R');clf
        for i = 1:3 % all, include Fovea, not include
            ax(i) = subplot(3,1,i);
            if i == 1
                similarity = RFsimilarity;
                
            elseif i == 2 % seperate fovea contatining and fovea excluded sites
                similarity = RFsimilarity(includefovea);
                
            else % seperate fovea contatining and fovea excluded sites
                similarity = RFsimilarity(~includefovea);
                title('fovea-excluded');
            end
            similarity(isnan(similarity)) = [];
            if i == 1
                histogram(similarity,0:5:100,'facecolor','w');
                title('ALL');
            elseif i == 2
                histogram(similarity,0:5:100,'DisplayStyle','stairs');
                title('fovea-containing');
            else
                histogram(similarity,0:5:100,'DisplayStyle','stairs');
                title('fovea-excluded');
            end
            xlabel('Similarity (%)');
            ylabel('Cases');
            ylimit = ylim;
            
            % 1. KS test, uniformity test
            [~,p_ks] = kstest(similarity,[similarity,unifcdf(similarity)]); % h=1,p<0.05 不是均匀分布
            % 2. modality test   p>0.05为符合该mode
            [~, p_mod_in] = modalityTestForYong(similarity,[1 2],1000);
            ylimit = ylim;
            str = [sprintf('n = %g',length(similarity)) newline ...
                sprintf('uniformity test (KS test), p = %.2g',p_ks) newline...
                sprintf('modality test (single peak), p = %.2g',p_mod_in(1))];
            text(5,ylimit(2)-ylimit(2)/10,str);
        end
        
        linkaxes([ax(1),ax(2)]);
        linkaxes([ax(1),ax(3)]);
        suptitle('Prefer OpticFlow Similarity between T and R');

        SetFigure(10);
        figN = figN+1;
        
        % RF ecc- RF size - RF similarity
        set(figure(figN),'Position',[600,100, 1200,800], 'Name', '');clf
        subplot(2,3,1)
        scatter(rf_ecc, rf_dia/2, mapminmax(RFsimilarity',20,100), RFsimilarity,'filled');
        hold on
%         scatter(rf_ecc(includefovea), rf_dia(includefovea)/2, mapminmax(RFsimilarity(includefovea)',1,100), RFsimilarity(includefovea),'filled');
        
        xlim([0 55]);
        ylim([0 55]);
        xlabel('RF eccencity');
        ylabel('RF equivalent radius');
        SameScale(1,1); % above this line: contain fovea; below line: exclude fovea
        text(30,40,'Fovea containing');
        text(40,35,'Fovea excluding');
        
        
        axis square

        color_map_RGB = colormap_set([0 0 0; 1 0 0]);
        colormap(color_map_RGB);
        colorbar
        
        hold on
        scatter(repmat(50,1,10), linspace(40,55,10), linspace(20,100,10), [1:10:100],'filled');
        text(50,38,'Color bar','horizontalalign','center');
        
%         % group 3 similarity groups
%         group_range = linspace(0,100,4);
%         for g = 1:3
%             group_this{g} = logical(RFsimilarity > group_range(g) & RFsimilarity <= group_range(g+1));
%         end
%         figure
%         for g = 1:3
%             subplot(1,3,g)
%             scatter(rf_ecc(group_this{g}), rf_dia(group_this{g}),'r');
%             xlim([0 50]);
%             ylim([10 110]);
%             xlabel('RF eccencity');
%             ylabel('RF diameter');
%             axis square
%             
%             hold on
%             % plot diameter/2 > ecc line (containing fovea)
%             plot([0 50],[0 125],'k-'); % above this line: contain fovea; below line: exclude fovea
%             
%              if g == 1
%                 title('RFsimilarity = 0~33 (different)','horizontalAlignment','center');
%             elseif g == 2
%                 title('RFsimilarity = 33~66 (medium)','horizontalAlignment','center');
%             else
%                 title('RFsimilarity = 66~100 (similar)','horizontalAlignment','center');
%             end
%         end

        % similarity vs RF eccencity
        subplot(2,3,2)
        plot(rf_ecc,RFsimilarity,'k.');
        xlabel('RF eccencity');
        ylabel('Similarity');
        [xmid, y, y_sem] = get_slideWin(rf_ecc,RFsimilarity,10,5,[0 50]);
        hold on
        errorbar(xmid,y,y_sem);
        [r,pp] = plot_corr_line(rf_ecc,RFsimilarity,'MethodOfCorr','spearman','FittingMethod',2,'LineWidth',2); % PCA fitting
        axis square
        xlim([0 50]);
        ylimit = ylim;
        xlimit = xlim;
        text(xlimit(1),ylimit(2),sprintf('r = %.2g, p = %.2g',r,pp),'color','r');
        
        % similarity vs RF size
        subplot(2,3,3)
        plot(rf_dia,RFsimilarity,'k.');
        xlabel('RF diameter');
        ylabel('Similarity');
        [xmid, y, y_sem] = get_slideWin(rf_dia,RFsimilarity,20,10,[10 110]);
        hold on
        errorbar(xmid,y,y_sem);
        [r,pp] = plot_corr_line(rf_dia,RFsimilarity,'MethodOfCorr','spearman','FittingMethod',2,'LineWidth',2); % PCA fitting
        axis square
        xlim([10 110]);
        ylimit = ylim;
        xlimit = xlim;
        text(xlimit(1),ylimit(2),sprintf('r = %.2g, p = %.2g',r,pp),'color','r');
        
        % similarity vs RF size & eccencity
        subplot(2,3,4)
        plot(rf_ecc./rf_dia,RFsimilarity,'k.');  % rf_ecc./rf_dia
        xlabel('RF ecc/diameter');
        ylabel('Similarity');
        % bin
        [xmid, y, y_sem] = get_slideWin(rf_ecc./rf_dia,RFsimilarity,0.2,0.1,[0 1.4]);
        hold on
        errorbar(xmid,y,y_sem);
        [r,pp] = plot_corr_line(rf_ecc./rf_dia,RFsimilarity,'MethodOfCorr','spearman','FittingMethod',2,'LineWidth',2); % PCA fitting
        plot_corr_line(rf_ecc./rf_dia,RFsimilarity,'MethodOfCorr','spearman','FittingMethod',1,'LineWidth',1); % linear fitting
        axis square
        xlim([0 1.4]);
        ylimit = ylim;
        xlimit = xlim;
        text(xlimit(1),ylimit(2),sprintf('r = %.2g, p = %.2g',r,pp),'color','r');
        
        % similarity vs RF angle
        subplot(2,3,5)
        plot(rf_ang,RFsimilarity,'k.');
        xlabel('RF angle');
        ylabel('Similarity');
        % bin
        [xmid, y, y_sem] = get_slideWin(rf_ang,RFsimilarity,20,10,[0 90]);
        hold on
        errorbar(xmid,y,y_sem);
        [r,pp] = plot_corr_line(rf_ang,RFsimilarity,'MethodOfCorr','spearman','FittingMethod',2,'LineWidth',2); % PCA fitting
        axis square
        xlim([0 90]);
        ylimit = ylim;
        xlimit = xlim;
        text(xlimit(1),ylimit(2),sprintf('r = %.2g, p = %.2g',r,pp),'color','r');
        SetFigure(15);
        figN = figN + 1;
        
        % 回归分析
        X = [];
        X = [ones(size(RFsimilarity)),rf_ang,rf_ecc,rf_dia,...
            rf_ang.*rf_ecc,rf_ang.*rf_dia,rf_ecc.*rf_dia,...
            rf_ang.*rf_ecc.*rf_dia];
        % remove NaN and normalized
        Y = RFsimilarity;
        remove_thid = any(isnan([X Y]),2);
        Y(remove_thid,:) = [];
        X(remove_thid,:) = [];
        X = zscore(X);
        X(:,1) = ones(size(X,1),1);
        [nnn,kkk] = size(X);
        kkk = kkk - 1; % 不包括常数项
        
        [B,BINT,R,RINT,STATS] = regress(Y,X);
        
        tbl = table(Y, X(:,2), X(:,3), X(:,4), 'Var', {'Similarity','RF_angle','RF_eccentricity','RF_diameter'});
        mdl = fitlm(tbl,'Similarity ~ RF_angle + RF_eccentricity + RF_diameter + RF_angle*RF_eccentricity + RF_angle*RF_diameter + RF_eccentricity*RF_diameter + RF_angle*RF_eccentricity*RF_diameter')
        [p,F,d] = coefTest(mdl)
        
        % 后3项目p>0.05，去掉
        mdl2 = fitlm(tbl,'Similarity ~ RF_angle + RF_eccentricity + RF_diameter + RF_angle*RF_eccentricity');
        % P>0.05: do not have a significant impact on similarity
        % P<0.05: the term is statistically significant.
        % 结果：都贡献
        
        % coefTest: F-test for the hypothesis that all regression coefficients (except for the intercept) are zero versus at least one differs from zero,
        % 针对整个model
        [p,F,d] = coefTest(mdl2);
        
        
        
    end


    function f2p9(debug) % 'Clustering'
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        
        % 参考 Ringbell_MSTd_cluster.m
        % 深度相隔100，200，300，400um， 允许误差20um
        % group_result(n).cluster_FILE
        
        % 1. cluster correlation
        % prepare array for data
        for d = 1:length(depth100) % 相隔100，200，300，400um
            for ind = 1:3
                cluster_r{ind,d} = nan(length(loadM_data),5); % large enough for most 4 time comparision in one depth, 同一个cite分配5个位置预存，例如dep=100，可以存[1000 1100],[1100 1200],[1200 1300],[1300 1400],[1400 1500] (实际上大部分情况只会用到前两个位置。)
                cluster_p{ind,d} = nan(length(loadM_data),5);
            end
            diff_spiral_index{d} = nan(length(loadM_data),5);
            diff_spiral_index_ref{d} = nan(length(loadM_data),5);
            diff_spiral_index_ref_plus{d} = nan(length(loadM_data),5);
        end
        
        for n = 1:length(loadM_data)
            if methods_of_select{1,1}(n)
                for d = 1:length(depth100) % 相隔100，200，300，400um
                    if ~isnan(group_result(n).depth_diff_plot{d}(1))
                        depth_diff_plot = group_result(n).depth_diff_plot{d};
                        depth_all = [group_result(n).cluster.depth];
                        for j = 1:size(depth_diff_plot,1) % j:多个相邻相同间隔(例如都是100)的site，例如[1000 1100],[1100 1200],[1200 1300],...
                            for ind = 1:3
                                if ind ~= 3
                                    resp1 = group_result(n).cluster(depth_all == depth_diff_plot(j,1)).resp(:,ind); % 相邻的第一个深度
                                    resp2 = group_result(n).cluster(depth_all == depth_diff_plot(j,2)).resp(:,ind); % 相邻的第二个深度
                                else % ******** how to compare 【Motion tpye cluster】 ?????????
                                    resp1 = group_result(n).cluster(depth_all == depth_diff_plot(j,1)).resp(:);
                                    resp2 = group_result(n).cluster(depth_all == depth_diff_plot(j,2)).resp(:);
                                end
                                
                                % 计算corr
                                [cluster_r{ind,d}(n,j), cluster_p{ind,d}(n,j)] = corr(resp1,resp2,'type','Pearson');
                                
                                % spiral index difference
                                spiral_index1 = group_result(n).cluster(depth_all == depth_diff_plot(j,1)).spiral_index; % ref
                                spiral_index2 = group_result(n).cluster(depth_all == depth_diff_plot(j,2)).spiral_index; % ref+100
                                diff_spiral_index{d}(n,j) = abs(spiral_index2 - spiral_index1);
                                
                                % for scatter plot
                                diff_spiral_index_ref{d}(n,j) = spiral_index1;
                                diff_spiral_index_ref_plus{d}(n,j) = spiral_index2;
                            end
                        end
                    end
                end
            else
                cluster_r{ind,d}(n,:) = nan(1,5);
                cluster_p{ind,d}(n,:) = nan(1,5);
                diff_spiral_index{ind,d}(n,:) = nan(1,5);
                diff_spiral_index_ref{d}(n,:) = nan(1,5);
                diff_spiral_index_ref_plus{d}(n,:) = nan(1,5);
            end
        end
        
        % 2. 计算上下100um深度下的prefered direction的stardand deviation
        pd_std = cell2mat({group_result.pd_std}'); % circular std
        pd_std2 = cell2mat({group_result.pd_std2}'); % normal std
        
        
        % reorganization and select
        for d = 1:4
            for ind = 1:3
                corr_coe_temp = cluster_r{ind,d}(:); % 相隔相同距离的合并
                corr_coe{ind,d} = corr_coe_temp(~isnan(corr_coe_temp)); % 去掉nan
                corr_coe_p_temp = cluster_p{ind,d}(:); 
                corr_coe_p{ind,d} = corr_coe_p_temp(~isnan(corr_coe_p_temp));
            end
            diff_spiral_index_temp = diff_spiral_index{d}(:);
            diff_si{d} = diff_spiral_index_temp(~isnan(diff_spiral_index_temp));
            
            % for scatter plot
            si_ref_temp = diff_spiral_index_ref{d}(:);
            si_ref_plus_temp = diff_spiral_index_ref_plus{d}(:);
            si_ref{d} = si_ref_temp(~isnan(si_ref_temp)); % si_ref 和 si_ref_plus 一一对应的
            si_ref_plus{d} = si_ref_plus_temp(~isnan(si_ref_plus_temp));
        end
        
        
        % ploting
        set(figure(figN),'Position',[350,100, 1200,800], 'Name', 'Cluster (Correlation)');
        
        edges = [-1:0.2:1];
        %         plot_depth_range = [1:3];
        plot_depth_range = [1:length(depth100)];
        ha = tight_subplot(max(plot_depth_range),4,[.1 .08],[.1 .07],[.07 .03],1);
        
        for ind = 1:3
            for d = plot_depth_range
                axes(ha(ind+(d-1)*4))
                nbins = hist(corr_coe{ind,d},edges);
                nbins_sig = hist(corr_coe{ind,d}(corr_coe_p{ind,d}<0.05),edges);
                hold on
                bar(edges,nbins,1,'facecolor','w','edgecolor','k');
                bar(edges,nbins_sig,1,'facecolor',c{ind},'edgecolor','k');
                xlim([-1.1 1.1]);
                ylimt = max(nbins);
                ylimt1 = ylimt + ylimt/3;
                ylim([0 ylimt1]);
                
                %  median correlation values and sign test，是否与0有显著区别
                median_cor{ind,d} = median(corr_coe{ind,d});
                p_sign_test{ind,d} = signtest(corr_coe{ind,d});
                mean_cor{ind,d} = mean(corr_coe{ind,d});
                sem_cor{ind,d} = std(corr_coe{ind,d})/sqrt(length(corr_coe{ind,d}));
                [~,p_ttest{ind,d}] = ttest(corr_coe{ind,d});
                
                % 箭头标记median值
                text(median_cor{ind,d},ylimt+ylimt/4,'\downarrow\bf','fontsize',15,'horizontalalignment','center');
                if p_sign_test{ind,d}<0.001
                    text(median_cor{ind,d},ylimt1,'***\bf','fontsize',20,'horizontalalignment','center');
                elseif (p_sign_test{ind,d}>0.001 && p_sign_test{ind,d}<0.01)
                    text(median_cor{ind,d},ylimt1,'**\bf','fontsize',20,'horizontalalignment','center');
                elseif (p_sign_test{ind,d}>0.01 && p_sign_test{ind,d}<0.05)
                    text(median_cor{ind,d},ylimt1,'*\bf','fontsize',20,'horizontalalignment','center'); 
                end
                text(-1, ylimt+ylimt/4, sprintf('median = %.2g',median_cor{ind,d}));
                text(-1, ylimt, sprintf('p = %.2g (sign test)',p_sign_test{ind,d}))
                text(-1, ylimt-ylimt/4, sprintf('%.2g±%.2g',mean_cor{ind,d},sem_cor{ind,d}));
                text(-1, ylimt-ylimt/2, sprintf('p = %.2g (t test)',p_ttest{ind,d}));

                ylabel('Cases');
                if ind == 2
                    delta = char(hex2dec('0394'));
                    titl = sprintf('depth = %.0f um',d*100);
                    s = strcat(delta, titl);
                    title(s);
                end
            end
        end
        
        % spiral index difference
        for d = plot_depth_range
            axes(ha(4+(d-1)*4))
            edges = [0:0.2:2];
            nbins = hist(diff_si{d},edges);
            bar(edges+0.1,nbins,0.8,'facecolor',c{ind},'edgecolor','k'); % 让每个bar中心对着原来的中心
            xlim([-0.1 1.6]);
            
            %  median
            median_si{d} = median(diff_si{d});
            p_si_sign_test{d} = signtest(diff_si{d},0.5,'tail','left'); % p<0.05 显著小于0.5
            
            % 箭头标记median值，是否显著小于0.5
            ylimt = ylim;
            ylimt = ylimt(2);
            hold on
%             plot([0.5 0.5],[0 ylimt],'k:');
            text(median_si{d},ylimt+ylimt/5,'\downarrow\bf','fontsize',15);
            if p_si_sign_test{d}<0.001
                text(median_si{d}-0.05,ylimt,'***\bf','fontsize',20);
            elseif (p_si_sign_test{d}>0.001 && p_sign_test{ind,d}<0.01)
                text(median_si{d}-0.05,ylimt,'**\bf','fontsize',20);
            elseif (p_si_sign_test{d}>0.01 && p_sign_test{ind,d}<0.05)
                text(median_si{d}-0.05,ylimt,'*\bf','fontsize',20);   
            end
            
            if d == max(plot_depth_range)
                xlabel('|Spiral Index|');
            end
        end
        SetFigure(13);figN = figN+1;
        
        set(figure(figN),'Position',[500,100, 800,800], 'Name', 'Spiral Index scatter plot');clf
        for d = plot_depth_range
            subplot(ceil(max(plot_depth_range)/2),2,d)
            
            scatter(si_ref{d},si_ref_plus{d},'MarkerEdgeColor','k','markerfacecolor',c{99});
            
            hold on
            plot([0 0],[-1 1],'k:');
            plot([-1 1],[0 0],'k:');
            SameScale(1)
            xlim([-1 1]);ylim([-1 1]);
            xlabel('Spiral Index (ref)');
            ylabel_temp = strcat('Spiral Index (ref+',num2str(plot_depth_range(d)*100),')');
            ylabel(ylabel_temp);
            
            % show the similar dot (same sign)
            patch([-1 0 0 -1],[-1 -1 0 0],[1,2,3,4],'edgecolor','none','facecolor',c{99},'facealpha',.2);
            patch([0 1 1 0],[0 0 1 1],[1,2,3,4],'edgecolor','none','facecolor',c{99},'facealpha',.2);
            
            % correlation
            [r,p] = corr(si_ref{d},si_ref_plus{d},'type','pearson');
            text(-0.9,0.9,sprintf('Pearson"s r = %1.3g',r));
            text(-0.9,0.75,sprintf('p = %.2g',p));
            
            delta = char(hex2dec('0394'));
            titl = sprintf('depth = %.0f um',d*100);
            s = strcat(delta, titl);
            title(s);
        end
        SetFigure(13);figN = figN+1;
        
%         % preferred direction SD
%         set(figure(figN),'Position',[500,100, 1300,550], 'Name', 'PD SD');
%         h = LinearCorrelation({
%             pd_std(methods_of_select{1,1},1);
%             },...
%             {
%             pd_std(methods_of_select{1,1},2);
%             },...
%             'CombinedIndex',[],'PlotCombinedOnly',0,...
%             'Xlabel','Preferred Translation SD (deg)','Ylabel','Prefferred Rotation SD(deg)',...
%             'FaceColors',{c{99}},'Markers',{'o'},...
%             'LineStyles',{'k-'},'MarkerSize',10,...
%             'figN',figN,'XHist',8,'YHist',8,...
%             'XHistStyle','grouped','YHistStyle','grouped','SameScale',1,...
%             'Method','Pearson','FittingMethod',2,...
%             'dot_leg',{'preferred direction SD'},...
%             'logx',1,'logy',1);
%         SetFigure(13);figN = figN + 1;
    end

    function f2p24(debug) % 'Clustering for different similarity sites' % cluster
        if debug  ; dbstack;   keyboard;      end
        methods_of_select = {select_typical, 'Typical cells'};
        
        % 1. cluster correlation average ±100um
        cluster_index = cell2mat({group_result.cluster_index}');
        cluster_index = cluster_index(methods_of_select{1,1},:);
        
        % 2. similarity
        try
            RFsimilarity = [group_result.RFsimilarity]';
        catch
            run_this = strcmpi({handles.function_handles_spi{2,2}{:,1}}','Save Prefer OpticFlow similarity (local planar motion)');
            feval(handles.function_handles_spi{2,2}{run_this,2},false); % first run "Spiral Tuning" - "Spiral Tuning" - "Save Prefer OpticFlow (local planar motion)"
            RFsimilarity = [group_result.RFsimilarity]';
        end
        RFsimilarity = RFsimilarity(methods_of_select{1,1});
        
        % group 3 similarity groups
        group_range = linspace(0,100,4);
        for g = 1:3
            group_this{g} = logical(RFsimilarity > group_range(g) & RFsimilarity <= group_range(g+1));
        end
        
        set(figure(figN),'Position',[350,100, 1300,550], 'Name', 'Clustering for different similarity sites');clf
        for g = 1:3
            subplot(1,3,g)
            cluster_this = cluster_index(group_this{g},:);
            select = sum(~isnan(cluster_this),2)==2; % select both not nan
            cluster_this = cluster_this(select,:);
            N_this = size(cluster_this,1);
            
%             % paired plot
%             plot([1 2],[cluster_this(:,2),cluster_this(:,1)],'ko-');
%             set(gca,'xtick',[1 2],'xticklabel',{'R cluster','T cluster'});
%             xlim([0.5 2.5]);
%             ylim([-1 1]);
%             if g == 1
%                 xlabel('RFsimilarity = 0~33 (different)');
%             elseif g == 2
%                 xlabel('RFsimilarity = 33~66 (medium)');
%             else
%                 xlabel('RFsimilarity = 66~100 (similar)');
%             end
%             ylabel('Cluster index (±100μm)');
%             
%             % median
%             medianPSE = median([cluster_this(:,2),cluster_this(:,1)]);
%             hold on
%             plot([1 2],medianPSE,'r.-','markersize',10,'linewidth',2);
%             
%             
%             % Wilcoxon Rank-Sum Test
%             p = ranksum(cluster_this(:,2),cluster_this(:,1));
%             text(1.5,1.1,sprintf('p=%.3g',p),'horizontalAlignment','center');
            
            
            % scatter plot
            plot(cluster_this(:,1),cluster_this(:,2),'k.','markersize',10);
            ylim([-1 1]); xlim([-1 1]);
            SameScale(1);
            [r,pp] = plot_corr_line(cluster_this(:,1),cluster_this(:,2),'MethodOfCorr','pearson','FittingMethod',2,'LineWidth',2);
            str = ['Pearson cor.' newline sprintf('r=%.3g',r) newline sprintf('p=%.3g',pp) newline sprintf('N = %g',N_this)];
            text(-0.9,-0.5,str);

            if g == 1
                title('RFsimilarity = 0~33 (different)');
            elseif g == 2
                title('RFsimilarity = 33~66 (medium)');
            else
                title('RFsimilarity = 66~100 (similar)');
            end
            xlabel('Translation cluster index (±100μm)');
            ylabel('Rotation cluster index (±100μm)');
        end
        SetFigure(15)

         [~,~,~,monkey] = SetTitle([],monkey_included_for_analysis);
        str =['Clustering for different similarity sites',monkey];
        suptitle(str); 
        figN = figN+1;
    end

    function f2p10(debug) % Preferred Direction around 0 degree (compare to CueS task)
        if debug  ; dbstack;   keyboard;      end
        
        % prefer axis distribution
        p_aroundF = cell2mat({group_result.p_aroundF}');
        pref_axis_num_temp = cell2mat({group_result.pref_axis_num}');
        pref_axis_num_temp(isnan(p_aroundF)) = nan;
        pref_axis_num_temp(sum(isnan(pref_axis_num_temp),2)==1,:) = nan;
        
        % combine with anova test
        for d = 1:2 % fine and coarse
            pref_axis_num_temp(p_aroundF(:,d)>0.05,d) = 0; % anova不显著的直接认为pref axis num=0
        end
        
        set(figure(figN),'Position',[350,100, 1300,550], 'Name', 'Prefer axis num around Forward motion');clf
        
        for d = 1:2 % fine and coarse
            for j = 1:4 % pref axis = 0 1 2 3
                if j == 1 % pref axis = 0
                    pref_axis_num(d,j) = sum(p_aroundF(:,d) > 0.05 | pref_axis_num_temp(:,d)==0);
                    pref_axis_select{d,j} = logical(p_aroundF(:,d) > 0.05 | pref_axis_num_temp(:,d)==0);
                else % % pref axis = 1,2,3
                    pref_axis_num(d,j) = sum(p_aroundF(:,d) < 0.05 & pref_axis_num_temp(:,d)==(j-1));
                    pref_axis_select{d,j} = logical(p_aroundF(:,d) < 0.05 & pref_axis_num_temp(:,d)==(j-1));
                end
            end
            subplot(1,2,d)
            %             bar(pref_axis_num(i,:));
            %             set(gca,'xticklabel',[0 1 2 3]);
            pie(pref_axis_num(d,:),[1 0 0 0]);
            
            if d == 1
                title('45 degree around Exp.');
            else
                title('90 degree around Exp.');
                legend('No pref axis','1 pref axis','2 pref axes','3 pref axes');
            end
        end
        suptitle('Prefer axis num around Forward motion');
        SetFigure(15);
        figN = figN+1;
        
        % 对比45°与90°的axes num区别，点对点,重叠的加粗线
        % 两边时0，1，2，3，全连接就是16根线
        set(figure(figN),'Position',[350,100, 550,550], 'Name', 'Prefer axis num between 45° and 90°');
        for d = 1:4
            for j = 1:4
                temp_num(d,j) = sum(pref_axis_num_temp(:,1) == (d-1) & pref_axis_num_temp(:,2) == (j-1));
                hold on
                try
                    if d==j
                        plot([1,2],[d-1,j-1],'-r','linewidth',temp_num(d,j)); % 45 and 90 is same num
                    else
                        plot([1,2],[d-1,j-1],'-k','linewidth',temp_num(d,j)); % difference num
                    end
                end
            end
        end
        set(gca,'ytick',[0 1 2 3]);
        yyaxis right;
        ylim([0 3]);
        set(gca,'ytick',[0 1 2 3]);
        
        set(gca,'xtick',[1 2],'xticklabel',{'45 degree','90 degree'});
        SetFigure(15);
        figN = figN+1;
        
        % rearragement,seperate fine and coarse
        for n = 1:length(loadM_data)
            pref_axis_fine{n,1} = group_result(n).pref_axis{1};
            pref_axis_coarse{n,1} = group_result(n).pref_axis{2};
        end
        
        % fine
        select_pref_axis{1,1} =  cell2mat({pref_axis_fine{pref_axis_select{1,2}}}'); % num = 1
        select_pref_axis{1,2} =  cell2mat({pref_axis_fine{pref_axis_select{1,3}}}'); % num = 2
        % coarse
        select_pref_axis{2,1} =  cell2mat({pref_axis_coarse{pref_axis_select{2,2},1}}'); % num = 1
        select_pref_axis{2,2} =  cell2mat({pref_axis_coarse{pref_axis_select{2,3},1}}'); % num = 2
        
        %         % check, 找到只fine，pref one axis 且为R（1）/ L（3）的unit    确实pref R的更多！
        %         select_pref_one = pref_axis_select{1,2};
        %          for n = 1:length(loadM_data)
        %              if pref_axis_fine{n,1} == 1
        %                  select_pref_R(n,1) = true;
        %              else
        %                  select_pref_R(n,1) = false;
        %              end
        %
        %              if pref_axis_fine{n,1} == 3
        %                  select_pref_L(n,1) = true;
        %              else
        %                  select_pref_L(n,1) = false;
        %              end
        %          end
        %
        %         sum(select_pref_R  & select_pref_one)
        %         sum(select_pref_L  & select_pref_one)
        
        % prefer axis的具体分布 only select 1 and 2 prefer axes
        for d = 1:2 % fine and coarse
            % prefer one axis distribution
            pref1_type{d} = hist(select_pref_axis{d,1},4); % R CW L CCW
            
            % prefer two axis distribution
            % R/CW CW/L L/CCW CCW/R L/R CW/CCW
            %  12  23    34    41    31   24
            pref2_type_temp = [1 2; 2 3; 3 4; 4 1; 3 1; 2 4]; % 6 types
            for j = 1:length(pref2_type_temp)
                pref2_type{d}(j) = sum( (select_pref_axis{d,2}(:,1) == pref2_type_temp(j,1) & select_pref_axis{d,2}(:,2) == pref2_type_temp(j,2)) | ...
                    (select_pref_axis{d,2}(:,2) == pref2_type_temp(j,1) & select_pref_axis{d,2}(:,1) == pref2_type_temp(j,2)));
            end
        end
        
        % bar plot:prefer one, prefer two axis, and pref opposite axis(L/R, CW/CCW)
        set(figure(figN),'Position',[350,100, 1000,500], 'Name', 'Prefer axis(seperate One and Two)');
        for d = 1:2 % fine and coarse
            subplot(1,2,d)
            bar([pref1_type{d} pref2_type{d}]);
            set(gca,'xticklabel',{'R','CW','L','CCW','R/CW','CW/L','L/CCW','CCW/R','L/R','CW/CCW'},'xTickLabelRotation',45);
            if d==1
                title('45 degree');
            else
                title('90 degree');
            end
        end
        suptitle('One or Two prefer axis');
        SetFigure(13);
        figN = figN+1;
        
        % compass plot
        % plot prefer one and two axes together
        az_axis = deg2rad([0:45:315]);
        % pool pref1_type, pref2_type together
        for d = 1:2 % fine and coarse
            pref_type{d}(1) = pref1_type{d}(1); % R
            pref_type{d}(2) = pref2_type{d}(1); % R/CW
            pref_type{d}(3) = pref1_type{d}(2); % CW
            pref_type{d}(4) = pref1_type{d}(2); % CW/L
            pref_type{d}(5) = pref1_type{d}(3); % L
            pref_type{d}(6) = pref1_type{d}(3); % L/CCW
            pref_type{d}(7) = pref1_type{d}(4); % CCW
            pref_type{d}(8) = pref1_type{d}(4); % CCW/R
        end
        
        set(figure(figN),'Position',[350,100, 1000,500], 'Name', ' Prefer axis');
        % plot
        for d = 1:2 % fine and coarse
            subplot(1,2,d)
            [x,y] = pol2cart(az_axis,pref_type{d});
            compass(x,y,'r');
            
            % add L/R and CW/CCW
            hold on
            [x,y] = pol2cart(deg2rad([0 180]),[pref2_type{d}(5) pref2_type{d}(5)]);
            compass(x,y,'k');
            
            [x,y] = pol2cart(deg2rad([90 270]),[pref2_type{d}(6) pref2_type{d}(6)]);
            compass(x,y,'k');
            
            if d==1
                title('45 degree');
            else
                title('90 degree');
            end
        end
        suptitle('Prefer axes (One or Two axes)');
        SetFigure(13);
        figN = figN+1;
    end

    function f2p11(debug) % 'Tuning Index'
        if debug  ; dbstack;   keyboard;      end
        
        % 1. d prime
        d_prime_fine = cell2mat({group_result.d_prime_fine}');
        d_prime_coarse = cell2mat({group_result.d_prime_coarse}');
        d_prime_maxDiff = cell2mat({group_result.d_prime_maxDiff}');
        d_prime_maxAxis = cell2mat({group_result.d_prime_maxAxis}');
        d_prime_maxFR = cell2mat({group_result.d_prime_maxFR}');
        
        % 2. DDI
        DDI_includeExpCon = cell2mat({group_result.DDI_includeExpCon}');
        DDI_excludeExpCon = cell2mat({group_result.DDI_excludeExpCon}');
        
        % 3.HTI
        HTI = cell2mat({group_result.HTI}');
        
        % 4. Neural sensitivity
        neu_sen_fine = cell2mat({group_result.neu_sen_fine}');
        neu_sen_coarse = cell2mat({group_result.neu_sen_coarse}');
        
        % slope at 0 degree
        ahead_slope = cell2mat({group_result.ahead_slope}');
        
        % % modulation index, DeAnglis and UKA, 2002, add by Lwh202101
        % direction/disparity tuning index, DeAngelis1, Newsome, 2004. tuning strength
        MDI = cell2mat({group_result.MDI}');
        
        % pair discrimination index
        pair_DI = cell2mat({group_result.pair_DI}');
        
        % non-pair discrimination index
        nonpair_DI = cell2mat({group_result.nonpair_DI}');
        nonpair_DI2 = cell2mat({group_result.nonpair_DI2}');
        
        % separable predictions
        SVM_correct_rate = cell2mat({group_result.SVM_correct_rate}');
        AuROC = cell2mat({group_result.AuROC}');

        % modulation index, DeAngelis and UKA, 2002
        modulation_index{1} = d_prime_fine; % d' = (Rright45 - Rleft45) / sqrt(1/2*(STDright^2+STDleft^2)) 
        modulation_index{2} = d_prime_coarse; % d' = (Rright90 - Rleft90) / sqrt(1/2*(STDright^2+STDleft^2)) 
        modulation_index{3} = d_prime_maxDiff; % d' = (Rmax - Rmin) / sqrt(1/2*(STDmax^2+STDmin^2)) 
        modulation_index{4} = d_prime_maxAxis; % d' = (Rmax180 - Rmin180) / sqrt(1/2*(STDmax180^2+STDmin180^2)) 
        modulation_index{5} = d_prime_maxFR; % d' = (Rmax - Rc-lat) / sqrt(1/2*(STDmax^2+STDc-lat^2)) 
        modulation_index{6} = DDI_includeExpCon; % DDI = delta / (delta + 2 * sqrt(SSE/(N-M)))
        modulation_index{7} = DDI_excludeExpCon; % DDI = delta / (delta + 2 * sqrt(SSE/(N-M)))
        modulation_index{8} = HTI; % HTI = |ΣRi| / Σ|Ri|   AMPvectorsum / sum(|mean-spon|)  (Ri = mean-spon)
        modulation_index{9} = neu_sen_fine; %  lg(|d'45|)
        modulation_index{10} = neu_sen_coarse; %  lg(|d'90|)
        modulation_index{11} = ahead_slope;
        modulation_index{12} = MDI; % modulation index = (Rmax-Rmin)/(Rmax-Spon)
        modulation_index{13} = pair_DI; % 1/n * Σ( (Rfar(i)-Rnear(i)) / |Rfar(i)-Rnear(i)| + std)
        modulation_index{14} = nonpair_DI; % |Rleft-Rright| / (|Rleft-Rright| + 2*sqrt(SSE/(N-2)))
        modulation_index{15} = nonpair_DI2; % |Rleft-Rright| / (|Rleft-Rright| + STDavg)
        modulation_index{16} = SVM_correct_rate;
        modulation_index{17} = AuROC;
        
        tl{1} = strcat('d''','(fine)');
        tl{2} = strcat('d''','(coarse)');
        tl{3} = strcat('d''','(maxDiff)');
        tl{4} = strcat('d''','(maxAxis)');
        tl{5} = strcat('d''','(maxFR)');
        tl{6} = 'DDI(inc)';
        tl{7} = 'DDI(exc)';
        tl{8} = 'HTI';
        tl{9} = 'Neural Sen.(fine)';
        tl{10} = 'Neural Sen.(coarse)';
        tl{11} = 'Zero slope';
        tl{12} = 'Modulation index'; 
        tl{13} = 'Pair DI';
        tl{14} = 'Non-pair DI';
        tl{15} = 'Non-pair DI2';
        tl{16} = 'SVM correct rate';
        tl{17} = 'AuROC';
        
        methods_of_select = {select_typical, 'Typical cells'};
        select = methods_of_select{1,1};
        
        set(figure(figN),'Position', [100,200 1700,600], 'Name', 'Tuning Index');clf;
        ha = tight_subplot(3,length(modulation_index),[.1 .02],[.1 .2],[.03 .01]);
        for m = 1:length(modulation_index)
            for d = 1:2 % h or r
                axes(ha(m+(d-1)*length(modulation_index)))
                if length(loop_num) == 1 && loop_num == 0 % only select bottom line
                    histogram(modulation_index{m}(select,d),'facecolor',c{99});
                elseif length(loop_num) == 1 % select specific cell type
                    histogram(modulation_index{m}(select,d),'facecolor',c{loop_num});
                else % select T+R+S cell
                    histogram(modulation_index{m}(select,d),'facecolor',c{99});
                end
                
                if d == 1
                    title(tl{m});
                end
                
                if m == 1
                    if d == 1
                        ylabel('Trans. tuning');
                    elseif d == 2
                        ylabel('Rot. tuning');
                    end
                end
                
                axis tight
                
                if m == 6 || m == 7 || m == 8 ||m == 13 || m == 14 || m == 15 || m == 16
                    xlim([0 1]);
                elseif m==17
                    xlim([0.5 1]);
                end
            end
            
            axes(ha(m+(3-1)*length(modulation_index)))
            TT = modulation_index{m}(select,1);
            RR = modulation_index{m}(select,2);
            h_r_index(:,m) = (abs(TT)-abs(RR))./(abs(TT)+abs(RR));
            if length(loop_num) == 1 % only select one cell type
                histogram(h_r_index(:,m),'facecolor',c{loop_num});
            else
                histogram(h_r_index(:,m),'facecolor',c{99});
            end
            
            
            if m == 1
                ylabel('(T-R)/(T+R)');
            end
            axis tight
        end
        figN = figN+1;
        
        %         figure
        %         for i = 1:sum(select)
        %             h = plot([1:length(modulation_index)],h_r_index(i,:),'.-');
        %             ylim([-1 1]);
        %
        %             if i == 1
        %                 set(gca, 'xtick', [1:length(modulation_index)]); %设置x坐标轴刻度数据点位置
        %                 set(gca,'xTickLabel',tl);%不显示y坐标轴的值
        %                 RotateTickLabel(gca,90);
        %                 hold on
        %                 plot([1 17],[0 0],'k:');
        %             end
        %         end
    end

    function f2p16(debug) % 'Gaussian Fit Result'
        if debug  ; dbstack;   keyboard;      end
        
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'}; % 默认T cell+R cell+ Spi cell
        
        % ---------------------- get data ----------------------------
        fitOK = cell2mat({group_result.fitOK}'); %  if output.exitflag > 0 && cf_{s}.sigma >= 0.3 && cf_{s}.sigma < 500; fitOK = 1;
        fitOK = fitOK(methods_of_select{1,1},:);
        
        GaussfitRsquare = cell2mat({group_result.GaussfitRsquare}');
        GaussfitRsquare = GaussfitRsquare(methods_of_select{1,1},:);
        
        % only select 同时fitOK=1
        select = sum(fitOK,2)==2;
        
        set(figure(figN),'Position',[610,70,500,300], 'Name', ' Gaussian Fit R square');
        h = LinearCorrelation({
            GaussfitRsquare(select,1);
            },...
            {
            GaussfitRsquare(select,2);
            },...
            'CombinedIndex',[],'PlotCombinedOnly',0,...
            'Xlabel','Translation Fitting','Ylabel','Rotation Fitting',...
            'FaceColors',{c{99}},'Markers',{'o'},...
            'LineStyles',{'k-'},'MarkerSize',10,...
            'figN',figN,'XHist',6,'YHist',6,...
            'XHistStyle','grouped','YHistStyle','grouped','SameScale',1,...
            'Method','Pearson','FittingMethod',2,...
            'dot_leg',{'R^2 (both fit OK)'}); figN = figN + 1;
        
        % add 不同时fitOK的unit
        axes(h.ax_raw);
        xlims = xlim; ylims = ylim;
        selectOnlyT = logical(fitOK(:,1)==1 & fitOK(:,2)==0);
        selectOnlyR = logical(fitOK(:,2)==1 & fitOK(:,1)==0);
        h1 = plot(GaussfitRsquare(selectOnlyT,1),ylims(1)*ones(size(GaussfitRsquare(selectOnlyT,1))),'ko','markerfacecolor',c{1},'markersize',10);
        h2 = plot(xlims(1)*ones(size(GaussfitRsquare(selectOnlyR,2))),GaussfitRsquare(selectOnlyR,2),'ko','markerfacecolor',c{2},'markersize',10);
        legend([h.group.dots h1 h2],'R^2 (both fit OK)','R^2 (only Trans. fit OK)','R^2 (only Rot. fit OK)');
        
        % compare GaussfitRsquare between T and R (all cell)
        % 大样本，非正态分布，还是用ttest吧
        [~,p_ttest] = ttest2(GaussfitRsquare(:,1),GaussfitRsquare(:,2));
        text(-0.1,0,sprintf('Non-paired ttest, p = %0.3g',p_ttest));
        
        % FitOK的比例
        fitOK_num = sum(fitOK); % T and R
        both_fitOK_num = sum(select);
        all_num = size(fitOK,1);
        
        axes('position',[0.72 0.1 0.2 0.4]);
        xlabe = categorical({'Translation Fitting','Rotation Fitting'});
        
        h1 = bar(xlabe,[both_fitOK_num fitOK_num(1)-both_fitOK_num,all_num-fitOK_num(1); both_fitOK_num fitOK_num(2)-both_fitOK_num,all_num-fitOK_num(2)],'stacked');
        set(h1(1),'FaceColor',c{99});
        set(h1(2),'FaceColor',c{1});
        set(h1(3),'FaceColor','w');
        
        hold on
        h2 = bar(xlabe,[0 0 0; both_fitOK_num fitOK_num(2)-both_fitOK_num,all_num-fitOK_num(2)],'stacked');
        set(h2(1),'FaceColor',c{99});
        set(h2(2),'FaceColor',c{2});
        set(h2(3),'FaceColor','w');
        legend('Both fit OK','One fit OK','Fit not good','location','bestoutside')
        title('Gaussian Fitting ');
        
        SetFigure(13);
        
        
        %% 比较vector preferred direction 与 Gaussian fit preferred direction
        prefGauss = cell2mat({group_result.prefGauss}');
        pref_direc = cell2mat({group_result.pref_direc}');
        
        set(figure(figN),'Position',[610,70,700,300], 'Name', ' GaussianPref and VectorSumPref');
        for s = 1:2 % t or r
            subplot(1,2,s)
            
            aa = pref_direc(methods_of_select{1,1},s);
            bb = prefGauss(methods_of_select{1,1},s);
            nan_this = any(isnan([aa bb]),2);
            aa = aa(~nan_this);
            bb = bb(~nan_this);

            plot(aa,bb,'k.');
            SameScale(1);
            xticks([-180:45:180]);
            yticks([-180:45:180]);
            axis equal
            xlabel('Vector-sum');
            ylabel('Gaussian fitting');
            if s == 1
                title('Translation preferred direction');
            else
                title('Rotation preferred direction');
            end
            
            [r,p] = corr(aa,bb);
            text(-180,180,sprintf('r = %.3g, p = %.3g',r,p));
        end
        SetFigure(13);
        % 基本在一条对角线上
        figN = figN+1;
    end


    function  f2p17(debug) % 'Gaussian Bandwidth in T and R'
        if debug  ; dbstack;   keyboard;      end
        
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'}; % 默认T cell+R cell+ Spi cell
        
        % ---------------------- get data ----------------------------
        % 只对fitOK的细胞
        fitOK = cell2mat({group_result.fitOK}'); %  if output.exitflag > 0 && cf_{s}.sigma >= 0.3 && cf_{s}.sigma < 500; fitOK = 1;
        fitOK = fitOK(methods_of_select{1,1},:);
        widGauss = cell2mat({group_result.widGauss}');
        widGauss = widGauss(methods_of_select{1,1},:);
        
        % only select 同时fitOK=1
        select = sum(fitOK,2)==2;
        
        set(figure(figN),'Position',[610,70,500,300], 'Name', ' Gaussian Bandwidth');
        h = LinearCorrelation({
            widGauss(select,1);
            },...
            {
            widGauss(select,2);
            },...
            'CombinedIndex',[],'PlotCombinedOnly',0,...
            'Xlabel','Translation Bandwidth','Ylabel','Rotation Bandwidth',...
            'FaceColors',{c{99}},'Markers',{'o'},...
            'LineStyles',{'k-'},'MarkerSize',10,...
            'figN',figN,'XHist',7,'YHist',7,...
            'XHistStyle','grouped','YHistStyle','grouped','SameScale',1,...
            'Method','Pearson','FittingMethod',2,...
            'dot_leg',{'R^2'});
        
        % compare Bandwidth between T and R
        % 大样本，非正态分布，还是用ttest吧
        axes(h.ax_raw);
        [~,p_ttest] = ttest2(widGauss(:,1),widGauss(:,2));
        text(200,350,sprintf('Non-paired ttest, p = %0.3g',p_ttest));
        
        SetFigure(13);
        figN = figN + 1;
    end

    function f2p18(debug) %'Bandwidth and DDI'
        if debug  ; dbstack;   keyboard;      end
        
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'}; % 默认T cell+R cell+ Spi cell
        
        % ---------------------- get data ----------------------------
        % 只对fitOK的细胞
        fitOK = cell2mat({group_result.fitOK}'); %  if output.exitflag > 0 && cf_{s}.sigma >= 0.3 && cf_{s}.sigma < 500; fitOK = 1;
        fitOK = fitOK(methods_of_select{1,1},:);
        
        widGauss = cell2mat({group_result.widGauss}');
        widGauss = widGauss(methods_of_select{1,1},:);
        
        DDI{1} = cell2mat({group_result.DDI_includeExpCon}');
        DDI{2} = cell2mat({group_result.DDI_excludeExpCon}');
        DDI{1} = DDI{1}(methods_of_select{1},:);
        DDI{2} = DDI{2}(methods_of_select{1},:);
        
        labe{1,1} = {'Translation Bandwidth','Translation DDI'};
        labe{1,2} = {'Rotation Bandwidth','Rotation DDI'};
        labe{2,1} = {'Translation Bandwidth','Rotation DDI'};
        labe{2,2} = {'Rotation Bandwidth','Translation DDI'};
        
        set(figure(figN),'Position',[610,70,800,800], 'Name', 'Gaussian Bandwidth and DDI');
        for hr = 1:2
            subplot(2,2,hr)
            plot(widGauss(fitOK(:,hr)==1,hr),DDI{1}(fitOK(:,hr)==1,hr),'ko'); % Translation bandwidth 对应 Translation DDI
            [r,pp] = plot_corr_line(widGauss(fitOK(:,hr)==1,hr),DDI{1}(fitOK(:,hr)==1,hr),'MethodOfCorr','pearson','FittingMethod',2,'LineWidth',2);
            xlimit = xlim;
            ylimit = ylim;
            text((xlimit(2)-xlimit(1))/2+xlimit(1),ylimit(2)-(ylimit(2)-ylimit(1))/20,sprintf('r=%0.3g, p=%0.3g',r,pp));
            xlabel(labe{1,hr}(1));
            ylabel(labe{1,hr}(2));
            
            subplot(2,2,hr+2)
            plot(widGauss(fitOK(:,hr)==1,hr),DDI{1}(fitOK(:,hr)==1,3-hr),'ko'); % Translation bandwidth 对应 Rotation DDI
            [r,pp] = plot_corr_line(widGauss(fitOK(:,hr)==1,hr),DDI{1}(fitOK(:,hr)==1,3-hr),'MethodOfCorr','Pearson','FittingMethod',2,'LineWidth',2);
            xlimit = xlim;
            ylimit = ylim;
            text((xlimit(2)-xlimit(1))/2+xlimit(1),ylimit(2)-(ylimit(2)-ylimit(1))/20,sprintf('r=%0.3g, p=%0.3g',r,pp));
            xlabel(labe{2,hr}(1));
            ylabel(labe{2,hr}(2));
        end
        figN = figN + 1;
    end

    function f2p21(debug) % 'AfterStim Spiral Tuning'
        if debug  ; dbstack;   keyboard;      end
        methods_of_select = {select_typical, 'Typical cells'; };
        
        AB_r = nan(length(loadM_data),2);
        AB_p = nan(length(loadM_data),2);
        spiral_index_AB = nan(length(loadM_data),2);
        modulation_site = 0;
        no_mod_site = 0;
        
        for n = 1:length(loadM_data)
            if ~isnan(group_result(n).AS_spit_resp(1)) && ~isnan(group_result(n).spit_resp(1)) && methods_of_select{1,1}(n)
                % only select after stim still have significant modulation unit
                if group_result(n).AS_p_stim < 0.05
                    % 1. person correlation same as cluster
                    for hr = 1:2
                        [AB_r(n,hr), AB_p(n,hr)] = corr(group_result(n).AS_spit_resp(hr,:)',group_result(n).spit_resp(hr,:)','type','Pearson');
                    end
                    
                    % spiral index difference
                    spiral_index_AB(n,:) = [group_result(n).spiral_index group_result(n).AS_spiral_index]; % before and after
                    
                    modulation_site = modulation_site + 1;
                else
                    no_mod_site = no_mod_site + 1;
                end
            end
        end
        
        % ploting
        set(figure(figN),'Position',[350,100, 1200,400], 'Name', 'Before/After Stim tuning Correlation');clf
        edges = [-1:0.2:1];
        for hr = 1:2
            subplot(1,3,hr)
            nbins = hist(AB_r(:,hr),edges);
            nbins_sig = hist(AB_r(AB_p(:,hr)<0.05,hr),edges);
            hold on
            bar(edges,nbins,1,'facecolor','w','edgecolor','k');
            bar(edges,nbins_sig,1,'facecolor',c{hr},'edgecolor','k');
            xlim([-1.1 1.1]);
            ylimt = max(nbins);
            ylimt1 = ylimt + ylimt/3;
            ylim([0 ylimt1]);
            xlabel('Correlation coefficient');
            axis square
            
            %  median correlation values and sign test，是否与0有显著区别
            median_cor = nanmedian(AB_r(:,hr));
            p_sign_test = signtest(AB_r(:,hr));
            
            % 箭头标记median值
            text(median_cor,ylimt+ylimt/4,'\downarrow\bf','fontsize',15,'horizontalalignment','center');
            if p_sign_test<0.001
                text(median_cor,ylimt1,'***\bf','fontsize',20,'horizontalalignment','center');
            elseif (p_sign_test>0.001 && p_sign_test<0.01)
                text(median_cor,ylimt1,'**\bf','fontsize',20,'horizontalalignment','center');
            elseif (p_sign_test>0.01 && p_sign_test<0.05)
                text(median_cor,ylimt1,'*\bf','fontsize',20,'horizontalalignment','center');
            end
            
            str = [sprintf('meadian = %.2g',median_cor) newline sprintf('p = %.2g',p_sign_test) newline  sprintf('N = %g',modulation_site)];
            text(-1, ylimt/2, str);
            
            ylabel('Cases');
            if hr == 1
                title('Transaltion plane');
            else
                title('Rotation plane');
            end
        end
        
        subplot(1,3,3)
        scatter(spiral_index_AB(:,1),spiral_index_AB(:,2),'MarkerEdgeColor','k','markerfacecolor',c{99});
        hold on
        plot([0 0],[-1 1],'k:');
        plot([-1 1],[0 0],'k:');
        SameScale(1)
        xlim([-1 1]);ylim([-1 1]);
        xlabel('Spiral Index Before Stim');
        ylabel('Spiral Index After Stim');
        
        % show the similar dot (same sign)
        patch([-1 0 0 -1],[-1 -1 0 0],[1,2,3,4],'edgecolor','none','facecolor',c{99},'facealpha',.2);
        patch([0 1 1 0],[0 0 1 1],[1,2,3,4],'edgecolor','none','facecolor',c{99},'facealpha',.2);
        
        % correlation
        [r,p] = corr(spiral_index_AB(:,1),spiral_index_AB(:,2),'type','pearson','row','complete');
        text(-0.9,0.9,sprintf('Pearson"s r = %1.3g',r));
        text(-0.9,0.75,sprintf('p = %.2g',p));
        title('Spiral Index');
        
        suptitle('Before-After Stim tuning Correlation')
        SetFigure(13);figN = figN+1;
    end

    function f2p22(debug) % 'ROC analysis and compare to d prime'
        if debug  ; dbstack;   keyboard;      end
        methods_of_select = {select_typical, 'Typical cells'; };
        
        d_prime = cell2mat({group_result.d_prime_coarse}'); % use d' on ±90°
        AuROC = cell2mat({group_result.AuROC}');
        d_prime = log(abs(d_prime(methods_of_select{1,1},:)));
        AuROC = AuROC(methods_of_select{1,1},:);
        
        % ploting
        set(figure(figN),'Position',[350,100, 1200,600], 'Name', 'ROC vs d prime');clf
        for hr = 1:2
            ax(hr) = subplot(1,2,hr);
            plot(d_prime(:,hr),AuROC(:,hr),'.','color',c{hr})
            axis square
            
            [r,p] = plot_corr_line(d_prime(:,hr),AuROC(:,hr),'MethodOfCorr','Pearson','FittingMethod',2,'LineWidth',2);
            str = ['Pearson correlation' newline sprintf('r=%2.2g, p=%2.2g',r,p)];
            text(-5,0.95,str);
            
            ylim([0.5 1])
            xlabel('log(|d prime|)');
            ylabel('ROC');
        end
        
        linkaxes([ax(1),ax(2)],'xy');
        suptitle('ROC vs d prime');
        SetFigure(13);figN = figN+1;
    end




%%
    function f3p1(debug) % 'Spiral Index'
        if debug  ; dbstack;   keyboard;      end
        methods_of_select = {select_typical, 'Typical cells'; };
        
        spiral_index = cell2mat({group_result.spiral_index}');
%         spiral_index = cell2mat({group_result.spiral_index_angle}');
%         spiral_index = cell2mat({group_result.SVM_index}');
        
        set(figure(figN),'Position',[100,70, 500,300], 'Name', ' Spiral Index');
        edges = linspace(-1,1,21);
        histogram(spiral_index(methods_of_select{1,1}),edges);
        [~,p] = ttest(spiral_index(methods_of_select{1,1}),0);
        text(-0.8,75,sprintf('p = %2.2g',p),'color','k');
        text(-0.7,60,sprintf('N = %g',length(spiral_index(methods_of_select{1,1}))));
        
        xlabel('Spiral index');
        ylabel('Cases');
        SetFigure(13);
        figN = figN+1;
        
        % 'Spiral Index among cell type'
        set(figure(figN),'Position',[610,70,500,900], 'Name', ' Spiral Index among cell type');
        
        % bar plot
        xbins_c = [];
        binsize = 0.2;
        xbins_c{1} = -1:binsize:1;
        xbins_c{2} = -1:binsize:1;
        xbins_c{3} = -1:binsize:1;
        
        rawx = [];
        for d = 4-loop_num % h cell, r cell, spiral cell
            hold on
            rawx{d} = spiral_index(([group_result.cell_type]' == d & methods_of_select{1,1}));
            
            subplot(3,1,d)
            %             histogram(temp,xbins_c{d},'facecolor',c{d},'normalization','probability');
            histogram(rawx{d},xbins_c{d},'facecolor',c{d});
            Nbin{d} = histcounts(rawx{d},xbins_c{d});
            Pbin{d} = histcounts(rawx{d},xbins_c{d},'normalization','probability');
           
            % mean
            mean_si(d) = nanmean(rawx{d});
            [~,p_ttest(d)] = ttest(rawx{d},0); % 与0比较
            
            % median
            median_si(d) = nanmedian(rawx{d});
            p_signtest(d) = signtest(rawx{d}); % 与0比较c
            
            
            ylimit = ylim;
            text(-1.1,ylimit(2),sprintf('median = %.2g, p = %.2g',median_si(d),p_signtest(d)),'color',c{d});
            hold on
            plot([median_si(d) median_si(d)],ylimit,'k-');
            xlim([-1.2 1.2]);
            xlabel('Spiral index');
            ylabel('Cases');
        end
%         legend('Spiral cell','Rotation cell','Translation cell');
        SetFigure(13);
        figN = figN+1;
        
        % group plot
        set(figure(figN),'Position',[800,70,500,400], 'Name', ' Spiral Index among cell type (group)');clf
        [h,x_centers_adj,hist_x] = plotGroupedBar(rawx,10,[-1 1],c,'count');
        xlim([-1.1 1.1])
        SetFigure(13);
        figN = figN+1;
        
        % stacked plot
        set(figure(figN),'Position',[1000,70,500,400], 'Name', ' Spiral Index among cell type (stacked)');clf
        bar(x_centers_adj(2,:),hist_x','stacked');
        xlim([-1.1 1.1]);
        hold on
        ylimit = ylim;
        for i = 1:3
            plot([mean_si(i) mean_si(i)],ylimit,'-','color',c{i});
        end
        SetFigure(13);
        figN = figN+1;  
        
        % 插值垂直画
        xedges_mid = xbins_c{1}+binsize/2;
        xedges_mid(end) = [];
        set(figure(figN),'Position',[1200,70,700,400], 'Name', ' Spiral Index among cell type (interpolation)');clf
        for i = 1:3 % h r spi
            % remove fist and last several zero
            select = [find(Pbin{i}~=0,1,'first'): find(Pbin{i}~=0,1,'last')];
            Pbin_adj = Pbin{i}(select);
            xedges_mid_adj = xedges_mid(select);
            
            % add first and last 0 for interpolation
            Pbin_adj = [0 Pbin_adj 0];
            xedges_mid_adj = [min(xedges_mid_adj)-binsize/2 xedges_mid_adj max(xedges_mid_adj)+binsize/2];
            xedges_mid_interp = min(xedges_mid_adj):0.01:max(xedges_mid_adj);
            Pbin_interp = interp1(xedges_mid_adj,Pbin_adj,xedges_mid_interp,'pchip');
            
            hold on
            plot(Pbin_interp+i,xedges_mid_interp,'-','color',c{i});
            plot(-Pbin_interp+i,xedges_mid_interp,'-','color',c{i}); % copy for other side
            
            % add cross = meadian
%             median_si = mean_si;
            plot([-0.1+i 0.1+i],[median_si(i) median_si(i)],'-','color',c{i}); % plot mean
            plot([i i],[median_si(i)-0.1 median_si(i)+0.1],'-','color',c{i});
            
            ylim([-1 1]);
            xlabel('Cases/N');
            ylabel('Spiral Index');
            
            text(i,-0.9,sprintf('meadian = %.2g', median_si(i)));
        end
        SetFigure(13);
        figN = figN+1;  
        
        
    end

    function f3p2(debug) % 'Firing rate/Vector-sum project on basic axis'
        if debug  ; dbstack;   keyboard;      end
        
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'}; % 默认T cell+R cell+ Spi cell
        
        fr_axis_ratio = cell2mat({group_result.fr_axis_ratio}');
        vec_axis_ratio = cell2mat({group_result.vec_axis_ratio}');
        d_prime_ratio = cell2mat({group_result.d_prime_ratio}');

        fr_axis_ratio = fr_axis_ratio(methods_of_select{1,1});
        vec_axis_ratio = vec_axis_ratio(methods_of_select{1,1});
        d_prime_ratio = d_prime_ratio(methods_of_select{1,1});
        
        % group together
        ratio_index{1} = d_prime_ratio;
        ratio_index{2} = fr_axis_ratio;
        ratio_index{3} = vec_axis_ratio;
        
        til{1} = 'd prime ratio';
        til{2} = 'FR on basic axix ratio';
        til{3} = 'Vector-sum projection on basic axix ratio';
        
        set(figure(figN),'Position',[100,70,1500,500], 'Name', 'Some Ratio for T-R');
        for i = 1:3
            subplot(1,3,i)
            histogram(log10(ratio_index{i}));
            str = ['Bias to R',repmat(32,1,20),'Bias to T' newline 'log10(',til{i}, ')'];
            xlabel(str)
            ylabel('Cases');
            axis tight
            xlimit = xlim;
            xlimit_max = max(abs(xlimit));
            xlim([-xlimit_max xlimit_max]);
            title(til{i});
        end
       
        SetFigure(13);
        figN = figN+1;
    end

    function f4p1(debug) % PCA for spiral tuning
        if debug  ; dbstack;   keyboard;      end
        
        methods_of_select = {select_typical, 'Typical cells'}; % 默认T cell+R cell+ Spi cell
        
        %%% true data fitting parameter
        fitOK = cell2mat({group_result.fitOK}'); %  if output.exitflag > 0 && cf_{s}.sigma >= 0.3 && cf_{s}.sigma < 500; fitOK = 1;
        for n = 1:length(group_result)
            fittingresult = group_result(n).fittingresult;
            if ~isempty(fittingresult)
                for hr = 1:2
                    a(n,hr) = fittingresult{hr}.a;
                    b(n,hr) = fittingresult{hr}.b;
                    pref(n,hr) = fittingresult{hr}.pref;
                    sigma(n,hr) = fittingresult{hr}.sigma;
                end
            else
                a(n,:) = nan(1,2);
                b(n,:) = nan(1,2);
                pref(n,:) = nan(1,2);
                sigma(n,:) = nan(1,2);
            end
        end
        
        sigma((sigma(:,1) > 10),1) = nan; % abnormal one！
        sigma((sigma(:,2) > 10),2) = nan; % abnormal one！
        select = logical([methods_of_select{1},methods_of_select{1}] & fitOK == true);
        fitparameter{1} = a;
        fitparameter{2} = b;
        fitparameter{3} = pref;
        fitparameter{4} = sigma;
        
        set(figure(figN),'Position',[100,600,500,500], 'Name', 'fitting parameters');
        for i = 1:4
            subplot(2,2,i)
            for hr = 1:2
                hold on
                histogram(fitparameter{i}(:,hr));
                if i == 1
                    title('a');
                elseif i == 2
                    title('b');
                elseif i == 3
                    title('pref');
                else
                    title('sigma');
                end
            end
        end
        figN = figN + 1;
        
        %%% group raw data
        T_response_all = [];
        response_all = [];
        response_trial = [];
        for n = 1:length(group_result)
            if methods_of_select{1}(n)
                % mean fr
                TR_component = repmat([-180:45:135]',2,1);
                spit_resp_temp = group_result(n).spit_resp; % mean firing rate: -180,  -135,  -90,  -45,  0,  45,  90,  135
                T_response = spit_resp_temp(1,:);
                R_response = spit_resp_temp(2,:);
                T_response_all = [T_response_all,T_response'];
                response_all_this = [T_response';R_response'];
                %             response_all_this([3,13]) = []; % remove the repeat 0 and 180 in rotation
                remove_this = logical(TR_component==0 | TR_component==-180);
                %             response_all_this(remove_this) = []; % remove the repeat 0 and 180 in both translation and rotation
                
                response_all = [response_all,response_all_this];
                
                %             % trial
                %             spit_resp_trial_temp = group_result(n).spit_resp_trial; % mean firing rate: -180,  -135,  -90,  -45,  0,  45,  90,  135
                %             if ~isnan(spit_resp_trial_temp{1}(1))
                %                 response_trial_this = [spit_resp_trial_temp{1}(:); spit_resp_trial_temp{2}(:)];
                %                 try
                %                 response_trial = [response_trial response_trial_this];
                %                 end
                %             end
            end
        end
        
        % remove nan
        remove_this = any(isnan(T_response_all));
        T_response_all(:,remove_this) = [];
        
        remove_this1 = any(isnan(response_all));
        response_all(:,remove_this1) = [];
        
        %         % Translation only
        %         % normalize
        %         T_response_all_zscore = zscore(T_response_all,[],2);
        %         [row,col]=size(T_response_all_zscore);
        %         c_mat = cov(T_response_all_zscore);
        %         [F,V]=eigs(c_mat);
        %         meanx = mean(T_response_all_zscore);
        %         tempx=repmat(meanx,row,1);
        %         score=(T_response_all_zscore-tempx)*F;
        %         pca=score(:,1:2);
        %         figure
        %         scatter(pca(:,1),pca(:,2));
        
        % Normalize
        response_all_zscore = zscore(response_all,[],2); % normalize
        [row,col]=size(response_all_zscore);
        
        % cell type
        cell_type = [group_result.cell_type]';
        cell_type = cell_type(methods_of_select{1});
        
        % PCA
        % 1. 分布计算
        c_mat = cov(response_all_zscore); % Covariancec 协方差矩阵
        [V,D] = eigs(c_mat); % V:特征向量组成的矩阵,D: 特征值组成的对角矩阵
        ev = diag(D); % 提取特征值
        ev = ev(:, end:-1:1); % eig计算出的特征值是升序的，这里手动倒序（W同理）
        W = V(:, end:-1:1); % 同上
        sum(W.*W, 1) % 可以验证每个特征向量各元素的平方和均为1
        Wr = W(:, 1:3); % 提取前3个主成分的特征向量
        Tr = response_all_zscore * Wr; % 新坐标空间的数据点
        
        meanx = mean(response_all_zscore);
        tempx = repmat(meanx,row,1);
        score1 = (response_all_zscore-tempx)*V;
        PCA_score = score1(:,1:3);
        
        % 2. 直接用函数pca
        [coeff, score, latent, tsquare] = pca(response_all_zscore);
        data_PCA =  response_all_zscore*coeff(:,1:10);
        latent1 = 100*latent./sum(latent); % 特征值
        cumsum(latent)./sum(latent)
        figure
        pareto(latent1); % pareto仅绘制累积分布的前95%, 图中的线表示的累计变量解释程度
        xlabel('PC');
        ylabel('Variance Explained (%)');
        response_pca = data_PCA(:,1:3);
        figure;
        biplot(coeff(:,1:3),'Scores',score(:,1:3));
        SetFigure(15);
        
        hr_condition = [-180,  -135,  -90,  -45,  0,  45,  90,  135];
        
        % 原始数据重构
        re_data = tempx + score(:,1:3)*coeff(:,1:3)';
        figure
        plot3(re_data(1,:),re_data(2,:),re_data(3,:),'.');
        
        %         %%% DCA, can not work!
        %         Xs{1} = response_all_zscore;
        %         Xs{2} = repmat([hr_condition zeros(1,8)]',1,size(response_all_zscore,2));
        %         Xs{3} = repmat([zeros(1,8) hr_condition]',1,size(response_all_zscore,2));
        %
        %         [U, dcovs] = dca(Xs);
        %
        %         x1 = U{1}(:,1)' * Xs{1};
        %         x2 = U{2}(:,1)' * Xs{2};
        %         x3 = U{3}(:,1)' * Xs{3};
        %         plot3(x1, x2, x3,'.k');
        %
        %         
        
        %%% plot PCA result
        size(PCA_score,1);
        for hr = 1:2
            pca_this = PCA_score((hr-1)*size(PCA_score,1)/2+1:hr*size(PCA_score,1)/2,:);
            pca_this(end+1,:) = pca_this(1,:); % repeat the first dot for T and R for a circle
            
            set(figure(figN),'Position',[100,70,500,500], 'Name', 'PCA');
            hold on
            % plot every condition for PCA
            plot3(pca_this(:,1),pca_this(:,2),pca_this(:,3),'o-','color',c{hr});
            xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
            
            % fitting a plane
            xyz0 = mean(pca_this,1);
            centerPlane = bsxfun(@minus, pca_this, xyz0);
            [U,S,V] = svd(centerPlane);
            aa(hr) = V(1,3);
            bb(hr) = V(2,3);
            cc(hr) = V(3,3);
            dd(hr) = -dot([aa(hr) bb(hr) cc(hr)],xyz0);
            
            % plot fitting plane
            xfit = min(pca_this(:,1))*1.2:0.1:max(pca_this(:,1))*1.2;
            yfit = min(pca_this(:,2))*1.2:0.1:max(pca_this(:,2))*1.2;
            [XFIT,YFIT] = meshgrid(xfit, yfit);
            ZFIT = -(dd(hr) + aa(hr) * XFIT + bb(hr) * YFIT) / cc(hr);
            surface(XFIT,YFIT,ZFIT,'facealpha',0.5,'facecolor',c{hr},'edgecolor','none');
            
            scale_f = 1;
            xlim([-15*scale_f 15*scale_f]);
            ylim([-15*scale_f 15*scale_f]);
            zlim([-15*scale_f 15*scale_f]);
            
            %             % annotation
            %             for ann = 1:length(hr_condition)
            %                 text(pca_this(ann,1),pca_this(ann,2),pca_this(ann,3),num2str(hr_condition(ann)),'color',c{hr});
            %             end
            
            %             % 投影in the same figure
            %             xlimit = xlim;
            %             ylimit = ylim;
            %             zlimit = zlim;
            %             plot3(pca_this(:,1),pca_this(:,2),zlimit(1).*ones(1,length(pca_this(:,1))),'o-','color',c{hr});
            %             plot3(pca_this(:,1),ylimit(2).*ones(1,length(pca_this(:,1))),pca_this(:,3),'o-','color',c{hr});
            %             plot3(xlimit(2).*ones(1,length(pca_this(:,1))), pca_this(:,2),pca_this(:,3),'o-','color',c{hr});
            
            xticks([0]); yticks([0]); zticks([0]);
            
            view(60, 10);
            SetFigure(15);
            grid on
            box on
            
            % 投影in the separate figure
            set(figure(figN+1),'Position',[650,70,500,500], 'Name', 'PCA projection');
            subplot(2,2,1); hold on
            plot(pca_this(:,1),pca_this(:,2),'o-','color',c{hr});
            xlabel('PC1'); ylabel('PC2'); axis square; xlim([-15 15]); ylim([-15 15]);box on
            subplot(2,2,2); hold on
            plot(pca_this(:,1),pca_this(:,3),'o-','color',c{hr});
            xlabel('PC1'); ylabel('PC3'); axis square; xlim([-15 15]); ylim([-15 15]);box on
            subplot(2,2,3); hold on
            plot(pca_this(:,2),pca_this(:,3),'o-','color',c{hr});
            xlabel('PC2'); ylabel('PC3'); axis square; xlim([-15 15]); ylim([-15 15]);box on
            SetFigure(15);
        end
        
        figN = figN + 2;
        
        % add cell type
        figure
        hold on
        for ct = 1:3
            plot3(coeff(cell_type==ct,1),coeff(cell_type==ct,2),coeff(cell_type==ct,3),'.','color',c{ct});
         end
        xlim([-0.25 0.25])
        ylim([-0.25 0.25])
        zlim([-0.25 0.25])
        axis square
        % plot fitting plane
        xfit = -0.25:0.01:0.25;
        yfit = -0.25:0.01:0.25;
        [XFIT,YFIT] = meshgrid(xfit, yfit);
        for hr = 1:2
            ZFIT = [];
            ZFIT = -(dd(hr) + aa(hr) * XFIT + bb(hr) * YFIT) / cc(hr);
            surface(XFIT,YFIT,ZFIT,'facealpha',0.5,'facecolor',c{hr},'edgecolor','none');
        end
        
        
        
        
        
        
        % PC1 contribution
        figure
        [sortPC1 sortIdx] = sort(abs(V(:,2)),'descend');
        hold on
        plot(abs(V(sortIdx,1)),'k.');
        plot(abs(V(sortIdx,2)),'.');
        plot(abs(V(sortIdx,3)),'b.');
        
        
%         %% fake neural data to PCA

    end


    function Show_individual_cell(~,~,h_line, select_for_this, couples)    % Show individual cell selected from the figure. @HH20150424
        
        if nargin < 5
            couples = 1; % HH20150515 Show coupled cells in the figure
        end
        
        h_marker = guidata(gcbo);
        
        allX = get(h_line,'xData');
        allY = get(h_line,'yData');
        n_single_group = length(allX)/couples;
        
        % ------- Recover the cell number -------
        
        if ismember('control',get(gcf,'currentModifier'))  % Select cell from file name and show in the figure. @HH20150527
            fileN = input('Which cell do you want from the figure?    ','s');
            available_cells = find(select_for_this);
            
            if fileN(1) == '#'  % HH20160905. Direct input the original cell #
                ori_cell_no = str2double(fileN(2:end));
                ind = sum(select_for_this(1:ori_cell_no));
            else
                
                nn = 1; iffind = [];
                while nn <= length(available_cells) && isempty(iffind) % Find the first match
                    iffind = strfind(group_result(available_cells(nn)).cellID{1}{1},fileN);
                    nn = nn + 1;
                end
                
                if ~isempty(iffind)
                    ind = nn - 1;
                else
                    fprintf('Are you kidding...?\n');
                    return
                end
                
                ori_cell_no = find(cumsum(select_for_this)==ind,1);  % 对应loadM_data or group_result中的顺序
            end
        else   % Select cell from figure
            pos = get(gca,'currentPoint'); posX = pos(1,1); posY = pos(1,2);
            [min_dis,ind] = min(abs(((posX-allX)/range(xlim)).^2+((posY-allY)/range(ylim)).^2));
            if min_dis > (range(xlim)^2+range(ylim)^2)/100 +inf ; return; end
            
            ind = mod(ind-1,n_single_group) + 1; % Deal with coupled cells
            ori_cell_no = find(cumsum(select_for_this)==ind,1);
        end
        
        % Plotting
        if ~isempty(h_marker) ; try delete(h_marker); catch ; end ;end
        
        all_inds = mod(1:length(allX),n_single_group) == mod(ind,n_single_group); % Deal with coupled cells
        h_marker = plot(allX(all_inds),allY(all_inds),'x','color','m','markersize',15,'linew',3);
        
        %         % Do plot
        %         Plot_HD([],[],ori_cell_no);
        
        % 暂时只show data
        loadM_data(ori_cell_no);
        path_temp = [];

        open_name_temp = group_result(ori_cell_no).SpiT_FILE;
        dot_loc = find(open_name_temp=='.'); % 找到.的位置
        open_name_temp = open_name_temp(1:dot_loc-1);
        
        c_loc = find(open_name_temp=='c'); % 找到c的位置
        monkey_this = str2double(open_name_temp(2:c_loc-1));

        if monkey_this == 7
            path_temp = 'Z:\Data\Tempo\Batch\Ringbell_SpiT';
        elseif monkey_this == 16
            path_temp = 'Z:\Data\Tempo\Batch\Arthas_SpiT';
        else
            keyboard
        end
                
        open_name = strcat(path_temp,'\',open_name_temp,'_5_SpiT_fig_10.bmp');
        figure(9999);clf;
        d = imread(open_name);
        d = imrotate(d,0);
        imshow(d,'InitialMagnification','fit');
        guidata(gcbo,h_marker);
    end

end