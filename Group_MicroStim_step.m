% 202204227 需要大改得地方，将2AFC
% task不在当作独立的一个任务，而是将2AFC-T和2AFC-R合并成一个新的任务，方便后续直接调取2-AFC数据！！






function function_handles_stim = Group_MicroStim(loadM_data,handles)
set(0,'defaultFigureRenderer','painters');
set(0,'defaultFigureRendererMode','manual');

%% Extract some data to global value for easier access
% 将常用的参数提取出来
global group_result;

% Add flexible monkey mask here (but I've still decided to choose monkey for analysis below). HH20150723
monkey_included_for_loading = [7 16];

% select the cell num
too_small_thred = [];
too_small_thred_File = [];

% microstim data
for n = 1:length(loadM_data)
    main_depth_loc = [];
    main_depth_loc = [(loadM_data(n).Data.Depth_code)]==0;
    select_stim_data{1} = loadM_data(n).Data(main_depth_loc).SpikeChan(1).centerEye_Protocol{8};  % stimCueS
    select_stim_data{2} = loadM_data(n).Data(main_depth_loc).SpikeChan(1).centerEye_Protocol{9};  % stimHonly
    select_stim_data{3} = loadM_data(n).Data(main_depth_loc).SpikeChan(1).centerEye_Protocol{10};  % stimRonly
    select_stim_data{4} = loadM_data(n).Data(main_depth_loc).SpikeChan(1).centerEye_Protocol{11};  % stimCo
    select_stim_data{5} = loadM_data(n).Data(main_depth_loc).SpikeChan(1).centerEye_Protocol{12};  % stimCoH
    select_stim_data{6} = loadM_data(n).Data(main_depth_loc).SpikeChan(1).centerEye_Protocol{13};  % stimCoR
    
    if sum(cellfun(@isempty,select_stim_data)) == 6  % all empty
        % 找出是否在别的深度也做了 stim task
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
            for nn1 = 1:length(depth_temp)
                select_stim_data{1} = loadM_data(n).Data(depth_ori == depth_temp(nn1)).SpikeChan(1).centerEye_Protocol{7};  % stimCueS
                select_stim_data{2} = loadM_data(n).Data(depth_ori == depth_temp(nn1)).SpikeChan(1).centerEye_Protocol{8};  % stimHonly
                select_stim_data{3} = loadM_data(n).Data(depth_ori == depth_temp(nn1)).SpikeChan(1).centerEye_Protocol{9};  % stimRonly
                select_stim_data{4} = loadM_data(n).Data(depth_ori == depth_temp(nn1)).SpikeChan(1).centerEye_Protocol{10};  % stimCo
                select_stim_data{5} = loadM_data(n).Data(depth_ori == depth_temp(nn1)).SpikeChan(1).centerEye_Protocol{11};  % stimCoH
                select_stim_data{6} = loadM_data(n).Data(depth_ori == depth_temp(nn1)).SpikeChan(1).centerEye_Protocol{12};  % stimCoR
                
                if sum(cellfun(@isempty,select_stim_data)) < 5  % not all empty
                    break % 找到有
                end
            end
        end
    end
    
    for t = 1:6 % 6 types of stimulation
        % stimCueS stimHonly stimRonly stimCo stimCoH stimCoR
        if isfield(select_stim_data{t},'FILE') && size(select_stim_data{t},2) && select_stim_data{t}.repetition >= 5
            group_result(n).stim_FILE{t} = select_stim_data{t}.FILE;
            group_result(n).Stim_repN(t) = select_stim_data{t}.repetition;
            group_result(n).coherence{t} = select_stim_data{t}.unique_coherence;
            group_result(n).Stim_unique_heading{t} = select_stim_data{t}.unique_degree{1}';
            group_result(n).Stim_unique_rotation{t} = select_stim_data{t}.unique_degree{2}';
            group_result(n).Stim_unique_mt{t} = select_stim_data{t}.unique_degree{3}';
            group_result(n).unique_switch_index{t} = select_stim_data{t}.unique_switch_index;
            group_result(n).Is_step{t} = select_stim_data{t}.Is_step;
            group_result(n).step_num{t} = select_stim_data{t}.step_num;
            group_result(n).stimAmp{t} = select_stim_data{t}.stimAmp;
            
            % change change index between nostim and stim:choice_change_index = sign * sqrt(chi2):
            % 对于MT(ind=3)时候，越负代表电刺激后（相对电刺激前）选择Rotation target变多了，越正代表电刺激后选择Translation target变多了，=0代表电刺激没有作用
            group_result(n).choice_change_index_diff{t} = select_stim_data{t}.choice_change_index_diff{:}; % 每一个choice_change_index 都有"Left"和"Right"两种choice的变化
            group_result(n).choice_change_index_in0{t} = select_stim_data{t}.choice_change_index_in0{:};
            group_result(n).choice_change_index_around0{t} = select_stim_data{t}.choice_change_index_around0{:};
            
            % chi-square test between nostim and stim
            group_result(n).chi_square_p_diff{t} = select_stim_data{t}.chi_square_p_diff{:};
            group_result(n).chi_square_p_in0{t} = select_stim_data{t}.chi_square_p_in0{:};
            group_result(n).chi_square_p_around0{t} = select_stim_data{t}.chi_square_p_around0{:};
            
            % psychometric curve condition list
            unique_condition_everysite_h{n,t} = (group_result(n).Stim_unique_heading{t});
            unique_condition_everysite_r{n,t} = (group_result(n).Stim_unique_rotation{t});
            unique_condition_everysite_co{n,t} = (group_result(n).coherence{t});
            unique_condition_everysite_mt{n,t} = (group_result(n).Stim_unique_mt{t});
            
            psy_perf_temp = squeeze(select_stim_data{t}.psy_perf);
            group_result(n).psy_perf{t} = psy_perf_temp(:,:,group_result(n).step_num{t}); % {n,k,ind,step} % 列：h r mt, 行:n=1电刺激前，n=2电刺激后
            
            psy_perf_true_temp = squeeze(select_stim_data{t}.psy_perf_true50); % PSE in 50% CR for 4 target
            group_result(n).psy_perf_true50{t} = psy_perf_true_temp(:,:,group_result(n).step_num{t}); % {n,k,ind,step} % 列：h r mt, 行:n=1电刺激前，n=2电刺激后
            % psy_perf_true25 in 'group_result(n).Stim_Bias_psy_true{t,ind}'

            group_result(n).correct_rate{t} = select_stim_data{t}.correct_rate;
            group_result(n).correct_rate_true{t} = select_stim_data{t}.correct_rate_true;
            
            % eyetrace
            group_result(n).LReye{t} = select_stim_data{t}.plot_LR;
            group_result(n).eye_data_VisualOn_X{t} = select_stim_data{t}.eye_data_VisualOn_X;
            group_result(n).eye_data_VisualOn_Y{t} = select_stim_data{t}.eye_data_VisualOn_Y;
            group_result(n).stimOnBin{t} = select_stim_data{t}.stimOnBin;
            group_result(n).stimOffBin{t} = select_stim_data{t}.stimOffBin;
            group_result(n).stimOnBin_trial{t} = select_stim_data{t}.stimOnBin_trial;
            group_result(n).stimOffBin_trial{t} = select_stim_data{t}.stimOffBin_trial;
            
            % ************************************** data seperated by control and stim **********************************************
            % ************************************************************************************************************************
            % choice array
            for nn = 1:2 % control and stim
                for ind = 1:2 % Heading or Rotation
                    group_result(n).choice_array{t}{ind,nn} = select_stim_data{t}.choice_array{nn,1,ind}; % eveyy cell, R-L-CCW-CW
                    group_result(n).choice_array2{t}{ind,nn} = select_stim_data{t}.choice_array2{nn,1,ind}; % every cell, R-CW-L-CCW
                    group_result(n).mt_right{t}{ind,nn} = select_stim_data{t}.mt_right_pro{nn,1,ind};
                    group_result(n).mt_wrong{t}{ind,nn} = select_stim_data{t}.mt_wrong_pro{nn,1,ind};
                    group_result(n).choose_all{t}{ind,nn} = select_stim_data{t}.choose_all{nn,1,ind};
                end
            end
            % calculate the choice change around zero: stim / crtl
            if t == 1|| t == 4 % only in 4 target task
                heading_num = length(group_result(n).Stim_unique_heading{t});
                rotation_num = length(group_result(n).Stim_unique_rotation{t});
                heading_around0 = [floor(heading_num/2):floor((heading_num+1)/2)+1];
                rotation_around0 = [floor(rotation_num/2),floor((rotation_num+1)/2)+1]; % 0°与heading part相同，不需要重复取
                
                ctrl_choose_trans = sum(sum(group_result(n).choice_array{t}{1,1}(1:2,heading_around0))) + sum(sum(group_result(n).choice_array{t}{2,1}(1:2,rotation_around0)));
                ctrl_choose_rot = sum(sum(group_result(n).choice_array{t}{1,1}(3:4,heading_around0))) + sum(sum(group_result(n).choice_array{t}{2,1}(3:4,rotation_around0)));
                stim_choose_trans = sum(sum(group_result(n).choice_array{t}{1,2}(1:2,heading_around0))) + sum(sum(group_result(n).choice_array{t}{2,2}(1:2,rotation_around0)));
                stim_choose_rot = sum(sum(group_result(n).choice_array{t}{1,2}(3:4,heading_around0))) + sum(sum(group_result(n).choice_array{t}{2,2}(3:4,rotation_around0)));
                group_result(n).choice_change_around0{t} = [stim_choose_trans/ctrl_choose_trans stim_choose_rot/ctrl_choose_rot];
            else
                group_result(n).choice_change_around0{t} = nan(1,2);
            end
            
            % ********************************************** data seperated by ind1,2,3 **********************************************
            % ************************************************************************************************************************
            % add step analysis  Lwh 20191211,  include all step (original data only the last step)
            % 每一列是step, 每一行是h, r, mt
            for ind = 1:3 % h, r, mt
                % *********************************** traditional psy curve ***********************************
                group_result(n).Stim_Bias_psy{t,ind} = cell2mat(squeeze(select_stim_data{t}.Bias_psy(:,:,ind,:))); % 每一行中第一小行是电刺激前的bias，第二小行是电刺激后的bias (cell2mat后自动清除了空置的step)
                group_result(n).Stim_Thresh_psy{t,ind} = cell2mat(squeeze(select_stim_data{t}.Thresh_psy(:,:,ind,:)));
                
                group_result(n).diffBias{t}(ind,:) =  group_result(n).Stim_Bias_psy{t,ind}(2,:) - group_result(n).Stim_Bias_psy{t,ind}(1,:); % after-before  >0: bias to left; <0: bias to right
                group_result(n).diffThresh{t}(ind,:) = group_result(n).Stim_Thresh_psy{t,ind}(2,:) - group_result(n).Stim_Thresh_psy{t,ind}(1,:);
                
                % normalize PSE shift: divided by the threshold under non-stim trials   (Xuefei Yu, 2018)
                if group_result(n).Stim_Thresh_psy{t,ind}(1,:)<0.1 % check! threshld is too small! 主要在rotation task中
                    group_result(n).nor_diffBias{t}(ind,:) = nan;
                    group_result(n).nor_PSE{t}(ind,:) = nan;
                else
                    group_result(n).nor_diffBias{t}(ind,:) = group_result(n).diffBias{t}(ind,:) ./ group_result(n).Stim_Thresh_psy{t,ind}(1);
                    group_result(n).nor_PSE{t}(ind,:) = group_result(n).Stim_Bias_psy{t,ind}(1) ./ group_result(n).Stim_Thresh_psy{t,ind}(1); % normalize PSE for control trial
                end

                group_result(n).p_mu{t}(ind,:) = cell2mat(squeeze(select_stim_data{t}.p_mu(:,ind,:)))';
                group_result(n).p_theta{t}(ind,:) = cell2mat(squeeze(select_stim_data{t}.p_theta(:,ind,:)))';
                
                % Psy_correct直接找final step的值，暂时不管前面的step
                group_result(n).Stim_Psy_correct{t,ind} = select_stim_data{t}.Psy_correct(:,:,ind,end); % 每一行中第一小行是电刺激前的，第二小行是电刺激后
                
                % *********************************** psy curve in 4 targets, use PSE in 25% CR and threhold in 50% ***********************************
                %                 group_result(n).Stim_Bias_psy_true50{t,ind} = cell2mat(squeeze(select_stim_data{t}.Bias_psy_true50(:,:,ind,:))); % 每一行中第一小行是电刺激前的bias，第二小行是电刺激后的bias (cell2mat后自动清除了空置的step)
                group_result(n).Stim_Bias_psy_true{t,ind} = cell2mat(squeeze(select_stim_data{t}.Bias_psy_true25(:,:,ind,:)));
                group_result(n).Stim_Thresh_psy_true{t,ind} = cell2mat(squeeze(select_stim_data{t}.Thresh_psy_true50(:,:,ind,:)));
                
                %                 group_result(n).diffBias_true50{t}(ind,:) =  group_result(n).Stim_Bias_psy_true50{t,ind}(2,:) - group_result(n).Stim_Bias_psy_true50{t,ind}(1,:); % after-before  >0: bias to left; <0: bias to right
                group_result(n).diffBias_true{t}(ind,:) =  group_result(n).Stim_Bias_psy_true{t,ind}(2,:) - group_result(n).Stim_Bias_psy_true{t,ind}(1,:); % after-before  >0: bias to left; <0: bias to right
                group_result(n).diffThresh_true{t}(ind,:) = group_result(n).Stim_Thresh_psy_true{t,ind}(2,:) - group_result(n).Stim_Thresh_psy_true{t,ind}(1,:);
                group_result(n).nor_diffBias_true{t}(ind,:) = group_result(n).diffBias_true{t}(ind,:) ./ group_result(n).Stim_Thresh_psy_true{t,ind}(1);
                group_result(n).p_mu_true{t}(ind,:) = cell2mat(squeeze(select_stim_data{t}.p_mu_true(:,ind,:)))';
                group_result(n).p_theta_true{t}(ind,:) = cell2mat(squeeze(select_stim_data{t}.p_theta_true(:,ind,:)))';
                group_result(n).Stim_Psy_correct_true{t,ind} = select_stim_data{t}.Psy_correct_true(:,:,ind,end); % 每一行中第一小行是电刺激前的，第二小行是电刺激后的
                % ***************************************************************************************************************************************
                
                % change change proportion: % 第一列L/CW/T, 第二列R/CCW/R,
                % 第一行heading task, 第二行rotation task,第三行motion type
                % stim - ctrl / ctrl
                group_result(n).choice_change_pro_diff{t}(ind,:) = select_stim_data{t}.choice_change_pro_diff{ind}; % (stimR-ctrlR)/ctrlR, (stimT-ctrlT)/ctrlT
                group_result(n).choice_change_pro_in0{t}(ind,:) = select_stim_data{t}.choice_change_pro_in0{ind};
                group_result(n).choice_change_pro_around0{t}(ind,:) = select_stim_data{t}.choice_change_pro_around0{ind};
                
                
                % time course of stim effect
                %  slide_positive_pro{n,k,ind}(slide_num) = sum(select_this_slide & choice_LEFT_2T{ind}) / sum(select_this_slide & (choice_LEFT_2T{ind}| choice_RIGHT_2T{ind}));
                group_result(n).slide_mid_trial{t,ind} = [select_stim_data{t}.slide_mid_trial{1,:,ind};select_stim_data{t}.slide_mid_trial{2,:,ind}]; % [control ; stim]
                group_result(n).slide_positive_pro{t,ind} = [select_stim_data{t}.slide_positive_pro{1,:,ind};select_stim_data{t}.slide_positive_pro{2,:,ind}]; % [control ; stim]
                group_result(n).slide_negative_pro{t,ind} = [select_stim_data{t}.slide_negative_pro{1,:,ind};select_stim_data{t}.slide_negative_pro{2,:,ind}]; % [control ; stim]
                group_result(n).slide_CR{t,ind} = [select_stim_data{t}.slide_CR{1,:,ind};select_stim_data{t}.slide_CR{2,:,ind}]; % [control ; stim]
            end

            % ************************************************************************************************************************
            % ---------------------------------------------------------------------------------------------------------------
            % remove the coarse task rotation(4 6) data in stars_type==1 (triangle stimulus)
            % stimCueS stimHonly stimRonly stimCo stimCoH stimCoR
            %             if group_result(n).stars_type == 0
            %                 keyboard
            %             end
            
            if group_result(n).stars_type == 1
                if t == 4 || t == 6
                    group_result(n).Stim_unique_rotation{t} = nan(size(group_result(n).Stim_unique_rotation{t}));
                    group_result(n).Stim_unique_mt{t} = nan(size(group_result(n).Stim_unique_mt{t}));
                    group_result(n).unique_switch_index{t} = 1;
                    
                    unique_condition_everysite_r{n,t} = group_result(n).Stim_unique_rotation{t};
                    unique_condition_everysite_co{n,t} = group_result(n).coherence{t};
                    unique_condition_everysite_mt{n,t} = group_result(n).Stim_unique_mt{t};
                    
                    group_result(n).choice_change_index_diff{t}([2,3]) = nan(1,2);
                    group_result(n).choice_change_index_in0{t}([2,3]) = nan(1,2);
                    group_result(n).choice_change_index_around0{t}([2,3]) = nan(1,2);
                    
                    group_result(n).chi_square_p_diff{t}([2,3]) = nan(1,2);
                    group_result(n).chi_square_p_in0{t}([2,3]) = nan(1,2);
                    group_result(n).chi_square_p_around0{t}([2,3]) = nan(1,2);
                    
                    group_result(n).mt_ctrl_correct_rate{t} = nan;
                    group_result(n).mt_stim_correct_rate{t} = nan;
                    group_result(n).TR_ctrl_correct_rate{t}(2) = nan; % only remove the R task CR
                    group_result(n).TR_stim_correct_rate{t}(2) = nan;
                    group_result(n).ctrl_CR{t}([2,3]) = nan(1,2);
                    group_result(n).stim_CR{t}([2,3]) = nan(1,2);
                    group_result(n).correct_rate{t}(:,2:3) = nan(2,2);
                    group_result(n).correct_rate_true{t}(:,2:3) = nan(2,2);
                    
                    % change to nan
                    for ind = 2:3 % r, mt
                        group_result(n).Stim_Psy_correct{t,ind} = {NaN;NaN};
                        group_result(n).Stim_Psy_correct_true{t,ind} = {NaN;NaN};
                        group_result(n).Stim_Bias_psy{t,ind} = nan(2,1);
                        group_result(n).Stim_Bias_psy_true{t,ind} = nan(2,1);
                        group_result(n).Stim_Thresh_psy{t,ind} = nan(2,1);
                        group_result(n).Stim_Thresh_psy_true{t,ind} = nan(2,1);
                        group_result(n).diffBias{t}(ind,:) = nan;
                        group_result(n).diffBias_true{t}(ind,:) = nan;
                        group_result(n).diffThresh{t}(ind,:) = nan;
                        group_result(n).diffThresh_true{t}(ind,:) = nan;
                        group_result(n).nor_diffBias{t}(ind,:) = nan;
                        group_result(n).nor_diffBias_true{t}(ind,:) = nan;
                        group_result(n).nor_PSE{t}(ind,:) = nan;
                        group_result(n).p_mu{t}(ind,:) = nan;
                        group_result(n).p_theta{t}(ind,:) = nan;
                        group_result(n).p_mu_true{t}(ind,:) = nan;
                        group_result(n).p_theta_true{t}(ind,:) = nan;
                        group_result(n).choice_change_pro_diff{t}(ind,:) = nan(1,2);
                        group_result(n).choice_change_pro_in0{t}(ind,:) = nan(1,2);
                        group_result(n).choice_change_pro_around0{t}(ind,:) = nan(1,2);
                        group_result(n).slide_mid_trial{t,ind} = nan(2,1);
                        group_result(n).slide_positive_pro{t,ind}= nan(2,1);
                        group_result(n).slide_negative_pro{t,ind} = nan(2,1);
                        group_result(n).slide_CR{t,ind} = nan(2,1);
                    end
                    
                    for nn = 1:2 % control and stim
                        ind1 = 2;
                        group_result(n).choice_array{t}{ind1,nn} = nan(size(group_result(n).choice_array{t}{ind1,nn}));
                        group_result(n).choice_array2{t}{ind1,nn} = nan(size(group_result(n).choice_array2{t}{ind1,nn}));
                        group_result(n).mt_right{t}{ind1,nn} = nan(size(group_result(n).mt_right{t}{ind1,nn}));
                        group_result(n).mt_wrong{t}{ind1,nn} = nan(size(group_result(n).mt_wrong{t}{ind1,nn}));
                        group_result(n).choose_all{t}{ind1,nn} = nan(size(group_result(n).choose_all{t}{ind1,nn}));
                    end
                end
            end
            % ---------------------------------------------------------------------------------------------------------------
            

        else
            group_result(n).stim_FILE{t} = nan;
            group_result(n).Stim_repN(t) = nan;
            group_result(n).coherence{t} = nan;
            group_result(n).Stim_unique_heading{t} = nan(1,7);
            group_result(n).Stim_unique_rotation{t} = nan(1,7);
            group_result(n).Stim_unique_mt{t} = nan(1,7);
            group_result(n).unique_switch_index{t} = nan;
            group_result(n).Is_step{t} = nan;
            group_result(n).step_num{t} = nan;
            group_result(n).stimAmp{t} = nan;
            
            group_result(n).diffBias{t} = NaN(3,1);
            group_result(n).diffThresh{t} = NaN(3,1);
            group_result(n).nor_diffBias{t} = NaN(3,1);
            group_result(n).nor_PSE{t} = NaN(3,1);
            group_result(n).p_mu{t} = NaN(3,1);
            group_result(n).p_theta{t} = NaN(3,1);
            
            group_result(n).diffBias_true{t} = NaN(3,1);
            group_result(n).diffThresh_true{t} = NaN(3,1);
            group_result(n).nor_diffBias_true{t} = NaN(3,1);
            group_result(n).p_mu_true{t} = NaN(3,1);
            group_result(n).p_theta_true{t} = NaN(3,1);
            
            % psychometric curve condition list
            unique_condition_everysite_h{n,t} = nan(1,7);
            unique_condition_everysite_r{n,t} = nan(1,7);
            unique_condition_everysite_co{n,t} = nan(1,7);
            unique_condition_everysite_mt{n,t} = nan(1,7);
            
            group_result(n).psy_perf{t} = mat2cell(nan(2,6),[1 1],[2 2 2]);
            group_result(n).psy_perf_true50{t} = mat2cell(nan(2,6),[1 1],[2 2 2]);
%             group_result(n).psy_perf_true25{t} = mat2cell(nan(2,6),[1 1],[2 2 2]);
            
            for ind = 1:3
                group_result(n).Stim_Psy_correct{t,ind} = {NaN;NaN};
                group_result(n).Stim_Bias_psy{t,ind} = nan(1,2);
                group_result(n).Stim_Thresh_psy{t,ind} = nan(2,1);
                group_result(n).slide_mid_trial{t,ind} = nan(2,1);
                group_result(n).slide_positive_pro{t,ind} = nan(2,1);
                group_result(n).slide_negative_pro{t,ind} = nan(2,1);
                group_result(n).slide_CR{t,ind} = nan(2,1);
                
                group_result(n).Stim_Psy_correct_true{t,ind} = {NaN;NaN};
                group_result(n).Stim_Bias_psy_true{t,ind} = nan(1,2);
                group_result(n).Stim_Thresh_psy_true{t,ind} = nan(2,1);
            end
            
            group_result(n).choice_array{t} = mat2cell(nan(8,14),[4 4],[7 7]);
            group_result(n).choice_array2{t} = mat2cell(nan(8,14),[4 4],[7 7]);
            group_result(n).mt_right_pro{t} = mat2cell(nan(2,14),[1 1],[7 7]);
            group_result(n).mt_wrong_pro{t} = mat2cell(nan(2,14),[1 1],[7 7]);
            group_result(n).choose_all{t} = mat2cell(nan(2,14),[1 1],[7 7]);
            group_result(n).mt_ctrl_correct_rate{t} = nan;
            group_result(n).mt_stim_correct_rate{t} = nan;
            
            group_result(n).choice_change_pro_diff{t} = nan(3,2);
            group_result(n).choice_change_pro_in0{t} = nan(3,2);
            group_result(n).choice_change_pro_around0{t} = nan(3,2);
            
            group_result(n).choice_change_index_diff{t} = nan(1,3); % 每一个choice_change_index 都有"Left"和"Right"两种choice的变化
            group_result(n).choice_change_index_in0{t} = nan(1,3);
            group_result(n).choice_change_index_around0{t} = nan(1,3);
            group_result(n).chi_square_p_diff{t} = nan(1,3);
            group_result(n).chi_square_p_in0{t} =nan(1,3);
            group_result(n).chi_square_p_around0{t} = nan(1,3);
            
            group_result(n).choice_change_around0{t} = nan(1,2);
            
            group_result(n).correct_rate{t} = nan(2,3);
            group_result(n).correct_rate_true{t} = nan(2,3);
            
            group_result(n).mt_ctrl_correct_rate{t} = nan;
            group_result(n).mt_stim_correct_rate{t} = nan;
            group_result(n).TR_ctrl_correct_rate{t} = [nan nan];
            group_result(n).TR_stim_correct_rate{t} = [nan nan];
            
            group_result(n).ctrl_CR{t} = nan(1,3);
            group_result(n).stim_CR{t} = nan(1,3); 

            group_result(n).LReye{t} = nan;
            group_result(n).eye_data_VisualOn_X{t} = nan;
            group_result(n).eye_data_VisualOn_Y{t} = nan;
            group_result(n).stimOnBin{t} = nan;
            group_result(n).stimOffBin{t} = nan;
            group_result(n).stimOnBin_trial{t} = nan;
            group_result(n).stimOffBin_trial{t} = nan;  
        end
    end
end


% step analysis
Is_step = {group_result.Is_step}';
Is_step = cellfun(@(x) cell2mat(x')',Is_step,'UniformOutput',0);
Is_step = cell2mat(Is_step);

step_num_all = {group_result.step_num}'; % the maximum step
step_num_all = cellfun(@(x) cell2mat(x')',step_num_all,'UniformOutput',0);
step_num_all = cell2mat(step_num_all);
step_num_max = max(max(step_num_all));
step_num_min = min(min(step_num_all));


% Reorganize data into matrices which are easier to plot and compare
% stimCueS stimHonly stimRonly stimCo stimCoH stimCoR

% 区分每一个step（每一行），如果该细胞没有跑step analysis，仍然会有很多step的数据（按照最大可以分的step来整理data的），不过前面的step都是nan，只有最后一个step才有值
% step_num_all = 1时，rep=5
% step_num_all = 2时，rep=6
% 。。。
for n = 1:length(loadM_data)
    for t = 1:6
        % 区分每一种stim task/ stim protocol
        if isnan(step_num_all(n,t))
            last_step = 1;
        else
            last_step = step_num_all(n,t);
        end
        % stimAmp
        stimAmp(n,t) = group_result(n).stimAmp{t}; % 每一个stim type的电刺激强度 （理论上同一个cell不同task的电刺激强度保持一致）
        
        % prepare for step analysis，每一行一个step，每一列h,r,spi
        diffBias_task_step_temp{t}{n} = nan(step_num_max,3);
        diffThresh_task_step_temp{t}{n} = nan(step_num_max,3);
        nor_diffBias_task_step_temp{t}{n} = nan(step_num_max,3);
        diffBias_p_task_step_temp{t}{n} = nan(step_num_max,3);
        diffThresh_p_task_step_temp{t}{n} = nan(step_num_max,3);
        
        if Is_step(n,t)==1
            diffBias_task_step_temp{t}{n}(1:last_step,:) =  group_result(n).diffBias{t}'; % 每一列h,r,mt，% after-before  >0: bias to left; <0: bias to right
            diffThresh_task_step_temp{t}{n}(1:last_step,:) =  group_result(n).diffThresh{t}';
            nor_diffBias_task_step_temp{t}{n}(1:last_step,:) = group_result(n).nor_diffBias{t}';
            diffBias_p_task_step_temp{t}{n}(1:last_step,:) = group_result(n).p_mu{t}';
            diffThresh_p_task_step_temp{t}{n}(1:last_step,:) = group_result(n).p_theta{t}';
        else
            diffBias_task_step_temp{t}{n}(1,:) =  group_result(n).diffBias{t}'; % 每一列h,r,mt
            diffThresh_task_step_temp{t}{n}(1,:) =  group_result(n).diffThresh{t}';
            nor_diffBias_task_step_temp{t}{n}(1,:) = group_result(n).nor_diffBias{t}';
            diffBias_p_task_step_temp{t}{n}(1,:) = group_result(n).p_mu{t}';
            diffThresh_p_task_step_temp{t}{n}(1,:) = group_result(n).p_theta{t}';
        end
        
        % the final step, original data
        % if Is_step ~= 1, that:diffBias_task = diffBias_task_step_temp
        if Is_step(n,t) == 1 % do the step analysis
            diffBias_task{t}(n,:) = diffBias_task_step_temp{t}{n}(last_step,:); % after-before  >0: bias to left; <0: bias to right
            diffThresh_task{t}(n,:) = diffThresh_task_step_temp{t}{n}(last_step,:);
            nor_diffBias_task{t}(n,:) = nor_diffBias_task_step_temp{t}{n}(last_step,:);
            diffBias_p_task{t}(n,:) = diffBias_p_task_step_temp{t}{n}(last_step,:);
            diffThresh_p_task{t}(n,:) =  diffThresh_p_task_step_temp{t}{n}(last_step,:);
        else % not do the step analysis
            diffBias_task{t}(n,:) = diffBias_task_step_temp{t}{n}(1,:);
            diffThresh_task{t}(n,:) = diffThresh_task_step_temp{t}{n}(1,:);
            nor_diffBias_task{t}(n,:) = nor_diffBias_task_step_temp{t}{n}(1,:);
            diffBias_p_task{t}(n,:) = diffBias_p_task_step_temp{t}{n}(1,:);
            diffThresh_p_task{t}(n,:) =  diffThresh_p_task_step_temp{t}{n}(1,:);
        end
        
        for ind = 1:3 % separate h r mt
            diffBias_task_step{t,ind}(n,:) = diffBias_task_step_temp{t}{n}(:,ind)';
            diffThresh_task_step{t,ind}(n,:) = diffThresh_task_step_temp{t}{n}(:,ind)';
            nor_diffBias_task_step{t,ind}(n,:) = nor_diffBias_task_step_temp{t}{n}(:,ind)';
            diffBias_p_task_step{t,ind}(n,:) = diffBias_p_task_step_temp{t}{n}(:,ind)';
            diffThresh_p_task_step{t,ind}(n,:) = diffThresh_p_task_step_temp{t}{n}(:,ind)';
        end
        
        % CR
        ctrl_CR{t}(n,:) =  group_result(n).correct_rate{t}(1,:);
        stim_CR{t}(n,:) =  group_result(n).correct_rate{t}(2,:);
        ctrl_CR_true{t}(n,:) =  group_result(n).correct_rate_true{t}(1,:);
        stim_CR_true{t}(n,:) =  group_result(n).correct_rate_true{t}(2,:);

        
        % normalize PSE for control trial
        nor_PSE{t}(n,:) = group_result(n).nor_PSE{t}';
        
        %  ind=3时第一列是R chocie变化，第二列是T choice变化
        for ind = 1:3 % separate h r mt 
            choice_change_pro_diff{t,ind}(n,:) = group_result(n).choice_change_pro_diff{t}(ind,:); % change change proportion: % 第一列L/CW/T, 第二列R/CCW/R,  第一行heading task, 第二行rotation task,第三行motion type
            choice_change_pro_in0{t,ind}(n,:) = group_result(n).choice_change_pro_in0{t}(ind,:);
            choice_change_pro_around0{t,ind}(n,:) = group_result(n).choice_change_pro_around0{t}(ind,:);
        end
        
        % chi-square test between nostim and stim
        chi_square_p_diff{t}(n,:) = group_result(n).chi_square_p_diff{t};
        chi_square_p_in0{t}(n,:) = group_result(n).chi_square_p_in0{t};
        chi_square_p_around0{t}(n,:) = group_result(n).chi_square_p_around0{t};
        
        % change change index between nostim and stim
        choice_change_index_diff{t}(n,:) = group_result(n).choice_change_index_diff{t};
        choice_change_index_in0{t}(n,:) = group_result(n).choice_change_index_in0{t};
        choice_change_index_around0{t}(n,:) = group_result(n).choice_change_index_around0{t};
        
        choice_change_around0{t}(n,:) = group_result(n).choice_change_around0{t};
        
        % fitting result
        for ind = 1:3
            psy_perf_ctrl{t,ind}(n,:) = group_result(n).psy_perf{t}{1,ind};  % 列：h r mt, 行:n=1电刺激前，n=2电刺激后
            psy_perf_stim{t,ind}(n,:) = group_result(n).psy_perf{t}{2,ind};  % 列：h r mt, 行:n=1电刺激前，n=2电刺激后
        end
        
        % general error level
        diffBias_true_task{t}(n,:) = group_result(n).diffBias_true{t}'; % 每一列h,r,mt，% after-before  >0: bias to left; <0: bias to right
        nor_diffBias_true_task{t}(n,:) = group_result(n).nor_diffBias_true{t}';
        diffThresh_true_task{t}(n,:) = group_result(n).diffThresh_true{t}';
        diffBias_p_true_task{t}(n,:) = group_result(n).p_mu_true{t}';
        diffThresh_p_true_task{t}(n,:) = group_result(n).p_theta_true{t}';
    end % end t
    
    if sum(~isnan(stimAmp(n,:)))>0
        unique_stimAmp_temp = unique(stimAmp(n,~isnan(stimAmp(n,:))));
        if length(unique_stimAmp_temp) == 1
            cell_unique_stimAmp(n,1) = unique_stimAmp_temp;
        else
            keyboard; %同一个cell用了多种电刺激强度
        end
    else
        cell_unique_stimAmp(n,1) = nan;
    end
end

% save origin bias shift : after-before  >0: bias to left; <0: bias to right  (与prefer方向无关，仅反应bias方向)
diffBias_task_origin = diffBias_task;
nor_diffBias_task_origin = nor_diffBias_task;
diffBias_task_true_origin = diffBias_true_task;
nor_diffBias_true_task_origin = nor_diffBias_true_task;

diffBias_mt_fine = diffBias_task{1}(:,3); % stimCueS
diffBias_mt_coarse = diffBias_task{4}(:,3); % stimCo

% positive: stim to pref direction
% pref right(code = 2): after - before < 0
% pref left(code = 1): after - before > 0
% Flip dPSE to positive if pref is right,
global_pref_code = cell2mat({group_result.global_pref_code}'); % come from global tuning, mt pref_direction come from spiral vector-sum
for t = 1:6
    diffBias_task{t}(global_pref_code == 2) = -diffBias_task{t}(global_pref_code == 2);
    nor_diffBias_task{t}(global_pref_code == 2) = -nor_diffBias_task{t}(global_pref_code == 2);
    diffBias_true_task{t}(global_pref_code == 2) = -diffBias_true_task{t}(global_pref_code == 2);
    nor_diffBias_true_task{t}(global_pref_code == 2) = -nor_diffBias_true_task{t}(global_pref_code == 2);
    
    % step analysis
    for ind = 1:3
        diffBias_task_step{t,ind}(global_pref_code(:,ind) == 2,:) = diffBias_task_step{t,ind}(global_pref_code(:,ind) == 2,:).* -1;
        nor_diffBias_task_step{t,ind}(global_pref_code(:,ind) == 2,:) = nor_diffBias_task_step{t,ind}(global_pref_code(:,ind) == 2,:).* -1;
    end
end

% replace mt bias shift not changed by preference
% shift < 0, bias to Translation
% shift > 0, bias to Rotation
diffBias_task_originalmt = diffBias_task;
% 只有{1} 和 {4} 有motion type task
diffBias_task_originalmt{1}(:,3) = diffBias_mt_fine;
diffBias_task_originalmt{4}(:,3) = diffBias_mt_coarse;


%% ***************************************************************************************
% *****************************************************************************************************

% 此时有3种diffBias：
% 1.diffBias_task_origin: after-before  >0: bias to left/CW; <0: bias to right/CCW
% 2.diffBias_task: positive: stim to pref direction
% 3.diffBias_task_originalmt: h and r:stim to pref direction; mt:after-before  >0: bias to Rotation; <0: bias to Heading

% combine fine and coarse task PSE shift
% normalize PSE shift: divided by the threshold under non-stim trials   (Xuefei Yu, 2018)
% ===============  Condition base  ====================
% stimCueS stimHonly stimRonly stimCo stimCoH stimCoR
%    1         2         3       4        5      6
% =====================================================

% *****************************************************************************************************
% *****************************************************************************************************



%% Cell Selection and Cell Counter
select_all = [];
select_typical = [];
select_no_typical = [];
loop_num = [];
selection_num = [];
t_criterion_txt = [];
global_pref_p = [];
monkey_included_for_analysis = [];
stimtask_type = cell(1,2);
cell_selection();

    function cell_selection(typical_cell_selection_num)  % Cell Selection and Cell
        if nargin < 1
            typical_cell_selection_num = 2; % Default,select T+R+Spi cell (all tuning cell)
        end
        
        area = {loadM_data.Area}';
        
        % properties of each cell, for selection
        global_pref_p = cell2mat({group_result.global_pref_p}'); % anova for all direction (h and r)
        
        RF = cell2mat({group_result.RF}');
        rf_size = RF(:,3).*RF(:,4);
        rf_dia = 2.* sqrt(rf_size / pi);
        rf_ecc = sqrt(RF(:,1).^2+RF(:,2).^2);
        
        rf_dia_limit = 45; % 区分大小RF
        rf_ecc_limit = 45/2;
        
        %         % 不同stim task的repetition不一样，没法直接限定rep。。。
        %         Stim_repN = cell2mat({group_result.Stim_repN}');
        
        % u-stim amplitude, 暂时一个细胞只有一个电刺激强度（不考虑一个细胞做了两种电刺激强度）
        % 暂时没想到从细胞层面选typical cell的方法，因为同一个细胞不同的protocol也可能有个不同的电刺激强度，故不再此筛选
        
        % loose standard: somehow have tuning (ranksum>spon & all conditino have tuning, anova<0.05)
        select_all_all_monkey = cellfun(@(x) strncmp(x,'MST',3),area) &  [group_result.p_stim]'<0.05 & [group_result.SpiT_repN]' >= 3 & global_pref_p(:,3)<0.05;
        
        % Extremely loose standards: 有反应就行，不一定有tuning
        % select_all_all_monkey = cellfun(@(x) strncmp(x,'MST',3),area) &  [group_result.p_stim]'<0.05 & [group_result.SpiT_repN]' >= 3;
        
        % restrictive standard： well-defined tuning cell
        % select_all_all_monkey = cellfun(@(x) strncmp(x,'MST',3),area) &  [group_result.p_stim]'<0.05 & [group_result.SpiT_repN]' >= 3 & (([group_result.cell_type] == 1)' | ([group_result.cell_type] == 2)'  | ([group_result.cell_type] == 3)' );
        
        typical_cell_selection_criteria = ...
            {%  Logic                                   Notes
            % Bottom-line
            'Bottom-line (all)', select_all_all_monkey; % Just all bottom-line cells
            
            % cell type base
            'T+R+S cell', select_all_all_monkey & (([group_result.cell_type] == 1)' | ([group_result.cell_type] == 2)'  | ([group_result.cell_type] == 3)' );
            'Translation Cell', select_all_all_monkey & ([group_result.cell_type] == 1)';
            'Rotation Cell', select_all_all_monkey & ([group_result.cell_type] == 2)';
            'Spiral Cell', select_all_all_monkey & ([group_result.cell_type] == 3)';
            
            % cell properties
            'Large RF cell', select_all_all_monkey & rf_dia > rf_dia_limit;
            'Small RF cell',select_all_all_monkey & rf_dia <= rf_dia_limit;
            'Center RF cell',select_all_all_monkey & rf_ecc <= rf_ecc_limit;
            'Peripheral RF cell',select_all_all_monkey & rf_ecc > rf_ecc_limit;
            'Small and Center RF cell',select_all_all_monkey & rf_ecc <= rf_ecc_limit & rf_dia <= rf_dia_limit;
            'Small and Peripheral RF cell',select_all_all_monkey & rf_ecc > rf_ecc_limit & rf_dia <= rf_dia_limit;
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
        
        % -------- Update task type selection. Lwh 202102 --------
        tartype_select = logical([get(findall(gcbf,'tag','two_target'),'value')  get(findall(gcbf,'tag','four_target'),'value')]);
        tasktype_select = logical([get(findall(gcbf,'tag','fine_task'),'value')  get(findall(gcbf,'tag','coarse_task'),'value')]);
        
        if ~tasktype_select(1) && ~tasktype_select(2)
            disp('********** Atleast choose one task type! **********');
            return
        elseif ~tartype_select(1) && ~tartype_select(2)
            disp('********** Atleast choose one target type! **********');
            return
        end
        
        if tartype_select(1) && ~tartype_select(2) % onbly 2-target
            stimtask_type{1} = [2 3]; % fine
            stimtask_type{2} = [5 6]; % coarse
        elseif ~tartype_select(1) && tartype_select(2) % onbly 4-target
            stimtask_type{1} = [1]; % fine
            stimtask_type{2} = [4]; % coarse
        elseif tartype_select(1) && tartype_select(2) % 2 and 4 target
            stimtask_type{1} = [1 2 3]; % fine
            stimtask_type{2} = [4 5 6]; % coarse
        end
        stimtask_type{1} = stimtask_type{1} .* tasktype_select(1); % {1} fine task
        stimtask_type{2} = stimtask_type{2} .* tasktype_select(2); % {2} coarse task
        
        for i = 1:2 % fine or coarse
            if sum(stimtask_type{i}==0)>0
                stimtask_type{i} = []; % remove 0
            end
        end
        
        %         stimtask_type{1}
        %         stimtask_type{2}
        %         disp('--------------------------');
    end


%% Final Preparation

% =========== Data for common use ============
function_handles_stim = [];

% ================ Miscellaneous ===================
set(0,'defaultAxesColorOrder',[41 89 204; 248 28 83; 14 153 46]/255);

colors = [41 89 204; 248 28 83; 14 153 46]/255;
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
figN = 300;

%% ====================================== Function Handles =============================================%
% Today I reconstruct the group analysis codes in a more logical and effcient way:
% 1. I divided each figure into small nested functions in Group_HD
% 2. Group_HD returns a cell called "function_handles" back to GROUP_GUI (with all the data available for each nested function)
% 3. Now I can choose which figure(s) to plot or which function to debug in GROUP_GUI
% 4. I feel so happy.
% @HH20150425;

function_handles_stim = {
    'Stim Effect in PSE shift',{
    'PSE shift (20uA + 40uA, DotStim)',@f1p1;
    'PSE shift (20uA or 40uA, DotStim)',@f1p2;
    'PSE shift (20uA + 40uA, TriangleStim)',@f1p3;
    'PSE shift (20uA or 40uA, TriangleStim)',@f1p4;
    'PSE shift (20uA + 40uA, Dot+Triangle)',@f1p5;
    'PSE shift (20uA or 40uA, Dot+Triangle)',@f1p6;
    'general error PSE shift (20+40uA)',@f1p9;
    '',@f1p101;
    '',@f1p101;
    'Stim effect and site location (Grid-xy)', @f1p7
    'PSE shift between different stim current or targets num', @f1p8
    }
    'Stim Effect in Threshold changed',{
    'Threshold changed (20uA + 40uA, DotStim)',@f2p1;
    'Threshold changed (20uA or 40uA, DotStim)',@f2p2;
    'Threshold changed (20uA + 40uA, TriangleStim)',@f2p3;
    'Threshold changed (20uA or 40uA, TriangleStim)',@f2p4;
    'Threshold changed (20uA + 40uA, Dot+Triangle)',@f2p5;
    'Threshold changed (20uA or 40uA, Dot+Triangle)',@f2p6;
    '',@f2p101;
    '',@f2p101;
    'Threshold changed between different stim current or targets num', @f2p7
    'Rotation Thr. changed between different targets num in same site', @f2p8
    'Stim induced CR change',@f2p9;
    'CR change between different stim current',@f2p10;
    'CR change between different target num',@f2p11;
    }
    
    'u-stim effect and some modulation',{
    'PSE shift affected by Modulation Index',@f3p1;
    %     'Motion type PSE shift and modulation index',@f3p2;  % move to Motion type list
    'Translation V.S. Rotation (PSE shift)',@f3p2;
    'Fine task V.S. Coarse task (PSE shift)',@f3p3;
    'Fine task V.S. Coarse task (Threshold change)',@f3p4;
    'PSE shift affect by repetition',@f3p5;
    'Time course of microstimulation effects',@f3p9;
    'PSE shift affect by maxFR',@f3p6;
    'PSE shift affect by cluster',@f3p7;
    'PSE shift affect by RF',@f3p8;
    'Regression for PSE shift',@f3p10;
    'PSE shit and Preferred direction',@f3p11;
    "PSE shift and d'", @f3p12;
    'PSE shift and some T-R ratio', @f3p13;
    'PSE shift and Preferred OpticFlow',@f3p14;
    'PSE shift and Preferred OpticFlow Fine vs Coarse',@f3p17;
    'PSE shift for fovea-containing RF sites',@f3p15;
    'PSE shift and CR',@f3p16;
    'PSE shift in normal psy and general error level',@f3p18;

    %         'PSE shift and CP'
    
    %     'RF effect in u-stim'
    }
    'Motion type',{
    'Motion type PSE shift',@f4p10;
    'Motion type PSE shift with Modulation index',@f4p1;
    'Motion type PSE shift with Spiral index',@f4p12;
    'Motion type PSE shift with direction preference',@f4p13;
    'Motion type PSE shift with T/R threshold change',@f4p17;
    'Choice changed across spiral index',@f4p7;
    'Choice changed across cell type',@f4p16;
    'Choice changed across direction preference',@f4p14;
    'Choice changed between control and stim and chi-square test',@f4p8;
    'Choice changed with some T-R ratio',@f4p15;
    'near 0° choice changed',@f4p9
    'Four choice changed in stimulation',@f4p5;
    'Translation and Rotatin choice changed in stimulation',@f4p6;
    'Around zero choice and d prime',@f4p2;
    'Four Type cell in Four type stim.',@f4p3; % pref R L CCW CW 四种细胞分别在R L CCW CW刺激下，电刺激的影响
    'a new Motion type psy curve (only preferred half-axis for T/R)',@f4p4;
    'show example site', @f4p18;
    '暂时不用的code',@f4p99;
    
    }
    
    'Tuning properties in StimCueS task',{
    'Preferred direction in stim task',@f5p1;
    'RF',@f5p2;
    'RF type and stimulation effect',@f5p3;
    'Some tuning index in 4-target task',@f5p4;
    'Cell type in stim task',@f5p5;
    }
    
    'Other',{
    'Psychometric function in Two-target task',@f6p1;
    'Psychometric function in Four-target task',@f6p2;
    'Comparison of 2T and 4T task performance',@f6p6;
    'Microstimulation site location',@f6p3;
    'Motion type behavioral performance',@f6p4;
    'CR of Translation and Rotation task',@f6p5;
    'Eyetrace during VisualOn period',@f6p7;
    'Flow-pattern μ distribution in 4AFC',@f6p8;
    }
    
    'NoShow',{@cell_selection};
    };



% ---------------------- get data ----------------------------
% ------------------ default parameter -----------------------
% ------------------------ Plot ------------------------------

%% ====================================== Function Definitions =============================================%
% ===============  Condition base  ====================
% stimCueS stimHonly stimRonly stimCo stimCoH stimCoR
%    1         2         3       4        5      6
% =====================================================

% stim_type:task类型，见上方
% effect_type : 1: PSE, 2: threshold

%% PSE shift
% % is_original_bias: 0: bias positive = prefer direction; 1: only original bias in MT;  2: original bias(after-before); 3： normalize PSE; 4 normalize original bias

    function f1p1(debug) % 20uA + 40uA, DotStim
        if debug  ; dbstack;   keyboard;      end
        effect_type = 1; % 1: PSE, 2: threshold
        is_original_bias = 0; % 0: bias positive = prefer direction; 1: only original bias in MT;  2: original bias(after-before);
        is_separate_stimAmp = 0; % 默认为0，combine 20 and 40 stim amplitude
        select_stars_type = [0]; % 0: dot stimulus; 1: triangle stimulus
        methods_of_select = {select_typical, 'Typical cells'};
        [diff_thistask, S_unit, NS_unit,unique_stimAmp_this] = get_stim_data_uA(stimtask_type,methods_of_select{1,1},is_separate_stimAmp,is_original_bias,select_stars_type);
        plot_stim_effect_distribution_uA(stimtask_type,monkey_included_for_analysis,effect_type,diff_thistask,S_unit,NS_unit,unique_stimAmp_this,is_separate_stimAmp,is_original_bias);
    end

    function f1p2(debug) % 20uA or 40uA, DotStim
        if debug  ; dbstack;   keyboard;      end
        effect_type = 1; % 1: PSE, 2: threshold
        is_original_bias = 0; % 0: bias positive = prefer direction; 1: only original bias in MT;  2: original bias(after-before);
        is_separate_stimAmp = 1; % 默认为0，combine 20 and 40 stim amplitude
        select_stars_type = [0]; % 0: dot stimulus; 1: triangle stimulus
        methods_of_select = {select_typical, 'Typical cells'};
        [diff_thistask, S_unit, NS_unit,unique_stimAmp_this] = get_stim_data_uA(stimtask_type,methods_of_select{1,1},is_separate_stimAmp,is_original_bias,select_stars_type);
        plot_stim_effect_distribution_uA(stimtask_type,monkey_included_for_analysis,effect_type,diff_thistask,S_unit,NS_unit,unique_stimAmp_this,is_separate_stimAmp,is_original_bias);
    end

    function f1p3(debug) % 20uA + 40uA, TriangleStim
        if debug  ; dbstack;   keyboard;      end
        effect_type = 1; % 1: PSE, 2: threshold
        is_original_bias = 0; % 0: bias positive = prefer direction; 1: only original bias in MT;  2: original bias(after-before);
        is_separate_stimAmp = 0; % 默认为0，combine 20 and 40 stim amplitude
        select_stars_type = [1]; % 0: dot stimulus; 1: triangle stimulus
        methods_of_select = {select_typical, 'Typical cells'};
        [diff_thistask, S_unit, NS_unit,unique_stimAmp_this] = get_stim_data_uA(stimtask_type,methods_of_select{1,1},is_separate_stimAmp,is_original_bias,select_stars_type);
        plot_stim_effect_distribution_uA(stimtask_type,monkey_included_for_analysis,effect_type,diff_thistask,S_unit,NS_unit,unique_stimAmp_this,is_separate_stimAmp,is_original_bias);
    end

    function f1p4(debug) % 20uA or 40uA, TriangleStim)
        if debug  ; dbstack;   keyboard;      end
        effect_type = 1; % 1: PSE, 2: threshold
        is_original_bias = 0; % 0: bias positive = prefer direction; 1: only original bias in MT;  2: original bias(after-before);
        is_separate_stimAmp = 1; % 默认为0，combine 20 and 40 stim amplitude
        select_stars_type = [1]; % 0: dot stimulus; 1: triangle stimulus
        methods_of_select = {select_typical, 'Typical cells'};
        [diff_thistask, S_unit, NS_unit,unique_stimAmp_this] = get_stim_data_uA(stimtask_type,methods_of_select{1,1},is_separate_stimAmp,is_original_bias,select_stars_type);
        plot_stim_effect_distribution_uA(stimtask_type,monkey_included_for_analysis,effect_type,diff_thistask,S_unit,NS_unit,unique_stimAmp_this,is_separate_stimAmp,is_original_bias);
    end

    function f1p5(debug) % 20uA + 40uA, Dot+Triangle
        if debug  ; dbstack;   keyboard;      end
        effect_type = 1; % 1: PSE, 2: threshold
        
        if isempty(stimtask_type{1})+isempty(stimtask_type{2}) == 1 % only fine or only coarse task
            is_original_bias = 0; % 0: bias positive = prefer direction; 1: only original bias in MT;  2: original bias(after-before); 3： normalize PSE
        else % combine fine and coarse, use normalized PSE shift
            is_original_bias = 3;
        end
        
        is_separate_stimAmp = 0; % 默认为0，combine 20 and 40 stim amplitude
        select_stars_type = [0 1]; % 0: dot stimulus; 1: triangle stimulus
        methods_of_select = {select_typical, 'Typical cells'};
        [diff_thistask, S_unit, NS_unit,unique_stimAmp_this] = get_stim_data_uA(stimtask_type,methods_of_select{1,1},is_separate_stimAmp,is_original_bias,select_stars_type);
        plot_stim_effect_distribution_uA(stimtask_type,monkey_included_for_analysis,effect_type,diff_thistask,S_unit,NS_unit,unique_stimAmp_this,is_separate_stimAmp,is_original_bias);
            
        % motion type 正方向：spiral_indexd定义的preference决定
        % 统计显著的cases数目，显著且positive的case 的数目
        figure(figN-1);
        ha = get(gcf,'children');
        subplot_this = [];
        % get subplot according title
        for i = 1:length(ha)
            if strcmp(ha(i).Title.String,'Translation')
                subplot_this(i) = 1;
            elseif strcmp(ha(i).Title.String,'Roll')
                subplot_this(i) = 2;
            elseif strcmp(ha(i).Title.String,'Flow pattern')
                subplot_this(i) = 3;
            end
        end
        
        hold on
        for hr = 1:3
            sig_num(hr) = sum(S_unit{effect_type}(:,hr));
            nosig_num(hr) = sum(NS_unit{effect_type}(:,hr));
            sig_pos_num(hr) = sum(S_unit{effect_type}(:,hr) & diff_thistask{effect_type}(:,hr)>0);
            
            ylimit = ha(subplot_this==hr).YLim;
            xlimit = ha(subplot_this==hr).XLim;
            str = [];
            str = [sprintf('Sig num = %g, Sig&Pos num = %g',sig_num(hr),sig_pos_num(hr)) newline ...
                sprintf('Sig = %.4g%%, Sig&Pos = %.4g%%',sig_num(hr)/(sig_num(hr)+nosig_num(hr))*100,sig_pos_num(hr)/(sig_num(hr))*100)];
            text(ha(subplot_this==hr),xlimit(1),ylimit(2)/2,str);
        end
    end

    function f1p9(debug) % 'general error PSE shift (20+40uA)'
        if debug  ; dbstack;   keyboard;      end
        effect_type = 1; % 1: PSE, 2: threshold
        
        if isempty(stimtask_type{1})+isempty(stimtask_type{2}) == 1 % only fine or only coarse task
            is_original_bias = 5; % 0: bias positive = prefer direction; 1: only original bias in MT;  2: original bias(after-before); 3： normalize PSE
        else % combine fine and coarse, use normalized PSE shift
            is_original_bias = 6;
        end
        
        is_separate_stimAmp = 0; % 默认为0，combine 20 and 40 stim amplitude
        select_stars_type = [0 1]; % 0: dot stimulus; 1: triangle stimulus
        methods_of_select = {select_typical, 'Typical cells'};
        [diff_thistask, S_unit, NS_unit,unique_stimAmp_this] = get_stim_data_uA(stimtask_type,methods_of_select{1,1},is_separate_stimAmp,is_original_bias,select_stars_type);
        plot_stim_effect_distribution_uA(stimtask_type,monkey_included_for_analysis,effect_type,diff_thistask,S_unit,NS_unit,unique_stimAmp_this,is_separate_stimAmp,is_original_bias);
            
        % motion type 正方向：spiral_indexd定义的preference决定
        % 统计显著的cases数目，显著且positive的case 的数目
        figure(figN-1);
        ha = get(gcf,'children');
        subplot_this = [];
        % get subplot according title
        for i = 1:length(ha)
            if strcmp(ha(i).Title.String,'Translation')
                subplot_this(i) = 1;
            elseif strcmp(ha(i).Title.String,'Roll')
                subplot_this(i) = 2;
            elseif strcmp(ha(i).Title.String,'Flow pattern')
                subplot_this(i) = 3;
            end
        end
        
        hold on
        for hr = 1:3
            sig_num(hr) = sum(S_unit{effect_type}(:,hr));
            nosig_num(hr) = sum(NS_unit{effect_type}(:,hr));
            sig_pos_num(hr) = sum(S_unit{effect_type}(:,hr) & diff_thistask{effect_type}(:,hr)>0);
            
            ylimit = ha(subplot_this==hr).YLim;
            xlimit = ha(subplot_this==hr).XLim;
            str = [];
            str = [sprintf('Sig num = %g, Sig&Pos num = %g',sig_num(hr),sig_pos_num(hr)) newline ...
                sprintf('Sig = %.4g%%, Sig&Pos = %.4g%%',sig_num(hr)/(sig_num(hr)+nosig_num(hr))*100,sig_pos_num(hr)/(sig_num(hr))*100)];
            text(ha(subplot_this==hr),xlimit(1),ylimit(2)/2,str);
        end
    end

    function f1p6(debug) % 20uA or 40uA, Dot+Triangle
        if debug  ; dbstack;   keyboard;      end
        effect_type = 1; % 1: PSE, 2: threshold
        is_original_bias = 0; % 0: bias positive = prefer direction; 1: only original bias in MT;  2: original bias(after-before);
        is_separate_stimAmp = 1; % 默认为0，combine 20 and 40 stim amplitude
        select_stars_type = [0 1]; % 0: dot stimulus; 1: triangle stimulus
        methods_of_select = {select_typical, 'Typical cells'};
        [diff_thistask, S_unit, NS_unit,unique_stimAmp_this] = get_stim_data_uA(stimtask_type,methods_of_select{1,1},is_separate_stimAmp,is_original_bias,select_stars_type);
        plot_stim_effect_distribution_uA(stimtask_type,monkey_included_for_analysis,effect_type,diff_thistask,S_unit,NS_unit,unique_stimAmp_this,is_separate_stimAmp,is_original_bias);
    end

    function f1p7(debug) % Stim effect (PSE shift significant v.s. not significant) and site location (Grid-xy)
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        effect_type = 1; % 1: PSE, 2: threshold
        select_stars_type = [0 1]; % 0: dot stimulus; 1: triangle stimulus
        is_separate_stimAmp = 0; % 默认为0，combine 20 and 40 stim amplitude
        is_original_bias = 0;  % 1: only original bias in MT;  2: original bias(after-before);   3: bias positive = prefer direction
        
        % ---------------------- get data ----------------------------
        % PSE shift significance
        [diff_thistask_temp, S_unit_temp, NS_unit_temp] = get_stim_data_uA(stimtask_type,methods_of_select{1,1},is_separate_stimAmp,is_original_bias,select_stars_type);
        diff_thistask = diff_thistask_temp{1,effect_type}; % diff_thistask_temp{amp,effect_type}
        S_unit = S_unit_temp{1,effect_type};
        NS_unit = NS_unit_temp{1,effect_type};
        
        % GridX, GridY
        GridX_temp = cell2mat({group_result.GridX}');
        GridY_temp = cell2mat({group_result.GridY}');
        
        % repmat GridX, GridY, methods_of_select to same size with diff_thistask
        task_type_num = sum(cellfun(@(x) length(x),stimtask_type));
        
        GridX = repmat(GridX_temp,task_type_num,1);
        GridY = repmat(GridY_temp,task_type_num,1);
        methods_of_select_this = repmat(methods_of_select{1,1},task_type_num,1);
        
        uniqueX = unique(GridX(methods_of_select_this));
        uniqueY = unique(GridY(methods_of_select_this));
        
        % always separate monkey
        if length(monkey_included_for_analysis)>1
            disp('********** Please choose ONE monkey! **********');
            return
        end
        
        % ------------------------ Plot ------------------------------
        %         % proportion
        %         set(figure(figN),'name','Stim effect and x-y position (Proportion)','pos',[150,50,600,800]); clf
        %         ha = tight_subplot(length(uniqueY), length(uniqueX),[.03 .03], [.05 .05], [.05 .05]);
        %         for xx = 1:length(uniqueX)
        %             for yy = 1:length(uniqueY)
        %                 axes(ha(xx+(yy-1)*length(uniqueX)))
        %                 % do not show the axis and label
        %                 ha(xx+(yy-1)*length(uniqueX)).YAxis.Visible = 'off';
        %                 ha(xx+(yy-1)*length(uniqueX)).XAxis.Visible = 'off';
        %
        %                 % show the x and y number
        %                 if xx == 1 && yy~=1
        %                     text(0,0.5,num2str(uniqueY(length(uniqueY)+1-yy)));
        %                 end
        %                 if yy == length(uniqueY)
        %                     text(1.2,-0.3,num2str(uniqueX(xx)));
        %                 end
        %
        %                 % x-y note
        %                 if xx == 1 && yy == 1
        %                     plot([1 2],[0.5 0.5],'k-');
        %                     hold on
        %                     plot([1.5 1.5],[.25 .75],'k-');
        %                     text(1,0.1,'Anterior');
        %                     text(1,0.9,'Posterior');
        %                     text(0,0.5,'Middle');
        %                     text(2,0.5,'Lateral');
        %                     ha(1).YAxis.Visible = 'off';
        %                     ha(1).XAxis.Visible = 'off';
        %                 end
        %
        %                 select = logical(GridX == uniqueX(xx) & GridY == uniqueY(length(uniqueY)+1-yy)); % Y 轴从大到小排列，原点(0,0)位于左下角
        %                 S_unit_num = sum(S_unit(select,:));
        %                 NS_unit_num = sum(NS_unit(select,:));
        %
        %                 if sum(S_unit_num + NS_unit_num)>0
        %                     for mt = 1:2 % h or r
        %                         S_unit_pref = logical(select & diff_thistask(:,mt)>0 & S_unit(:,mt)); % significant bias to prefer direction
        %                         S_unit_pref_num = sum(S_unit_pref);
        %                         sig_pro(mt) = S_unit_num(mt) / (S_unit_num(mt) + NS_unit_num(mt));
        %                         sig_pref_pro(mt) = S_unit_pref_num / (S_unit_num(mt) + NS_unit_num(mt));
        %                         hold on
        %                         bar(mt,1,'facecolor',[.69 .69 .69]); % all
        %                         bar(mt,sig_pro(mt),'facecolor','k'); % significant but not prefer direction
        %                         bar(mt,sig_pref_pro(mt),'facecolor',c{mt}); % significant and prefer direction
        %                     end
        %                     set(gca,'ytick', []) ;
        %                     set(gca,'xtick', []) ;
        %                     text(1.2,-0.07,sprintf('N≈%0.3g',S_unit_num(1)+NS_unit_num(1)),'fontsize',8.5); % translation 和 rotation 数量可能不同，但是差别应该不大
        %                 end
        %
        %                 xlim([0.5 2.5]);
        %                 ylim([0 1]);
        %                 box off
        %
        %                 if xx == length(uniqueX) && yy == length(uniqueY)
        %                     hold on
        %                     bar(1,1,'facecolor',[.69 .69 .69]); % all
        %                     bar(1,2/3,'facecolor','k'); % significant but not prefer direction
        %                     bar(1,1/3,'facecolor',c{1}); % significant and prefer direction
        %                     ylim([0 1000]);
        %                     legend('not significant proportion','sig bias to non-prefer proportion','sig bias to prefer proportion','Location','northoutside')
        %                 end
        %             end
        %         end
        %         figN = figN+1;
        
        
        % ------------------------ Plot ------------------------------
        % exact num
        set(figure(figN),'name','Stim effect and x-y position (Number)','pos',[800,50,600,800]); clf
        ha = tight_subplot(length(uniqueY), length(uniqueX),[.03 .03], [.05 .05], [.05 .05]);
        max_num = 0;
        
        % set ylimit
        for xx = 1:length(uniqueX)
            for yy = 1:length(uniqueY)
                select = logical(GridX == uniqueX(xx) & GridY == uniqueY(length(uniqueY)+1-yy)); % Y 轴从大到小排列，原点(0,0)位于左下角
                S_unit_num = sum(S_unit(select,:));
                NS_unit_num = sum(NS_unit(select,:));
                all_num = S_unit_num+NS_unit_num;
                max_num = max(max(all_num),max_num);
            end
        end
        
        % plot
        for xx = 1:length(uniqueX)
            for yy = 1:length(uniqueY)
                axes(ha(xx+(yy-1)*length(uniqueX)));
                % do not show the axis and label
                %                 ha(xx+(yy-1)*length(uniqueX)).YAxis.Visible = 'off';
                %                 ha(xx+(yy-1)*length(uniqueX)).XAxis.Visible = 'off';
                
                select = logical(GridX == uniqueX(xx) & GridY == uniqueY(length(uniqueY)+1-yy)); % Y 轴从大到小排列，原点(0,0)位于左下角
                S_unit_num = sum(S_unit(select,:));
                NS_unit_num = sum(NS_unit(select,:));
                all_num = S_unit_num+NS_unit_num;
                
                if sum(S_unit_num + NS_unit_num)>0
                    for mt = 1:2 % h or r
                        S_unit_pref = logical(select & diff_thistask(:,mt)>0 & S_unit(:,mt)); % significant bias to prefer direction
                        S_unit_pref_num = sum(S_unit_pref);
                        hold on
                        bb = bar(mt,all_num(mt),'facecolor',[.69 .69 .69]); % significant + non-sig
                        bar(mt,S_unit_num(mt),'facecolor','k'); % significant
                        bar(mt,S_unit_pref_num,'facecolor',c{mt}); % significant and prefer direction
                        
                        text(mt,all_num(mt)+max_num/10,num2str(all_num(mt)),'HorizontalAlignment','center'); % show all num
                    end
                    set(gca,'ytick', []) ;
                    set(gca,'xtick', []) ;
                end
                
                xlim([0.5 2.5]);
                ylim([0 max_num])
                
                % show the x and y number
                if xx == 1
                    text(0.4,max_num/2,num2str(uniqueY(length(uniqueY)+1-yy)),'HorizontalAlignment','right'); % ylabel
                end
                if yy == length(uniqueY)
                    text(1.5,-max_num/12,num2str(uniqueX(xx)),'HorizontalAlignment','center'); % xlabel
                end
                box off
            end
        end
        set(findall(gcf,'type','axes'),'ylim',[0 max_num]);
        figN = figN+1;
    end

    function f1p8(debug) % PSE shift between different stim current or targets num
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        effect_type = 1; % 1: PSE, 2: threshold
        select_stars_type = [0 1]; % 0: dot stimulus; 1: triangle stimulus

        % ---------------------- get data ----------------------------
        stimtask_type_this{1,1} = [2,3]; stimtask_type_this{1,2} = []; % fine and 2-target
        stimtask_type_this{2,1} = []; stimtask_type_this{2,2} = [5,6]; % coarse and 2-target
        stimtask_type_this{3,1} = [1]; stimtask_type_this{3,2} = []; % fine and 4-target
        stimtask_type_this{4,1} = []; stimtask_type_this{4,2} = [4]; % coarse and 4-target
        
        % data for 20 uA and 40 uA
        is_separate_stimAmp = 1; % 默认为0，combine 20 and 40 stim amplitude
%         is_original_bias = 3;  % normlize PSE shift (divided by the threshold under non-stim trials
        is_original_bias = 0;
        for i = 1:4 % fine2T, coarse2T, fine4T, coarse4T
            diff_thistask = [];
            [diff_thistask] = get_stim_data_uA({stimtask_type_this{i,:}},methods_of_select{1,1},is_separate_stimAmp,is_original_bias,select_stars_type);
            PSE_20uA{i,1} = diff_thistask{1,effect_type}(~isnan(diff_thistask{1,effect_type}(:,1)),1); % translation % remove NaN
            PSE_20uA{i,2} = diff_thistask{1,effect_type}(~isnan(diff_thistask{1,effect_type}(:,2)),2); % rotation
            
            PSE_40uA{i,1} = diff_thistask{2,effect_type}(~isnan(diff_thistask{2,effect_type}(:,1)),1); % translation
            PSE_40uA{i,2} = diff_thistask{2,effect_type}(~isnan(diff_thistask{2,effect_type}(:,2)),2); % rotation
        end
        
        % combine 2 and 4 task: PSE_20uA_all{fine/coarse, T/R}
        PSE_20uA_all{1,1} = [PSE_20uA{1,1}; PSE_20uA{3,1}]; % fine translation
        PSE_20uA_all{1,2} = [PSE_20uA{1,2}; PSE_20uA{3,2}]; % fine rotation
        PSE_20uA_all{2,1} = [PSE_20uA{2,1}; PSE_20uA{4,1}]; % coarse translation
        PSE_20uA_all{2,2} = [PSE_20uA{2,2}; PSE_20uA{4,2}]; % coarse rotation
        
        % combine 2 and 4 task
        PSE_40uA_all{1,1} = [PSE_40uA{1,1}; PSE_40uA{3,1}]; % fine translation
        PSE_40uA_all{1,2} = [PSE_40uA{1,2}; PSE_40uA{3,2}]; % fine rotation
        PSE_40uA_all{2,1} = [PSE_40uA{2,1}; PSE_40uA{4,1}]; % coarse translation
        PSE_40uA_all{2,2} = [PSE_40uA{2,2}; PSE_40uA{4,2}]; % coarse rotation
        
        
        % data for 2T and 4T, 
        is_separate_stimAmp = 0;
        is_original_bias = 0; % 正值为bias to prefer direction
        data_2T = []; data_4T = [];
        for fc = 1:2 % fine or coarse   % fine2T, coarse2T, fine4T, coarse4T
            diff_thistask1 = [];
            [diff_thistask1] = get_stim_data_uA({stimtask_type_this{fc,:}},methods_of_select{1,1},is_separate_stimAmp,is_original_bias,select_stars_type);
            data_2T{fc} = diff_thistask1{1,effect_type}(:,1:2);
            
            diff_thistask2 = [];
            [diff_thistask2] = get_stim_data_uA({stimtask_type_this{fc+2,:}},methods_of_select{1,1},is_separate_stimAmp,is_original_bias,select_stars_type);
            data_4T{fc} = diff_thistask2{1,effect_type}(:,1:2);
        end

        % ------------------------ Plot ------------------------------
        % 1. stim current comparison
        set(figure(figN),'name','stim current comparison','pos',[100,50,600,800]); clf
        for i = 1:4 % fine2T, coarse2T, fine4T, coarse4T
            for hr = 1:2
                subplot(4,2,hr+(i-1)*2)
                hold on

                [binedge, ~, xtick_this, ~] = properHistBin([PSE_20uA{i,hr};PSE_40uA{i,hr}],'small');
                
                h1 = histogram(PSE_20uA{i,hr},binedge,'facecolor',c{hr});
                h2 = histogram(PSE_40uA{i,hr},binedge,'facecolor',c{hr+10});

                h1.Normalization = 'probability';
                h2.Normalization = 'probability';
                
                xticks(xtick_this);
                xlim([xtick_this(1) xtick_this(end)])
                
                ymax = max([h1.Values h2.Values]);
                [ylimit_max, ytick_this] = properYtick(ymax);
                yticks(ytick_this)
                ylim([0 ylimit_max])
                
                % mean
                plot(mean(PSE_20uA{i,hr}),ylimit_max,'v','color',c{hr});
                plot(mean(PSE_40uA{i,hr}),ylimit_max,'v','color',c{hr+10});
                
                % non-paired ttest
                [~,p_ttest(i,hr)] = ttest2(PSE_20uA{i,hr},PSE_40uA{i,hr}); % nan 为该task没有同时做过两种current stim
                text(xtick_this(1),ylimit_max-ylimit_max/15,sprintf('t-test, p=%0.3g',p_ttest(i,hr)));
                str = [sprintf('n20 = %.3g',length(PSE_20uA{i,hr})) newline sprintf('n40 = %.3g',length(PSE_40uA{i,hr}))];
                text(xtick_this(1),ylimit_max/2,str);
                
                if i == 1 && hr == 1
                    ylabel('Cases/N');
                    legend('20uA','40uA');
                end
                xlabel('normalized PSE shift');
                
                if hr == 1
                    if i == 1
                        title('Fine 2-Tar');
                    elseif i == 2
                        title('Coarse 2-Tar');
                    elseif i ==3
                        title('Fine 4-Tar');
                    else
                        title('Coarse 4-Tar');
                    end
                end
            end
        end
        SetFigure(13)
        suptitle('PSE shift between different stim current');
        figN = figN+1;

        % combine 2 and 4 task
        set(figure(figN),'name','stim current comparison','pos',[200,100,800,800]); clf
        for fc = 1:2
            for hr = 1:2
                
                subplot(2,2,hr+(fc-1)*2)
                hold on
                all_20 = PSE_20uA_all{fc,hr};
                all_40 = PSE_40uA_all{fc,hr};
                
                xmax = ceil(max(abs([all_20; all_40]))/5)*5; % 取5的倍数
                xbins = linspace(-xmax,xmax,13);
                
                % remove nan
                all_20 = all_20(~isnan(all_20));
                all_40 = all_40(~isnan(all_40));
                
                h1 = histogram(all_20,xbins,'facecolor',c{hr});
                h2 = histogram(all_40,xbins,'facecolor',c{hr+10});
                h1.Normalization = 'probability';
                h2.Normalization = 'probability';
                
                ymax = max([h1.Values h2.Values]);
                ylim([0 ymax+ymax/2]);
                xlim([min(xbins)-(xbins(2)-xbins(1)) max(xbins)+(xbins(2)-xbins(1))]);
                
                % median
                plot(median(all_20),ymax+ymax/4,'v','color',c{hr});
                plot(median(all_40),ymax+ymax/4,'v','color',c{hr+10});
                
                % sign test with 0
                p_sign_20 = signtest(all_20,0);
                p_sign_40 = signtest(all_40,0);
                str = [sprintf('20uA: p = %.2g (sign test)',p_sign_20) newline sprintf('40uA: p = %.2g',p_sign_40) newline ...
                    sprintf('n_2_0 = %g', length(all_20)) newline sprintf('n_4_0 = %g', length(all_40))];
                text(xbins(1),ymax,str);
                
                % non-paired t test
                p_ttest = [];
                [~,p_ttest] = ttest2(all_20,all_40); % nan 为该task没有同时做过两种current stim
                text(xbins(1),ymax+ymax/2.5,sprintf('t test, p=%2.2g',p_ttest));
                
                % KW test for unequal sample size, nonparametic
                if ~isempty(all_20) && ~isempty(all_40)
                    X = [all_20 ones(size(all_20)); all_40 2.*ones(size(all_40))];
                    [p] = KWtest(X,0.05);
                    text(xbins(1),ymax+ymax/3,sprintf('KW test, p=%2.2g',p));
                else
                    text(xbins(1),ymax+ymax/3,'KW test, p=NaN');
                end
                
                axis square
                ylabel('Cases/N');
                xlabel('PSE shift');
                
                if fc == 1 && hr == 1
                    legend('20uA','40uA');
                end
                
                if fc == 1
                    title('fine task');
                else
                    title('coarse task');
                end
                
                set(gca,'xticklabelmode','manual');
            end
        end
        
        SetFigure(10);
        [~,~,~,monkey] = SetTitle(stimtask_type,monkey_included_for_analysis);
        str = ['PSE shift (2AFC+4AFC) between different stim current ', monkey];
        suptitle(str);
        figN = figN + 1;
        
        % 2. task type 
        % (only select 40uA stim current)
        % Update, combine 20uA and 40uA
        set(figure(figN),'name','target num comparison','pos',[1000,50,800,800]); clf
        for fc = 1:2 % fine or coarse   % fine2T, coarse2T, fine4T, coarse4T
            for hr = 1:2 % T or R
                subplot(2,2,hr+2*(fc-1))
                hold on
                
                data_2T_this = data_2T{fc}(:,hr);
                data_4T_this = data_4T{fc}(:,hr);
                data_2T_this(isnan(data_2T_this)) = [];
                data_4T_this(isnan(data_4T_this)) = [];

                xbins = linspace(-max(abs([data_2T_this;data_4T_this])),max(abs([data_2T_this;data_4T_this])),11);
%                 h1 = histogram(h1_data,xbins,'facecolor',c{hr}); % 2T
%                 h2 = histogram(h2_data,xbins,'facecolor',c{hr+10});  % 4T
                
                h1 = histogram(data_2T_this,xbins,'edgecolor',c{hr+10},'DisplayStyle','stairs','linewidth',2); % 2T
                h2 = histogram(data_4T_this,xbins,'edgecolor',c{hr},'DisplayStyle','stairs','linewidth',2);  % 4T
                
                h1.Normalization = 'probability';
                h2.Normalization = 'probability';
                
                ymax = max([h1.Values h2.Values]);
                ymax = 0.6;
                ylim([0 ymax+ymax/2]);
                xlim([min(xbins)-(xbins(2)-xbins(1)) max(xbins)+(xbins(2)-xbins(1))]);
                xlabel('PSE shift');
                
                % mean
                plot(nanmean(data_2T_this),ymax+ymax/4,'v','color',c{hr+10});
                plot(nanmean(data_4T_this),ymax+ymax/4,'v','color',c{hr});
                
                % non-paired ttest
                [~,p_ttest2(fc,hr)] = ttest2(data_2T_this,data_4T_this);
                
                % anotation
                str = [sprintf('t-test, p=%0.3g',p_ttest2(fc,hr)) newline ...
                    sprintf('2-T N = %0.3g',length(data_2T_this)) newline sprintf('4-T N = %0.3g',length(data_4T_this))];
                text(xbins(1),ymax+ymax/3,str);
                
                legend('2-Target','4-Target');
                if fc == 1 && hr == 1
                    ylabel('Cases/N');
                end
                
                if fc == 1
                    title('Fine task');
                else
                    title('Coarse task');
                end
            end
        end
        SetFigure(13)
        suptitle('PSE shift between different target number');
        figN = figN+1;
    end


%% Threshold changedd
    function f2p1(debug) % 20uA + 40uA, DotStim
        if debug  ; dbstack;   keyboard;      end
        effect_type = 2; % 1: PSE, 2: threshold
        is_original_bias = 0; % 0: bias positive = prefer direction; 1: only original bias in MT;  2: original bias(after-before);
        is_separate_stimAmp = 0; % 默认为0，combine 20 and 40 stim amplitude
        select_stars_type = [0]; % 0: dot stimulus; 1: triangle stimulus
        methods_of_select = {select_typical, 'Typical cells'};
        [diff_thistask, S_unit, NS_unit,unique_stimAmp_this] = get_stim_data_uA(stimtask_type,methods_of_select{1,1},is_separate_stimAmp,is_original_bias,select_stars_type);
        plot_stim_effect_distribution_uA(stimtask_type,monkey_included_for_analysis,effect_type,diff_thistask,S_unit,NS_unit,unique_stimAmp_this,is_separate_stimAmp,is_original_bias);
    end

    function f2p2(debug) % 20uA or 40uA, DotStim
        if debug  ; dbstack;   keyboard;      end
        effect_type = 2; % 1: PSE, 2: threshold
        is_original_bias = 0; % 0: bias positive = prefer direction; 1: only original bias in MT;  2: original bias(after-before);
        is_separate_stimAmp = 1; % 默认为0，combine 20 and 40 stim amplitude
        select_stars_type = [0]; % 0: dot stimulus; 1: triangle stimulus
        methods_of_select = {select_typical, 'Typical cells'};
        [diff_thistask, S_unit, NS_unit,unique_stimAmp_this] = get_stim_data_uA(stimtask_type,methods_of_select{1,1},is_separate_stimAmp,is_original_bias,select_stars_type);
        plot_stim_effect_distribution_uA(stimtask_type,monkey_included_for_analysis,effect_type,diff_thistask,S_unit,NS_unit,unique_stimAmp_this,is_separate_stimAmp,is_original_bias);
    end

    function f2p3(debug) % 20uA + 40uA, TriangleStim
        if debug  ; dbstack;   keyboard;      end
        effect_type = 2; % 1: PSE, 2: threshold
        is_original_bias = 0; % 0: bias positive = prefer direction; 1: only original bias in MT;  2: original bias(after-before);
        is_separate_stimAmp = 0; % 默认为0，combine 20 and 40 stim amplitude
        select_stars_type = [1]; % 0: dot stimulus; 1: triangle stimulus
        methods_of_select = {select_typical, 'Typical cells'};
        [diff_thistask, S_unit, NS_unit,unique_stimAmp_this] = get_stim_data_uA(stimtask_type,methods_of_select{1,1},is_separate_stimAmp,is_original_bias,select_stars_type);
        plot_stim_effect_distribution_uA(stimtask_type,monkey_included_for_analysis,effect_type,diff_thistask,S_unit,NS_unit,unique_stimAmp_this,is_separate_stimAmp,is_original_bias);
    end

    function f2p4(debug) % 20uA or 40uA, TriangleStim)
        if debug  ; dbstack;   keyboard;      end
        effect_type = 2; % 1: PSE, 2: threshold
        is_original_bias = 0; % 0: bias positive = prefer direction; 1: only original bias in MT;  2: original bias(after-before);
        is_separate_stimAmp = 1; % 默认为0，combine 20 and 40 stim amplitude
        select_stars_type = [1]; % 0: dot stimulus; 1: triangle stimulus
        methods_of_select = {select_typical, 'Typical cells'};
        [diff_thistask, S_unit, NS_unit,unique_stimAmp_this] = get_stim_data_uA(stimtask_type,methods_of_select{1,1},is_separate_stimAmp,is_original_bias,select_stars_type);
        plot_stim_effect_distribution_uA(stimtask_type,monkey_included_for_analysis,effect_type,diff_thistask,S_unit,NS_unit,unique_stimAmp_this,is_separate_stimAmp,is_original_bias);
    end

    function f2p5(debug) % 20uA + 40uA, Dot+Triangle
        if debug  ; dbstack;   keyboard;      end
        effect_type = 2; % 1: PSE, 2: threshold
        is_original_bias = 0; % 0: bias positive = prefer direction; 1: only original bias in MT;  2: original bias(after-before);
        is_separate_stimAmp = 0; % 默认为0，combine 20 and 40 stim amplitude
        select_stars_type = [0 1]; % 0: dot stimulus; 1: triangle stimulus
        methods_of_select = {select_typical, 'Typical cells'};
        [diff_thistask, S_unit, NS_unit,unique_stimAmp_this] = get_stim_data_uA(stimtask_type,methods_of_select{1,1},is_separate_stimAmp,is_original_bias,select_stars_type);
        plot_stim_effect_distribution_uA(stimtask_type,monkey_included_for_analysis,effect_type,diff_thistask,S_unit,NS_unit,unique_stimAmp_this,is_separate_stimAmp,is_original_bias);
        
        figure(figN-1);
        ha = get(gcf,'children');
        subplot_this = [];
        % get subplot
        for i = 1:length(ha)
            if strcmp(ha(i).Title.String,'Translation')
                subplot_this(i) = 1;
            elseif strcmp(ha(i).Title.String,'Roll')
                subplot_this(i) = 2;
            elseif strcmp(ha(i).Title.String,'Flow pattern')
                subplot_this(i) = 3;
            end
        end
        
        hold on
        for hr = 1:3
            sig_num(hr) = sum(S_unit{effect_type}(:,hr));
            nosig_num(hr) = sum(NS_unit{effect_type}(:,hr));
            sig_pos_num(hr) = sum(S_unit{effect_type}(:,hr) & diff_thistask{effect_type}(:,hr)>0);
            
            ylimit = ha(subplot_this==hr).YLim;
            xlimit = ha(subplot_this==hr).XLim;
            str = [];
            str = [sprintf('Sig num = %g, Sig&Pos num = %g',sig_num(hr),sig_pos_num(hr)) newline ...
                sprintf('Sig = %.4g%%, Sig&Pos = %.4g%%',sig_num(hr)/(sig_num(hr)+nosig_num(hr))*100,sig_pos_num(hr)/(sig_num(hr))*100)];
            text(ha(subplot_this==hr),xlimit(1),ylimit(2)/2,str);
        end    
    end

    function f2p6(debug) % 20uA or 40uA, Dot+Triangle
        if debug  ; dbstack;   keyboard;      end
        effect_type = 2; % 1: PSE, 2: threshold
        is_original_bias = 0; % 0: bias positive = prefer direction; 1: only original bias in MT;  2: original bias(after-before);
        is_separate_stimAmp = 0; % 默认为0，combine 20 and 40 stim amplitude
        select_stars_type = [0 1]; % 0: dot stimulus; 1: triangle stimulus
        methods_of_select = {select_typical, 'Typical cells'};
        [diff_thistask, S_unit, NS_unit,unique_stimAmp_this] = get_stim_data_uA(stimtask_type,methods_of_select{1,1},is_separate_stimAmp,is_original_bias,select_stars_type);
        plot_stim_effect_distribution_uA(stimtask_type,monkey_included_for_analysis,effect_type,diff_thistask,S_unit,NS_unit,unique_stimAmp_this,is_separate_stimAmp,is_original_bias);
    end

    function f2p7(debug) % Threshold changed between different stim current or targets num
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        effect_type = 2; % 1: PSE, 2: threshold
        select_stars_type = [0 1]; % 0: dot stimulus; 1: triangle stimulus
        is_separate_stimAmp = 1; % 默认为0，combine 20 and 40 stim amplitude
        is_original_bias = 0;  % 1: only original bias in MT;  2: original bias(after-before);   3: bias positive = prefer direction
        
        % ---------------------- get data ----------------------------
        stimtask_type_this{1,1} = [2,3]; stimtask_type_this{1,2} = []; % fine and 2-target
        stimtask_type_this{2,1} = []; stimtask_type_this{2,2} = [5,6]; % coarse and 2-target
        stimtask_type_this{3,1} = [1]; stimtask_type_this{3,2} = []; % fine and 4-target
        stimtask_type_this{4,1} = []; stimtask_type_this{4,2} = [4]; % coarse and 4-target
        
        for i = 1:4 % fine2T, coarse2T, fine4T, coarse4T
            diff_thistask = [];
            [diff_thistask] = get_stim_data_uA({stimtask_type_this{i,:}},methods_of_select{1,1},is_separate_stimAmp,is_original_bias,select_stars_type);
            Thre_20uA{i,1} = diff_thistask{1,effect_type}(~isnan(diff_thistask{1,effect_type}(:,1)),1); % translation % remove NaN
            Thre_20uA{i,2} = diff_thistask{1,effect_type}(~isnan(diff_thistask{1,effect_type}(:,2)),2); % rotation
            
            Thre_40uA{i,1} = diff_thistask{2,effect_type}(~isnan(diff_thistask{2,effect_type}(:,1)),1); % translation
            Thre_40uA{i,2} = diff_thistask{2,effect_type}(~isnan(diff_thistask{2,effect_type}(:,2)),2); % rotation
        end
        
        
        % ------------------------ Plot ------------------------------
        % 1. stim current comparison
        set(figure(figN),'name','stim current comparison','pos',[100,50,600,800]); clf
        for i = 1:4 % fine2T, coarse2T, fine4T, coarse4T
            for hr = 1:2
                subplot(4,2,hr+(i-1)*2)
                hold on
                
                xbins = linspace(min([Thre_20uA{i,hr};Thre_40uA{i,hr}]),max([Thre_20uA{i,hr};Thre_40uA{i,hr}]),9);
                h1 = histogram(Thre_20uA{i,hr},xbins,'facecolor',c{hr});
                h2 = histogram(Thre_40uA{i,hr},xbins,'facecolor',c{hr+10});
                
                h1.Normalization = 'probability';
                h2.Normalization = 'probability';
                
                ymax = max([h1.Values h2.Values]);
                ylim([0 ymax+ymax/2]);
                xlim([min(xbins)-(xbins(2)-xbins(1)) max(xbins)+(xbins(2)-xbins(1))]);
                
                % mean
                plot(mean(Thre_20uA{i,hr}),ymax+ymax/4,'v','color',c{hr});
                plot(mean(Thre_40uA{i,hr}),ymax+ymax/4,'v','color',c{hr+10});
                
                % non-paired ttest
                [~,p_ttest(i,hr)] = ttest2(Thre_20uA{i,hr},Thre_40uA{i,hr}); % nan 为该task没有同时做过两种current stim
                text(xbins(1),ymax+ymax/3,sprintf('t-test, p=%0.3g',p_ttest(i,hr)));
                
                if i == 1 && hr == 1
                    ylabel('Cases/N');
                    legend('20uA','40uA');
                end
                
                if hr == 1
                    if i == 1
                        text(min(xbins)-3*(xbins(2)-xbins(1)),ymax/2,'Fine 2-Tar','rotation',90);
                    elseif i == 2
                        text(min(xbins)-3*(xbins(2)-xbins(1)),ymax/2,'Coarse 2-Tar','rotation',90);
                    elseif i ==3
                        text(min(xbins)-3*(xbins(2)-xbins(1)),ymax/2,'Fine 4-Tar','rotation',90);
                    else
                        text(min(xbins)-3*(xbins(2)-xbins(1)),ymax/2,'Coarse 4-Tar','rotation',90);
                    end
                end
            end
        end
        SetFigure(13)
        suptitle('Threshold change between different stim current');
        figN = figN+1;
        
        % 合并2T和4T，只比较20uA和40uA
        set(figure(figN),'name','stim current comparison','pos',[300,50,600,600]); clf
        for i = 1:2 % fine2T+fine4T, coarse2T+coarse4T
            for hr = 1:2
                subplot(2,2,hr+(i-1)*2)
                hold on
                
                if i == 1 % fine
                    all_20 = [Thre_20uA{1,hr};Thre_20uA{3,hr}];
                    all_40 = [Thre_40uA{1,hr};Thre_40uA{3,hr}];
                elseif i == 2 % coarse
                    all_20 = [Thre_20uA{2,hr};Thre_20uA{4,hr}];
                    all_40 = [Thre_40uA{2,hr};Thre_40uA{4,hr}];
                end
                
                xbins = linspace(min([all_20;all_40]),max([all_20;all_40]),9);
                h1 = histogram(all_20,xbins,'facecolor',c{hr});
                h2 = histogram(all_40,xbins,'facecolor',c{hr+10});
                
                h1.Normalization = 'probability';
                h2.Normalization = 'probability';
                
                ymax = max([h1.Values h2.Values]);
                ylim([0 ymax+ymax/2]);
                xlim([min(xbins)-(xbins(2)-xbins(1)) max(xbins)+(xbins(2)-xbins(1))]);
                
                % mean
                plot(mean(all_20),ymax+ymax/4,'v','color',c{hr});
                plot(mean(all_40),ymax+ymax/4,'v','color',c{hr+10});
                
                % non-paired ttest
                [~,p_ttest2(i,hr)] = ttest2(all_20,all_40); % nan 为该task没有同时做过两种current stim
                text(xbins(1),ymax+ymax/3,sprintf('t-test, p=%0.3g',p_ttest2(i,hr)));
                
                if i == 1 && hr == 1
                    
                    legend('20uA','40uA');
                end
                ylabel('Cases/N');
                xlabel('Threshold change');
                
                if hr == 1
                    if i == 1
                        title('Fine task')
                    elseif i == 2
                        title('Coarse task')
                    end
                end
            end
        end
        SetFigure(13)
        suptitle('Threshold change between different stim current');
        figN = figN+1;
        

        % 2. task type 
        % (only select 40uA stim current)
        % update, combine 20uA and 40uA
        
        set(figure(figN),'name','target num comparison','pos',[1000,50,800,800]); clf
        for fc = 1:2 % fine or coarse
            for hr = 1:2 % T or R
                subplot(2,2,hr+2*(fc-1))
                hold on
                
                data_2T = [Thre_20uA{fc,hr}; Thre_40uA{fc,hr}]; % 2-target
                data_4T = [Thre_20uA{fc+2,hr}; Thre_40uA{fc+2,hr}]; % 4-target
                data_2T(isnan(data_2T)) = [];
                data_4T(isnan(data_4T)) = [];
                
                xbins = linspace(min([data_2T;data_4T]),max([data_2T;data_4T]),9);
                h1 = histogram(data_2T,xbins,'edgecolor',c{hr+10},'DisplayStyle','stairs','linewidth',2); % 2T
                h2 = histogram(data_4T,xbins,'edgecolor',c{hr},'DisplayStyle','stairs','linewidth',2);  % 4T
                h1.Normalization = 'probability';
                h2.Normalization = 'probability';
                
                ymax = max([h1.Values h2.Values]);
                ylim([0 ymax+ymax/2]);
                xlim([min(xbins)-(xbins(2)-xbins(1)) max(xbins)+(xbins(2)-xbins(1))]);
                xlabel('threshold change');
                
                % mean
                plot(nanmean(data_2T),ymax+ymax/4,'v','color',c{hr+10});
                plot(nanmean(data_4T),ymax+ymax/4,'v','color',c{hr});
                
                % non-paired ttest
                [~,p_ttest2(fc,hr)] = ttest2(data_2T,data_4T);
                
                % anotation
                str = [sprintf('t-test, p=%0.3g',p_ttest2(fc,hr)) newline ...
                    sprintf('2-T N = %0.3g',length(data_2T)) newline sprintf('4-T N = %0.3g',length(data_4T))];
                text(xbins(1),ymax+ymax/3,str);
                
                legend('2-Target','4-Target');
                if fc == 1 && hr == 1
                    ylabel('Cases/N');
                end
                
                if fc == 1
                    title('Fine task');
                else
                    title('Coarse task');
                end
                
            end
        end
        SetFigure(13)
        suptitle('Threshold change between different target number');
        figN = figN+1;
    end

    function f2p8(debug) % Rotation Thr. changed between different targets num in same site
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        effect_type = 2; % 1: PSE, 2: threshold
        select_stars_type = [0 1]; % 0: dot stimulus; 1: triangle stimulus
        is_separate_stimAmp = 0; % 默认为0，combine 20 and 40 stim amplitude
        is_original_bias = 0;  % 1: only original bias in MT;  2: original bias(after-before);   3: bias positive = prefer direction
        
        % ---------------------- get data ----------------------------
        % all for rotation task threshold change
        Thr_change_temp{1}(:,1) = diffThresh_task{3}(:,2); % fine2T
        Thr_change_temp{1}(:,2) = diffThresh_task{1}(:,2); % fine4T
        Thr_change_temp{2}(:,1) = diffThresh_task{6}(:,2); % coarse2T
        Thr_change_temp{2}(:,2) = diffThresh_task{4}(:,2); % coarse4T
        
        % remove nan and save site with both 2 and 4 target task
        Thr_change = cellfun(@(x) x(sum(~isnan(x),2)==2,:),Thr_change_temp,'UniformOutput',0);
        
        % ------------------------ Plot ------------------------------
        % scatter plot
        set(figure(figN),'name','Rotation threshold change between different target num in same site','pos',[1000,50,600,400]); clf
        for fc = 1:2 % fine or coarse
            subplot(1,2,fc)
            scatter(Thr_change{fc}(:,1),Thr_change{fc}(:,2),'markerfacecolor',c{2});
            SameScale(1);
            xlabel('2-Target task');
            ylabel('4-Target task');
            if fc == 1
                title('Fine task Thr. change');
            else
                title('Coarse task Thr. change');
            end
            hold on
            xlimit = xlim;
            ylimit = ylim;
            plot(xlimit,[0 0],'k:');
            plot([0 0],ylimit,'k:');
            text(xlimit(1)+(xlimit(2)-xlimit(2)/1.1),ylimit(2)-(ylimit(2)-ylimit(2)/1.1),sprintf('N = %g',size(Thr_change{1},1)));
        end
        SetFigure(13)
        suptitle('Rotation threshold change between different target num');
    end

    function f2p9(debug) % 'Stim induced CR change'
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        
        % ---------------------- get data ----------------------------
        % 合并所选的task tpye
        for fc = 1:length(stimtask_type)  % fine or coarse 分开
            ctrl_cr{fc} = [];
            stim_cr{fc} = [];
            ctrl_cr_true{fc} = [];
            stim_cr_true{fc} = [];
            
            if ~isempty(stimtask_type{fc})
                for i = 1:length(stimtask_type{fc}) % 具体task, 需要combine一起，所以最后的index没有i
                    select = logical(methods_of_select{1,1}); % select for stimAmp, stimtask_type_this and stars_type
                    
                    temp1 = ctrl_CR{stimtask_type{fc}(i)};
                    temp2 = stim_CR{stimtask_type{fc}(i)};
                    temp3 = ctrl_CR_true{stimtask_type{fc}(i)};
                    temp4 = stim_CR_true{stimtask_type{fc}(i)};
                    
                    temp1(~select,:) = nan; % remove undesired unit
                    temp2(~select,:) = nan; % remove undesired unit
                    temp3(~select,:) = nan; % remove undesired unit
                    temp4(~select,:) = nan; % remove undesired unit
                    
                    ctrl_cr{fc} = [ctrl_cr{fc}; temp1]; % combine selected stimtask_type_this
                    stim_cr{fc} = [stim_cr{fc}; temp2]; % combine selected stimtask_type_this
                    ctrl_cr_true{fc} = [ctrl_cr_true{fc}; temp3]; % combine selected stimtask_type_this
                    stim_cr_true{fc} = [stim_cr_true{fc}; temp4]; % combine selected stimtask_type_this
                end
            end
        end
        
        % CR select (ctrl) > 65%
        mt_performance_cri = [65 100];
        for fc = 1:2
            perform_select = logical(ctrl_cr{fc}(:,3) >= mt_performance_cri(1) & ctrl_cr{fc}(:,3) <= mt_performance_cri(2) ); % only select ctrl CR > 65%
            ctrl_cr{fc}(~perform_select,3) = nan;
            ctrl_cr_true{fc}(~perform_select,3) = nan;
            stim_cr{fc}(~perform_select,3) = nan;
            stim_cr_true{fc}(~perform_select,3) = nan;
        end
        
        % ------------------------ Plot ------------------------------
        %%% paired plot
        for i = 1:2 % 1 for normal CR, 2 for true CR (in 4AFC)
            if i == 1
                set(figure(figN),'Position',[100,50 1500,700], 'Name', 'Induced CR change');clf
                ctrl_this = ctrl_cr;
                stim_this = stim_cr;
            else
                set(figure(figN),'Position',[200,100 1500,700], 'Name', 'Induced true CR change');clf
                ctrl_this = ctrl_cr_true;
                stim_this = stim_cr_true;
            end

            for fc = 1:2
                for hr = 1:3
                    subplot(2,3,hr+(fc-1)*3)
                    num_this = sum(~any(isnan([ctrl_this{fc}(:,hr) stim_this{fc}(:,hr)]),2));
                    plot([1 2],[ctrl_this{fc}(:,hr) stim_this{fc}(:,hr)],'ko-');
                    xlim([0.5 2.5]);
                    
                    if i == 1
                        ylim([40 100]);
                    else
                        ylim([0 100]);
                    end
                    set(gca,'xtick',[1 2],'xticklabel',{'Nonstim.','Microstim.'});
                    axis square
                    
                    % mean
                    hold on
                    mean_this = nanmean([ctrl_this{fc}(:,hr) stim_this{fc}(:,hr)]);
                    plot([1 2],mean_this,'o-','color',c{hr});
                    ctrl_sem_cr = nanstd(ctrl_this{fc}(:,hr)) / sqrt(sum(~isnan(ctrl_this{fc}(:,hr))));
                    stim_sem_cr = nanstd(stim_this{fc}(:,hr)) / sqrt(sum(~isnan(stim_this{fc}(:,hr))));
                    
                    
                    %                     errorbar([1 2],nanmean([ctrl_this{fc}(:,hr) stim_this{fc}(:,hr)]),[ctrl_sem_cr stim_sem_cr],'o-','color',c{hr});
                    
                    % decrease num
                    dec_num = sum((ctrl_this{fc}(:,hr)-stim_this{fc}(:,hr)) > 0);
                    dec_pro = dec_num/num_this*100;
                    
                    % decrease amplitude based on non-stim
                    dec_amp = mean_this(1)-mean_this(2);
                    dec_amp_pro = (mean_this(1)-mean_this(2))/mean_this(1) * 100;
                    
                    % ttest
                    [~,p_t{fc,hr}] = ttest(ctrl_this{fc}(:,hr),stim_this{fc}(:,hr));
                    % signtest
                    p_s{fc,hr} = signtest(ctrl_this{fc}(:,hr),stim_this{fc}(:,hr));
                    
                    str = [];
                    str = [sprintf('ttest, p = %.2g',p_t{fc,hr}) newline sprintf('signtest, p = %.2g',p_s{fc,hr}) newline sprintf('all N = %g',num_this) newline sprintf('Decrease num = %g (%.3g%%)',dec_num, dec_pro) ...
                        newline sprintf('Decrease amp. = %2.2g (%2.2g%%)',dec_amp,dec_amp_pro)];
                    text(0.7,60,str,'color','r');
                    
                    if fc == 1
                        fine_or_coarse = 'fine';
                    else
                        fine_or_coarse = 'coarse';
                    end
                    
                    if hr == 1
                        h_or_r = ' Translation';
                    elseif hr == 2
                        h_or_r = ' Rotation';
                    else
                        h_or_r = ' Motion type (CR>65%)';
                    end
                    
                    str = [fine_or_coarse,h_or_r, ' task'];
                    title(str)
                end
            end
            
            SetFigure(15);
            [~,~,~,monkey] = SetTitle(stimtask_type,monkey_included_for_analysis);
            if i == 1
                str = ['Stim induced CR change', monkey];
            else
                str = ['Stim induced true CR change', monkey];
            end
            suptitle(str);
            figN = figN + 1;
        end
    end

    function f2p10(debug) % 'CR change between different stim current'
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        effect_type = 2; % 1: PSE, 2: threshold
        select_stars_type = [0 1]; % 0: dot stimulus; 1: triangle stimulus
        is_original_bias = 0;  % 1: only original bias in MT;  2: original bias(after-before);   3: bias positive = prefer direction
        
        % ---------------------- get data ----------------------------
        ctrl_CR_all = [];
        stim_CR_all = [];
        stimAmp_all = [];
        
        for fc = 1:length(stimtask_type) % fine or coarse 分开
            ctrl_CR_all{fc} = [];
            stim_CR_all{fc} = [];
            stimAmp_all{fc} = [];
            
            if ~isempty(stimtask_type{fc})
                for i = 1:length(stimtask_type{fc}) % 具体task, 需要combine一起，所以最后的index没有i
                    select = logical(methods_of_select{1,1}); % select for stimAmp, stimtask_type_this and stars_type
                    temp1 = []; temp2 = []; temp3 = [];
                    temp1 = ctrl_CR{stimtask_type{fc}(i)};
                    temp2 = stim_CR{stimtask_type{fc}(i)};
                    temp3 = stimAmp(:,stimtask_type{fc}(i));

                    temp1(~select,:) = nan; % remove undesired unit
                    temp2(~select,:) = nan; % remove undesired unit
                    temp3(~select,:) = nan; % remove undesired unit
                    
                    ctrl_CR_all{fc} = [ctrl_CR_all{fc}; temp1]; % combine selected stimtask_type_this
                    stim_CR_all{fc} = [stim_CR_all{fc}; temp2]; % combine selected stimtask_type_this
                    stimAmp_all{fc} = [stimAmp_all{fc}; temp3]; % combine selected stimtask_type_this
                end
                CR_ratio{fc} = stim_CR_all{fc}./ctrl_CR_all{fc};
            end
        end

        
        % ------------------------ Plot ------------------------------
        % hist plot
        set(figure(figN),'Position', [50,50,1500,1000], 'Name', 'CR change and stim current');clf
        xbins = linspace(0.5,1.5,13);
        
        for fc = 1:2 % fine and coarse
            for hr = 1:3
                all_20 = CR_ratio{fc}(stimAmp_all{fc}==20,hr);
                all_40 = CR_ratio{fc}(stimAmp_all{fc}==40,hr);
                
                % remove nan
                all_20 = all_20(~isnan(all_20));
                all_40 = all_40(~isnan(all_40));
                
                subplot(2,3,hr+(fc-1)*3)
                if ~isnan([all_20;all_40])
                    hold on
                    h1 = histogram(all_20,xbins,'facecolor',c{hr});
                    h2 = histogram(all_40,xbins,'facecolor',c{hr+10});
                    %                 h1.Normalization = 'probability';
                    %                 h2.Normalization = 'probability';
                    ymax = max([h1.Values h2.Values]);
                    ylim([0 ymax+ymax/2]);
                    xlim([min(xbins)-(xbins(2)-xbins(1)) max(xbins)+(xbins(2)-xbins(1))]);
                    % median
                    plot(median(all_20),ymax+ymax/4,'v','color',c{hr});
                    plot(median(all_40),ymax+ymax/4,'v','color',c{hr+10});
                    
                    % sign test with 1
                    p_sign_20 = signtest(all_20,1);
                    p_sign_40 = signtest(all_40,1);
                    str = [sprintf('20uA: p = %.2g (sign test)',p_sign_20) newline sprintf('40uA: p = %.2g',p_sign_40) newline ...
                        sprintf('n_20 = %g', length(all_20)) newline sprintf('n_40 = %g', length(all_40))];
                    text(xbins(1),ymax,str);
                    
                    %                 % non-paired ttest between 20 and 40 uA
                    %                 [~,p_ttest2] = ttest2(all_20,all_40); % nan 为该task没有同时做过两种current stim
                    %                 text(xbins(1),ymax+ymax/3,sprintf('non-paired t-test, p=%0.3g',p_ttest2));
                    
                    % KW test for unequal sample size, nonparametic
                    if ~isempty(all_20) && ~isempty(all_40)
                        X = [all_20 ones(size(all_20)); all_40 2.*ones(size(all_40))];
                        [p] = KWtest(X,0.05);
                        text(xbins(1),ymax+ymax/3,sprintf('KW test, p=%0.2g',p));
                    else
                        text(xbins(1),ymax+ymax/3,'KW test, p=NaN');
                    end
                    
                    
                    axis square
                    ylabel('Cases/N');
                    xlabel('CR ratio');
                    
                    if fc == 1 && hr == 1
                        legend('20uA','40uA');
                    end
                    
                    if fc == 1
                        title('fine task');
                    else
                        title('coarse task');
                    end
                end
            end
        end
        SetFigure(15);
        [~,~,~,monkey] = SetTitle(stimtask_type,monkey_included_for_analysis);
        str = ['CR change across stimAmp ', monkey];
        suptitle(str);
        figN = figN + 1;
        
        % bar plot
        set(figure(figN),'Position', [50,50,800,800], 'Name', 'CR change and stim current');clf
        for fc = 1:2 % fine and coarse
            for hr = 1:3
                all_20 = CR_ratio{fc}(stimAmp_all{fc}==20,hr);
                all_40 = CR_ratio{fc}(stimAmp_all{fc}==40,hr);
                % remove nan
                all_20 = all_20(~isnan(all_20));
                all_40 = all_40(~isnan(all_40));
                
                subplot(2,3,hr+(fc-1)*3)
                hold on
                plot([0.5 2.5],[1 1],'k--');
                
                % bar
                bar([1 2],[mean(all_20) mean(all_40)],0.6,'facecolor','none','edgecolor','k','linewidth',1.5);
                
                % scatter
                plot(rand(size(all_20))*0.4-0.2+1,all_20,'.','color',c{hr});
                plot(rand(size(all_40))*0.4-0.2+2,all_40,'.','color',c{hr});
%                 sem_20 = std(all_20)/sqrt(length(all_20));
%                 sem_40 = std(all_40)/sqrt(length(all_40));
                % mean±std
                errorbar([1 2],[mean(all_20) mean(all_40)],[std(all_20) std(all_40)],'color','k','LineStyle','none','linewidth',1.5);
                ylim([0.5 1.5]);
                xlim([0.5 2.5]);
                xticks([1 2]); xticklabels({'20uA','40uA'});
                xlabel('Electrical current');
                
                % t test with 1
                [~,p_20] = ttest(all_20,1);
                [~,p_40] = ttest(all_40,1);
                if ~isnan(p_20)
                    str1 = [sprintf('%2.2g±%2.2g',mean(all_20),std(all_20)) newline sprintf('p = %2.2g (t test)',p_20) newline sprintf('n_2_0 = %g', length(all_20))];
                    text(1,0.6,str1,'horizontalalignment','center');
                end
                
                if ~isnan(p_40)
                    str2 = [sprintf('%2.2g±%2.2g',mean(all_40),std(all_40)) newline sprintf('p = %2.2g (t test)',p_40) newline sprintf('n_4_0 = %g', length(all_40))];
                    text(2,0.6,str2,'horizontalalignment','center');
                end
                
                % non-paired ttest between 20 and 40 uA
                [~,p_ttest2] = ttest2(all_20,all_40); % nan 为该task没有同时做过两种current stim
                str = ['non-paired t-test' newline sprintf('p=%0.2g',p_ttest2)];
                text(0.6,1.4,str);
                
                if hr == 1
                    ylabel('CR ratio');
                end
                
                if fc == 1 && hr == 2
                    title('fine task');
                elseif fc == 2  && hr == 2
                    title('coarse task');
                end
            end
        end
        SetFigure(10);
        [~,~,~,monkey] = SetTitle(stimtask_type,monkey_included_for_analysis);
        str = ['CR change across stimAmp ', monkey];
        suptitle(str);
        figN = figN + 1;
    end

    function f2p11(debug) % 'CR change between different target num'
         if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        effect_type = 2; % 1: PSE, 2: threshold
        select_stars_type = [0 1]; % 0: dot stimulus; 1: triangle stimulus
        is_original_bias = 0;  % 1: only original bias in MT;  2: original bias(after-before);   3: bias positive = prefer direction
        
        
        % ---------------------- get data ----------------------------
        stimtask_type_this{1} = [2,3];  % fine and 2-target
        stimtask_type_this{2} = [5,6]; % coarse and 2-target
        stimtask_type_this{3} = [1]; % fine and 4-target
        stimtask_type_this{4} = [4]; % coarse and 4-target
        
        % data for 2T and 4T
        CR_ratio_2T = cell(1,2);
        CR_ratio_4T = cell(1,2);
        
        
        for fc = 1:2 % fine or coarse   % fine2T, coarse2T, fine4T, coarse4T
            select = logical(methods_of_select{1,1}); % select for stimAmp, stimtask_type_this and stars_type
            
            % 2T
            for i = 1:length(stimtask_type_this{fc})
                temp = [];
                temp = stim_CR{stimtask_type_this{fc}(i)} ./ ctrl_CR{stimtask_type_this{fc}(i)}; % CR ratio = stim / ctrl
                temp(~select,:) = nan; % remove undesired unit
                CR_ratio_2T{fc} = [CR_ratio_2T{fc}; temp];
            end
            
            
            % 4T
            for i = 1:length(stimtask_type_this{fc+2})
                temp = [];
                temp = stim_CR{stimtask_type_this{fc+2}(i)} ./ ctrl_CR{stimtask_type_this{fc+2}(i)}; % CR ratio = stim / ctrl
                temp(~select,:) = nan; % remove undesired unit
                CR_ratio_4T{fc} = [CR_ratio_4T{fc}; temp];
            end
        end
        
        % ------------------------ Plot ------------------------------
        % hist plot
        set(figure(figN),'name','CR change between different target num','pos',[50,50,800,800]); clf
        for fc = 1:2 % fine or coarse   % fine2T, coarse2T, fine4T, coarse4T
            for hr = 1:2 % T, R, MT
                subplot(2,2,hr+2*(fc-1))
                hold on
                
                data_2T_this = CR_ratio_2T{fc}(:,hr);
                data_4T_this = CR_ratio_4T{fc}(:,hr);
                data_2T_this(isnan(data_2T_this)) = [];
                data_4T_this(isnan(data_4T_this)) = [];
                
                xbins = linspace(0.5,1.5,9);
                %                 h1 = histogram(h1_data,xbins,'facecolor',c{hr}); % 2T
                %                 h2 = histogram(h2_data,xbins,'facecolor',c{hr+10});  % 4T
                
                h1 = histogram(data_2T_this,xbins,'edgecolor',c{hr+10},'DisplayStyle','stairs','linewidth',2); % 2T
                h2 = histogram(data_4T_this,xbins,'edgecolor',c{hr},'DisplayStyle','stairs','linewidth',2);  % 4T
                
                h1.Normalization = 'probability';
                h2.Normalization = 'probability';
                
                ymax = max([h1.Values h2.Values]);
                ymax = 0.6;
                ylim([0 ymax+ymax/2]);
                xlim([min(xbins)-(xbins(2)-xbins(1)) max(xbins)+(xbins(2)-xbins(1))]);
                xlabel('CR ratio');
                
                % mean
                plot(nanmean(data_2T_this),ymax+ymax/4,'v','color',c{hr+10});
                plot(nanmean(data_4T_this),ymax+ymax/4,'v','color',c{hr});
                
                % non-paired ttest
                [~,p_ttest2(fc,hr)] = ttest2(data_2T_this,data_4T_this);
                
                % anotation
                str = [sprintf('t-test, p=%0.3g',p_ttest2(fc,hr)) newline ...
                    sprintf('2-T N = %0.3g',length(data_2T_this)) newline sprintf('4-T N = %0.3g',length(data_4T_this))];
                text(xbins(1),ymax+ymax/3,str);
                
                legend('2-Target','4-Target');
                if fc == 1 && hr == 1
                    ylabel('Cases/N');
                end
                
                if fc == 1
                    title('Fine task');
                else
                    title('Coarse task');
                end
            end
        end
        SetFigure(13)
        suptitle('CR change between different target number');
        figN = figN+1;
        
    end


%% microstimulation effect and modulation
    function f3p1(debug) % 'PSE shift affected by Modulation Index''
        if debug  ; dbstack;   keyboard;      end
        
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        effect_type = 1; % 1: PSE, 2: threshold
        select_stars_type = [0 1]; % 0: dot stimulus; 1: triangle stimulus
        is_separate_stimAmp = 0;
        
        % First part: as more modulation index as possible, use original PSE shift
        % Second part: only select the d' and neural sensitivity, use PSE shift adapted to prefer direction
        
        % ---------------------- get data ----------------------------
        
%         is_original_bias = 2; % original bias: after-before  >0: bias to left/CW; <0: bias to right/CCW
%         [PSE_shift_ori_temp, S_unit_ori_temp, NS_unit_ori_temp, ~] = get_stim_data_uA(stimtask_type,methods_of_select{1,1},is_separate_stimAmp,is_original_bias,select_stars_type);
%         PSE_shift_ori = PSE_shift_ori_temp{effect_type};
%         S_unit_ori = S_unit_ori_temp{effect_type};
%         NS_unit_ori = NS_unit_ori_temp{effect_type};
        
        is_original_bias = 4;  % original bias and normalized (after-before  >0: bias to left/CW; <0: bias to right/CCW)
        [PSE_shift_temp, S_unit_temp, NS_unit_temp, ~] = get_stim_data_uA(stimtask_type,methods_of_select{1,1},is_separate_stimAmp,is_original_bias,select_stars_type);
        PSE_shift_ori_nor = PSE_shift_temp{effect_type};
        S_unit_ori_nor = S_unit_temp{effect_type};
        NS_unit_ori_nor = NS_unit_temp{effect_type};
        
        % find the outlier, may be different from PSE_shift_nor
        for i = 1:2
            reasonable_range = [nanmean(PSE_shift_ori_nor(:,i)) - 3*nanstd(PSE_shift_ori_nor(:,i)), nanmean(PSE_shift_ori_nor(:,i)) + 3*nanstd(PSE_shift_ori_nor(:,i))];
            outlier(:,i) = logical(PSE_shift_ori_nor(:,i)<reasonable_range(1) | PSE_shift_ori_nor(:,i)>reasonable_range(2));
        end

   
        % Modulation Index
        % some index have directionality (positive or negative, represente to L/R CW/CCW)：
        modulation_index{1} = cell2mat({group_result.d_prime_fine}'); % d' = (Rright45 - Rleft45) / sqrt(1/2*(STDright^2+STDleft^2))
        modulation_index{2} = cell2mat({group_result.d_prime_coarse}'); % d' = (Rright90 - Rleft90) / sqrt(1/2*(STDright^2+STDleft^2))
        modulation_index{3} = cell2mat({group_result.ahead_slope}');
        % some index do not have directionality
        modulation_index{4} = cell2mat({group_result.d_prime_maxDiff}'); % d' = (Rmax - Rmin) / sqrt(1/2*(STDmax^2+STDmin^2))
        modulation_index{5} = cell2mat({group_result.d_prime_maxAxis}'); % d' = (Rmax180 - Rmin180) / sqrt(1/2*(STDmax180^2+STDmin180^2))
        modulation_index{6} = cell2mat({group_result.d_prime_maxFR}'); % d' = (Rmax - Rc-lat) / sqrt(1/2*(STDmax^2+STDc-lat^2))
        modulation_index{7} = cell2mat({group_result.DDI_includeExpCon}'); % DDI = delta / (delta + 2 * sqrt(SSE/(N-M)))
        modulation_index{8} = cell2mat({group_result.DDI_excludeExpCon}'); % DDI = delta / (delta + 2 * sqrt(SSE/(N-M)))
        modulation_index{9} = cell2mat({group_result.HTI}'); % HTI = |ΣRi| / Σ|Ri|   AMPvectorsum / sum(|mean-spon|)  (Ri = mean-spon)
        modulation_index{10} = cell2mat({group_result.neu_sen_fine}'); %  lg(|d'45|)
        modulation_index{11} = cell2mat({group_result.neu_sen_coarse}'); %  lg(|d'90|)
        modulation_index{12} = cell2mat({group_result.MDI}'); % modulation index = (Rmax-Rmin)/(Rmax-Spon)
        modulation_index{13} = cell2mat({group_result.pair_DI}'); % 1/n * Σ( (Rfar(i)-Rnear(i)) / |Rfar(i)-Rnear(i)| + std)
        modulation_index{14} = cell2mat({group_result.nonpair_DI}'); % |Rleft-Rright| / (|Rleft-Rright| + 2*sqrt(SSE/(N-2)))
        modulation_index{15} = cell2mat({group_result.nonpair_DI2}'); % |Rleft-Rright| / (|Rleft-Rright| + STDavg)
        modulation_index{16} = cell2mat({group_result.SVM_correct_rate}');
        modulation_index{17} = cell2mat({group_result.AuROC}');
        
        tl{1} = strcat('d''','(fine)');
        tl{2} = strcat('d''','(coarse)');
        tl{3} = 'Zero slope';
        tl{4} = strcat('d''','(maxDiff)');
        tl{5} = strcat('d''','(maxAxis)');
        tl{6} = strcat('d''','(maxFR)');
        tl{7} = 'DDI(inc)';
        tl{8} = 'DDI(exc)';
        tl{9} = 'HTI';
        tl{10} = 'Neural Sen.(fine)';
        tl{11} = 'Neural Sen.(coarse)';
        tl{12} = 'Modulation index';
        tl{13} = 'Pair DI';
        tl{14} = 'Non-pair DI';
        tl{15} = 'Non-pair DI2';
        tl{16} = 'SVM correct rate';
        tl{17} = 'AuROC';
        
        modulation_index = cellfun(@(x) repmat(x,length(cell2mat(stimtask_type)),1),modulation_index,'UniformOutput',0);
        
        % ------------------------ Plot ------------------------------
        set(figure(figN),'Position', [50,540 2500,400], 'Name', 'PSE shift(origin) and Modulation Index');clf
        ha = tight_subplot(2,length(modulation_index),[.1 .02],[.15 .15],[.04 .01],1);
        % 具有方向性的modulation index(1-3)，plot with normalized original bias.散点图理想情况是反对角线
        % 绝对值越大，数值越大，discriminability越强，负号表示pref Left（CW）, 正号表示pref Right(CCW)
        % remove outlier to show and correlation
        for m = 1:3
            for hr = 1:2
                axes(ha(m+(hr-1)*length(modulation_index)));
                hold on
                plot(modulation_index{m}(NS_unit_ori_nor(:,hr) & ~outlier(:,hr),hr),PSE_shift_ori_nor(NS_unit_ori_nor(:,hr) & ~outlier(:,hr),hr),'ko','markersize',4);
                plot(modulation_index{m}(S_unit_ori_nor(:,hr) & ~outlier(:,hr),hr),PSE_shift_ori_nor(S_unit_ori_nor(:,hr) & ~outlier(:,hr),hr),'ko','markerfacecolor',c{hr},'markersize',4);
                xlimit = xlim; ylimit = ylim;
                xlimit = max(abs(xlimit));
                ylimit = max(abs(ylimit));
                xlim([-xlimit xlimit]);
                ylim([-ylimit ylimit]);
                xlimit = xlim;
                ylimit = ylim;
                plot(xlimit,[0 0],'k:'); plot(xlimit,[0 0],'k:');
                plot([0 0],ylimit,'k:'); plot(xlimit,[0 0],'k:');
                
                % Pearson correlation
                [r,p] = plot_corr_line(modulation_index{m}(~outlier(:,hr),hr),PSE_shift_ori_nor(~outlier(:,hr),hr),'MethodOfCorr','Pearson','FittingMethod',2);
                if p<0.05
                    box on
                    ax = gca;
                    ax.XColor = 'red';
                    ax.YColor = 'red';
                end
                
                % annotation
                text(xlimit(1)+(xlimit(2)-xlimit(1))/20,ylimit(2)-(ylimit(2)-ylimit(1))/20,sprintf('r=%0.3g,p=%0.3g',r,p),'fontsize',8.5);
                if hr==1; title(tl{m});end
                if hr==2
                    str2 = strcat('Pref R/CCW',repmat(32,1,5),'Pref L/CW');
                    text(0,ylimit(1)-(ylimit(2)-ylimit(1))/4,str2,'fontsize',8.5,'horizontalalignment','center');
                end
                if m==1
                    str = strcat('Bias R/CCW',repmat(32,1,5),'Bias L/CW');
                    text(xlimit(1)-(xlimit(2)-xlimit(1))/3.5,0,str,'rotation',90,'horizontalalignment','center');
                    if hr == 1
                        text(xlimit(1)-(xlimit(2)-xlimit(1))/2,0,'T Task','rotation',90,'fontsize',13,'horizontalalignment','center');
                    else
                        text(xlimit(1)-(xlimit(2)-xlimit(1))/2,0,'R Task','rotation',90,'fontsize',13,'horizontalalignment','center');
                    end
                end
            end
        end
        
        % 不具有方向性的modulation index(4-11)，只有大小，plot with abs(normalized original PSE)，理想情况是正对角线
        % remove outlier to show and correlation
        for m = 4:length(modulation_index) % 4:11  %后面几个modulation
            for hr = 1:2
                axes(ha(m+(hr-1)*length(modulation_index)))
                %                 % PSE shift 取log10
                %                 plot(modulation_index{n}(select_L,i),log10(abs(PSE_shift_ori(:,i))),'ko');
                plot(modulation_index{m}(~outlier(:,hr),hr),abs(PSE_shift_ori_nor(~outlier(:,hr),hr)),'ko','markersize',4);
                
                % correlation
                [r1,p1] = plot_corr_line(modulation_index{m}(~outlier(:,hr),hr),abs(PSE_shift_ori_nor(~outlier(:,hr),hr)),'MethodOfCorr','Pearson','FittingMethod',2);
                 if p1<0.05
                    box on
                    ax = gca;
                    ax.XColor = 'red';
                    ax.YColor = 'red';
                end
                
                % annotation
                xlimit = xlim;
                ylimit = ylim;
                text(xlimit(1)+(xlimit(2)-xlimit(1))/20,ylimit(2)-(ylimit(2)-ylimit(1))/20,sprintf('r=%0.3g,p=%0.3g',r1,p1),'fontsize',8.5);
                if hr==1 ; title(tl{m}); end
            end
        end
        
        str = ['Correlation between normalized PSE shift and tuning modulation index' newline '(Remvoe outlier)'];
        suptitle(str);
        figN = figN + 1;
    end

    function f3p2(debug) %'Comparison of Translation and rotation PSE shift'
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        effect_type = 1; % 1: PSE, 2: threshold
        select_stars_type = [0 1]; % 0: dot stimulus; 1: triangle stimulus
        is_separate_stimAmp = 0;
        %         is_original_bias = 0;  % 正值：bias to prefer direction
        is_original_bias = 3; % normlize PSE shift (divided by the threshold under non-stim trials
        
        % ---------------------- get data ----------------------------
        [PSE_shift_temp, S_unit_temp, ~, ~] = get_stim_data_uA(stimtask_type,methods_of_select{1,1},is_separate_stimAmp,is_original_bias,select_stars_type);
        PSE_shift = PSE_shift_temp{effect_type}(:,1:2); % only save h and r, not mt
        S_unit = S_unit_temp{effect_type}(:,1:2)*1;
        
        cell_type = [group_result.cell_type]';
        cell_type = repmat(cell_type,length(cell2mat(stimtask_type)),1);
        
        spiral_index = [group_result.spiral_index]';
        spiral_index = repmat(spiral_index,length(cell2mat(stimtask_type)),1);
        
        % 在two-target task中，因为函数get_stim_data_uA会叠加heading和rotation task，而不是将其合并，所以需要重新调整PSE等排布
        stimtask_type_this = cell2mat(stimtask_type);
        block_size = length(loadM_data);
        if sum(stimtask_type_this==3)>0 % 2 target fine task
            block_loc_this = find(stimtask_type_this==3);
            PSE_shift(block_size*(block_loc_this-2)+1:block_size*(block_loc_this-1),2) = PSE_shift(block_size*(block_loc_this-1)+1:block_size*block_loc_this,2); % 将rotation补充到heading那里
            PSE_shift(block_size*(block_loc_this-1)+1:block_size*block_loc_this,:) = ones(block_size,2)*9999; % 将原来rotation位置变为0，后面去掉
            S_unit(block_size*(block_loc_this-2)+1:block_size*(block_loc_this-1),2) = S_unit(block_size*(block_loc_this-1)+1:block_size*block_loc_this,2); % 将rotation补充到heading那里
            S_unit(block_size*(block_loc_this-1)+1:block_size*block_loc_this,:) = ones(block_size,2)*9999; % 将原来rotation位置变为0，后面去掉
            cell_type(block_size*(block_loc_this-1)+1:block_size*block_loc_this) = ones(block_size,1)*9999; % 将原来rotation位置变为0，后面去掉
            spiral_index(block_size*(block_loc_this-1)+1:block_size*block_loc_this) = ones(block_size,1)*9999; % 将原来rotation位置变为0，后面去掉
        end
        
        if sum(stimtask_type_this==6)>0 % 2 target coarse task
            block_loc_this = find(stimtask_type_this==6);
            PSE_shift(block_size*(block_loc_this-2)+1:block_size*(block_loc_this-1),2) = PSE_shift(block_size*(block_loc_this-1)+1:block_size*block_loc_this,2); % 将rotation补充到heading那里
            PSE_shift(block_size*(block_loc_this-1)+1:block_size*block_loc_this,:) = ones(block_size,2)*9999; % 原来的填充9999
            S_unit(block_size*(block_loc_this-2)+1:block_size*(block_loc_this-1),2) = S_unit(block_size*(block_loc_this-1)+1:block_size*block_loc_this,2); % 将rotation补充到heading那里
            S_unit(block_size*(block_loc_this-1)+1:block_size*block_loc_this,:) = ones(block_size,2)*9999;
            cell_type(block_size*(block_loc_this-1)+1:block_size*block_loc_this) = ones(block_size,1)*9999;
            spiral_index(block_size*(block_loc_this-1)+1:block_size*block_loc_this) = ones(block_size,1)*9999;
        end
        
        % 去掉0
        PSE_shift(PSE_shift(:,1)==9999,:) = [];
        S_unit(S_unit(:,1)==9999,:) = [];
        cell_type(cell_type(:,1)==9999,:) = [];
        spiral_index(spiral_index(:,1)==9999,:) = [];
        
        % remove the outlier to calculate the correlation
        for hr = 1:2
            normal_range = [nanmean(PSE_shift(:,hr))-3*nanstd(PSE_shift(:,hr)),nanmean(PSE_shift(:,hr))+3*nanstd(PSE_shift(:,hr))];
            outdot(:,hr) = logical(PSE_shift(:,hr)<normal_range(1) | PSE_shift(:,hr)>normal_range(2));
        end
        not_outdot = sum(outdot==0,2)==2;
        
        % find which do heading and rotation simultaneously
        do_both = sum(~isnan(PSE_shift),2)==2;
        site_num = sum(do_both);
        select_h_sig = logical(do_both & S_unit(:,1)==1 & S_unit(:,2)~=1);
        select_r_sig = logical(do_both & S_unit(:,1)~=1 & S_unit(:,2)==1);
        select_both_sig = logical(do_both & S_unit(:,1)==1 & S_unit(:,2)==1);
        select_not_sig = logical(do_both & S_unit(:,1)~=1 & S_unit(:,2)~=1);
        
        site_without_outlier = logical(do_both & not_outdot);
        site_num_without_outlier = sum(do_both & not_outdot);
        
        % ------------------------ Plot ------------------------------
        set(figure(figN),'name','Comparison of Translation and Rotation PSE shift','pos',[200,100,800,800]); clf
        hold on
        h1 = plot(PSE_shift(select_not_sig,1),PSE_shift(select_not_sig,2),'o','color',c{99},'markersize',7);
        h2 = plot(PSE_shift(select_h_sig,1),PSE_shift(select_h_sig,2),'+','color',c{1},'markersize',10,'linewidth',1.5);
        h3 = plot(PSE_shift(select_r_sig,1),PSE_shift(select_r_sig,2),'x','color',c{2},'markersize',10,'linewidth',1.5);
        h4 = plot(PSE_shift(select_both_sig,1),PSE_shift(select_both_sig,2),'o','color','k','markerfacecolor','k','markersize',7);

        [r,p] = plot_corr_line(PSE_shift(do_both,1),PSE_shift(do_both,2),'MethodOfCorr','Pearson','FittingMethod',1,'LineWidth',2,'LineStyles',{'-','color',c{99}}); % with outlier
        [r2,p2,h6] = plot_corr_line(PSE_shift(do_both & not_outdot,1),PSE_shift(do_both & not_outdot,2),'MethodOfCorr','Pearson','FittingMethod',1,'LineWidth',2); % without outlier
        
        axis tight
        SameScale(1);
        xlimit = xlim;
        ylimit = ylim;
        plot([0 0],ylimit,'k:');
        plot(xlimit,[0 0],'k:');
        text(xlimit(1)+(xlimit(2)-xlimit(1))/20,ylimit(2)-(ylimit(2)-ylimit(1))/20,sprintf('Pearson correlation with outlier: r=%0.3g, p=%0.3g',r,p),'fontsize',15);
        text(xlimit(1)+(xlimit(2)-xlimit(1))/20,ylimit(2)-(ylimit(2)-ylimit(1))/8,sprintf('Pearson correlation without outlier: r=%0.3g, p=%0.3g',r2,p2),'fontsize',15);
        text(xlimit(1)+(xlimit(2)-xlimit(1))/20,ylimit(2)-(ylimit(2)-ylimit(1))/2,sprintf('N= %g',site_num),'fontsize',15);
        
        xlabel('Normalized Translation PSE shift');
        ylabel('Normalized Rotation PSE shift');
        legend('Not sig','Trans. sig','Rot. sig','Both sig','Fitting with outlier','Fitting without outlier','location','best');
        SetFigure(15);
        
        til = SetTitle(stimtask_type,monkey_included_for_analysis);
        til2 = strcat('Translation V.S. Rotation in',32,til);
        title(til2)
        figN = figN+1;
        
        % 小图：without outlier
        set(figure(figN),'name','Comparison of Translation and Rotation PSE shift (without outlier)','pos',[800,100,500,500]); clf
        hold on
        plot(PSE_shift(select_not_sig & not_outdot,1),PSE_shift(select_not_sig & not_outdot,2),'o','color',c{99},'markersize',4);
        plot(PSE_shift(select_h_sig & not_outdot,1),PSE_shift(select_h_sig & not_outdot,2),'+','color',c{1},'markersize',7,'linewidth',1.2);
        plot(PSE_shift(select_r_sig & not_outdot,1),PSE_shift(select_r_sig & not_outdot,2),'x','color',c{2},'markersize',7,'linewidth',1.2);
        plot(PSE_shift(select_both_sig & not_outdot,1),PSE_shift(select_both_sig & not_outdot,2),'o','color','k','markerfacecolor','k','markersize',4);
        copyobj(h6,gca);
        axis tight
        SameScale(1);
        xlimit = xlim;
%         ylim(ylimit_ori*1.1);
        ylimit = ylim;
        plot([0 0],ylimit,'k:');
        plot(xlimit,[0 0],'k:');
        set(gca,'ytickmode','auto')
        SetFigure(15);
        box on
        %         set(gca,'ticklength',[0 0]);
        text(xlimit(1)+(xlimit(2)-xlimit(1))/20,ylimit(2)-(ylimit(2)-ylimit(1))/20,sprintf('N= %g',site_num_without_outlier),'fontsize',15);
        figN = figN+1;
        
        set(figure(figN),'name','Comparison of Translation and Rotation PSE shift (out 10 with log scale)','pos',[800,100,1000,500]); clf
        temp1 = PSE_shift(select_not_sig,:);
        temp2 = PSE_shift(select_h_sig,:);
        temp3 = PSE_shift(select_r_sig,:);
        temp4 = PSE_shift(select_both_sig,:);
        ha(1) = subplot(1,2,1); % x in 10, linnear scale
        hold on
        plot(temp1(temp1(:,1)<=10,1),temp1(temp1(:,1)<=10,2),'o','color',c{99},'markersize',4);
        plot(temp2(temp2(:,1)<=10,1),temp2(temp2(:,1)<=10,2),'+','color',c{1},'markersize',7,'linewidth',1.2);
        plot(temp3(temp3(:,1)<=10,1),temp3(temp3(:,1)<=10,2),'x','color',c{2},'markersize',7,'linewidth',1.2);
        plot(temp4(temp4(:,1)<=10,1),temp4(temp4(:,1)<=10,2),'o','color','k','markerfacecolor','k','markersize',4);
        axis tight
        SameScale(1);
        xlimit = xlim;
        ylimit = ylim;
        plot([0 0],ylimit,'k-');
        plot(xlimit,[0 0],'k-');
        
        [r,p] = plot_corr_line(PSE_shift(do_both,1),PSE_shift(do_both,2),'MethodOfCorr','Pearson','FittingMethod',1,'LineWidth',1,'LineStyles',{'-','color',c{99}}); % with outlier
        text(xlimit(1)+(xlimit(2)-xlimit(1))/20,ylimit(2)-(ylimit(2)-ylimit(1))/20,sprintf('Pearson correlation with outlier: r=%0.3g, p=%0.3g',r,p),'fontsize',15);
        
        
        
        ha(2) = subplot(1,2,2); % x out 10, log scale
        hold on
        plot(temp1(temp1(:,1)>10,1),temp1(temp1(:,1)>10,2),'o','color',c{99},'markersize',4);
        plot(temp2(temp2(:,1)>10,1),temp2(temp2(:,1)>10,2),'+','color',c{1},'markersize',7,'linewidth',1.2);
        plot(temp3(temp3(:,1)>10,1),temp3(temp3(:,1)>10,2),'x','color',c{2},'markersize',7,'linewidth',1.2);
        plot(temp4(temp4(:,1)>10,1),temp4(temp4(:,1)>10,2),'o','color','k','markerfacecolor','k','markersize',4);
        set(gca,'xscale','log');
%         semilogx
        
        linkaxes([ha(1),ha(2)],'y');
        xlim([1 100]);
        SetFigure(15);
        figN = figN+1;
        
        % 小图：cell type
        xx = PSE_shift(site_without_outlier,1);
        yy = PSE_shift(site_without_outlier,2);
        cell_type = cell_type(site_without_outlier);
        set(figure(figN),'name','Comparison of Translation and Rotation PSE shift (without outlier, cell types)','pos',[1000,300,1000,500]); clf
        subplot(1,2,1)
        hold on
        for ct = 1:3
            plot(xx(cell_type == ct),yy(cell_type == ct),'o','color',c{ct},'markersize',7);
        end
        axis tight
        SameScale(1);
        xlimit = xlim;
        ylimit = ylim;
        plot([0 0],ylimit,'k:');
        plot(xlimit,[0 0],'k:');
        set(gca,'ytickmode','auto');
        title('Across cell type');
        SetFigure(15);
        
        % across spiral_index
        max_color_value = 20;
        jet_color = colormap(jet(max_color_value));
        spiral_index = spiral_index(site_without_outlier);
        color_index = ceil(spiral_index*max_color_value/2 + max_color_value/2); % 将spiral index（-1~1）对应到max_color_value（1-20）个color上去
        
        % 修改画图顺序：spiral index=0先画，spiral index=±1后画
        [~,plot_id] = sort(abs(spiral_index));
        
        % 越Rotation，spiral index（-1）越小，color_index越小（1），颜色越蓝；
        % 越Heading，spiral index（1）越大，color index越大（max_color_value20），颜色越红
        subplot(1,2,2)
        for i = 1:length(plot_id)
            selected_color = jet_color(color_index(plot_id(i)),:);
            plot(xx(plot_id(i)),yy(plot_id(i)),'o','color',selected_color,'markerfacecolor',selected_color,'markersize',7);
            hold on
        end
        colorbar('Ticks',[0,1],'TickLabels',{'Rotation','Translation'});
        axis tight
        SameScale(1);
        xlimit = xlim;
        ylimit = ylim;
        plot([0 0],ylimit,'k:');
        plot(xlimit,[0 0],'k:');
        set(gca,'ytickmode','auto');
        title('Across spiral index');
        SetFigure(15);

        figN = figN+1;
    end

    function f3p14(debug) % 'T V.S. R V.S. Preferred OpticFlow(PSE shift)'
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        effect_type = 1; % 1: PSE, 2: threshold
        select_stars_type = [0 1]; % 0: dot stimulus; 1: triangle stimulus
        is_separate_stimAmp = 0;
        %         is_original_bias = 0;  % 正值：bias to prefer direction
        is_original_bias = 3; % normlize PSE shift (divided by the threshold under non-stim trials
        
        % ---------------------- get data ----------------------------
        % prefer optic flow similarity
        try
            RFsimilarity = [group_result.RFsimilarity]';
        catch
            run_this = strcmpi({handles.function_handles_spi{2,2}{:,1}}','Save Prefer OpticFlow similarity (local planar motion)');
            feval(handles.function_handles_spi{2,2}{run_this,2},false); % first run "Spiral Tuning" - "Spiral Tuning" - "Save Prefer OpticFlow (local planar motion)"
            RFsimilarity = [group_result.RFsimilarity]';
        end

        [PSE_shift_temp, S_unit_temp, ~, ~] = get_stim_data_uA(stimtask_type,methods_of_select{1,1},is_separate_stimAmp,is_original_bias,select_stars_type);
        PSE_shift = PSE_shift_temp{effect_type}(:,1:2); % only save h and r, not mt
        S_unit = S_unit_temp{effect_type}(:,1:2)*1;
        
        % 在two-target task中，因为函数get_stim_data_uA会叠加heading和rotation task，而不是将其合并，所以需要重新调整PSE等排布
        stimtask_type_this = cell2mat(stimtask_type);
        block_size = length(loadM_data);
        if sum(stimtask_type_this==3)>0 % 2 target fine task
            block_loc_this = find(stimtask_type_this==3);
            PSE_shift(block_size*(block_loc_this-2)+1:block_size*(block_loc_this-1),2) = PSE_shift(block_size*(block_loc_this-1)+1:block_size*block_loc_this,2); % 将rotation补充到heading那里
            PSE_shift(block_size*(block_loc_this-1)+1:block_size*block_loc_this,:) = ones(block_size,2)*9999; % 将原来rotation位置变为999，后面去掉
            S_unit(block_size*(block_loc_this-2)+1:block_size*(block_loc_this-1),2) = S_unit(block_size*(block_loc_this-1)+1:block_size*block_loc_this,2); % 将rotation补充到heading那里
        end
        
        if sum(stimtask_type_this==6)>0 % 2 target coarse task
            block_loc_this = find(stimtask_type_this==6);
            PSE_shift(block_size*(block_loc_this-2)+1:block_size*(block_loc_this-1),2) = PSE_shift(block_size*(block_loc_this-1)+1:block_size*block_loc_this,2); % 将rotation补充到heading那里
            PSE_shift(block_size*(block_loc_this-1)+1:block_size*block_loc_this,:) = ones(block_size,2)*9999; % 原来的填充9999
            S_unit(block_size*(block_loc_this-2)+1:block_size*(block_loc_this-1),2) = S_unit(block_size*(block_loc_this-1)+1:block_size*block_loc_this,2); % 将rotation补充到heading那里
        end
        
        RFsimilarity = repmat(RFsimilarity,[length(stimtask_type_this),1]);

        % 去掉999
        RFsimilarity(PSE_shift(:,1)==9999,:) = [];
        S_unit(PSE_shift(:,1)==9999,:) = [];
        PSE_shift(PSE_shift(:,1)==9999,:) = [];
        
        % 1. group 3 similarity groups
        % 2. 连续分布
        group_range = linspace(0,100,4);
        for g = 1:3
            group_this{g} = logical(RFsimilarity > group_range(g) & RFsimilarity <= group_range(g+1));
        end
        
        % remove the outlier to calculate the correlation
        for hr = 1:2
            normal_range = [nanmean(PSE_shift(:,hr))-3*nanstd(PSE_shift(:,hr)),nanmean(PSE_shift(:,hr))+3*nanstd(PSE_shift(:,hr))];
            outdot(:,hr) = logical(PSE_shift(:,hr)<normal_range(1) | PSE_shift(:,hr)>normal_range(2));
        end
        outdot = false(size(outdot)); % do not remove outlier
        not_outdot = sum(outdot==0,2)==2;
        
        % find which do heading and rotation simultaneously
        do_both = sum(~isnan(PSE_shift),2)==2;
        site_num = sum(do_both);
        select_h_sig = logical(do_both & S_unit(:,1)==1 & S_unit(:,2)~=1);
        select_r_sig = logical(do_both & S_unit(:,1)~=1 & S_unit(:,2)==1);
        select_both_sig = logical(do_both & S_unit(:,1)==1 & S_unit(:,2)==1);
        select_not_sig = logical(do_both & S_unit(:,1)~=1 & S_unit(:,2)~=1);
        
        site_num_without_outlier = sum(do_both & not_outdot);
        
        % ------------------------ Plot ------------------------------
        %%% 1. 3 groups, paired plot
        set(figure(figN),'name','grouped OpticFlow Similarity','pos',[200,100,1300,500]); clf
        for g = 1:3
            group_this{g} = logical(group_this{g} & do_both & not_outdot);
            ax(g) = subplot(1,3,g);
            num_this = size([PSE_shift(group_this{g},2),PSE_shift(group_this{g},1)],1);
            plot([1 2],[PSE_shift(group_this{g},2),PSE_shift(group_this{g},1)],'ko-');
            set(gca,'xtick',[1 2],'xticklabel',{'RΔPSE','TΔPSE'});
            xlim([0.5 2.5]);
            if g == 1
            xlabel('RFsimilarity = 0~33 (different)');
            elseif g == 2
                xlabel('RFsimilarity = 33~66 (medium)');
            else
                xlabel('RFsimilarity = 66~100 (similar)');
            end
            ylabel('Normalized PSE shift');
            
            hold on
            plot([0.5 2.5],[0 0],'r-');
            
            % mean+sem
            meanPSE = mean([PSE_shift(group_this{g},2),PSE_shift(group_this{g},1)]);
            semPSE = std([PSE_shift(group_this{g},2),PSE_shift(group_this{g},1)]) / sqrt(length(PSE_shift(group_this{g},2)));
            % plot
            errorbar([1 2],meanPSE,semPSE,'r.-','markersize',10);
            
            text(0.8,15,sprintf('n = %g',num_this));
            
            % sign test for T and R PSE shift > 0 ？
            R_p = signtest(PSE_shift(group_this{g},2));
            T_p = signtest(PSE_shift(group_this{g},1));
            text(0.8,14,'Sign test with 0');
            text(1,max(PSE_shift(group_this{g},2))+1,sprintf('p = %.3g',R_p),'horizontalAlignment','center');
            text(2,max(PSE_shift(group_this{g},1))+1,sprintf('p = %.3g',T_p),'horizontalAlignment','center');
            ylim([-8 15]);
        end
        %         linkaxes([ax(1),ax(2)],'xy'); % 同步多个坐标区的范围
        %         linkaxes([ax(1),ax(3)],'xy'); % 同步多个坐标区的范围
        % 第三个图有个较大的异常值，保留，人为设置y轴
        suptitle(SetTitle(stimtask_type,monkey_included_for_analysis))
        SetFigure(15)
        figN = figN + 1;
        
        %%% 3 groups, scatter plot
        set(figure(figN),'name','T/R PSE shift and grouped Prefer OpticFlow Similarity','pos',[200,150,800,600]); clf
        color_value = nan(size(RFsimilarity));
        color_temp = [0 50 100];
        for g = 1:3
            color_value(group_this{g}) = color_temp(g);
        end
        
        scatter(PSE_shift(do_both & not_outdot,1),PSE_shift(do_both & not_outdot,2),70,color_value(do_both & not_outdot,:),'filled');
        color_this = [0 0 1; 0 1 0; 1 0 0];
        color_map_RGB = colormap_set([0 0 1; 0 1 0; 1 0 0]);
        colormap(color_map_RGB);
        xlabel('Normalized Translation PSE shift');
        ylabel('Normalized Rotation PSE shift');
        cb = colorbar('Ticks',[0 100],'Ticklabels',{'Difference','Same'});
        cb.Label.String = 'Preferred OpticFlow Similarity';
        ylim([-7.5 10]);
        xlim([-7.5 25]);
        xticks([-5:5:25]);
        yticks([-5:5:10]);
%         axis tight
%         SameScale(1);

        xlimit = xlim;
        ylimit = ylim;
        hold on
        plot([0 0],ylimit,'k:');
        plot(xlimit,[0 0],'k:');
        % text(xlimit(1)+(xlimit(2)-xlimit(1))/20,ylimit(2)-(ylimit(2)-ylimit(1))/2,sprintf('N= %g',site_num),'fontsize',15);
        text(xlimit(1)+(xlimit(2)-xlimit(1))/20,ylimit(2)-(ylimit(2)-ylimit(1))/20,sprintf('N= %g',site_num_without_outlier),'fontsize',15);
        
        for g = 1:3
            hold on
            TPSEshift = PSE_shift(do_both & not_outdot & group_this{g},1);
            RPSEshift = PSE_shift(do_both & not_outdot & group_this{g},2);
            if g == 1
                error_ellipse(cov(TPSEshift,RPSEshift),'conf',0.95,'mu',[mean(TPSEshift) mean(RPSEshift)],'style','b'); % confidence ellipses for data,0-33
            elseif g == 2
                error_ellipse(cov(TPSEshift,RPSEshift),'conf',0.95,'mu',[mean(TPSEshift) mean(RPSEshift)],'style','g'); % confidence ellipses for data,33-66
            elseif g == 3
                error_ellipse(cov(TPSEshift,RPSEshift),'conf',0.95,'mu',[mean(TPSEshift) mean(RPSEshift)],'style','r'); % confidence ellipses for data,66-100
            end
        end
        title('95% confidence ellipses');
        SetFigure(15);
        figN = figN+1;
        
        % difference between groups:
        % ktest
        [~, pValue12] = kstest_2s_2d(PSE_shift(group_this{1},:), PSE_shift(group_this{2},:)); % different vs medium
        [~, pValue23] = kstest_2s_2d(PSE_shift(group_this{2},:), PSE_shift(group_this{3},:)); % medium vs similar
        [~, pValue13] = kstest_2s_2d(PSE_shift(group_this{1},:), PSE_shift(group_this{3},:)); % different vs similar
        str = ['2d KStest' newline sprintf('different vs medium, p = %.2g',pValue12) ...
            newline sprintf('medium vs similar, p = %.2g',pValue23) ...
            newline sprintf('different vs similar, p = %.2g',pValue13)];
        text(xlimit(1)+(xlimit(2)-xlimit(1))/20,6, str);
        
        
        %         % Permanova
        %         % interaction: ANCOVA
        %         X = [PSE_shift(group_this{1},1);PSE_shift(group_this{2},1)];
        %         Y = [PSE_shift(group_this{1},2);PSE_shift(group_this{2},2)];
        %         G = [ones(size(PSE_shift(group_this{1},1)));2*ones(size(PSE_shift(group_this{2},1)))];
        %         % [h,a,c,s] = aoctool(x,y, group,alpha,xname,yname,gname,displayopt,model)
        %         [h,atab,ctab,stats] = aoctool(X,Y,G,0.05,'d prime','PSE shift','Stimulus type','off');
        %         p_TR = atab{3,6}; % across RFsimilarity, T and R PSE shift (covariate effect)
        %         p_T_similar = atab{4,6}; % interaction between T PSE shift and RFsimilarity
        %         p_R_similar = atab{2,6}; % interaction between R PSE shift and RFsimilarity
        
        
        %%% 2. continue distribution, scatter plot for T and R
        set(figure(figN),'name','T/R PSE shift and Prefer OpticFlow Similarity','pos',[200,200,800,600]); clf
        color_value = RFsimilarity;
        scatter(PSE_shift(do_both & not_outdot,1),PSE_shift(do_both & not_outdot,2),70,color_value(do_both & not_outdot,:),'filled');
        color_map_RGB = colormap_set([0 0 0; 1 0 0]);
        colormap(color_map_RGB);
        
        xlabel('Normalized Translation PSE shift');
        ylabel('Normalized Rotation PSE shift');
        axis square
        cb = colorbar('Ticks',[0 100],'Ticklabels',{'Difference','Same'});
        cb.Label.String = 'Preferred OpticFlow Similarity';
        %         axis tight
        %         SameScale(1);
        xlimit = xlim;
        ylimit = ylim;
        hold on
        plot([0 0],ylimit,'k:');
        plot(xlimit,[0 0],'k:');
        % text(xlimit(1)+(xlimit(2)-xlimit(1))/20,ylimit(2)-(ylimit(2)-ylimit(1))/2,sprintf('N= %g',site_num),'fontsize',15);
        text(xlimit(1)+(xlimit(2)-xlimit(1))/20,ylimit(2)-(ylimit(2)-ylimit(1))/20,sprintf('N= %g',site_num_without_outlier),'fontsize',15);
        
        SetFigure(15);
        figN = figN+1;
        
        %%% 2D histogram
        XX = PSE_shift(do_both & not_outdot,1);
        YY = PSE_shift(do_both & not_outdot,2);
        CC = RFsimilarity(do_both & not_outdot);
        binnum = 21;
        xedges = linspace(xlimit(1),xlimit(2),binnum);
        yedges = linspace(ylimit(1),ylimit(2),binnum);
        
        meanSimilarity = nan(binnum-1,binnum-1);
        for i = 1:length(yedges)-1 % 第i行
            for j = 1:length(xedges)-1 % 第j列
                select = logical(YY>yedges(i) & YY<=yedges(i+1) & XX>xedges(j) & XX<=xedges(j+1));
                if sum(select) > 0
                    meanSimilarity(i,j) = mean(CC(select));
                end
            end
        end
        set(figure(figN),'name','T and R PSE shift and Prefer OpticFlow Similarity','pos',[800,250,800,600]); clf
        %         xmid = xedges(1:end-1) + (xedges(2) - xedges(1))/2;
        %         ymid = yedges(1:end-1) + (yedges(2) -  yedges(1))/2;
        %         imagesc(xmid, ymid, meanSimilarity);
        %         set(gca,'YDir','normal'); % flip up-down, make sure smaller y bin in the bottom, smaller x bin in the left
        %         axis square
        %         colorbar
        %         SetFigure(15);
        %         figN = figN+1;
        
        % step window
        win = [2 1];
        step = [1 0.5];
        [xmid_slide, ymid_slide, z, z_sem] = get_slideWin2D(XX,YY,CC,win,step,[-10 10],[-5 5]);
        [aa,bb] = size(z);
        z_plot = z;
        z_plot(aa+1,:) = z_plot(aa,:); % 伪造最后一行
        z_plot(:,bb+1) = z_plot(:,bb); % 伪造最后一列
        
        x_plot = unique(xmid_slide) - step(1)/2; % 调整坐标为每个格子的左下角网格交点
        x_plot(end+1) = x_plot(end) + step(1);
        
        y_plot = unique(ymid_slide) - step(2)/2; % 调整坐标为每个格子的左下角网格交点
        y_plot(end+1) = y_plot(end) + step(2);
        [x_plot, y_plot] = meshgrid(x_plot,y_plot);
        
        s = pcolor(x_plot, y_plot,z_plot); % 每个格子中心对应原来xmid，ymid
        axis square
        xlabel('Normalized Translation PSE shift');
        ylabel('Normalized Rotation PSE shift');
        %         color_map_RGB = colormap_set([0 0 0; 1 0 0]);
        %         colormap(color_map_RGB);
        
        s.FaceColor = 'interp'; % 跨格子插值
        cb = colorbar('Ticks',[0 100],'Ticklabels',{'Difference','Same'});
        cb.Label.String = 'Preferred OpticFlow Similarity';
        cb.Limits = [0 100];
        SetFigure(15);
        figN = figN+1;
        
        %         % save figer， 然后用ai打开（不能直接拖）
        %         xticks([])
        %         yticks([])
        %         axis off
        %         Fig = figure(figN);
        %         cd
        %         saveas(Fig,'309','pdf')
        
        %%% 3. continue distribution, separate for T and R
        set(figure(figN),'name','T or R PSE shift and Prefer OpticFlow Similarity','pos',[800,250,1200,600]); clf
        for hr = 1:2
            subplot(3,2,[hr hr+2])
            hold on
            plot(RFsimilarity(not_outdot),PSE_shift(not_outdot,hr),'ko');
            plot(RFsimilarity(not_outdot & S_unit(:,hr)),PSE_shift((S_unit(:,hr) & not_outdot),hr),'o','color','k','markerfacecolor',c{hr});
            axis tight
            ylimit = ylim;
            ylimit = ylimit.*1.1;
            ylimit = [-max(abs(ylimit)) max(abs(ylimit))];
            ylim(ylimit);
            xlimit = [0 100];
            plot(xlimit,[0 0],'g-');
            xlim(xlimit);
            if is_original_bias == 3
                ylabel('Normalized Induced PSE shift (°)');
            elseif is_original_bias == 0
                ylabel('Induced PSE shift (°)');
            end
            xlabel('Similarity (%)');
            
            % step mean
            ax(hr) = subplot(3,2,hr+4);
            win = 20;
            step = 5;
            [xmid, y, y_sem, y_step_data] = get_slideWin(RFsimilarity(not_outdot),PSE_shift(not_outdot,hr),win,step,[-10 110]);
            shadedErrorBar(xmid,y,y_sem,'k');  
%             errorbar(xmid,y,y_sem);
            axis tight
            xlim([0 100]);
            ylimit = ylim;
            ylimit = [-max(abs(ylimit))*1.1 max(abs(ylimit))*1.1];
            ylim(ylimit);
            hold on
            plot([0 100],[0 0],'g-');
            
            % t test for mean PSE
            for i = 1:length(xmid)
                [~,p(i)] = ttest(y_step_data{i});
                
                % plot significance
                if p(i)<0.05
                    plot([xmid(i)-step/2, xmid(i)+step/2],[ylimit(2) ylimit(2)],'r-','linewidth',2);
                end
            end
        end
        suptitle(SetTitle(stimtask_type,monkey_included_for_analysis));
        SetFigure(15);
        figN = figN+1;
    end

    function f3p17(debug) % 'PSE shift and Preferred OpticFlow Fine vs Coarse'
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        effect_type = 1; % 1: PSE, 2: threshold
        select_stars_type = [0 1]; % 0: dot stimulus; 1: triangle stimulus
        is_separate_stimAmp = 0;
        is_original_bias = 3; % normlize PSE shift (divided by the threshold under non-stim trials
        
        % ---------------------- get data ----------------------------
        stimtask_type_fine{1} = stimtask_type{1};
        stimtask_type_fine{2} = [];
        stimtask_type_coarse{1} = [];
        stimtask_type_coarse{2} = stimtask_type{2};
        
        [PSE_shift_temp1, S_unit_temp1] = get_stim_data_uA(stimtask_type_fine,methods_of_select{1,1},is_separate_stimAmp,is_original_bias,select_stars_type);
        [PSE_shift_temp2, S_unit_temp2] = get_stim_data_uA(stimtask_type_coarse,methods_of_select{1,1},is_separate_stimAmp,is_original_bias,select_stars_type);
        PSE_shift{1} = PSE_shift_temp1{effect_type}(:,1:2); % only save h and r, not mt
        PSE_shift{2} = PSE_shift_temp2{effect_type}(:,1:2); % only save h and r, not mt
        S_unit{1} = S_unit_temp1{effect_type}(:,1:2)*1;
        S_unit{2} = S_unit_temp2{effect_type}(:,1:2)*1;
        
        block_size = length(loadM_data);
        for fc = 1:2
            if fc == 1
                if sum(stimtask_type_fine{fc}==3)>0 % 2 target fine task
                    block_loc_this = find(stimtask_type_fine{fc}==3);
                    PSE_shift{fc}(block_size*(block_loc_this-2)+1:block_size*(block_loc_this-1),2) = PSE_shift{fc}(block_size*(block_loc_this-1)+1:block_size*block_loc_this,2); % 将rotation补充到heading那里
                    PSE_shift{fc}(block_size*(block_loc_this-1)+1:block_size*block_loc_this,:) = ones(block_size,2)*9999; % 将原来rotation位置变为999，后面去掉
                    S_unit{fc}(block_size*(block_loc_this-2)+1:block_size*(block_loc_this-1),2) = S_unit{fc}(block_size*(block_loc_this-1)+1:block_size*block_loc_this,2); % 将rotation补充到heading那里
                end
            elseif fc == 2
                if sum(stimtask_type_coarse{fc}==6)>0 % 2 target coarse task
                    block_loc_this = find(stimtask_type_coarse{fc}==6);
                    PSE_shift{fc}(block_size*(block_loc_this-2)+1:block_size*(block_loc_this-1),2) = PSE_shift{fc}(block_size*(block_loc_this-1)+1:block_size*block_loc_this,2); % 将rotation补充到heading那里
                    PSE_shift{fc}(block_size*(block_loc_this-1)+1:block_size*block_loc_this,:) = ones(block_size,2)*9999; % 将原来rotation位置变为999，后面去掉
                    S_unit{fc}(block_size*(block_loc_this-2)+1:block_size*(block_loc_this-1),2) = S_unit{fc}(block_size*(block_loc_this-1)+1:block_size*block_loc_this,2); % 将rotation补充到heading那里
                end
            end
        end
        
        % prefer optic flow similarity
        try
            RFsimilarity_temp = [group_result.RFsimilarity]';
        catch
            run_this = strcmpi({handles.function_handles_spi{2,2}{:,1}}','Save Prefer OpticFlow similarity (local planar motion)');
            feval(handles.function_handles_spi{2,2}{run_this,2},false); % first run "Spiral Tuning" - "Spiral Tuning" - "Save Prefer OpticFlow (local planar motion)"
            RFsimilarity_temp = [group_result.RFsimilarity]';
        end
        RFsimilarity{1} = repmat(RFsimilarity_temp,[length(stimtask_type_fine{1}),1]);
        RFsimilarity{2} = repmat(RFsimilarity_temp,[length(stimtask_type_coarse{2}),1]);
        
        % 去掉999
        for fc = 1:2
            RFsimilarity{fc}(PSE_shift{fc}(:,1)==9999,:) = [];
            S_unit{fc}(PSE_shift{fc}(:,1)==9999,:) = [];
            PSE_shift{fc}(PSE_shift{fc}(:,1)==9999,:) = [];
            
            % group 3 similarity groups
            group_range = linspace(0,100,4);
            for g = 1:3
                group_this{fc,g} = logical(RFsimilarity{fc} > group_range(g) & RFsimilarity{fc} <= group_range(g+1));
            end
        end

        % ------------------------ Plot ------------------------------
        set(figure(figN),'name','T/R PSE shift and grouped Prefer OpticFlow Similarity (Fine vs Coarse)','pos',[100,100,2000,500]); clf
        for g = 1:3
            hold on
            for fc = 1:2
                
                % find which do heading and rotation simultaneously
                not_do_both = any(isnan(PSE_shift{fc}),2);
                select = logical(group_this{fc,g} & ~not_do_both);
                
                plot([(g-1)*4+(fc-1)*2+1 (g-1)*4+(fc-1)*2+2],[PSE_shift{fc}(select,2),PSE_shift{fc}(select,1)],'ko-'); % R vs T
                
                num_this = size(PSE_shift{fc}(select,:),1);
                text((g-1)*4+(fc-1)*2+1.5,12,sprintf('n = %g',num_this),'horizontalAlignment','center');
                
                % mean+sem
                medianPSE = median(PSE_shift{fc}(select,:));
                %                 semPSE = std(PSE_shift{fc}(select,:)) ./ sqrt(num_this);
                medianPSE = fliplr(medianPSE); % rotation first
                %                 semPSE = fliplr(semPSE); % rotation first
                
                %                 % plot
                %                 errorbar([(g-1)*4+(fc-1)*2+1 (g-1)*4+(fc-1)*2+2],medianPSE,semPSE,'r.-','markersize',10);
                plot([(g-1)*4+(fc-1)*2+1 (g-1)*4+(fc-1)*2+2],medianPSE,'ro-','markersize',10);
                
                % ***** sign test ****** for T and R PSE shift > 0 ？
                R_p = signtest(PSE_shift{fc}(select,2));
                T_p = signtest(PSE_shift{fc}(select,1));
                text((g-1)*4+(fc-1)*2+1,max(PSE_shift{fc}(select,2))+1,sprintf('p = %.3g',R_p),'horizontalAlignment','center');
                text((g-1)*4+(fc-1)*2+2,max(PSE_shift{fc}(select,1))+1,sprintf('p = %.3g',T_p),'horizontalAlignment','center');
            end
            
            set(gca,'xtick',[1:12],'xticklabel',repmat({'RΔPSE','TΔPSE'},1,6));
            xlim([0.5 12.5]);
            ylabel('Normalized PSE shift');
            
            if g == 1
                text((g-1)*4+2.5,-10,'RFsimilarity = 0~33 (different)','horizontalAlignment','center');
            elseif g == 2
                text((g-1)*4+2.5,-10,'RFsimilarity = 33~66 (medium)','horizontalAlignment','center');
            else
                text((g-1)*4+2.5,-10,'RFsimilarity = 66~100 (similar)','horizontalAlignment','center');
            end
            
            plot([0.5 12.5],[0 0],'r-');
            ylim([-8 15]);
            
            % difference between fine and coarse:
            % kStest
            [~, pValue] = kstest_2s_2d(PSE_shift{1}(select,:), PSE_shift{2}(select,:)); % different vs medium
            text((g-1)*4+2.5,14,sprintf('2D KS test, p=% 2.3g',pValue),'horizontalAlignment','center');
        end
        title(SetTitle(stimtask_type,monkey_included_for_analysis));
        SetFigure(15)
        figN = figN + 1;
    end


    function f3p15(debug) % 'PSE shift for fovea-containing RF sites'
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        effect_type = 1; % 1: PSE, 2: threshold
        select_stars_type = [0 1]; % 0: dot stimulus; 1: triangle stimulus
        is_separate_stimAmp = 0;
        is_original_bias = 3; % normlize PSE shift (divided by the threshold under non-stim trials
        
        % make sure select fine and coarse task simultaneously
        both_fc = sum(cellfun(@(x) isempty(x),stimtask_type));
        if both_fc==1
            disp('********** Please choose BOTH fine and coarse task ! **********');
            return
        end
        
        % ---------------------- get data ----------------------------
        % seperate fine and coarse task
        stimtask_type_fine{1} = stimtask_type{1};
        stimtask_type_fine{2} = [];
        stimtask_type_coarse{1} = [];
        stimtask_type_coarse{2} = stimtask_type{2};
        
        [PSE_shift_temp1, S_unit_temp1, NS_unit_temp1] = get_stim_data_uA(stimtask_type_fine,methods_of_select{1,1},is_separate_stimAmp,is_original_bias,select_stars_type);
        [PSE_shift_temp2, S_unit_temp2, NS_unit_temp2] = get_stim_data_uA(stimtask_type_coarse,methods_of_select{1,1},is_separate_stimAmp,is_original_bias,select_stars_type);
        PSE{1} = PSE_shift_temp1{effect_type}(:,1:2); % only save h and r, not mt
        PSE{2} = PSE_shift_temp2{effect_type}(:,1:2);
        S_unit{1} = S_unit_temp1{effect_type}(:,1:2);
        S_unit{2} = S_unit_temp2{effect_type}(:,1:2);
        NS_unit{1} = NS_unit_temp1{effect_type}(:,1:2);
        NS_unit{2} = NS_unit_temp2{effect_type}(:,1:2);
        
        % 在two-target task中，因为函数get_stim_data_uA会叠加heading和rotation task，而不是将其合并，所以需要重新调整PSE等排布
        stimtask_type_this = cell2mat(stimtask_type_fine);
        block_size = length(loadM_data);
        if sum(stimtask_type_this==3)>0 % 2 target fine task
            block_loc_this = find(stimtask_type_this==3);
            PSE{1}(block_size*(block_loc_this-2)+1:block_size*(block_loc_this-1),2) = PSE{1}(block_size*(block_loc_this-1)+1:block_size*block_loc_this,2); % 将rotation补充到heading那里
            PSE{1}(block_size*(block_loc_this-1)+1:block_size*block_loc_this,:) = ones(block_size,2)*9999; % 将原来rotation位置变为999，后面去掉
            S_unit{1}(block_size*(block_loc_this-2)+1:block_size*(block_loc_this-1),2) = S_unit{1}(block_size*(block_loc_this-1)+1:block_size*block_loc_this,2); % 将rotation补充到heading那里
            NS_unit{1}(block_size*(block_loc_this-2)+1:block_size*(block_loc_this-1),2) = NS_unit{1}(block_size*(block_loc_this-1)+1:block_size*block_loc_this,2); % 将rotation补充到heading那里
        end
        
        stimtask_type_this = cell2mat(stimtask_type_coarse);
        if sum(stimtask_type_this==6)>0 % 2 target coarse task
            block_loc_this = find(stimtask_type_this==6);
            PSE{2}(block_size*(block_loc_this-2)+1:block_size*(block_loc_this-1),2) = PSE{2}(block_size*(block_loc_this-1)+1:block_size*block_loc_this,2); % 将rotation补充到heading那里
            PSE{2}(block_size*(block_loc_this-1)+1:block_size*block_loc_this,:) = ones(block_size,2)*9999; % 原来的填充9999
            S_unit{2}(block_size*(block_loc_this-2)+1:block_size*(block_loc_this-1),2) = S_unit{2}(block_size*(block_loc_this-1)+1:block_size*block_loc_this,2); % 将rotation补充到heading那里
            NS_unit{2}(block_size*(block_loc_this-2)+1:block_size*(block_loc_this-1),2) = NS_unit{2}(block_size*(block_loc_this-1)+1:block_size*block_loc_this,2); % 将rotation补充到heading那里
        end
        
        for fc = 1:2
            select_PSE9999{fc} = PSE{fc}(:,1)==9999;
            S_unit{fc}(select_PSE9999{fc},:) = [];
            NS_unit{fc}(select_PSE9999{fc},:) = [];
            PSE{fc}(select_PSE9999{fc},:) = [];
        end
        
        % RF size: 转化为圆的半径比较
        % RF eccentricity: rf中心距离（0，0）距离
        rf =  cell2mat({group_result.RF}'); % X,Y,W,H  (X,Y为RF中心点坐标)
        rf_area = rf(:,3) .* rf(:,4);
        rf_dia = 2.* sqrt(rf_area / pi); % 等效直径
        rf_ecc = sqrt(rf(:,1).^2 + rf(:,2).^2); % RF中心点距离原点的距离
        
        % 1. Square RF, included fovea
        criti_value = 0; % in visual degree
        %         criti_value_range = [0:0.1:40];
        
        includefovea_temp = nan(length(rf),1);
        for d = 1:length(rf)
            if ~isnan(rf(d,1))
                if (rf(d,1)-rf(d,3)/2) < -criti_value && (rf(d,1)+rf(d,3)/2) > criti_value && (rf(d,2)-rf(d,4)/2) < -criti_value && (rf(d,2)+rf(d,4)/2) > criti_value
                    includefovea_temp(d) = 1;
                else
                    includefovea_temp(d) = 0;
                end
            end
        end
                
        % 2. Circle RF, included fovea
        %         for i = 1:length(criti_value_range)
        
        includefovea{1} = repmat(includefovea_temp,[length(stimtask_type_fine{1}),1]); % copy for fine and coarse
        includefovea{2} = repmat(includefovea_temp,[length(stimtask_type_coarse{2}),1]); % copy for fine and coarse
        
        
        %             criti_value = criti_value_range(i);
        %             includefovea2_temp = logical(rf_dia/2 > (rf_ecc + criti_value));
        includefovea2_temp = logical((rf_dia/2)./rf_ecc > criti_value);
        includefovea2{1} = repmat(includefovea2_temp,[length(stimtask_type_fine{1}),1]); % copy for fine and coarse
        includefovea2{2} = repmat(includefovea2_temp,[length(stimtask_type_coarse{2}),1]); % copy for fine and coarse
        
        % 3. select RF size
        size_criti_value = 0; % do not select
        large_RF_temp = logical(rf_dia >= size_criti_value);
        large_RF{1} = repmat(large_RF_temp,[length(stimtask_type_fine{1}),1]); % copy for fine and coarse
        large_RF{2} = repmat(large_RF_temp,[length(stimtask_type_coarse{2}),1]); % copy for fine and coarse
        
        % prefer optic flow similarity
        try
            RFsimilarity = [group_result.RFsimilarity]';
        catch
            run_this = strcmpi({handles.function_handles_spi{2,2}{:,1}}','Save Prefer OpticFlow similarity (local planar motion)');
            feval(handles.function_handles_spi{2,2}{run_this,2},false); % first run "Spiral Tuning" - "Spiral Tuning" - "Save Prefer OpticFlow (local planar motion)"
            RFsimilarity = [group_result.RFsimilarity]';
        end
        
        % find included fovea sites' similarity
        RFsimilarity_temp{1} = repmat(RFsimilarity,[length(stimtask_type_fine{1}),1]); % copy for fine and coarse
        RFsimilarity_temp{2} = repmat(RFsimilarity,[length(stimtask_type_coarse{2}),1]); % copy for fine and coarse
        
        % 去掉999
            for fc = 1:2
                includefovea{fc}(select_PSE9999{fc},:) = [];
                includefovea2{fc}(select_PSE9999{fc},:) = [];
                RFsimilarity_temp{fc}(select_PSE9999{fc},:) = [];
                large_RF{fc}(select_PSE9999{fc},:) = [];
                
                select = logical(includefovea2{fc}==1 & large_RF{fc});
                p_signtest = signtest(PSE{fc}(select,2)); % 自动省略nan
                [~,p_ttest] = ttest(PSE{fc}(select,2)); % 自动省略nan
%                 
            end
%         end
        
%         figure
%         for fc = 1:2
%             subplot(2,2,fc)
%             plot(criti_value_range,p_signtest,'k-');
%             hold on
%             plot([criti_value_range(1) criti_value_range(end)],[0.05 0.05],'k:');
%             
%             subplot(2,2,fc+2)
%             plot(criti_value_range,p_ttest,'k-');
%             hold on
%             plot([criti_value_range(1) criti_value_range(end)],[0.05 0.05],'k:');
%         end
% %         
        
        % ------------------------ Plot ------------------------------
        set(figure(figN),'name','R PSE shift for fovea-containing RF sites','pos',[100,100,1300,700]); clf
        % 只画Rotation
        for fc = 1:2 % fine and coarse
            for i = 1:2 % include or not include fovea
                subplot(2,4,fc+(i-1)*2)
                if i == 1
                    select = logical(includefovea{fc}==1 & large_RF{fc});
                elseif i == 2
                    select = logical(includefovea{fc}==0 & large_RF{fc});
                end
                
                % num
                N_this = sum(~isnan(PSE{fc}(select,2)));
                
                [binedge, binmid, xtick_this, bin_size] = properHistBin(PSE{fc}(select,2),'small');
                histS = histcounts(PSE{fc}(S_unit{fc}(:,2) & select,2),binedge);
                histNS = histcounts(PSE{fc}(NS_unit{fc}(:,2) & select,2),binedge);
                hbars = bar(binmid,[histS' histNS'],1,'stacked','LineWidth',1.5);
                set(hbars,'EdgeColor','none','FaceColor',c{2});
                set(hbars(2),'EdgeColor','none','FaceColor','y'); % for segement
                xlim([xtick_this(1) xtick_this(end)]);
                xticks(xtick_this);
                ymax = max([histS+histNS]);
                [ylimit_max, ytick_this] = properYtick(ymax);
                yticks(ytick_this);
                ylim([0 ylimit_max])
                
                % median
                PSE_median = nanmedian(PSE{fc}(select,2));
                %ttest for significant
                p_signtest = signtest(PSE{fc}(select,2)); % 自动省略nan
                [~,p_ttest] = ttest(PSE{fc}(select,2)); % 自动省略nan
                
                % add Asterisks marker
                if p_signtest<0.001
                    asterisk = '***';
                elseif p_signtest<0.01
                    asterisk = '**';
                elseif p_signtest<0.05
                    asterisk = '*';
                else
                    asterisk = 'n.s.';
                end
                
                %annotation
                hold on
                plot(PSE_median,ylimit_max,'vk','linewidth',1,'markerfacecolor',c{2});  %三角标记
                % plot median
                text(xtick_this(1),ylimit_max-ylimit_max/10,sprintf('%2.2g (p(signtest) = %2.2g)',PSE_median, p_signtest),'fontsize',10,'color',c{2});
                text(xtick_this(1),ylimit_max-ylimit_max/6,sprintf('p(t test) = %2.2g',p_ttest),'fontsize',10,'color',c{2});
                % add Asterisks marker
                text(PSE_median,ylimit_max+ylimit_max/20,asterisk,'fontsize',15,'color',c{2},'HorizontalAlign','center');
                ylabel('Cases');
                text(xtick_this(1),ylimit_max,sprintf('n = %g',N_this));
                
                if i == 1 && fc == 1
                    title('Included fovea fine R');
                elseif i == 1 && fc == 2
                    title('Included fovea coarse R');
                elseif i == 2 && fc == 1
                    title('Not included fovea fine R');
                elseif i == 2 && fc == 2
                    title('Not included fovea coarse R');
                end
                xlabel('Normalized PSE');
                axis square
                
                % save this sites similarity
                 % num
                N_this = sum(~isnan(PSE{fc}(select,2)));
                
                select_this = logical(~isnan(PSE{fc}(select,2)));
                RFsimilarity_this = RFsimilarity_temp{fc}(select_this);
                RFsimilarity_this = RFsimilarity_this(~isnan(RFsimilarity_this));
                N_RFsimilarity = length(RFsimilarity_this);
                
                subplot(2,4,fc+(i-1)*2+4)
                histogram(RFsimilarity_this,[0:10:100]);
                ylimit = ylim;
                text(5,ylimit(2),sprintf('N = %.3g',N_RFsimilarity));
            end
        end
        suptitle(SetTitle(stimtask_type,monkey_included_for_analysis))
        SetFigure(10);
        figN = figN+1;
        
        set(figure(figN),'name','R PSE shift for fovea-containing RF sites, fine vs coarse','pos',[200,100,700,700]); clf
        % do both
        PSE_fine = PSE{1}(includefovea{1}==1 & large_RF{1},2); % only r PSE, contain fovea
        PSE_coarse = PSE{2}(includefovea{2}==1 & large_RF{2},2); % only r PSE, contain fovea
        
        do_both = ~any(isnan([PSE_fine,PSE_coarse]),2);
        plot(PSE_fine(do_both),PSE_coarse(do_both),'ko');
        axis tight
        SameScale(1);
        [r,p] = plot_corr_line(PSE_fine(do_both),PSE_coarse(do_both),'MethodOfCorr','Pearson','FittingMethod',1,'LineWidth',2);
        text(-4,6,sprintf('Pearson corr. r=%.3g, p=%.3g',r,p));
        xlabel('Normalized fine roll PSE');
        ylabel('Normalized coarse roll PSE');
        text(-4,5,sprintf('N = %g',sum(do_both)));
        
        % sign test for fine or coarse
        p_fine = signtest(PSE_fine(do_both));
        p_coarse = signtest(PSE_coarse(do_both));
        
        title('Fovea-containing RF sits, fine-roll vs coarse-roll');
        SetFigure(14);
        figN = figN+1;

    end


    function f3p3(debug) % 'Comparison of Fine and Coarse task in PSE shift'
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        is_separate_stimAmp = 0;
        is_original_bias = 3; % normalize PSE shift, divided by the threshold under non-stim trials
        select_stars_type = [0 1];
        effect_type = 1;
        
        % make sure select fine and coarse task simultaneously
        both_fc = sum(cellfun(@(x) isempty(x),stimtask_type));
        if both_fc==1
            disp('********** Please choose BOTH fine and coarse task ! **********');
            return
        end
        
        % ---------------------- get data ----------------------------
        % seperate fine and coarse task
        stimtask_type_fine{1} = stimtask_type{1};
        stimtask_type_fine{2} = [];
        stimtask_type_coarse{1} = [];
        stimtask_type_coarse{2} = stimtask_type{2};
        
        [diff_thistask_fine, S_unit_fine] = get_stim_data_uA(stimtask_type_fine,methods_of_select{1,1},is_separate_stimAmp,is_original_bias,select_stars_type);
        [diff_thistask_coarse, S_unit_coarse] = get_stim_data_uA(stimtask_type_coarse,methods_of_select{1,1},is_separate_stimAmp,is_original_bias,select_stars_type);
        fine_PSE = diff_thistask_fine{effect_type};
        coarse_PSE = diff_thistask_coarse{effect_type};
        fine_sig = S_unit_fine{effect_type};
        coarse_sig = S_unit_coarse{effect_type};
        
        % 选取共同做了fine and coarse task
        do_fine = ~isnan(fine_PSE);
        do_coarse = ~isnan(coarse_PSE);
        for i = 1:3 % h r mt
            do_both(:,i) = logical(do_fine(:,i) & do_coarse(:,i));
        end
        site_num = sum(do_both);
        
        % outlier
        for i = 1:2 % h r
            f_temp = fine_PSE(:,i);
            c_temp = coarse_PSE(:,i);
            f_reasnonal_range = [nanmean(f_temp)-3*nanstd(f_temp) nanmean(f_temp)+3*nanstd(f_temp)];
            c_reasnonal_range = [nanmean(c_temp)-3*nanstd(c_temp) nanmean(c_temp)+3*nanstd(c_temp)];
            f_outlier(:,i) = logical(f_temp<f_reasnonal_range(1) | f_temp>f_reasnonal_range(2));
            c_outlier(:,i) = logical(c_temp<c_reasnonal_range(1) | c_temp>c_reasnonal_range(2));
            
            outlier(:,i) = logical(f_outlier(:,i) | c_outlier(:,i));
        end
        outlier = zeros(size(outlier));

        % ------------------------ Plot ------------------------------
        set(figure(figN),'name','Comparison of Fine and Coarse task','pos',[100,100,1500,700]); clf
        for i = 1:2 % h r
            subplot(1,2,i)
            if site_num(i)~=0
                %                 plot(fine_PSE(do_both(:,i),i),coarse_PSE(do_both(:,i),i),'o','color',c{i});
                hold on
                % actual remove num
                remove_num = sum(do_both(:,i) & outlier(:,i));
                
                % remove outlier
                do_both(:,i) = logical(do_both(:,i) & ~outlier(:,i));
                % acutal site num
                site_num = sum(do_both);
                
                select_fine_sig = logical(do_both(:,i) & fine_sig(:,i) & ~coarse_sig(:,i));
                select_coarse_sig = logical(do_both(:,i) & ~fine_sig(:,i) & coarse_sig(:,i));
                select_both_sig = logical(do_both(:,i) & fine_sig(:,i) & coarse_sig(:,i));
                select_not_sig = logical(do_both(:,i) & ~fine_sig(:,i) & ~coarse_sig(:,i));
                
                plot(fine_PSE(select_not_sig,i),coarse_PSE(select_not_sig,i),'o','color',c{99},'markersize',7);
                plot(fine_PSE(select_fine_sig,i),coarse_PSE(select_fine_sig,i),'+','color','k','markersize',10,'linewidth',1.5);
                plot(fine_PSE(select_coarse_sig,i),coarse_PSE(select_coarse_sig,i),'x','color','k','markersize',10,'linewidth',1.5);
                plot(fine_PSE(select_both_sig,i),coarse_PSE(select_both_sig,i),'o','color','k','markerfacecolor',c{i},'markersize',7,'linewidth',1);
                [r,ppp,h1] = plot_corr_line(fine_PSE(do_both(:,i),i),coarse_PSE(do_both(:,i),i),'MethodOfCorr','Pearson','FittingMethod',1,'LineWidth',2); % 有outlier，用最小二乘法fit
                
                axis tight
                SameScale
%                 h6 = SameScale(1); % same scale, diagonal line, keep tick not change
%                 set(get(get(h6,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                
                xlimit = xlim;
                ylimit = ylim;
                text(xlimit(1)+(xlimit(2)-xlimit(1))/20,ylimit(2)-(ylimit(2)-ylimit(1))/20,sprintf('Pearson correlation: r=%0.3g, p=%0.3g',r,ppp),'fontsize',15);
                text(xlimit(1)+(xlimit(2)-xlimit(1))/20,ylimit(2)-(ylimit(2)-ylimit(1))/8,sprintf('N = %g',site_num(i)),'fontsize',15);
                
                text(xlimit(1)+(xlimit(2)-xlimit(1))/20,ylimit(2)-(ylimit(2)-ylimit(1))/4,sprintf('Remoce outlier = %g',remove_num),'fontsize',15);
                
                xlabel('normalize Fine PSE shift');
                ylabel('normalize Coarse PSE shift'); 
            end
            if i == 1
                legend('Not sig','Fine sig','Coarse sig','Both sig','Fitting','location','best');
            end
            
            if i == 1
                title('Translatioin');
            elseif i ==2
                title('Rotation');
            end
                
        end
        SetFigure(15);
        suptitle(SetTitle(stimtask_type,monkey_included_for_analysis))
        figN = figN+1;
    end

    function f3p4(debug) % 'Fine task V.S. Coarse task (Threshold change)'
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        is_separate_stimAmp = 0;
        is_original_bias = 0; 
        select_stars_type = [0 1];
        effect_type = 2;
        
        % make sure select fine and coarse task simultaneously
        both_fc = sum(cellfun(@(x) isempty(x),stimtask_type));
        if both_fc==1
            disp('********** Please choose BOTH fine and coarse task ! **********');
            return
        end
        
        % ---------------------- get data ----------------------------
        % seperate fine and coarse task
        stimtask_type_fine{1} = stimtask_type{1};
        stimtask_type_fine{2} = [];
        stimtask_type_coarse{1} = [];
        stimtask_type_coarse{2} = stimtask_type{2};
        
        [diff_thistask_fine, S_unit_fine] = get_stim_data_uA(stimtask_type_fine,methods_of_select{1,1},is_separate_stimAmp,is_original_bias,select_stars_type);
        [diff_thistask_coarse, S_unit_coarse] = get_stim_data_uA(stimtask_type_coarse,methods_of_select{1,1},is_separate_stimAmp,is_original_bias,select_stars_type);
        fine_thr = diff_thistask_fine{effect_type};
        coarse_thr = diff_thistask_coarse{effect_type};
        fine_sig = S_unit_fine{effect_type};
        coarse_sig = S_unit_coarse{effect_type};
        
        % 选取共同做了fine and coarse task
        do_fine = ~isnan(fine_thr);
        do_coarse = ~isnan(coarse_thr);
        for i = 1:3 % h r mt
            do_both(:,i) = logical(do_fine(:,i) & do_coarse(:,i));
        end
        site_num = sum(do_both);
        
        % outlier
        for i = 1:2 % h r
            f_temp = fine_thr(:,i);
            c_temp = coarse_thr(:,i);
            f_reasnonal_range = [nanmean(f_temp)-3*nanstd(f_temp) nanmean(f_temp)+3*nanstd(f_temp)];
            c_reasnonal_range = [nanmean(c_temp)-3*nanstd(c_temp) nanmean(c_temp)+3*nanstd(c_temp)];
            f_outlier(:,i) = logical(f_temp<f_reasnonal_range(1) | f_temp>f_reasnonal_range(2));
            c_outlier(:,i) = logical(c_temp<c_reasnonal_range(1) | c_temp>c_reasnonal_range(2));
            
            outlier(:,i) = logical(f_outlier(:,i) | c_outlier(:,i));
            outlier(:,i) = false(size(outlier(:,i))); % not select
        end
        
        % ------------------------ Plot ------------------------------
        set(figure(figN),'name','Comparison of Fine and Coarse task','pos',[100,100,1500,700]); clf
        for i = 1:2 % h r
            subplot(1,2,i)
            if site_num(i)~=0
                %                 plot(fine_PSE(do_both(:,i),i),coarse_PSE(do_both(:,i),i),'o','color',c{i});
                hold on
                % actual remove num
                remove_num = sum(do_both(:,i) & outlier(:,i));
                
                % remove outlier
                do_both(:,i) = logical(do_both(:,i) & ~outlier(:,i));
                % acutal site num
                site_num = sum(do_both);
                
                select_fine_sig = logical(do_both(:,i) & fine_sig(:,i) & ~coarse_sig(:,i));
                select_coarse_sig = logical(do_both(:,i) & ~fine_sig(:,i) & coarse_sig(:,i));
                select_both_sig = logical(do_both(:,i) & fine_sig(:,i) & coarse_sig(:,i));
                select_not_sig = logical(do_both(:,i) & ~fine_sig(:,i) & ~coarse_sig(:,i));
                
                plot(fine_thr(select_not_sig,i),coarse_thr(select_not_sig,i),'o','color',c{99},'markersize',7);
                plot(fine_thr(select_fine_sig,i),coarse_thr(select_fine_sig,i),'+','color','k','markersize',10,'linewidth',1.5);
                plot(fine_thr(select_coarse_sig,i),coarse_thr(select_coarse_sig,i),'x','color','k','markersize',10,'linewidth',1.5);
                plot(fine_thr(select_both_sig,i),coarse_thr(select_both_sig,i),'o','color','k','markerfacecolor',c{i},'markersize',7,'linewidth',1);
                [r,ppp,h1] = plot_corr_line(fine_thr(do_both(:,i),i),coarse_thr(do_both(:,i),i),'MethodOfCorr','Pearson','FittingMethod',2,'LineWidth',2);
                
                axis tight
%                 h6 = SameScale(1); % same scale, diagonal line, keep tick not change
                
                xlimit = xlim;
                ylimit = ylim;
                text(xlimit(1)+(xlimit(2)-xlimit(1))/20,ylimit(2)-(ylimit(2)-ylimit(1))/20,sprintf('Pearson correlation: r=%0.3g, p=%0.3g',r,ppp),'fontsize',15);
                text(xlimit(1)+(xlimit(2)-xlimit(1))/20,ylimit(2)-(ylimit(2)-ylimit(1))/8,sprintf('N = %g',site_num(i)),'fontsize',15);
                
                text(xlimit(1)+(xlimit(2)-xlimit(1))/20,ylimit(2)-(ylimit(2)-ylimit(1))/4,sprintf('Remoce outlier = %g',remove_num),'fontsize',15);
                
                xlabel('Fine threshold change');
                ylabel('Coarse threshold change');
                
                xlim(xlimit*1.1);
                ylim(ylim*1.1);
            end
            if i == 1
                legend('Not sig','Fine sig','Coarse sig','Both sig','Fitting','location','best');
            end
            
            if i == 1
                title('Translatioin');
            elseif i ==2
                title('Rotation');
            end
        end
        SetFigure(15);
        suptitle(SetTitle(stimtask_type,monkey_included_for_analysis))
        figN = figN+1;
        
    end

    function f3p5(debug) %'PSE shift affect by repetition'
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        stimtask_type_this{1} = [1 2 3]; % fine
        stimtask_type_this{2} = [4 5 6]; % coarse
        
        % ---------------------- get data ----------------------------
        for fc = 1:2 % find or coarse
            for ind = 1:3 % h r mt
                diffBias_step_temp = [];
                diffBias_step_pref_temp = [];
                diffBias_step_non_temp = [];
                
                diffBias_step_p_temp = [];
                diffBias_step_p_pref_temp = [];
                diffBias_step_p_non_temp = [];
                
                % 将不同stim protocol拼接起来
                for i = 1:length(stimtask_type_this{fc})
                    % 1. bias
                    % 暂时是使用原来的bias，而不用除了threshold的 nor_diffBias_task_step, 取绝对值
                    diffBias_step_temp = [diffBias_step_temp; diffBias_task_step{stimtask_type_this{fc}(i),ind}(methods_of_select{1,1},:)];
                    % 根据最后一个bias区分bias“左”或“右”
                    bias_pref = logical(diffBias_task{stimtask_type_this{fc}(i)}(:,ind)>0);  % bias to preferred direction
                    bias_non = logical(diffBias_task{stimtask_type_this{fc}(i)}(:,ind)<0);
                    
                    diffBias_step_pref_temp = [diffBias_step_pref_temp; diffBias_task_step{stimtask_type_this{fc}(i),ind}(methods_of_select{1,1} & bias_pref,:)];
                    diffBias_step_non_temp = [diffBias_step_non_temp; diffBias_task_step{stimtask_type_this{fc}(i),ind}(methods_of_select{1,1} & bias_non,:)];
                    
                    % 2. P value
                    diffBias_step_p_temp = [diffBias_step_p_temp; diffBias_p_task_step{stimtask_type_this{fc}(i),ind}(methods_of_select{1,1},:)];
                    % P value 区分pref， 根据tuning
                    diffBias_step_p_pref_temp = [diffBias_step_p_pref_temp; diffBias_p_task_step{stimtask_type_this{fc}(i),ind}(methods_of_select{1,1} & bias_pref,:)];
                    diffBias_step_p_non_temp = [diffBias_step_p_non_temp; diffBias_p_task_step{stimtask_type_this{fc}(i),ind}(methods_of_select{1,1} & bias_non,:)];
                end
                
                % 1. bias
                diffBias_step_abs = abs(diffBias_step_temp);
                % 当使用nor_bias的时候可能会出现异常值（除以threshold时发生），去掉异常值
                diffBias_step_abs(diffBias_step_abs>10) = nan;
                % normalize
                % normc(M)：对m的每一列进行归一化（0，1）
                diffBias_step_abs_nor = normr(diffBias_step_abs);
                diffBias_step_abs_nor(isnan(diffBias_step_abs)) = nan;
                diffBias_step_pref_nor = normr(diffBias_step_pref_temp);
                diffBias_step_pref_nor(isnan(diffBias_step_pref_temp)) = nan;
                diffBias_step_non_nor = normr(diffBias_step_non_temp);
                diffBias_step_non_nor(isnan(diffBias_step_non_temp)) = nan;
                % 取mean+sem
                diffBias_step_nor_mean{fc,ind} = nanmean(diffBias_step_abs_nor);
                diffBias_step_nor_sem{fc,ind} = nanstd(diffBias_step_abs_nor)./sqrt(sum(~isnan(diffBias_step_abs_nor)));
                diffBias_step_pref_nor_mean{fc,ind} = nanmean(diffBias_step_pref_nor);
                diffBias_step_pref_nor_sem{fc,ind} = nanstd(diffBias_step_pref_nor)./sqrt(length(diffBias_step_pref_nor));
                diffBias_step_non_nor_mean{fc,ind} = nanmean(diffBias_step_non_nor);
                diffBias_step_non_nor_sem{fc,ind} = nanstd(diffBias_step_pref_nor)./sqrt(length(diffBias_step_non_nor));
                
                % 2. p value, 计算<0.05的数目
                diffBias_step_p_num{fc,ind} = sum(~isnan(diffBias_step_p_temp));
                diffBias_step_p_pref_num{fc,ind} = sum(~isnan(diffBias_step_p_pref_temp));
                diffBias_step_p_non_num{fc,ind} = sum(~isnan(diffBias_step_p_non_temp));
                diffBias_step_p_sig{fc,ind} = sum(diffBias_step_p_temp < 0.05);
                
                diffBias_step_p_pref_sig{fc,ind} = sum(diffBias_step_p_pref_temp < 0.05)   ;  % bias to pref and sig
                diffBias_step_p_pref_nonsig{fc,ind} = sum(diffBias_step_p_pref_temp >= 0.05); % bias to pref but not sig
                diffBias_step_p_non_sig{fc,ind} = sum(diffBias_step_p_non_temp < 0.05);       % bias to non-pref and sig
                diffBias_step_p_non_nonsig{fc,ind} = sum(diffBias_step_p_non_temp >= 0.05);   % bias to non-pref and not sig
            end
        end
        
        % ------------------------ Plot ------------------------------
        % 1. Bias shift
        set(figure(figN),'name','PSE shift affected by stim repetition','pos',[200,50,1000,900]); clf
        ha = tight_subplot(4,3,[0.07 0.05],[0.08,0.05],[0.08 0.05]);
        for fc = 1:2 % find and coarse
            for ind = 1:3 % h r mt
                % 1.1. seperate bias to pref and non-pref
                axes(ha(ind+(fc-1)*3))
                hold on
                set(errorbar([1:step_num_max]+4, diffBias_step_pref_nor_mean{fc,ind},  diffBias_step_pref_nor_sem{fc,ind}),'color',c{ind},'linestyle','-','linewidth',1.5);
                set(errorbar([1:step_num_max]+4, diffBias_step_non_nor_mean{fc,ind},  diffBias_step_non_nor_sem{fc,ind}),'color',c{ind},'linestyle',':','linewidth',1.5);
                set(gca,'xticklabelmode','auto','yticklabelmode','auto');
                xlimit = xlim;
                ylimit = ylim;
                ylimit = max(abs(ylimit));
                ylim([-ylimit ylimit]);
                plot(xlimit,[0 0],'k:');
                if fc == 1 && ind == 1
                    ylabel('Fine PSE shift');
                elseif fc == 2 && ind == 1
                    ylabel('Coarse PSE shift');
                end
                if fc == 1 && ind == 1
                    legend('Bias Pref','Bias Null','location','best');
                end
                
                % 1.2. normalize abs all PSE shift
                axes(ha(ind+(fc-1)*3+6))
                shadedErrorBar([1:step_num_max]+4,diffBias_step_nor_mean{fc,ind},diffBias_step_nor_sem{fc,ind},{'color',c{ind},'linewidth',1.5},0.5)
                % set(errorbar([1:step_num_max]+4, diffBias_step_nor_mean{fc,ind},  diffBias_step_nor_sem{fc,ind}),'color',c{ind},'linewidth',1.5);
                if fc == 1 && ind == 1
                    ylabel('| Fine PSE shift |');
                elseif fc == 2 && ind == 1
                    ylabel('| Coarse PSE shift |');
                end
                if fc == 2
                    xlabel('First N repetition');
                end
            end
        end
        SetFigure(11);figN = figN+1;
        
        % 2. P value
        set(figure(figN),'name','P for PSE shift affected by stim repetition','pos',[50,50,1000,900]); clf
        ha2 = tight_subplot(4,3,[0.07 0.05],[0.08,0.05],[0.08 0.05]);
        for fc = 1:2 % find and coarse
            for ind = 1:3 % h r mt
                % 1.1. seperate p value to pref and non-pref
                axes(ha2(ind+(fc-1)*3))
                hold on
                %                 diffBias_step_p_pref_sig{fc,ind} = sum(diffBias_step_p_pref_temp < 0.05)   ;  % bias to pref and sig
                %                 diffBias_step_p_pref_nonsig{fc,ind} = sum(diffBias_step_p_pref_temp >= 0.05); % bias to pref but not sig
                %                 diffBias_step_p_non_sig{fc,ind} = sum(diffBias_step_p_non_temp < 0.05);       % bias to non-pref and sig
                %                 diffBias_step_p_non_nonsig{fc,ind} = sum(diffBias_step_p_non_temp >= 0.05);   % bias to non-pref and not sig
                
                bar([diffBias_step_p_pref_sig{fc,ind}',diffBias_step_p_pref_nonsig{fc,ind}'],1,'stack');
                bar([-diffBias_step_p_non_sig{fc,ind}',-diffBias_step_p_non_nonsig{fc,ind}'],1,'stack');
                set(gca,'xticklabelmode','auto','yticklabelmode','auto');
                
                % 1.2. normalize abs all PSE shift
                axes(ha2(ind+(fc-1)*3+6))
                hold on
                bar([diffBias_step_p_pref_sig{fc,ind}'./ diffBias_step_p_num{fc,ind}',diffBias_step_p_pref_nonsig{fc,ind}'./ diffBias_step_p_num{fc,ind}'],1,'stack');
                bar([-diffBias_step_p_non_sig{fc,ind}'./ diffBias_step_p_num{fc,ind}',-diffBias_step_p_non_nonsig{fc,ind}'./ diffBias_step_p_num{fc,ind}'],1,'stack');
                set(gca,'xticklabelmode','auto','yticklabelmode','auto');
            end
        end
        SetFigure(11);figN = figN+1;
    end

    function f3p9(debug) % 'Time course of microstimulation effects' Salzman,  Murasugi, Britten, Newsome,1992 Figure 10
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        select_stars_type = [0 1]; % 0: dot stimulus; 1: triangle stimulus
        
        % ---------------------- get data ----------------------------
        % select for stars_type: 0, 1, [0 1](combine)
        stars_type = cell2mat({group_result.stars_type})';
        select1 = false(size(stars_type)); % select for stars_type
        for st = 1:length(select_stars_type)
            select1_temp = logical(stars_type == select_stars_type(st));
            select1 = logical(select1_temp | select1);
        end
        
        ctrlPSE = arrayfun(@(x) nan(1000,100), [1:3;1:3],'UniformOutput',0);  % large enough
        stimPSE = arrayfun(@(x) nan(1000,100), [1:3;1:3],'UniformOutput',0);
        ctrlCR = arrayfun(@(x) nan(1000,100), [1:3],'UniformOutput',0);  % large enough
        stimCR = arrayfun(@(x) nan(1000,100), [1:3],'UniformOutput',0);
        
        nn = ones(2,3);
        for n = 1:length(loadM_data)
            if (methods_of_select{1,1}(n) && select1(n))
                for fc = 1:length(stimtask_type) % fine or coarse，需要combine一起，所以最后的index没有fr
                    if ~isempty(stimtask_type{fc})
                        for i = 1:length(stimtask_type{fc}) % 具体task, 需要combine一起，所以最后的index没有i
                            for ind = 1:3
                                % PSE shift
                                if ~isnan(group_result(n).slide_positive_pro{stimtask_type{fc}(i),ind}(1,1)) && group_result(n).p_mu{stimtask_type{fc}(i)}(ind)<0.05 % select significant stim effect site
                                    this_length1 = length(group_result(n).slide_positive_pro{stimtask_type{fc}(i),ind}(1,:));
                                    this_length2 = length(group_result(n).slide_positive_pro{stimtask_type{fc}(i),ind}(2,:));
                                    
                                    %  slide_positive_pro{n,k,ind}(slide_num) = sum(select_this_slide & choice_LEFT_2T{ind}) / sum(select_this_slide & (choice_LEFT_2T{ind}| choice_RIGHT_2T{ind}));
                                    ctrlPSE{1,ind}(nn(1,ind),1:this_length1) = group_result(n).slide_positive_pro{stimtask_type{fc}(i),ind}(1,:);
                                    stimPSE{1,ind}(nn(1,ind),1:this_length2) = group_result(n).slide_positive_pro{stimtask_type{fc}(i),ind}(2,:);
                                    nn(1,ind) = nn(1,ind)+1;
                                end
                                
                                if ~isnan(group_result(n).slide_positive_pro{stimtask_type{fc}(i),ind}(1,1)) && group_result(n).p_mu{stimtask_type{fc}(i)}(ind)>0.05 % select not significant stim effect site
                                    this_length1 = length(group_result(n).slide_positive_pro{stimtask_type{fc}(i),ind}(1,:));
                                    this_length2 = length(group_result(n).slide_positive_pro{stimtask_type{fc}(i),ind}(2,:));
                                    
                                    ctrlPSE{2,ind}(nn(2,ind),1:this_length1) = group_result(n).slide_positive_pro{stimtask_type{fc}(i),ind}(1,:);
                                    stimPSE{2,ind}(nn(2,ind),1:this_length2) = group_result(n).slide_positive_pro{stimtask_type{fc}(i),ind}(2,:);
                                    nn(2,ind) = nn(2,ind)+1;
                                end
                                
                                % CR
                                this_length = length(group_result(n).slide_CR{stimtask_type{fc}(i),ind}(1,:));
                                ctrlCR{ind}(nn(1,ind),1:this_length) = group_result(n).slide_CR{stimtask_type{fc}(i),ind}(1,:);
                                stimCR{ind}(nn(1,ind),1:this_length) = group_result(n).slide_CR{stimtask_type{fc}(i),ind}(2,:);
                            end
                        end
                    end
                end
            end
        end
        
        slide_window = 30; % set by CueSwitching_performance
        slide_step = 10;
        first_window_mid = (1+slide_window)/2;
        for ind = 1:3
            for sig = 1:2
                ctrlPSE{sig,ind}(all(isnan(ctrlPSE{sig,ind}),2),:) = []; % Trim the large matrix by removing unnecessary NaNs. HH20150704
                ctrlPSE{sig,ind}(:,all(isnan(ctrlPSE{sig,ind}),1)) = []; % Trim the large matrix by removing unnecessary NaNs. HH20150704
                stimPSE{sig,ind}(all(isnan(stimPSE{sig,ind}),2),:) = []; % Trim the large matrix by removing unnecessary NaNs. HH20150704
                stimPSE{sig,ind}(:,all(isnan(stimPSE{sig,ind}),1)) = []; % Trim the large matrix by removing unnecessary NaNs. HH20150704
                
                % mean and sem PSE shift
                mean_ctrlPSE{sig,ind} = nanmean(ctrlPSE{sig,ind},1);
                mean_stimPSE{sig,ind} = nanmean(stimPSE{sig,ind},1);
                
                sem_ctrlPSE{sig,ind} = nanstd(ctrlPSE{sig,ind},1) ./ sum(~isnan(ctrlPSE{sig,ind}),1);
                sem_stimPSE{sig,ind} = nanstd(stimPSE{sig,ind},1) ./ sum(~isnan(stimPSE{sig,ind}),1);
                
                % 横坐标
                step_num = size(ctrlPSE{sig,ind},2);
                PSEslide_mid_trial{sig,ind} = (first_window_mid:slide_step:1000); % long enough
                PSEslide_mid_trial{sig,ind}(step_num+1:end) = [];
            end
            
            % mean CR
            ctrlCR{ind}(all(isnan(ctrlCR{ind}),2),:) = []; % Trim the large matrix by removing unnecessary NaNs. HH20150704
            ctrlCR{ind}(:,all(isnan(ctrlCR{ind}),1)) = []; % Trim the large matrix by removing unnecessary NaNs. HH20150704
            stimCR{ind}(all(isnan(stimCR{ind}),2),:) = []; % Trim the large matrix by removing unnecessary NaNs. HH20150704
            stimCR{ind}(:,all(isnan(stimCR{ind}),1)) = []; % Trim the large matrix by removing unnecessary NaNs. HH20150704
            
            mean_ctrlCR{ind} = nanmean(ctrlCR{ind},1);
            mean_stimCR{ind} = nanmean(stimCR{ind},1);
            
            step_num = size(ctrlCR{ind},2);
            CRslide_mid_trial{ind} = (first_window_mid:slide_step:1000); % long enough
            CRslide_mid_trial{ind}(step_num+1:end) = [];

        end
        
        % ------------------------ Plot ------------------------------
        for sig = 1:2
            if sig == 1
                set(figure(figN),'name','Time course of induced PSE shift(sig site)','pos',[50,500,1000,300]); clf
            else
                set(figure(figN+1),'name','Time course of induced PSE shift(not sig site)','pos',[50,50,1000,300]); clf
            end
            for ind = 1:3
                subplot(1,3,ind)
                if ~isnan(nanmean(ctrlPSE{sig,ind}))
                    hold on
                    plot(PSEslide_mid_trial{sig,ind},mean_ctrlPSE{sig,ind},'ko');
                    plot(PSEslide_mid_trial{sig,ind},mean_stimPSE{sig,ind},'*','color',c{ind});
                    %                     errorbar(slide_mid_trial{sig,ind},mean_control{sig,ind},sem_control{sig,ind},'ko'); % only plot positive pproportion
                    %                     errorbar(slide_mid_trial{sig,ind},mean_stim{sig,ind}, sem_stim{sig,ind},'*','color',c{ind}); % only plot positive pproportion
                    plot([0 max(PSEslide_mid_trial{sig,ind})],[0.5 0.5],'k:');
                    
                    % fitting
                    plot_corr_line(PSEslide_mid_trial{sig,ind},mean_ctrlPSE{sig,ind},'MethodOfCorr','pearson','FittingMethod',1,'LineStyles','k-','LineWidth',2);
                    plot_corr_line(PSEslide_mid_trial{sig,ind},mean_stimPSE{sig,ind},'MethodOfCorr','pearson','FittingMethod',1,'LineStyles','y-','LineWidth',2);
                end
                xlabel('Trial Number');
                ylabel('Proportion of Stim choice'); % 电刺激、control条件下，最终电刺激bias方向的选择的比例
                ylim([0 1]);
                if ind == 1
                    title("Translation");
                elseif ind == 2
                    title('Rotation');
                else
                    title("Motion type");
                end
            end
            SetFigure(12);
            til = SetTitle(stimtask_type,monkey_included_for_analysis);
            if sig == 1
                til2 = strcat(til,' (Significant site)');
            else
                til2 = strcat(til,' (Not significant site)');
            end
            suptitle(til2);
        end
        figN = figN + 2;
        
        set(figure(figN),'name','Time course of CR','pos',[1000,50,1000,300]); clf
        for ind = 1:3
            subplot(1,3,ind)
            if ~isnan(nanmean(ctrlPSE{sig,ind}))
                hold on
                plot(CRslide_mid_trial{ind},mean_ctrlCR{ind},'ko');
                plot(CRslide_mid_trial{ind},mean_stimCR{ind},'*','color',c{ind});
                plot([0 max(CRslide_mid_trial{ind})],[0.5 0.5],'k:');
                
                % fitting
                plot_corr_line(CRslide_mid_trial{ind},mean_ctrlCR{ind},'MethodOfCorr','pearson','FittingMethod',1,'LineStyles','k-','LineWidth',2);
                plot_corr_line(CRslide_mid_trial{ind},mean_stimCR{ind},'MethodOfCorr','pearson','FittingMethod',1,'LineStyles','y-','LineWidth',2);
            end
            xlabel('Trial Number');
            ylabel('CR'); % 电刺激、control条件下，最终电刺激bias方向的选择的比例
            ylim([0 1]);
            if ind == 1
                title("Translation");
            elseif ind == 2
                title('Rotation');
            else
                title("Motion type");
            end
        end
        SetFigure(12);
        til = SetTitle(stimtask_type,monkey_included_for_analysis);
        til2 = strcat('CR for ',til);
        suptitle(til2);
        figN = figN + 2;
    end


    function f3p6(debug) %'PSE shift affect by maxFR'
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        effect_type = 1; % 1: PSE, 2: threshold
        is_separate_stimAmp = 0;
        is_original_bias = 0;  % 正值：bias to prefer direction
        select_stars_type = [0 1]; % 0: dot stimulus; 1: triangle stimulus
        
        % ---------------------- get data ----------------------------
        % 为了区分不同task对star type的利用，分开heading和rotation
        [diff_temp, S_unit_temp, NS_unit_temp, ~] = get_stim_data_uA(stimtask_type,methods_of_select{1,1},is_separate_stimAmp,is_original_bias,select_stars_type);
        PSE_shift = diff_temp{effect_type};
        S_unit = S_unit_temp{effect_type};
        NS_unit = NS_unit_temp{effect_type};
        
        % maxFR
        maxFR = cell2mat({group_result.maxFR}');
        maxFR = repmat(maxFR,length(cell2mat(stimtask_type)),1);
        
        % ------------------------ Plot ------------------------------
        set(figure(figN),'name','PSE shift affected by maxFR','pos',[50,50,1000,500]); clf
        for hr = 1:2
            subplot(1,2,hr)
            plot(maxFR(:,hr),PSE_shift(:,hr),'ko','markersize',8);
            hold on
            plot(maxFR(S_unit(:,hr)==1,hr),PSE_shift(S_unit(:,hr)==1,hr),'ko','markerfacecolor',c{hr},'markersize',8);
            [r,p] = plot_corr_line(maxFR(:,hr),PSE_shift(:,hr),'MethodOfCorr','pearson','FittingMethod',2,'LineStyles','k-','LineWidth',2);
            xlimit = xlim;
            ylimit = ylim;
            text(xlimit(1)+(xlimit(2)-xlimit(1))/20,ylimit(2)-(ylimit(2)-ylimit(1))/20,sprintf('r=%0.3g,p=%0.3g',r,p),'fontsize',15);
            xlabel('maxFR (Hz)');
            ylabel('PSE shift');
        end
        SetFigure(13);
        suptitle(SetTitle(stimtask_type,monkey_included_for_analysis))
        figN = figN+1;
    end

    function f3p7(debug) % 'PSE shift affect by cluster'
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        effect_type = 1; % 1: PSE, 2: threshold
        is_separate_stimAmp = 40;
        select_stars_type = [0 1]; % 0: dot stimulus; 1: triangle stimulus
        is_original_bias = 0;  % 正值：bias to prefer direction
        
        % ---------------------- get data ----------------------------
        [PSE_shift_temp, S_unit_temp, NS_unit_temp, ~] = get_stim_data_uA(stimtask_type,methods_of_select{1,1},is_separate_stimAmp,is_original_bias,select_stars_type);
        PSE_shift = PSE_shift_temp{effect_type};
        S_unit = S_unit_temp{effect_type};
        NS_unit = NS_unit_temp{effect_type};
        
        % cluster
        % 1. SD for Pd
        pd_std = cell2mat({group_result.pd_std}'); % circular std
        pd_std2 = cell2mat({group_result.pd_std2}'); % simple std
        pd_std = repmat(pd_std,length(cell2mat(stimtask_type)),1);
        pd_std2 = repmat(pd_std2,length(cell2mat(stimtask_type)),1);
        
        % 取对数
        pd_std = log(pd_std);
        
        [~,p] = ttest(pd_std(:,1),pd_std(:,2));
        
        % 2. correlation coeficient in 100 um
        
        
        % ------------------------ Plot ------------------------------
        set(figure(figN),'name','PSE shift affected by cluster','pos',[50,50,1000,600]); clf
        for hr = 1:2
            subplot(1,2,hr)
            plot(pd_std(:,hr),PSE_shift(:,hr),'ko');
            hold on
            plot(pd_std(S_unit(:,hr)==1,hr),PSE_shift(S_unit(:,hr)==1,hr),'ko','markerfacecolor',c{hr},'markersize',8);
            [r,p] = plot_corr_line(pd_std(:,hr),PSE_shift(:,hr),'MethodOfCorr','pearson','FittingMethod',2,'LineStyles','k-','LineWidth',2);
            xlimit = xlim;
            ylimit = ylim;
            ylim([-max(abs(ylimit)) max(abs(ylimit))]);
            plot([xlimit(1) xlimit(2)],[0 0],'k:');
            text(xlimit(1)+(xlimit(2)-xlimit(1))/20,ylimit(2)-(ylimit(2)-ylimit(1))/20,sprintf('r=%0.3g,p=%0.3g',r,p),'fontsize',15);
            axis square
            set(gca,'xtick',[-2:2:6]);
%             set(gca,'xtick',[-2:2:6],'xticklabel',{'10^{-2}','10^{0}','10^{2}','10^{4}','10^{6}'});
            xlabel('SD of preferred direction (log)');
            ylabel('PSE shift');
            if hr == 1
                title('Translation');
            else
                title('Rotation');
            end
        end
        SetFigure(15)
        
        til = SetTitle(stimtask_type,monkey_included_for_analysis);
        til2 = strcat('PSE shift affected by cluster in',32,til);
        suptitle(til2);
        figN = figN+1;
    end

    function f3p8(debug) % 'PSE shift affect by RF'
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        effect_type = 1; % 1: PSE, 2: threshold
        is_separate_stimAmp = 0;
        is_original_bias = 0;  % 正值：bias to prefer direction
        select_stars_type = [0 1];
        
        % ---------------------- get data ----------------------------
        [PSE_shift_temp, S_unit_temp, NS_unit_temp, ~] = get_stim_data_uA(stimtask_type,methods_of_select{1,1},is_separate_stimAmp,is_original_bias,select_stars_type);
        PSE_shift = PSE_shift_temp{effect_type};
        S_unit = S_unit_temp{effect_type};
        NS_unit = NS_unit_temp{effect_type};
        
        % RF size, RF eccentricity
        RF = cell2mat({group_result.RF}');
        rf_size = RF(:,3).*RF(:,4);
        rf_dia = 2.* sqrt(rf_size / pi);
        rf_ecc = sqrt(RF(:,1).^2+RF(:,2).^2);
        
        rf_dia = repmat(rf_dia,length(cell2mat(stimtask_type)),1);
        rf_ecc = repmat(rf_ecc,length(cell2mat(stimtask_type)),1);
        
        % ------------------------ Plot ------------------------------
        set(figure(figN),'name','PSE shift affected by RF diameter','pos',[50,50,700,350]); clf
        for hr = 1:2
            subplot(1,2,hr)
            plot(rf_dia,PSE_shift(:,hr),'ko');
            hold on
            plot(rf_dia(S_unit(:,hr)==1),PSE_shift(S_unit(:,hr)==1,hr),'ko','markerfacecolor',c{hr});
            [r,p] = plot_corr_line(rf_dia,PSE_shift(:,hr),'MethodOfCorr','pearson','FittingMethod',2,'LineStyles','k-','LineWidth',2);
            xlimit = xlim;
            ylimit = ylim;
            text(xlimit(1)+(xlimit(2)-xlimit(1))/20,ylimit(2)-(ylimit(2)-ylimit(1))/20,sprintf('r=%0.3g,p=%0.3g',r,p),'fontsize',15);
            xlabel('RF diameter');
            ylabel('PSE shift');
            set(gca,'xtick',[0:20:100]);
            axis square;
        end
        SetFigure(13);figN = figN+1;
        
        set(figure(figN),'name','PSE shift affected by RF eccentricity','pos',[750,50,700,350]); clf
        for hr = 1:2
            subplot(1,2,hr)
            plot(rf_ecc,PSE_shift(:,hr),'ko');
            hold on
            plot(rf_ecc(S_unit(:,hr)==1),PSE_shift(S_unit(:,hr)==1,hr),'ko','markerfacecolor',c{hr});
            [r,p] = plot_corr_line(rf_ecc,PSE_shift(:,hr),'MethodOfCorr','pearson','FittingMethod',2,'LineStyles','k-','LineWidth',2);
            xlimit = xlim;
            ylimit = ylim;
            text(xlimit(1)+(xlimit(2)-xlimit(1))/20,ylimit(2)-(ylimit(2)-ylimit(1))/20,sprintf('r=%0.3g,p=%0.3g',r,p),'fontsize',15);
            xlabel('RF eccentricity');
            ylabel('PSE shift');
            %             set(gca,'xtick',[0:10:50]);
            axis square;
        end
        SetFigure(13);figN = figN+1;
        
        set(figure(figN),'name','PSE shift affected by RF','pos',[750,50,700,350]); clf
        for hr = 1:2
            subplot(1,2,hr)
            plot3(rf_dia,rf_ecc,PSE_shift(:,hr),'ko');
            hold on
            plot3(rf_dia(S_unit(:,hr)==1),rf_ecc(S_unit(:,hr)==1),PSE_shift(S_unit(:,hr)==1,hr),'ko','markerfacecolor',c{hr});
            xlabel('RF diameter');
            ylabel('RF eccentricity');
            zlabel('PSE shift');
            %             set(gca,'xtick',[0:10:50]);
        end
        SetFigure(13);figN = figN+1;
    end

    function f3p10(debug) % regression
        if debug  ; dbstack;   keyboard;      end
        methods_of_select = {select_typical, 'Typical cells'};
        effect_type = 1; % 1: PSE, 2: threshold
        is_separate_stimAmp = 0;
        is_original_bias = 2;  % 0: bias positive = prefer direction; 1: only original bias in MT;  2: original bias;
        select_stars_type = [0 1];
        
        % ---------------------- get data ----------------------------
        % PSE (only translation and rotation PSE shift)
        % 区分fine and coarse，2 or 4-Target合并
        % stimCueS stimHonly stimRonly stimCo stimCoH stimCoR
        disp('******   combine 2- and 4-Target and separate fine and coarse task   ******');
        
        select = repmat(methods_of_select{1,1},2,1);
        
        % original PSE shift (after-before)
        diffBias2Tfine_ori = [diffBias_task_origin{2}(:,1) diffBias_task_origin{3}(:,2)];
        diffBias2Tcoarse_ori = [diffBias_task_origin{5}(:,1) diffBias_task_origin{6}(:,2)];

        diffBias4Tfine_ori = diffBias_task_origin{1}(:,[1:2]);
        diffBias4Tcoarse_ori = diffBias_task_origin{4}(:,[1:2]);
        
        diffBias_fine_ori = [diffBias2Tfine_ori;diffBias4Tfine_ori];
        diffBias_coarse_ori = [diffBias2Tcoarse_ori;diffBias4Tcoarse_ori];
        diffBias_fine_ori = diffBias_fine_ori(select,:);
        diffBias_coarse_ori = diffBias_coarse_ori(select,:);
        
        % normal PSE shift (pref-null)
        diffBias2Tfine = [diffBias_task{2}(:,1) diffBias_task{3}(:,2)];
        diffBias2Tcoarse = [diffBias_task{5}(:,1) diffBias_task{6}(:,2)];

        diffBias4Tfine = diffBias_task{1}(:,[1:2]);
        diffBias4Tcoarse = diffBias_task{4}(:,[1:2]);
        
        diffBias_fine = [diffBias2Tfine_ori;diffBias4Tfine];
        diffBias_coarse = [diffBias2Tcoarse_ori;diffBias4Tcoarse];
        diffBias_fine = diffBias_fine(select,:);
        diffBias_coarse = diffBias_coarse(select,:);
        
        % RF size, RF eccentricity
        RF = cell2mat({group_result.RF}');
        rf_size = RF(:,3).*RF(:,4);
        rf_dia_temp = 2.* sqrt(rf_size / pi);
        rf_ecc_temp = sqrt(RF(:,1).^2+RF(:,2).^2);
        rf_dia = repmat(rf_dia_temp,2,1);
        rf_ecc = repmat(rf_ecc_temp,2,1);
        rf_dia = rf_dia(select,:);
        rf_ecc = rf_ecc(select,:);
        
        % recording X/Y location
        GridX_temp = cell2mat({group_result.GridX}');
        GridY_temp = cell2mat({group_result.GridY}');
        GridX = repmat(GridX_temp,2,1);
        GridY = repmat(GridY_temp,2,1);
        GridX = GridX(select,:);
        GridY = GridY(select,:);
        
        % direction preference
        pref_direc_temp = cell2mat({group_result.pref_direc}');
        pref_direc = repmat(pref_direc_temp,2,1);
        pref_direc = pref_direc(select,:);
        % 将prefer directino转为1维:-90 -45/-135 0/180/-180 45/135 90 (距离±90的角度大小，且正负号表示原来的左或者右)
        pref_LR_diff = nan(size(pref_direc(:,[1 2])));
        for i=1:2
            pref_LR_diff(pref_direc(:,i)<0,i) = -90+abs(pref_direc(pref_direc(:,i)<0,i)-(-90));
            pref_LR_diff(pref_direc(:,i)>=0,i) = 90-abs(pref_direc(pref_direc(:,i)>0,i)-(90));
        end

        % cluster index: Gu:cluster index可定义为每个位点的tuning和临近两个位点的pearson correlation 的平均值（上下各一个位点）
        cluster_index = cell2mat({group_result.cluster_index}');
        cluster_index = repmat(cluster_index,2,1);
        cluster_index = cluster_index(select,:);
        
        % 暂时用90°的d'
        dprime = cell2mat({group_result.d_prime_coarse}'); % d' = (Rright90 - Rleft90) / sqrt(1/2*(STDright^2+STDleft^2))
        dprime = repmat(dprime,2,1);
        dprime = dprime(select,:);
        dprime = abs(dprime);
        
        % ------------------------ Plot ------------------------------
        %%% 1. regression with d', cluster, PD
        % combine 2T and 4T, separate T and R, fine and coarse
        diffBias{1} = diffBias_fine; % normal PSE shift, positive-Pref, negative-Null
        diffBias{2} = diffBias_coarse; % normal PSE shift
        
        % 因为PSE shift和prefer direction是非线性关系，调整PD为距离±90°的距离
        pref_direc_90 = abs(abs(pref_direc)-90);
        
        row_this = 0; % for subplot
        set(figure(figN),'name','PSE shift correlation with d prime/cluster/PD','pos',[50,50,1000,1000]); clf
        for hr = 1:2
            for fc = 1:2
                y = diffBias{fc}(:,hr); 
%                 X = [ones(size(y)) dprime(:,hr) cluster_index(:,hr) pref_direc(:,hr)]; 
                X = [dprime(:,hr) cluster_index(:,hr) pref_direc_90(:,hr)]; % d', cluster index, prefer direction
                ylabletext{1} = '|d prime|';
                ylabletext{2} = 'cluster index';
                ylabletext{3} = 'preferred direction';
                
                % 1. linear regression
                [B,BINT,R,RINT,STATS] = regress(y,[ones(size(y)) X]); % 常数项，d', cluster index, prefer direction
                R2 = STATS(1);
                
                % 2. partial correlation
                [rr,pp] = partialcorr([y X],'rows','complete');
                partial_corr_coef = rr(1,2:4);
                partial_corr_p = pp(1,2:4);
                
                
                % 3. correlation
                for i = 1:size(X,2) % first column is 1
                    subplot(4,size(X,2),row_this+i)
                    plot(X(:,i),y,'k.');
                    [r,p] = plot_corr_line(X(:,i),y,'MethodOfCorr','Pearson');
                    xlimit = xlim;
                    ylimit = ylim;
                    str = [sprintf('r = %.2g',r) newline sprintf('p = %.2g',p)];
                    text(xlimit(1),ylimit(2),str);
                    xlabel(ylabletext{i});
                    
                    if hr == 1 && fc == 1
                        ylabel('Trans. fine PSE shift');
                    elseif hr == 1 && fc == 2
                        ylabel('Trans. coarse PSE shift');
                    elseif hr == 2 && fc == 1
                        ylabel('Rot. fine PSE shift');
                    elseif hr == 2 && fc == 2
                        ylabel('Rot. coarse PSE shift');
                    end
                    
                    str = ['Partial correlation' newline sprintf('Corr. coef. = %.2g',partial_corr_coef(i)) ...
                        newline sprintf('p = %.2g',partial_corr_p(i))];
                    text(xlimit(2)/2,ylimit(2),str);
                end
                
                row_this = row_this + size(X,2);
            end
        end
        SetFigure(10);
        

        %%% 2. regression
        figure
        subplot(2,2,1)
        hist(pref_direc(:,1));title('Translation prefer direction');
        subplot(2,2,3)
        hist(pref_direc(:,2));title('Rotation prefer direction');
        subplot(2,2,2)
        hist(pref_LR_diff(:,1));title('Translation |prefer direction-90|');
        subplot(2,2,4)
        hist(pref_LR_diff(:,2));title('Rotation |prefer direction-90|');
        
        % PSE shift - pref_LR_diff correlation
        figure
        subplot(1,2,1) % fine
        plot(pref_LR_diff(:,1),diffBias_fine(:,1),'.','color',c{1});
        [r,p] = plot_corr_line(pref_LR_diff(:,1),diffBias_fine(:,1),'MethodOfCorr','Pearson','FittingMethod',2,'LineStyles','r-');
        hold on
        plot(pref_LR_diff(:,2),diffBias_fine(:,2),'.','color',c{2});
        [r,p] = plot_corr_line(pref_LR_diff(:,2),diffBias_fine(:,2),'MethodOfCorr','Pearson','FittingMethod',2,'LineStyles','b-');
        xlabel('|prefer direction-90|');
        ylabel('PSE shift');
        
        subplot(1,2,2) % coarse
        plot(pref_LR_diff(:,1),diffBias_coarse(:,1),'.','color',c{1});
        [r,p] = plot_corr_line(pref_LR_diff(:,1),diffBias_coarse(:,1),'MethodOfCorr','Pearson','FittingMethod',2,'LineStyles','r-');
        hold on
        plot(pref_LR_diff(:,2),diffBias_coarse(:,2),'.','color',c{2});
        [r,p] = plot_corr_line(pref_LR_diff(:,2),diffBias_coarse(:,2),'MethodOfCorr','Pearson','FittingMethod',2,'LineStyles','b-');
        xlabel('|prefer direction-90|');
        ylabel('PSE shift');
        
        
        % regression
        y = (diffBias_coarse(:,1));
        X = [ones(size(y)) pref_LR_diff(:,1) pref_LR_diff(:,2) rf_dia rf_ecc rf_dia.*rf_ecc GridX GridY GridX.*GridY];
        % 常数项，HeadingPref, RotationPref, RFSize，RFecc，RFsize*RFecc，GridX, GridT, GridX*Y
        figure
        for i = 1:size(X,2)-1
            subplot(1,size(X,2)-1,i)
            plot(X(:,i+1),y,'k.');
            plot_corr_line(X(:,i+1),y,'MethodOfCorr','Pearson','FittingMethod',2,'LineStyles','k-');
        end

        % linear regression
        [B,BINT,R,RINT,STATS] = regress(y,X);
        
        X = [ones(size(y)) pref_LR_diff(:,2) rf_dia rf_ecc rf_dia.*rf_ecc GridX GridY GridX.*GridY];
        % linear regression
        [B1,BINT1,R1,RINT1,STATS1] = regress(y,X);
        
        % lasso regression
        [b, fitinfo] = lasso(X, y,'CV',10,'Alpha',1);
        axTrace = lassoPlot(b,fitinfo); %交叉验证训练轨迹
        axCV = lassoPlot(b,fitinfo,'PlotType','CV');
        lam = fitinfo.IndexMinMSE;  % 最小MSE对应lambda
        mat = b(:,lam);             % 最优lambda对应的稀疏系数
        [row, ] = find(b(:,lam)~=0);% 非零系数对应的行
        Xla = X(:, row');
    end

    function f3p11(debug)  % PSE shit and Preferred direction
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        effect_type = 1; % 1: PSE, 2: threshold
        is_separate_stimAmp = 0;
%         is_original_bias = 0;  % 正值：bias to prefer direction
        is_original_bias = 3; % normlize PSE shift (divided by the threshold under non-stim trials
        select_stars_type = [0 1];
        
        % ---------------------- get data ----------------------------
        [PSE_shift_temp, S_unit_temp, NS_unit_temp, ~] = get_stim_data_uA(stimtask_type,methods_of_select{1,1},is_separate_stimAmp,is_original_bias,select_stars_type);
        PSE_shift = PSE_shift_temp{effect_type};
        S_unit = S_unit_temp{effect_type};
        NS_unit = NS_unit_temp{effect_type};
        
        % preferred direction
        pref_direc_temp = cell2mat({group_result.pref_direc}');
        pref_direc = repmat(pref_direc_temp,length(cell2mat(stimtask_type)),1); % 0 = exp
        
        % ------------------------ Plot ------------------------------
        set(figure(figN),'name','PSE shift affected by preferred direction','pos',[50,50,1200,600]); clf
        for hr = 1:2
            % scatter plot
            subplot(3,2,[hr+2 hr+4]);
            reasonable_range = [nanmean(PSE_shift(:,hr)) - 3*nanstd(PSE_shift(:,hr)) nanmean(PSE_shift(:,hr)) + 3*nanstd(PSE_shift(:,hr))];  % mean±3*std
            outdot = logical(PSE_shift(:,hr)<reasonable_range(1) | PSE_shift(:,hr)>reasonable_range(2));
            outdot = false(size(outdot)); % show all
            
            hold on
            if hr == 1
                plot(pref_direc(~outdot,hr),PSE_shift(~outdot,hr),'o','color','k');
                plot(pref_direc((S_unit(:,hr) & ~outdot),hr),PSE_shift((S_unit(:,hr) & ~outdot),hr),'o','color','k','markerfacecolor',c{hr});
            else
                plot(pref_direc(~outdot,hr),PSE_shift(~outdot,hr),'^','color','k');
                plot(pref_direc((S_unit(:,hr) & ~outdot),hr),PSE_shift((S_unit(:,hr) & ~outdot),hr),'^','color','k','markerfacecolor',c{hr});
            end
            
            axis tight
            ylimit = ylim;
            ylimit = ylimit.*1.1;
            ylimit = [-max(abs(ylimit)) max(abs(ylimit))];
            ylim(ylimit);
            
            if sum(outdot)>0
                text(-170,ylimit(2),sprintf('%g dots out of 3*std did not show',sum(outdot)));
%                 figure
%                 plot(pref_direc(outdot,hr),PSE_shift(outdot,hr),'o','color',c{hr});
%                 hold on
%                 plot(pref_direc((S_unit(:,hr) & outdot),hr),PSE_shift((S_unit(:,hr) & outdot),hr),'o','color',c{hr},'markerfacecolor',c{hr});
%                 xlimit = [-180,180];
%                 xlim(xlimit);
%                 xticks([-180:45:180]);
%                 ylimit1 = ylim;
%                 plot([-90 -90],ylimit1,'k:');
%                 plot([0 0],ylimit1,'k:');
%                 plot([90 90],ylimit1,'k:');
            end
            
            plot([-180 180],[0 0],'g-');
            plot([-90 -90],ylimit,'g-');
            plot([0 0],ylimit,'g-');
            plot([90 90],ylimit,'g-');
            xlimit = [-180,180];
            xlim(xlimit);
            xticks([-180:45:180]);
            if hr == 2
                xticklabels_temp = xticklabels;
                xticklabels_temp{[-180:45:180]==-90} = 'CW';
                xticklabels_temp{([-180:45:180]==90)} = 'CCW';
                xticklabels(xticklabels_temp);
            end
            if is_original_bias == 3
                ylabel('Normalized Induced PSE shift (°)');
            elseif is_original_bias == 0
                ylabel('Induced PSE shift (°)');
            end
            xlabel('Direction preference');
            
            % step mean
            ax(hr) = subplot(3,2,hr);
            win = 30;
            step = 10;
            
            if length(cell2mat(stimtask_type))==1
                win = 45;
                step = 10;
            end
            
            % step mean
            [xmid, y, y_sem, y_step_data] = get_slideWin(pref_direc(~outdot,hr),PSE_shift(~outdot,hr),win,step,[-180-win/2 180+win/2]);
            if hr == 1
                shadedErrorBar(xmid,y,y_sem,'r');
            else
                shadedErrorBar(xmid,y,y_sem,'b');
            end
            %             errorbar(xmid,y,y_sem);

            hold on
            axis tight
            ylimit = ylim;
            plot([-180 180],[0 0],'g-');
            plot([-90 -90],ylimit,'g-');
            plot([0 0],ylimit,'g-');
            plot([90 90],ylimit,'g-');
            xlim(xlimit);
            xticks([-180:45:180]);
            if hr == 2
                xticklabels_temp = xticklabels;
                xticklabels_temp{[-180:45:180]==-90} = 'CW';
                xticklabels_temp{([-180:45:180]==90)} = 'CCW';
                xticklabels(xticklabels_temp);
            end
 
             % t test for mean PSE
            for i = 1:length(xmid)
                [~,p(i)] = ttest(y_step_data{i});
                
                % plot significance
                if p(i)<0.05
                    plot([xmid(i)-step/2, xmid(i)+step/2],[ylimit(2) ylimit(2)],'-','color',c{hr},'linewidth',2);
                end
            end
            
%             % significant marker
%             plot(xmid(p<0.05),y(p<0.05)+(ylimit(2)-ylimit(1))/10,'*','color',c{hr});
%             ylim([ylimit(1)*1.1 ylimit(2)*1.1]);
            
        end
        
%         linkaxes([ax(1),ax(2)],'xy'); % 同步多个坐标区的范围     
        
        SetFigure(11)
        [til] = SetTitle(stimtask_type,monkey_included_for_analysis);
        str=['Correlation between microstimulation effect and direction preference' newline til];
        suptitle(str)
        figN = figN+1;
        
        % average reuslt grouped by four relative preferred direction
        % ---------------------- get data ----------------------------
        select = [];
        for hr = 1:2
            select{1} = logical(abs(pref_direc(:,hr))<=45 & abs(pref_direc(:,hr))>=0); % 0-45
            select{2} = logical(abs(pref_direc(:,hr))<=90 & abs(pref_direc(:,hr))>45); % 45-90
            select{3} = logical(abs(pref_direc(:,hr))<=135 & abs(pref_direc(:,hr))>90); % 90-135
            select{4} = logical(abs(pref_direc(:,hr))<=180 & abs(pref_direc(:,hr))>135); % 135-180
            
            for i = 1:4
                PSE_shift_group{hr,i} = PSE_shift(select{i},hr); % 0-45
                mean_PSE_group(hr,i) = nanmean(PSE_shift_group{hr,i});
                sem_PSE_group(hr,i) = nanstd(PSE_shift_group{hr,i}) / sqrt(sum(select{i}));
                
                % significant
                [~,pp(hr,i)] = ttest(PSE_shift_group{hr,i},0);
                if pp(hr,i) < 0.001
                    sig_mark{hr,i} = '***';
                elseif pp(hr,i)>=0.001 && pp(hr,i)<0.01
                    sig_mark{hr,i} = '**';
                elseif pp(hr,i)>=0.01 && pp(hr,i)<0.05
                    sig_mark{hr,i} = '*';
                else
                    sig_mark{hr,i} = [];
                end
            end
        end

        % ------------------------ Plot ------------------------------
        set(figure(figN),'name','PSE shift grouped by four relative preferred direction','pos',[500,50,1200,600]); clf
        for hr = 1:2
            ax(hr) = subplot(1,2,hr);
            barwitherr(sem_PSE_group(hr,:), mean_PSE_group(hr,:),'facecolor',c{hr});
            hold on
            axis tight
            for i = 1:4
                ylimit_this = mean_PSE_group(hr,i)+sem_PSE_group(hr,i)+0.1;
                text(i+0.2,ylimit_this,sprintf('p = %2.3g',pp(hr,i)));
                if ~isempty(sig_mark{hr,i})
                    
                    text(i,ylimit_this,sig_mark{hr,i},'HorizontalAlignment','center');
                end
            end
            xticklabels({'0-45','45-90','90-135','135-180'});
            xlabel('|Preferred direction from forward| (°)');
            if is_original_bias == 3
                ylabel('Normalized Induced PSE shift');
            elseif is_original_bias == 0
                ylabel('Induced PSE shift (°)');
            end
            ylimit = ylim;
            ylim(ylimit.*1.1);
        end
%         linkaxes([ax(1),ax(2)],'xy'); % 同步多个坐标区的范围, same y lim T and R    
        SetFigure(11)
        [til] = SetTitle(stimtask_type,monkey_included_for_analysis);
        str=['PSE shift grouped by four relative preferred direction' newline til];
        suptitle(str)
        figN = figN+1;
    end

    function f3p12(debug) % PSE shift and d'  Gu, 2012 Fig.5
        if debug  ; dbstack;   keyboard;      end
        
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        effect_type = 1; % 1: PSE, 2: threshold
        select_stars_type = [0 1]; % 0: dot stimulus; 1: triangle stimulus
        is_separate_stimAmp = 0;
        is_original_bias = 4; % original bias and normalized (after-before  >0: bias to left/CW; <0: bias to right/CCW) (后面取绝对值，所以没关系)
        
        [PSE_shift_temp, S_unit_temp, NS_unit_temp, ~] = get_stim_data_uA(stimtask_type,methods_of_select{1,1},is_separate_stimAmp,is_original_bias,select_stars_type);
        PSE_shift = PSE_shift_temp{effect_type};
        PSE_shift_abs = abs(PSE_shift);
        S_unit = S_unit_temp{effect_type};
        NS_unit = NS_unit_temp{effect_type};
        
        % find the outlier, may be different from PSE_shift_nor
        for hr = 1:2
            reasonable_range = [nanmean(PSE_shift(:,hr)) - 3*nanstd(PSE_shift(:,hr)), nanmean(PSE_shift(:,hr)) + 3*nanstd(PSE_shift(:,hr))];
            outlier(:,hr) = logical(PSE_shift(:,hr)<reasonable_range(1) | PSE_shift(:,hr)>reasonable_range(2));
        end
        outlier = false(size(outlier));
        
        % 暂时用90°的d‘
        dprime_all = cell2mat({group_result.d_prime_coarse}'); % d' = (Rright90 - Rleft90) / sqrt(1/2*(STDright^2+STDleft^2))
        dprime = repmat(dprime_all,length(cell2mat(stimtask_type)),1);
        dprime_abs = abs(dprime);
        
%         % ROC
%         AuROC = cell2mat({group_result.AuROC}');
%         AuROC = repmat(AuROC,length(cell2mat(stimtask_type)),1);
%         dprime_abs = AuROC;
        
        % save the plot d'
        for hr = 1:2
            select_d{hr} = logical(S_unit(:,hr) | NS_unit(:,hr));
        end
        
        % difference for T and R d' in stim task (non-paired)
        [~,p_ttest_d_stim] = ttest2(dprime_abs(select_d{1},1),dprime_abs(select_d{2},2));
        % difference for d' in tuning task (non-paired)  不取绝对值差别不显著，取绝对值差别显著202108
        [~,p_ttest_d_tun] = ttest2(abs(dprime_all(methods_of_select{1,1},1)),abs(dprime_all(methods_of_select{1,1},2)));
        
        
        % whether variations in d' across T and R account for differences
        % in the efficacy of microstimulation across T and R: ANCOVA
        x_temp = [dprime_abs(:,1);dprime_abs(:,2)];
        y_temp = [PSE_shift_abs(:,1);PSE_shift_abs(:,2)];
        g_temp = [ones(size(dprime_abs(:,1)));2*ones(size(dprime_abs(:,2)))];
        
        % remove nan
        x = x_temp(~isnan(y_temp));
        y = y_temp(~isnan(y_temp));
        g = g_temp(~isnan(y_temp));
        
        % [h,a,c,s] = aoctool(x,y, group,alpha,xname,yname,gname,displayopt,model)  
        [h,atab,ctab,stats] = aoctool(x,y,g,0.05,'d prime','PSE shift','Stimulus type','off');  
        d_p = atab{3,6}; % across stimulus tyoe, d' and PSE shift (covariate effect)
        d_stimulus_p = atab{4,6}; % interaction between d' and stimulus type
        stimulus_p = atab{2,6}; % interaction between stimulus type and PSE shift 
        
        % ------------------------ Plot ------------------------------
        % plot without log scale
        % remove the outlier of PSE shift
        % do not remove the outlier
        set(figure(figN),'Position', [200,200 1000,500], 'Name', 'Relationship between PSE shift and d prime (normal scale)');clf
        sy{1} = 'o';
        sy{2} = '^';
        for hr = 1:2
            subplot(1,2,hr)
            hold on
            plot(dprime_abs(NS_unit(:,hr) & ~outlier(:,hr),hr),PSE_shift_abs(NS_unit(:,hr) & ~outlier(:,hr),hr),sy{hr},'markersize',8);
            plot(dprime_abs(S_unit(:,hr) & ~outlier(:,hr),hr),PSE_shift_abs(S_unit(:,hr) & ~outlier(:,hr),hr),sy{hr},'markersize',8,'color','k','markerfacecolor',c{hr});
            xlabel('|d prime|');
            ylabel('|Induced PSE shift (°)|');
            axis tight
            xlimit = xlim;
            ylimit = ylim;
            xlim([0 xlimit(2)]);
            ylim([0 ylimit(2)]);
            [r,p] = plot_corr_line(dprime_abs(~outlier(:,hr),hr),PSE_shift_abs(~outlier(:,hr),hr),'MethodOfCorr','Pearson','FittingMethod',2,'LineStyles',{'-','color',c{hr}});
            text(xlimit(1)+(xlimit(2)-xlimit(1))/20,ylimit(2)-(ylimit(2)-ylimit(1))/20,sprintf('r=%0.3g,p=%0.3g',r,p));
            text(xlimit(1)+(xlimit(2)-xlimit(1))/20,ylimit(2),sprintf('Remove outlier (>3*sigma): %g sites',sum(outlier(:,hr))));
            axis square
            xlimit = xlim;
            ylimit = ylim;
            xlim([-xlimit(2)/20 xlimit(2)*1.1]);
            ylim([-ylimit(2)/20 ylimit(2)*1.1]);
        end
        SetFigure(13);
        [til] = SetTitle(stimtask_type,monkey_included_for_analysis);
        str = ['Relationship between PSE shift and d prime' newline til];
        suptitle(str);
        figN = figN+1;
        
        
        % plot with log10 scale, do not remove the outlier
        plot_together = 1; % plot T and R together; =0: separate
        % scatter plot
        set(figure(figN),'Position', [50,50 1000,1200], 'Name', 'Relationship between PSE shift and d prime');clf
        for hr = 1:2
            if plot_together
                ax(hr) = subplot(4,4,[1 2 5 6]); % plot together
            else
                ax(hr) = subplot(4,4,[hr*2-1 hr*2 hr*2+3 hr*2+4]); % plot seperate
            end
            
            loglog(dprime_abs(NS_unit(:,hr),hr),PSE_shift_abs(NS_unit(:,hr),hr),sy{hr},'markersize',8,'color','k');
            hold on
            loglog(dprime_abs(S_unit(:,hr),hr),PSE_shift_abs(S_unit(:,hr),hr), sy{hr},'markersize',8,'color','k','markerfacecolor',c{hr});
            
            xlabel('|d prime|');
            ylabel('|Induced PSE shift (°)|');
            axis tight
            axis square
            set(gca,'xtickmode','manual');
            
            xx = dprime_abs((NS_unit(:,hr)|S_unit(:,hr)),hr);
            yy = PSE_shift_abs((NS_unit(:,hr)|S_unit(:,hr)),hr);
            nan_this = any(isnan([xx,yy]),2);
            xx = xx(~nan_this);
            yy = yy(~nan_this);
            
            ax1(hr) = subplot(4,4,[3 4 7 8]);
            plot(log10(xx),log10(yy),sy{hr},'color','k','markersize',8);
            axis tight
            axis square
            hold on
            error_ellipse(cov(log10(xx)',log10(yy)'),'conf',0.95,'mu',[mean(log10(xx)) mean(log10(yy))]); % confidence ellipses for data
            title('95% confidence ellipses');
            
        end
        linkaxes([ax(1) ax(2)]); % 同步坐标轴
        linkaxes([ax1(1) ax1(2)]); % 同步坐标轴

        xlimit_temp = ax(1).XLim;
        ylimit_temp = ax(1).YLim;
        xlimit = [xlimit_temp(1)/2 xlimit_temp(2)*2.4];
        ylimit = [ylimit_temp(1)/2 ylimit_temp(2)*2.4];
        set(ax(1),'XLim', [min([xlimit ylimit]) max([xlimit ylimit])]);
        set(ax(1),'YLim', [min([xlimit ylimit]) max([xlimit ylimit])]);
        
        % same set in log scale
        xylimit = [min([xlimit ylimit]) max([xlimit ylimit])];
        xylimitlog = log10(xylimit);
        set(ax1(1),'XLim', xylimitlog);
        set(ax1(1),'YLim', xylimitlog);
        
        
        % show tuning and ANCOVA result
        str = ['T-test for |d prime| tuning task:' newline ...
            sprintf('p = %.3g', p_ttest_d_tun) newline ...
            'ANCOVA result:' newline sprintf('d prime vs PSE shift: p = %.3g',d_p) newline ...
            sprintf('stimulus type vs PSE shift: p = %.3g',stimulus_p) newline...
            sprintf('d prime vs stumulus type: p = %.3g',d_stimulus_p)];
        text(-3,1.5,str);
        
        
        % mamrginal distributions
        % d' distributions
        nbins = 16;
        d_edges = logspace(log10(xlimit(1)),log10(xlimit(2)),nbins);
        % find the midpoint of edges
        for i = 1:nbins-1
            temp = logspace(log10(d_edges(i)),log10(d_edges(i+1)),3);
            x_this(i) = temp(2);
        end
        
        for hr = 1:2
            if plot_together
                ax(hr+2) = subplot(4,4,[9 10]);
            else
                ax(hr+2) = subplot(4,4,[hr*2+7 hr*2+8]);
            end
            y_this = histcounts(dprime_abs(select_d{hr},hr),d_edges, 'Normalization', 'probability');
            % plot 折线图
            plot(x_this,y_this,'-','color',c{hr});
            set(gca,'XScale','log');
            axis tight
            xlim(xlimit);
            
            % geometric mean
            g_mean = geomean(dprime_abs(select_d{hr},hr));
            hold on
            plot(g_mean,max(y_this),'v','color',c{hr});
            text(0.01,0.1+0.1*hr,sprintf('geomean = %.2g',g_mean),'color',c{hr});
        end
        linkaxes([ax(3) ax(4)],'y'); % 同步x坐标轴
        linkaxes([ax(1) ax(2) ax(3) ax(4)],'x'); % 同步x坐标轴
        
        % show ttest for d'
        text(10^-1,max(y_this),sprintf('t test, p = %0.3g',p_ttest_d_stim));
        title('mamrginal distribution of d prime')
        
        
        % PSE distribution
        nbins = 16;
        d_edges = logspace(log10(ylimit(1)),log10(ylimit(2)),nbins);
        % find the midpoint of edges
        for i = 1:nbins-1
            temp = logspace(log10(d_edges(i)),log10(d_edges(i+1)),3);
            x_this(i) = temp(2);
        end
        
        for hr = 1:2
            if plot_together
                ax(hr+4) = subplot(4,4,[13 14]);
            else
                ax(hr+4) = subplot(4,4,[hr*2+11 hr*2+12]);
            end
            y_this = histcounts(PSE_shift_abs(:,hr),d_edges, 'Normalization', 'probability');
            
            % plot 折线图
            plot(x_this,y_this,'-','color',c{hr});
            set(gca,'XScale','log');
            axis tight
            xlim(ylimit);
            
            % geometric mean
            g_mean = geomean(PSE_shift_abs(~isnan(PSE_shift_abs(:,hr)),hr));
            hold on
            plot(g_mean,max(y_this),'v','color',c{hr});
            text(0.01,0.01+0.01*hr,sprintf('geomean = %.2g',g_mean),'color',c{hr});
        end
        linkaxes([ax(5) ax(6)]); % 同步y坐标轴
        linkaxes([ax(1) ax(2) ax(5) ax(6)],'x'); % 同步x坐标轴
        title('mamrginal distribution of PSE shift')
        
        SetFigure(13);
        [til] = SetTitle(stimtask_type,monkey_included_for_analysis);
        str = ['Relationship between PSE shift and d prime' newline til];
        suptitle(str);
        figN = figN+1;
    end

    function f3p13(debug) %'PSE shift and some T-R ratio'
        if debug  ; dbstack;   keyboard;      end
        
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        effect_type = 1; % 1: PSE, 2: threshold
        select_stars_type = [0 1]; % 0: dot stimulus; 1: triangle stimulus
        is_separate_stimAmp = 0;
        is_original_bias = 4; % original bias and normalized (after-before  >0: bias to left/CW; <0: bias to right/CCW) (后面取绝对值，所以没关系)
        
         % ---------------------- get data ----------------------------
        [PSE_shift_temp, S_unit_temp, NS_unit_temp, ~] = get_stim_data_uA(stimtask_type,methods_of_select{1,1},is_separate_stimAmp,is_original_bias,select_stars_type);
        PSE_shift = PSE_shift_temp{effect_type};
        PSE_shift_abs = abs(PSE_shift);
        S_unit = S_unit_temp{effect_type};
        NS_unit = NS_unit_temp{effect_type};
        
        keyboard
    end

    function f3p16(debug) %'PSE shift and CR'
        if debug  ; dbstack;   keyboard;      end
        
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        effect_type = 1; % 1: PSE, 2: threshold
        select_stars_type = [0 1]; % 0: dot stimulus; 1: triangle stimulus
        is_separate_stimAmp = 0;
        is_original_bias = 3; % normalize PSE shift
        
        % ---------------------- get data ----------------------------
        % PSE shicft
        [PSE_shift_temp, S_unit_temp, NS_unit_temp, ~] = get_stim_data_uA(stimtask_type,methods_of_select{1,1},is_separate_stimAmp,is_original_bias,select_stars_type);
        PSE_shift = PSE_shift_temp{effect_type};
        S_unit = S_unit_temp{effect_type};
        NS_unit = NS_unit_temp{effect_type};
        
        % CR, CCI (参考函数get_stim_data_uA)
        ctrl_CR_all = [];
        stim_CR_all = [];
        
        CCI_all = [];
        chi_square_all = [];
        for fc = 1:length(stimtask_type) % fine or coarse，需要combine一起，所以最后的index没有fr
            if ~isempty(stimtask_type{fc})
                for i = 1:length(stimtask_type{fc}) % 具体task, 需要combine一起，所以最后的index没有i
                    select = logical(methods_of_select{1,1}); % select for stimAmp, stimtask_type_this and stars_type
                    temp1 = []; temp2 = []; temp3 = []; temp4 = [];
                    temp1 = ctrl_CR{stimtask_type{fc}(i)};
                    temp2 = stim_CR{stimtask_type{fc}(i)};
                    temp3 = choice_change_pro_around0{stimtask_type{fc}(i),3}(:,1); % only use translation CCI: (:,3)
                    temp4 = chi_square_p_around0{stimtask_type{fc}(i)}(:,3);
                    
                    temp1(~select,:) = nan; % remove undesired unit
                    temp2(~select,:) = nan; % remove undesired unit
                    temp3(~select,:) = nan; % remove undesired unit
                    temp4(~select,:) = nan; % remove undesired unit
                    
                    ctrl_CR_all = [ctrl_CR_all; temp1]; % combine selected stimtask_type_this
                    stim_CR_all = [stim_CR_all; temp2]; % combine selected stimtask_type_this
                    CCI_all = [CCI_all; temp3];
                    chi_square_all = [chi_square_all; temp4]; % p value
                end
            end
        end
        CR_ratio = stim_CR_all./ctrl_CR_all;
        
        x_this{1} = [CR_ratio CR_ratio(:,3)]; % add one column for CCI
        x_this{2} = [ctrl_CR_all ctrl_CR_all(:,3)];
        x_this{3} = [stim_CR_all stim_CR_all(:,3)];
        
        y_this = [PSE_shift CCI_all];

        % ------------------------ Plot ------------------------------
        set(figure(figN),'Position', [200,200 1500,1000], 'Name', 'Relationship between PSE shift and CR');clf
        for hr = 1:4 % t,r,mt, CCI
            for i = 1:3 % ratio, ctrl, stim CR
                subplot(3,4,hr+(i-1)*4)
                plot(x_this{i}(:,hr),y_this(:,hr),'k.');
                ylimit = ylim;
                ylimit_max = max(ylimit);
                ylim([-ylimit_max ylimit_max]);
                hold on
                if i == 1
                    xlim([0.5 1.5]);
                    plot([1 1],[-ylimit_max ylimit_max],'k-');
                    plot([0.5 1.5],[0 0],'k-');
                    xlabel('CR ratio (Stim/Ctrl)');
                else
                    xlim([40 100]);
                    plot([50 50],[-ylimit_max ylimit_max],'k-');
                    plot([40 100],[0 0],'k-');
                    if i == 2
                        xlabel('control CR');
                    elseif i == 3
                        xlabel('stim. CR');
                    end
                end
                xlimit = xlim;
                
                
                num_this = sum(sum([~isnan(x_this{i}(:,hr)) ~isnan(y_this(:,hr))],2)==2);
                text(xlimit(2)-xlimit(2)/10, ylimit_max, sprintf('n = %g',num_this));
                
                % spearman correlation
                [r,p] = plot_corr_line(x_this{i}(:,hr),y_this(:,hr),'MethodOfCorr','Spearman','FittingMethod',1,'LineStyles','k-');
                str = [sprintf('r=%.2g',r) newline sprintf('p=%.2g',p)];
                text(xlimit(1),ylimit_max-ylimit_max/10,str);
                
                axis square;
                if hr == 1
                    ylabel('PSE shift');
                    title('Translation PSE');
                elseif hr == 2
                    ylabel('PSE shift');
                    title('Roll PSE');
                elseif hr == 3
                    ylabel('PSE shift');
                    title('Flow-pattern PSE');
                elseif hr == 4
                    ylabel('CCI');
                    title('Flow-pattern CCI');
                end
            end
        end
        SetFigure(13);
        [til] = SetTitle(stimtask_type,monkey_included_for_analysis);
        str = ['Relationship between PSE shift and CR' newline til];
        suptitle(str);
        figN = figN+1;

    end

    function f3p18(debug) % 'PSE shift in normal psy and general error level'
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        effect_type = 1; % 1: PSE, 2: threshold
        select_stars_type = [0 1]; % 0: dot stimulus; 1: triangle stimulus
        is_separate_stimAmp = 0;
        
        % ---------------------- get data ----------------------------
        % PSE shift
        is_original_bias = 0; % normal psy
        PSE_shift_temp = []; S_unit_temp = []; NS_unit_temp = [];
        [PSE_shift_temp, S_unit_temp, NS_unit_temp, ~] = get_stim_data_uA(stimtask_type,methods_of_select{1,1},is_separate_stimAmp,is_original_bias,select_stars_type);
        PSE_shift = PSE_shift_temp{effect_type};
        S_unit = S_unit_temp{effect_type};
        NS_unit = NS_unit_temp{effect_type};
        
        is_original_bias = 5; % general error PSE
        PSE_shift_temp = []; S_unit_temp = []; NS_unit_temp = [];
        [PSE_shift_temp, S_unit_temp, NS_unit_temp, ~] = get_stim_data_uA(stimtask_type,methods_of_select{1,1},is_separate_stimAmp,is_original_bias,select_stars_type);
        PSE_shift_true = PSE_shift_temp{effect_type};
        S_unit_true = S_unit_temp{effect_type};
        NS_unit_true = NS_unit_temp{effect_type};

        % ------------------------ Plot ------------------------------
        set(figure(figN),'Position', [200,200 800,700], 'Name', 'Normal PSE vs General error level PSE');clf
        for hr = 1:2
            subplot(1,2,hr)
            plot(PSE_shift(:,hr),PSE_shift_true(:,hr),'ko');
            xlabel('Normal PSE');
            ylabel('General error level PSE');
            if hr == 1
                title('Translation');
            else
                title('Roll');
            end
            
            xtick_this = properXtick([PSE_shift(:,hr);PSE_shift_true(:,hr)]);
            xlim([xtick_this(1) xtick_this(end)]);
            ylim([xtick_this(1) xtick_this(end)]);
            xticks(xtick_this);
            yticks(xtick_this);
            
            hold on
            plot([xtick_this(1) xtick_this(end)],[0 0],'k-');
            plot([0 0], [xtick_this(1) xtick_this(end)],'k-');
            SameScale(1,1)
            
            % correlation
            [r,p] = corr(PSE_shift(:,hr), PSE_shift_true(:,hr),'type','Pearson','rows','complete');
            
            str = [sprintf('n = %g',sum(~any(isnan([PSE_shift(:,hr) PSE_shift_true(:,hr)]),2))) newline 'Pearson' newline ...
                sprintf('r = %2.2g (p = %2.2g)',r,p)];
            text(xtick_this(1),xtick_this(end),str);
            
            
        end
        suptitle(SetTitle(stimtask_type,monkey_included_for_analysis));
        SetFigure(15)
        figN = figN+1;
    end

%% Motion type
    function f4p1(debug) % 'Motion type PSE shift with Modulation index',
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        % 默认只选4 target task，这样才有motion type
        % *****************************************
        stimtask_type_this{1} = [1];
        stimtask_type_this{2} = [4];
        select_stars_type = [0 1]; % 0: dot stimulus; 1: triangle stimulus
        is_separate_stimAmp = 0; % combine 20 and 40 stim amplitude
        is_original_bias = 2; % original bias: after-before  >0: bias to left/CW; <0: bias to right/CCW
        effect_type = 1; % bias
        
        % ---------------------- get data ----------------------------
        for fc = 1:2
            [diff_temp, S_unit_temp, NS_unit_temp, ~] = get_stim_data_uA(stimtask_type_this(fc),methods_of_select{1,1},is_separate_stimAmp,is_original_bias,select_stars_type);
            mt_bias{fc} = diff_temp{effect_type}(:,3); % >0: bias to H; <0 bias to R
            S_unit_ori{fc} = S_unit_temp{effect_type}(:,3);
            NS_unit_ori{fc} = NS_unit_temp{effect_type}(:,3);
        end
        
        % Modulatio index: <0: pref R; >0 pref H
        modulation_index{1} = [group_result.ahead_slope_index]'; % 斜率绝对值
        modulation_index{2} = [group_result.d_prime_index_fine]';
        modulation_index{3} = [group_result.d_prime_index_coarse]';
        modulation_index{4} = [group_result.d_prime_maxDiff_index]'; % 差值最大的两个点,在Group_SpiT中去掉>6std的异常值
        modulation_index{5} = [group_result.d_prime_maxAxis_index]'; % 差值最大的一条轴（180相对两点）,在Group_SpiT中去掉>6std的异常值
        modulation_index{6} = [group_result.d_prime_maxFR_index]'; % FR最大点及其对侧点,在Group_SpiT中去掉>6std的异常值
        modulation_index{7} = [group_result.DDI_inc_index]';
        modulation_index{8} = [group_result.DDI_exc_index]';
        modulation_index{9} = [group_result.HTI_index]';
        modulation_index{10} = [group_result.neu_sen_index_fine]';
        modulation_index{11} = [group_result.neu_sen_index_coarse]';
        modulation_index{12} = [group_result.spiral_index]';
        
        d_prime_coarse = cell2mat({group_result.d_prime_coarse}');
        modulation_index{13} = abs(d_prime_coarse(:,1)); % heading d prime coarse
        modulation_index{14} = abs(d_prime_coarse(:,2)); % rotation d prime coarse
        
        neu_sen_coarse = cell2mat({group_result.neu_sen_coarse}');
        modulation_index{15} = neu_sen_coarse(:,1); % heading neural sensitivity
        modulation_index{16} = neu_sen_coarse(:,2); % rotation neural sensitivity
        
        % onlyt select 'Typical cells' (默认, all cell type: t cell, r cell and spiral cell)
        methods_of_select = repmat(methods_of_select{1,1},length(stimtask_type_this{1}),1);
        mt_bias = cellfun(@(x) x(methods_of_select,:),mt_bias,'uniformoutput',0);
        S_unit_ori = cellfun(@(x) x(methods_of_select,:),S_unit_ori,'uniformoutput',0);
        NS_unit_ori = cellfun(@(x) x(methods_of_select,:),NS_unit_ori,'uniformoutput',0);
        modulation_index = cellfun(@(x) x(methods_of_select,:),modulation_index,'uniformoutput',0);
        
        % ------------------------ Plot ------------------------------
        tl{1} = 'Zero slope index';
        tl{2} = strcat('d''','index(fine)');
        tl{3} = strcat('d''','index(coarse)');
        tl{4} = strcat('d''','index(maxDiff)');
        tl{5} = strcat('d''','index(maxAxis)');
        tl{6} = strcat('d''','index(maxFR)');
        tl{7} = 'DDI index(inc)';
        tl{8} = 'DDI index(exc)';
        tl{9} = 'HTI index';
        tl{10} = 'Neural Sens.Index(fine)';
        tl{11} = 'Neural Sens.Index(coarse)';
        tl{12} = 'Spiral Index';
        tl{13} = strcat('|Heading','d''','(coarse)|');
        tl{14} = strcat('|Rotation','d''','(coarse)|');
        tl{15} = 'Heading neural sens.';
        tl{16} = 'Rotation neural sens.';
        
        for fc = 1:2 % fine or coarse
            if fc == 1
                set(figure(figN-1+fc),'name','Fine Motion type PSE shift and Modulation index','pos',[150,150,1700,850]); clf
            else
                set(figure(figN-1+fc),'name','Coarse Motion type PSE shift and Modulation index','pos',[50,50,1700,850]); clf
            end
            ha = tight_subplot(3, 6,[.09 .04], [.07 .1], [.08 .05]);
            for nn1 = 1:length(modulation_index)
                axes(ha(nn1))
                % mt bias: % >0: bias to R; <0 bias to H
                % modulatio index: <0: pref R; >0 pref H (除了后4个index)
                % 理想情况是反对角线
                try
                    plot(modulation_index{nn1}(S_unit_ori{fc}),mt_bias{fc}(S_unit_ori{fc}),'ko','markerfacecolor',c{3});
                catch
                    keyboard
                end
                hold on
                plot(modulation_index{nn1}(NS_unit_ori{fc}),mt_bias{fc}(NS_unit_ori{fc}),'ko');
                
                % Pearson correlation
                [r,p] = plot_corr_line(modulation_index{nn1},mt_bias{fc},'MethodOfCorr','Pearson','FittingMethod',2,'LineStyles','k-');
                [r_sp,p_sp] = plot_corr_line(modulation_index{nn1},mt_bias{fc},'MethodOfCorr','Spearman','FittingMethod',2,'LineStyles','k-');
                
                % annotation
                xlimit = xlim; ylimit = ylim;
                ylimit = max(abs(ylimit));
                ylim([-ylimit ylimit]);
                ylimit = ylim;
                if nn1 < 13
                    plot([0 0],ylimit,'k:'); plot(xlimit,[0 0],'k:');
                end
                text(xlimit(1)+(xlimit(2)-xlimit(1))/20,ylimit(2)-(ylimit(2)-ylimit(1))/20,sprintf('r=%0.3g,p=%0.3g',r,p),'fontsize',8.5);
                text(xlimit(1)+(xlimit(2)-xlimit(1))/20,ylimit(2)-(ylimit(2)-ylimit(1))/10,sprintf('r=%0.3g,p=%0.3g',r_sp,p_sp),'fontsize',8.5);
                
                if nn1 == 1
                    text(xlimit(1)+(xlimit(2)-xlimit(1))/20,ylimit(1)-(ylimit(2)-ylimit(1))/7,'Pref R','fontsize',10);
                    text(xlimit(2)-(xlimit(2)-xlimit(1))/4,ylimit(1)-(ylimit(2)-ylimit(1))/7,'Pref H','fontsize',10);
                    str = strcat('Bias H',repmat(32,1,15),'Bias R');
                    text(xlimit(1)-(xlimit(2)-xlimit(1))/4.5,ylimit(1)+(ylimit(2)-ylimit(1))/8,str,'rotation',90);
                end
                
                if nn1 == 13 || nn1 == 15
                    xlabel('More Heading');
                elseif nn1 == 14 || nn1 == 16
                    xlabel('More Rotation');
                end
                
                title(tl{nn1});
            end
            SetFigure(10)
            figN = figN + 2;
        end
    end

    function f4p12(debug) %'Motion type PSE shift with Spiral index'
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        select_stars_type = [0 1]; % 0: dot stimulus; 1: triangle stimulus
        is_separate_stimAmp = 0; % combine 20 and 40 stim amplitude
        
        effect_type = 1; % bias
        
        % ---------------------- get data ----------------------------
        %         % must select 4-target task
        %         if sum(stimtask_type{1}==1)+sum(stimtask_type{2}==4)==0 % fine 4T和coarse 4T都没选
        %             disp('********** Please choose 4-Target task! **********');
        %             return
        %         end
        %         stimtask_type_this = stimtask_type;
        %         stimtask_type_this{1}(stimtask_type_this{1}~=1) = []; % only left 4-Target task
        %         stimtask_type_this{2}(stimtask_type_this{2}~=4) = [];
        
        % auto select
        stimtask_type_this{1} = 1; % 4 targets fine task
        stimtask_type_this{2} = 4; % 4 targets coarse task
        
        % 0: bias positive = prefer direction; 1: only original bias in MT;  2: original bias(after-before); 3： normalize PSE; 4 normalize original bias
        is_original_bias = 2; % original bias: after-before  >0: bias to left/CW/R; <0: bias to right/CCW/T
        [diff_temp, S_unit_temp, NS_unit_temp, ~] = get_stim_data_uA(stimtask_type_this,methods_of_select{1,1},is_separate_stimAmp,is_original_bias,select_stars_type);
        mt_bias = diff_temp{effect_type}(:,3); % >0: bias to R; <0 bias to T
        S_unit = S_unit_temp{effect_type}(:,3);
        NS_unit = NS_unit_temp{effect_type}(:,3);
        
        spiral_index = [group_result.spiral_index]';
        %         spiral_index = repmat(spiral_index,length(cell2mat(stimtask_type_this)),1);
        %         select = repmat(methods_of_select{1,1},length(cell2mat(stimtask_type_this)),1);
        
        % behavioral performance
        do_stimtask_type = cellfun(@(x) ~isempty(x),stimtask_type_this);
        mt_cr = [];
        for fc = 1:length(stimtask_type_this)
            mt_cr(:,fc) = ctrl_CR{stimtask_type_this{fc}}(:,3); % control task CR
        end

        % ------------------------ Plot ------------------------------
        % select with behavioral performance
%         mt_performance_cri{1} = [60 100]; % only select performance > 65% and < 100% correct rate
        mt_performance_cri{1} = [65 100];
%         mt_performance_cri{3} = [70 100];
%         mt_performance_cri{4} = [0 100];
        
        for ps = 1:length(mt_performance_cri)
            set(figure(figN),'name','Motion type PSE shift and Spiral index','pos',[100+ps*100,150,600,1200]); clf
            for fc = 1:2 % fine or coarse
                if fc == 1
                    PSEshift = mt_bias(1:length(loadM_data));
                    NS_this = NS_unit(1:length(loadM_data));
                    S_this = S_unit(1:length(loadM_data));
                elseif fc == 2
                    PSEshift = mt_bias(length(loadM_data)+1:end);
                    NS_this = NS_unit(length(loadM_data)+1:end);
                    S_this = S_unit(length(loadM_data)+1:end);
                end
                perform_select = logical(mt_cr(:,fc)>=mt_performance_cri{ps}(1) & mt_cr(:,fc)<=mt_performance_cri{ps}(2)); % select control task performance
                select = logical(methods_of_select{1,1} & perform_select);
                
                subplot(2,1,fc)
                plot(spiral_index(NS_this & select),PSEshift(NS_this & select),'ko');
                hold on
                plot(spiral_index(S_this & select),PSEshift(S_this & select),'ko','markerfacecolor',c{3});
                xlim([-1 1]);
                ylim([-1 1]);
                axis square
                
                ylimit = ylim;
                [r_p,p_p] = plot_corr_line(spiral_index(select),PSEshift(select),'MethodOfCorr','Pearson','FittingMethod',2,'LineStyles',{'-','color','k'});
                [r_s,p_s] = plot_corr_line(spiral_index(select),PSEshift(select),'MethodOfCorr','Spearman','FittingMethod',2,'LineStyles',{'-','color',c{99}});
                text(-0.9,ylimit(2)/1.1,sprintf('Pearson correlation: r=%.3g, p=%.3g',r_p,p_p));
                text(-0.9,ylimit(2)/1.2,sprintf('Spearman correlation: r=%.3g, p=%.3g',r_s,p_s));
                
                xstr = ['More Rotation',repmat(32,1,20),'More Translation' newline 'Spiral Index'];
                xlabel(xstr);
                ystr = ['Motion type PSE shift' newline 'Bias to Translation',repmat(32,1,20),'Bias to Rotation'  ];
                ylabel(ystr);
                
                text(-0.9,ylimit(2)-ylimit(2)/3,sprintf('CR = [%g-%g]%%, N = %g',mt_performance_cri{ps}(1),mt_performance_cri{ps}(2), sum(select)));

                if fc == 1
                    title('Fine');
                elseif fc == 2
                    title('Coarse'); 
                end
                
            end
            suptitle(SetTitle(stimtask_type,monkey_included_for_analysis));
            SetFigure(10);
            figN = figN+1;
        end
    end

    function f4p13(debug) % 'Motion type PSE shift with direction preference'
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        select_stars_type = [0 1]; % 0: dot stimulus; 1: triangle stimulus
        is_separate_stimAmp = 0; % combine 20 and 40 stim amplitude
        
        effect_type = 1; % bias
        
        % ---------------------- get data ----------------------------
        % must select 4-target task
        if sum(stimtask_type{1}==1)+sum(stimtask_type{2}==4)==0 % fine 4T和coarse 4T都没选
            disp('********** Please choose 4-Target task! **********');
            return
        end
        stimtask_type_this = stimtask_type;
        stimtask_type_this{1}(stimtask_type_this{1}~=1) = []; % only left 4-Target task
        stimtask_type_this{2}(stimtask_type_this{2}~=4) = [];
        
        if length(cell2mat(stimtask_type_this))==1 % only select fine or coarse, 不使用normalize PSE shift
            is_original_bias = 2; % original bias: after-before  >0: bias to left/CW; <0: bias to right/CCW
        else
            is_original_bias = 4; % normalize original PSE
        end
        [diff_temp, S_unit_temp, NS_unit_temp, ~] = get_stim_data_uA(stimtask_type_this,methods_of_select{1,1},is_separate_stimAmp,is_original_bias,select_stars_type);
        mt_bias = diff_temp{effect_type}(:,3); % >0: bias to T; <0 bias to R
        S_unit = S_unit_temp{effect_type}(:,3);
        NS_unit = NS_unit_temp{effect_type}(:,3);
        
        % preferred direction
        pref_direc_temp = cell2mat({group_result.pref_direc}');
        pref_direc = repmat(pref_direc_temp,length(cell2mat(stimtask_type_this)),1); % 0 = exp
        
        % ------------------------ Plot ------------------------------
        set(figure(figN),'name','MT PSE shift affected by preferred direction','pos',[50,50,1200,600]); clf
        for hr = 1:2
            % scatter plot
            subplot(3,2,[hr+2 hr+4]);
            reasonable_range = [nanmean(mt_bias) - 3*nanstd(mt_bias) nanmean(mt_bias) + 3*nanstd(mt_bias)];  % mean±3*std
            outdot = logical(mt_bias<reasonable_range(1) | mt_bias>reasonable_range(2));
            outdot = false(size(outdot)); % show all
            
            plot(pref_direc(~outdot,hr),mt_bias(~outdot),'o','color',c{hr});
            hold on
            plot(pref_direc((S_unit & ~outdot),hr),mt_bias(S_unit & ~outdot),'o','color',c{hr},'markerfacecolor',c{hr});
            axis tight
            ylimit = ylim;
            ylimit = ylimit.*1.1;
            
            if sum(outdot)>0
                text(-170,ylimit(2),sprintf('%g dots out of 3*std did not show',sum(outdot)));
            end
            plot([-180 180],[0 0],'k:');
            plot([-90 -90],ylimit,'k:');
            plot([0 0],ylimit,'k:');
            plot([90 90],ylimit,'k:');
            xlimit = [-180,180];
            xlim(xlimit);
            xticks([-180:45:180]);
            if hr == 2
                xticklabels_temp = xticklabels;
                xticklabels_temp{[-180:45:180]==-90} = 'CW';
                xticklabels_temp{([-180:45:180]==90)} = 'CCW';
                xticklabels(xticklabels_temp);
            end
            if is_original_bias == 2
                ylabel_temp = 'Normalized Induced MT PSE shift (°)';
            elseif is_original_bias == 4
                ylabel_temp = 'Induced MT PSE shift (°)';
            end
            str = [ylabel_temp newline 'Bias to Rot    Bias to Tran'];
            ylabel(str);
            xlabel('Direction preference');
            
            % step mean
            ax(hr) = subplot(3,2,hr);
            bin = 30;
            step = 10;
            binL = xlimit(1);
            binR = binL+bin;
            n = 1;
            while binR<=xlimit(2)
                select = logical(pref_direc(:,hr)>=binL & pref_direc(:,hr)<binR);
                mean_pref(n) = (binL+binR)/2;
                if sum(select)~=0 && sum(~isnan(mt_bias(select)))>0 % select中有prefer该方向的细胞，且做了电刺激
                    mean_PSE_shift(n) = nanmean(mt_bias(select));
                    sem_PSE_shit(n) = nanmean(mt_bias(select))/sqrt(sum(select));
                    
                    % significant
                    [~,p(n)] = ttest(mt_bias(select),0);
                else
                    mean_PSE_shift(n) = nan;
                    sem_PSE_shit(n) = nan;
                    p(n) = nan;
                end
                binL = binL+step;
                binR = binR+step;
                n = n+1;
            end
            
            % fill the nan
            mean_PSE_shift_nan = interp1(mean_pref(~isnan(mean_PSE_shift)),mean_PSE_shift(~isnan(mean_PSE_shift)),mean_pref);
            hold on
            plot(mean_pref,mean_PSE_shift_nan,':','color',c{hr});
            errorbar(mean_pref,mean_PSE_shift,sem_PSE_shit,'ko-','markerfacecolor',c{hr},'markeredgecolor','none');
            
            hold on
            ylimit = ylim;
            plot([-180 180],[0 0],'k:');
            plot([-90 -90],ylimit,'k:');
            plot([0 0],ylimit,'k:');
            plot([90 90],ylimit,'k:');
            xlim(xlimit);
            xticks([-180:45:180]);
            if hr == 2
                xticklabels_temp = xticklabels;
                xticklabels_temp{[-180:45:180]==-90} = 'CW';
                xticklabels_temp{([-180:45:180]==90)} = 'CCW';
                xticklabels(xticklabels_temp);
            end
            str = ['Bias to Rot    Bias to Tran'];
            ylabel(str);
            
            % significant marker
            plot(mean_pref(p<0.05),mean_PSE_shift(p<0.05)+(ylimit(2)-ylimit(1))/10,'*','color',c{hr});
            ylim([ylimit(1)*1.1 ylimit(2)*1.1]);
        end
        linkaxes([ax(1),ax(2)],'xy'); % 同步多个坐标区的范围
        SetFigure(11)
        [til] = SetTitle(stimtask_type_this,monkey_included_for_analysis);
        str=['Correlation betweenM MT PSE shift and direction preference' newline til];
        suptitle(str)
        figN = figN+1;
        
        % average reuslt grouped by four relative preferred direction
        % ---------------------- get data ----------------------------
        select = [];
        for hr = 1:2
            select{1} = logical(abs(pref_direc(:,hr))<=45 & abs(pref_direc(:,hr))>=0); % 0-45
            select{2} = logical(abs(pref_direc(:,hr))<=90 & abs(pref_direc(:,hr))>45); % 45-90
            select{3} = logical(abs(pref_direc(:,hr))<=135 & abs(pref_direc(:,hr))>90); % 90-135
            select{4} = logical(abs(pref_direc(:,hr))<=180 & abs(pref_direc(:,hr))>135); % 135-180
% 
%             select{1} = logical(abs(pref_direc(:,hr))<45); % -45-45
%             select{2} = logical((pref_direc(:,hr))<=-45 & (pref_direc(:,hr))>=-135); % -45~-135
%             select{3} = logical((pref_direc(:,hr))<=135 & (pref_direc(:,hr))>=45); % 45~135
%             select{4} = logical(abs(pref_direc(:,hr))>135); % 135-180,-135~-180
            
            for i = 1:4
                PSE_shift_group{hr,i} = mt_bias(select{i}); % 0-45
                mean_PSE_group(hr,i) = nanmean(PSE_shift_group{hr,i});
                sem_PSE_group(hr,i) = nanstd(PSE_shift_group{hr,i}) / sqrt(sum(select{i}));
                
                % significant
                [~,p] = ttest(PSE_shift_group{i},0);
                if p < 0.001
                    sig_mark{hr,i} = '***';
                elseif p>=0.001 && p<0.01
                    sig_mark{hr,i} = '**';
                elseif p>=0.01 && p<0.05
                    sig_mark{hr,i} = '*';
                else
                    sig_mark{hr,i} = [];
                end
            end
        end
        
        % ------------------------ Plot ------------------------------
        set(figure(figN),'name','MT PSE shift grouped by four relative preferred direction','pos',[500,50,1200,600]); clf
        for hr = 1:2
            ax(hr) = subplot(1,2,hr);
            barwitherr(sem_PSE_group(hr,:), mean_PSE_group(hr,:),'facecolor',c{hr});
            hold on
            axis tight
            for i = 1:4
                if ~isempty(sig_mark{hr,i})
                    ylimit_this = mean_PSE_group(hr,i)+sem_PSE_group(hr,i)+mean_PSE_group(hr,i)/10;
                    text(i,ylimit_this,sig_mark{hr,i},'HorizontalAlignment','center');
                end
            end
            xticklabels({'0-45','45-90','90-135','135-180'});
            %             xticklabels({'Exp','CW/Left','CCW/Right','Con'});
            xlabel('|Preferred direction from forward| (°)');
            str = ['Bias to Rot    Bias to Tran'];
            ylabel(str);
            ylimit = ylim;
            ylim(ylimit.*1.1);
        end
        linkaxes([ax(1),ax(2)],'xy'); % 同步多个坐标区的范围, same y lim T and R
        SetFigure(11)
        [til] = SetTitle(stimtask_type,monkey_included_for_analysis);
        str=['PSE shift grouped by four relative preferred direction' newline til];
        suptitle(str)
        figN = figN+1;
    end

    function f4p17(debug) % 'Motion type PSE shift with T/R threshold change'
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        select_stars_type = [0 1]; % 0: dot stimulus; 1: triangle stimulus
        is_separate_stimAmp = 0; % combine 20 and 40 stim amplitude

        % ---------------------- get data ----------------------------
        % must select 4-target task
        if sum(stimtask_type{1}==1)+sum(stimtask_type{2}==4)==0 % fine 4T和coarse 4T都没选
            disp('********** Please choose 4-Target task! **********');
            return
        end
        stimtask_type_this = stimtask_type;
        stimtask_type_this{1}(stimtask_type_this{1}~=1) = []; % only left 4-Target task
        stimtask_type_this{2}(stimtask_type_this{2}~=4) = [];
        
        % is_original_bias: 0: bias positive = prefer direction; 1: only original bias in MT;  2: original bias(after-before); 3： normalize PSE; 4 normalize original bias
        if length(cell2mat(stimtask_type_this))==1 % only select fine or coarse, 不使用normalize PSE shift
            is_original_bias = 2; 
        else
            is_original_bias = 4; 
        end
        [diff_temp, S_unit_temp, NS_unit_temp, ~] = get_stim_data_uA(stimtask_type_this,methods_of_select{1,1},is_separate_stimAmp,is_original_bias,select_stars_type);
        mt_bias = diff_temp{1}(:,3); % >0: bias to T; <0 bias to R
        S_unit = S_unit_temp{1}(:,3);
        NS_unit = NS_unit_temp{1}(:,3);
        
        T_threshold = diff_temp{2}(:,1);
        TS_unit = S_unit_temp{2}(:,1);
        TNS_unit = NS_unit_temp{2}(:,1);
        
        R_threshold = diff_temp{2}(:,2);
        RS_unit = S_unit_temp{2}(:,2);
        RNS_unit = NS_unit_temp{2}(:,2);
        
        
        % ------------------------ Plot ------------------------------
        set(figure(figN),'name','Motion type PSE shift with T/R threshold change','pos',[500,50,1200,600]); clf
        select1 = logical(T_threshold>0 & TS_unit);
        select2 = logical(R_threshold>0 & RS_unit);
        select3 = logical(select1 & select2);
        select4 = logical(~select1 & ~select2);
        hold on
%         grouped = [ones(size(mt_bias(select1)));ones(size(mt_bias(select2)))*2;ones(size(mt_bias(select3)))*3;ones(size(mt_bias(select4)))*4];
%         boxplot([mt_bias(select1); mt_bias(select2); mt_bias(select3); mt_bias(select4)],grouped)
        
        plot(1,mt_bias(select1),'ko');
        plot(1,mt_bias(select1 & S_unit),'ko','markerfacecolor','k');
        plot(2,mt_bias(select2),'ko');
        plot(2,mt_bias(select2 & S_unit),'ko','markerfacecolor','k');
        plot(3,mt_bias(select3),'ko');
        plot(3,mt_bias(select3 & S_unit),'ko','markerfacecolor','k');
        plot(4,mt_bias(select4),'ko');
        plot(4,mt_bias(select4 & S_unit),'ko','markerfacecolor','k');
        figN = figN+1;
    end


    function f4p14(debug) %'Choice changed across direction preference'
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        select_stars_type = [0 1]; % 0: dot stimulus; 1: triangle stimulus
        is_separate_stimAmp = 0; % combine 20 and 40 stim amplitude
        effect_type = 1; % bias
        
        % stimtask_type, always choose 4-target fine and coarse task
        stimtask_type_this{1} = 1; % 4 targets fine task
        stimtask_type_this{2} = 4; % 4 targets coarse task
        
        % ---------------------- get data ----------------------------
        select = logical(methods_of_select{1,1});
        
        % choice change
        % choice_change_pro为正表示该motion type choice在电刺激后比电刺激前选择增增加了，两列分别是T和R，proportion = (stim-ctrl)/ctrl   详细见CueSwitch_performance.m
        for fc = 1:2
            ind = 3;
            choice_change_pro_diff_4T{fc} = choice_change_pro_diff{stimtask_type_this{fc},ind}; % 两列分别是T和R
            choice_change_pro_in0_4T{fc} = choice_change_pro_in0{stimtask_type_this{fc},ind}; % 两列分别是T和R
            choice_change_pro_around0_4T{fc} = choice_change_pro_around0{stimtask_type_this{fc},ind}; % 两列分别是T和R
            
            chi_square_p_diff_4T{fc} = chi_square_p_diff{stimtask_type_this{fc}}(:,ind) ;
            chi_square_p_in0_4T{fc} = chi_square_p_in0{stimtask_type_this{fc}}(:,ind);
            chi_square_p_around0_4T{fc} = chi_square_p_around0{stimtask_type_this{fc}}(:,ind);
        end
        
        % select and group together
        choice_change_proportion{1} = cellfun(@(x) x(select,:),choice_change_pro_diff_4T,'uniformoutput',0);
        choice_change_proportion{2} = cellfun(@(x) x(select,:),choice_change_pro_in0_4T,'uniformoutput',0);
        choice_change_proportion{3} = cellfun(@(x) x(select,:),choice_change_pro_around0_4T,'uniformoutput',0);
        
        chi_square_p{1} = cellfun(@(x) x(select),chi_square_p_diff_4T,'uniformoutput',0);
        chi_square_p{2} = cellfun(@(x) x(select),chi_square_p_in0_4T,'uniformoutput',0);
        chi_square_p{3} = cellfun(@(x) x(select),chi_square_p_around0_4T,'uniformoutput',0);
        
        % preferred direction
        pref_direc_temp = cell2mat({group_result.pref_direc}');
        pref_direc = pref_direc_temp(select,:);
        
        % ------------------------ Plot ------------------------------
        for j = 1:3 % difficult trial, 0 degree, around 0 degree
            set(figure(figN),'name','Translation and Rotation choice changed across preferred direction','pos',[600*j-600,50,630,950]); clf
            ha = tight_subplot(3,2,[.1 .1],[.07 0.07],[.2 .05],1);
            for i = 1:2 % R or T (第一列是R chocie变化，第二列是T choice变化）
                for fc = 1:3 % 1:only fine, 2:only coatse, 3: combine
                    % scatter plot
                    axes(ha(i+(fc-1)*2))
                    if fc~=3
                        PD = pref_direc(:,3-i); % first column is R
                        choice_change = choice_change_proportion{j}{fc}(:,i); % first column is R
                        sig_site = logical(chi_square_p{j}{fc}<0.05);
                    else % combine fine and coarse
                        PD = [pref_direc(:,3-i);pref_direc(:,3-i)]; % first column is R
                        choice_change = [choice_change_proportion{j}{1}(:,i);choice_change_proportion{j}{2}(:,i)];
                        sig_site = logical([chi_square_p{j}{1}; chi_square_p{j}{2}]<0.05);
                    end
                    hold on
                    plot(PD,choice_change,'o','color',c{3-i});
                    plot(PD(sig_site),choice_change(sig_site),'o','color',c{3-i},'markerfacecolor',c{3-i});
                    plot([-1 1],[0 0],'k:');
                    if i==1
                        ylabel('Rot. Choice changed');
                    else
                        ylabel('Tran. Choice changed');
                    end
                    
                    if i == 1
                        xlabel('Rotation preference');
                    else
                        xlabel('Translation preference');
                    end
                    
                    axis tight
                    ylimit = ylim;
                    ylimit = ylimit.*1.1;
                    plot([-180 180],[0 0],'k:');
                    plot([-90 -90],ylimit,'k:');
                    plot([0 0],ylimit,'k:');
                    plot([90 90],ylimit,'k:');
                    xlimit = [-180,180];
                    xlim(xlimit);
                    xticks([-180:45:180]);
                    if i == 1
                        xticklabels_temp = xticklabels;
                        xticklabels_temp{[-180:45:180]==-90} = 'CW';
                        xticklabels_temp{([-180:45:180]==90)} = 'CCW';
                        xticklabels(xticklabels_temp);
                    end
                end
            end
            if j == 1 % difficult trial, 0 degree, around 0 degree
                suptitle('Difficult trial');
            elseif j == 2
                suptitle('0 degree');
            else
                suptitle('Around 0 degree');
            end
            SetFigure(13);
            figN = figN + 1;
        end
    end



    function f4p2(debug)  % around zero choice and d'
        if debug  ; dbstack;   keyboard;      end
        
        
        %         % ------------------ default parameter -----------------------
        %         methods_of_select = {select_typical, 'Typical cells'};
        %
        %         % 每次请只选Fine task 或者 coarse task，懒得写合并一起的，默认4 target，不需要选
        %         tasktype_select = logical([get(findall(gcbf,'tag','fine_task'),'value')  get(findall(gcbf,'tag','coarse_task'),'value')]);
        %         if sum(tasktype_select) == 0
        %             disp('********** Atleast choose one task type! **********');
        %             return
        %         end
        %
        %         if sum(tasktype_select) == 2
        %             disp('********** Please choose one task type! **********');
        %             return
        %         end
        %
        %         stimtask_type_this = cell2mat(stimtask_type); % 4 targets fine task
        %
        %         % ---------------------- get data ----------------------------
        %         az = [0:90:270]; % R-CW-L-CCW
        %         d_prime_single_fine = cell2mat({group_result.d_prime_single_fine}'); % R-L-CCW-CW
        %         d_prime_single_coarse = cell2mat({group_result.d_prime_single_coarse}'); % R-L-CCW-CW
        %
        %         d_prime_vector = (d_prime_single_coarse(:,[1 4 2 3])); % R-CW-L-CCW
        %
        %         for n = 1:length(d_prime_vector)
        %             % vector-sum
        %             if sum(isnan(d_prime_vector(n,:)))<1
        %                 [d_azi(n,1), ~, d_amp(n,1)] = vectorsumSpiral(d_prime_vector(n,:),az); % in deg
        %             else
        %                 d_azi(n,1) = nan;
        %                 d_amp(n,1) = nan;
        %             end
        %         end
        %
        %         % group_choice_array_0{cs}(n,i)
        %         % cs:Control = 1; Stim = 2;
        %         % i: % R L CCW CW
        %         for n = 1:length(loadM_data)
        %             if methods_of_select{1,1}(n)==1
        %                 for cs = 1:2 % control or stim
        %                     for i = 1:4 % R L CCW CW
        %                         % select near zero degree choice
        %                         heading_part_num = size(group_result(n).choice_array2{stimtask_type_this}{1,cs},2);
        %                         rotation_part_num = size(group_result(n).choice_array2{stimtask_type_this}{2,cs},2);
        %
        %                         % have zero degree?
        %                         if rem(heading_part_num,2) == 0 % 偶数，no zero degree
        %                             selectH = [heading_part_num/2,heading_part_num/2+1]; % [1 2 3 4 5 6]中的[3 4]
        %                         else
        %                             selectH = [floor(heading_part_num/2):floor(heading_part_num/2)+2]; % [1 2 3 4 5 6 7]中的[3 4 5]
        %                         end
        %                         if rem(rotation_part_num,2) == 0 % 偶数，no zero degree
        %                             selectR = [rotation_part_num/2,rotation_part_num/2+1]; % [1 2 3 4 5 6]中的[3 4]
        %                         else
        %                             selectR = [floor(rotation_part_num/2),floor(rotation_part_num/2)+2]; % [1 2 3 4 5 6 7]中的[3 5]，中间0度【4】与heading part重复了
        %                         end
        %
        %                         group_choice_array_0{cs}(n,i) = sum(group_result(n).choice_array2{stimtask_type_this}{1,cs}(i,selectH)) + ...
        %                             sum(group_result(n).choice_array2{stimtask_type_this}{2,cs}(i,selectR)); % i:R-CW-L-CCW
        %                     end
        %                 end
        %             else
        %                 group_choice_array_0{1}(n,:) =  nan(1,4);
        %                 group_choice_array_0{2}(n,:) =  nan(1,4);
        %             end
        %         end
        %
        %         % 电刺激前后差值的vector-sum
        %         for n = 1:length(group_choice_array_0{1})
        %             choice_change_ori(n,:) = group_choice_array_0{2}(n,:) - group_choice_array_0{1}(n,:); % after-befor, R-CW-L-CCW
        %             %             choice_change_nor(n,:) = mapminmax(choice_change_ori(n,:),0,1); % 差值归一化到0-1之间
        %             choice_change_nor(n,:) = choice_change_ori(n,:); % 差值归一化到0-1之间
        %
        %             %             [~,~,ic] = unique(abs(choice_change_ori(n,:)));
        %             %             choice_change_nor2(n,:) = mapminmax(ic',1,4); % 1
        %             %             是变化最小（或者是负值变化最大），4是变化最大
        %
        %             % vector-sum
        %             if sum(isnan(choice_change_nor(n,:)))<1
        %                 [change_azi(n,1), ~, change_amp(n,1)] = vectorsumSpiral(choice_change_nor(n,:),az);
        %             else
        %                 change_azi(n,1) = nan;
        %                 change_amp(n,1) = nan;
        %             end
        %         end
        %
        %         % ------------------------ Plot ------------------------------
        %         % 分别对4个方向进行plot
        %         set(figure(figN),'name','every d and choice change','pos',[100,100,700,700]); clf
        %         tl = [];
        %         tl{1} = 'R';
        %         tl{2} = 'CW';
        %         tl{3} = 'L';
        %         tl{4} = 'CCW';
        %         for i = 1:4 % R CW L CCW
        %             subplot(2,2,i)
        %             plot(d_prime_single_coarse(:,i),choice_change_nor(:,i),'ko');
        %             hold on
        %             [r,p] = plot_corr_line(d_prime_single_coarse(:,i),choice_change_nor(:,i),'MethodOfCorr','Pearson');
        %             xlimit = xlim;
        %             ylimit = ylim;
        %             text(xlimit(1)+(xlimit(2)-xlimit(1))/20,ylimit(2)-(ylimit(2)-ylimit(1))/20,sprintf('r = %0.3g, p = %0.3g',r_square,p));
        %             title(tl{i})
        %             xlabel('d prime with spon.');
        %             ylabel('Choice changed');
        %         end
        %         suptitle('Choice changed around zero degree with d prime');
        %
        %
        %         % 比较d' vector-sum结果和choice change vector-sum结果
        %         set(figure(figN+1),'name','d and choice change vertor-sum','pos',[50,100,400,400]); clf
        %         plot(d_azi,change_azi,'ko'); % 散点图不好表示，因为横纵坐标都是圆
        %         xlim([0 360]);
        %         ylim([0 360]);
        %         xlabel('d vertor-sum');
        %         ylabel('Choice change vertor-sum');
        %
        %         % Pearson correlation
        %         hold on
        %         [r,p] = plot_corr_line(d_azi,change_azi,'MethodOfCorr','Pearson');
        %
        %         set(figure(figN+2),'name','d and choice change vertor-sum','pos',[500,100,500,300]); clf
        %         % 两者做差
        %         diff_dc = abs(change_azi - d_azi);
        %         diff_dc(diff_dc>180) = 360 - diff_dc(diff_dc>180);
        %         edges = [0:30:180];
        %         histogram(diff_dc,edges);
        %         ylabel('Cases');
        %         xlabel('|vectorSum(d) - vectorSum(choice change)|');
        %
        %         meadian_d = nanmedian(diff_dc);
        %         std_d = nanstd(diff_dc);
        %
        %         [p,h] = signtest(diff_dc);
        %
        %         figN = figN + 3;
    end


    function f4p3(debug) % pref R L CCW CW 四种细胞分别在R L CCW CW刺激下，电刺激的影响
        if debug  ; dbstack;   keyboard;      end
        %         stimtask_type_this = 4; % 4 targets coarse task
        if sum(stimtask_type{1} == 1)+sum(stimtask_type{2} == 4) == 1
            if sum(stimtask_type{1} == 1)>0
                stimtask_type_this = 1; % 4 targets fine task
            elseif sum(stimtask_type{2} == 4)>0
                stimtask_type_this = 4; % 4 targets coarse task
            end
        else
            disp('*********** PLEASE choose fine OR coarse task************');
            return
        end
        
        methods_of_select = {select_typical, 'Typical cells'};
        
        for n = 1:length(loadM_data)
            if methods_of_select{1,1}(n)
                for cs = 1:2 % control or stim
                    %                 group_result(n).choice_array{stimtask_type_this}{1,cs} % Heading visual stim
                    %                 group_result(n).choice_array{stimtask_type_this}{2,cs} % Rotation visual stim
                    
                    % visual stim: R L CCW CW
                    heading_part_num = size(group_result(n).choice_array{stimtask_type_this}{1,cs},2);
                    rotation_part_num = size(group_result(n).choice_array{stimtask_type_this}{2,cs},2);
                    
                    % have zero degree?
                    if rem(heading_part_num,2) == 0 % 偶数，no zero degree e.g.heading = [1 2 3 4 5 6]
                        visual_L = [1:heading_part_num/2];% [1 2 3]
                        visual_R = [heading_part_num/2+1:heading_part_num];% [4 5 6]
                    else % have zero,e.g. heading = [1 2 3 4 5 6 7]
                        visual_L = [1:ceil(heading_part_num/2)];% [1 2 3 4]
                        visual_R = [ceil(heading_part_num/2):heading_part_num];% [4 5 6 7]
                    end
                    if rem(rotation_part_num,2) == 0 % 偶数，no zero degree e.g.rotation = [1 2 3 4 5 6] （横坐标是rotation angle,从CW-0-CCW）
                        visual_CW = [1:rotation_part_num/2];% [1 2 3]
                        visual_CCW = [rotation_part_num/2+1:rotation_part_num];% [4 5 6]
                    else % have zero,e.g. rotation = [1 2 3 4 5 6 7]
                        visual_CW = [1:ceil(rotation_part_num/2)];% [1 2 3 4]
                        visual_CCW = [ceil(rotation_part_num/2):rotation_part_num];% [4 5 6 7]
                    end
                    
                    % choose around zeor degree
                    visual_L_choose{n}(cs,:) = sum(group_result(n).choice_array{stimtask_type_this}{1,cs}(:,visual_L),2); % 第一列chooseR，第二列chooseL，CCW，CW； 第一行是刺激前，第二行是刺激后 ()
                    visual_R_choose{n}(cs,:) = sum(group_result(n).choice_array{stimtask_type_this}{1,cs}(:,visual_R),2);
                    visual_CW_choose{n}(cs,:) = sum(group_result(n).choice_array{stimtask_type_this}{2,cs}(:,visual_CW),2);
                    visual_CCW_choose{n}(cs,:) = sum(group_result(n).choice_array{stimtask_type_this}{2,cs}(:,visual_CCW),2);
                end
                
                % 行：R-L-CCW-CW（实际给的刺激）;  列:R-L-CCW-CW(选择的target)
                visual_choose_change{n}(1,:) = visual_R_choose{n}(2,:) - visual_R_choose{n}(1,:); % 第一行：visual R
                visual_choose_change{n}(2,:) = visual_L_choose{n}(2,:) - visual_L_choose{n}(1,:); % visual L
                visual_choose_change{n}(3,:) = visual_CCW_choose{n}(2,:) - visual_CCW_choose{n}(1,:); % visual CCW
                visual_choose_change{n}(4,:) = visual_CW_choose{n}(2,:) - visual_CW_choose{n}(1,:); % visual CW
            else
                visual_choose_change{n} = nan(4,4);
            end
        end
        
        % 对细胞prefer direction进行分类
        pref_direc = [];
        pref_direc = cell2mat({group_result.pref_direc}');
        dir_range = [-180,180];
        % 分为4个象限的细胞
        pref_dir_group{1} = logical(pref_direc(:,1)<0 & pref_direc(:,2)>0); % 左上，CCW，L
        pref_dir_group{2} = logical(pref_direc(:,1)>0 & pref_direc(:,2)>0); % 右上，CCW,R
        pref_dir_group{3} = logical(pref_direc(:,1)<0 & pref_direc(:,2)<0); % 左下，CW，L
        pref_dir_group{4} = logical(pref_direc(:,1)>0 & pref_direc(:,2)<0); % 右下，CW，R
        
        % pref direction在不同象限的细胞在不同视觉刺激下不同选择的变化
        for i = 1:4 % pref direction group
            cnum = size(visual_choose_change(pref_dir_group{i}),2);
            temp = visual_choose_change(pref_dir_group{i});
            for n = 1:cnum
                temp2 = temp{n};
                temp3(n,:,:) = temp2;
            end
            pref_dir_visual_choose_change{i} = squeeze(nanmean(temp3));
        end
        
        set(figure(figN),'name','Choice change for foure type cell in four type sim.','pos',[50,50,600,600]); clf
        for i = 1:4
            subplot(2,2,i)
            Z = pref_dir_visual_choose_change{i};
            % 伪造最后一行和最后一列
            [last_x,last_y] = size(Z);
            Z(last_x+1,:) = zeros(1,last_y);
            Z(:,last_y+1) = zeros(last_x+1,1);
            h = pcolor(Z);
            axis equal tight
            climit = max(abs(caxis));
            caxis([-climit,climit]); % 统一colorbar
            colorbar
            
            set(gca,'xtick',1.5:1:4.5);
            set(gca,'ytick',1.5:1:4.5);
            set(gca,'xticklabel',{'R' 'L' 'CCW' 'CW'});
            set(gca,'yticklabel',{'CW' 'CCW' 'L' 'R'});
            xlabel('Choose Changed');
            ylabel('Self-motion stim.');
            
            if i == 1
                title('Prefer L+CCW');
            elseif i == 2
                title('Prefer R+CCW');
            elseif i == 3
                title('Prefer L+CW');
            elseif i == 4
                title('Prefer R+CW');
            end
        end
        
        standard_color = [0 0 1; 1 1 1; 1 0 0];
        color_map_RGB = colormap_set(standard_color);
        colormap(color_map_RGB);
        
        [til] = SetTitle(stimtask_type,monkey_included_for_analysis);
        suptitle(til);
        
        SetFigure(10);
        figN = figN + 1;

    end

    function f4p4(debug) % 'a new Motion type psy curve (only preferred half-axis for T/R)'
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        stimtask_type_this = 1; % 4 targets fine task
        
        % ---------------------- get data ----------------------------
        % 做了4-target task
        select_4T_task = false(length(loadM_data),1);
        for n = 1:length(loadM_data)
            if ischar(group_result(n).stim_FILE{1}) % have done the fine 4 target task
                select_4T_task(n) = true;
                
                % normalize the x axis (from rotation to translation)
                unique_heading = group_result(n).Stim_unique_heading{stimtask_type_this};
                unique_rotation = group_result(n).Stim_unique_rotation{stimtask_type_this};
                nor_heading_temp = unique_heading(ceil(length(unique_heading)/2)+1:end);
                nor_rotation_temp = unique_rotation(ceil(length(unique_rotation)/2)+1:end);
                nor_heading = normr(nor_heading_temp);
                nor_rotation = normr(nor_rotation_temp);
                if sum(unique_heading==0) > 0
                    unique_mt{n} = [flipud(-nor_rotation),0,nor_heading];
                else
                    unique_mt{n} = [flipud(-nor_rotation),nor_heading];
                end
                
                % 根据pref拼接motion type curve，只选取各自prefer方向的那一段
                %                 global_pref_code(n,:)
            end
        end
    end

    function f4p18(debug) % 'show example cell'
        if debug  ; dbstack;   keyboard;      end
        
        file_select = [];
        for fc = 1:2
            if fc == 1
                stimtask = 1;
            else
                stimtask = 4;
            end
            
            % find MT significant and positive stim site and chi-square significant
            selectMT= logical(diffBias_p_task{stimtask}(:,3)<0.05 & diffBias_task{stimtask}(:,3)>0 & chi_square_p_around0{stimtask}(:,3)<0.05);
            findMT = find(diffBias_p_task{stimtask}(:,3)<0.05 & diffBias_task{stimtask}(:,3)>0  & chi_square_p_around0{stimtask}(:,3)<0.05);
            
            % sort by PSE shift: large to small
            [~,idx] = sort(diffBias_task{stimtask}(selectMT,3),'descend');
            
            % show file name
            for i = 1:length(findMT)
                file_select{i,fc} = group_result(findMT(i)).stim_FILE{stimtask};
            end
            
            save_this{fc} = {file_select{idx,fc}}';
        end
        
        save_this{1}
        save_this{2}
        keyboard
    end




    function f4p99(debug)
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        stimtask_type_this = 1; % 4 targets fine task
        
        % ---------------------- get data ----------------------------
        % group all L/R choice across Heading and Rotation part; all CW/CCW
        % choice across Heading and Rotation part, seperate from control and
        % stim
        
        % group_choice_array(n,cc,cs)
        % cc;Heading choice = 1; Rotation choice = 2
        % cs:Control = 1; Stim = 2;
        for n = 1:length(loadM_data)
            for cs = 1:2 % control or stim
                for i = 1:4 % R L CW CCW
                    group_choice_array{cs}(n,i) = sum(group_result(n).choice_array{stimtask_type_this}{1,cs}(i,:)) + ... % Heading part choice
                        sum(group_result(n).choice_array{stimtask_type_this}{2,cs}(i,:)); % Rotation part choice  % R-L-CCW-CW
                    
                    % select near zero degree choice
                    heading_part_num = size(group_result(n).choice_array{stimtask_type_this}{1,cs},2);
                    rotation_part_num = size(group_result(n).choice_array{stimtask_type_this}{2,cs},2);
                    
                    % have zero degree?
                    if rem(heading_part_num,2) == 0 % 偶数，no zero degree
                        selectH = [heading_part_num/2,heading_part_num/2+1]; % [1 2 3 4 5 6]中的[3 4]
                    else
                        selectH = [floor(heading_part_num/2):floor(heading_part_num/2)+2]; % [1 2 3 4 5 6 7]中的[3 4 5]
                    end
                    if rem(rotation_part_num,2) == 0 % 偶数，no zero degree
                        selectR = [rotation_part_num/2,rotation_part_num/2+1]; % [1 2 3 4 5 6]中的[3 4]
                    else
                        selectR = [floor(rotation_part_num/2),floor(rotation_part_num/2)+2]; % [1 2 3 4 5 6 7]中的[3 5]，中间0度【4】与heading part重复了
                    end
                    
                    group_choice_array_0{cs}(n,i) = sum(group_result(n).choice_array{stimtask_type_this}{1,cs}(i,selectH)) + ...
                        sum(group_result(n).choice_array{stimtask_type_this}{2,cs}(i,selectR)); % i:R-L-CCW-CW
                end
            end
        end
        
        % 不区分 heading target与rotation target
        for cs = 1:2
            group_choice_array_HR{cs}(:,1) = group_choice_array{cs}(:,1) + group_choice_array{cs}(:,2);
            group_choice_array_HR{cs}(:,2)= group_choice_array{cs}(:,3) + group_choice_array{cs}(:,4);
            
            group_choice_array_HR_0{cs}(:,1) = group_choice_array_0{cs}(:,1) + group_choice_array_0{cs}(:,2);
            group_choice_array_HR_0{cs}(:,2)= group_choice_array_0{cs}(:,3) + group_choice_array_0{cs}(:,4);
        end
        
        % remove nan
        group_choice_array{1}(any(isnan(group_choice_array{1}),2),:)=[];
        group_choice_array{2}(any(isnan(group_choice_array{2}),2),:)=[];
        group_choice_array_HR{1}(any(isnan(group_choice_array_HR{1}),2),:)=[];
        group_choice_array_HR{2}(any(isnan(group_choice_array_HR{2}),2),:)=[];
        group_choice_array_0{1}(any(isnan(group_choice_array_0{1}),2),:)=[];
        group_choice_array_0{2}(any(isnan(group_choice_array_0{2}),2),:)=[];
        group_choice_array_HR_0{1}(any(isnan(group_choice_array_HR_0{1}),2),:)=[];
        group_choice_array_HR_0{2}(any(isnan(group_choice_array_HR_0{2}),2),:)=[];
        
        
        
        % control 没有电刺激时候，猴子对4个选项有没有偏好性 （偏好Heading！）
        az = deg2rad([0:90:270]); % R-CW-L-CCW
        figure
        %         for i = 1:length(group_choice_array_0{1})
        %             polarplot([az az(1)],[group_choice_array_0{1}(i,1),group_choice_array_0{1}(i,3),group_choice_array_0{1}(i,2),group_choice_array_0{1}(i,4),group_choice_array_0{1}(i,1)]);
        %             hold on
        %         end
        sum_around_0 = sum(group_choice_array_0{1});
        polarplot([az az(1)],[sum_around_0(1) sum_around_0(3) sum_around_0(2) sum_around_0(4) sum_around_0(1)],'linewidth',2);
        
        % 做了4-target task的unit的tuning
        select_4T_task = [];
        nn1 = 1;
        for n = 1:length(loadM_data)
            if ischar(group_result(n).stim_FILE{1}) % have done the fine 4 target task
                select_4T_task(nn1,1) = n;
                nn1 = nn1+1;
            end
        end
        
        FourT_cell_type = cell2mat({group_result(select_4T_task).cell_type}); % heading =1;rotation = 2; spiral = 3
        FourT_spiral_index = cell2mat({group_result(select_4T_task).spiral_index})';
        
        figure
        hist(FourT_cell_type,[1:5]);
        
        % 同时做了heading only与rotation only task的unit
        select_2T_task = [];
        nn1 = 1;
        for n = 1:length(loadM_data)
            if ischar(group_result(n).stim_FILE{2}) && ischar(group_result(n).stim_FILE{3}) % have done the both fine 2 target task
                select_2T_task(nn1,1) = n;
                nn1 = nn1+1;
            end
        end
        TwoT_cell_type = cell2mat({group_result(select_2T_task).cell_type})'; % heading =1;rotation = 2; spiral = 3
        TwoT_spiral_index = cell2mat({group_result(select_2T_task).spiral_index})';
        
        %% 2 target task, compare Heading and rotation PSE shift
        % across cell type
        TwoT_H_PSE_shift_temp = nor_diffBias_task{2}(:,1);
        TwoT_R_PSE_shift_temp = nor_diffBias_task{3}(:,2);
        TwoT_H_PSE_shift = TwoT_H_PSE_shift_temp(select_2T_task);
        TwoT_R_PSE_shift = TwoT_R_PSE_shift_temp(select_2T_task);
        
        figure(205);clf
        set(205,'Name', 'Compare Heading and Rotation PSE shift in 2-target task', 'color','w');
        % across cell type
        subplot(1,2,1)
        col{1} = 'r';
        col{2} = 'b';
        col{3} = 'g';
        for ct = 1:3 % heading cell, r cell, spiral cell
            select = logical(TwoT_cell_type == ct);
            plot(TwoT_H_PSE_shift(select),TwoT_R_PSE_shift(select),'o','color',col{ct});
            hold on
        end
        legend('Heading cell','Rotation cell','Spiral cell');
        xlim([-4 10]);
        ylim([-4 10]);
        plot([-4 10],[-4 10],'k:');
        xlabel('normlize Heading PSE shift');
        ylabel('normlize Rotation PSE shift');
        
        % across spiral_index
        max_color_value = 20;
        jet_color = colormap(jet(max_color_value));
        color_index = ceil(TwoT_spiral_index*max_color_value/2 + max_color_value/2); % 将spiral index（-1~1）对应到max_color_value（1-20）个color上去
        % 越Rotation，spiral index（-1）越小，color_index越小（1），颜色越蓝；
        % 越Heading，spiral index（1）越大，color index越大（max_color_value20），颜色越红
        
        subplot(1,2,2)
        for i = 1:length(select_2T_task)
            try
                selected_color = jet_color(color_index(i),:);
                plot(TwoT_H_PSE_shift(i),TwoT_R_PSE_shift(i),'o','color',selected_color,'markerfacecolor',selected_color);
            catch
                plot(TwoT_H_PSE_shift(i),TwoT_R_PSE_shift(i),'ko'); % m7c1124 没有siral_index (该文件丢失)
            end
            hold on
        end
        colorbar('Ticks',[0,1],'TickLabels',{'Rotation','Heading'});
        xlim([-4 10]);
        ylim([-4 10]);
        plot([-4 10],[-4 10],'k:');
        xlabel('normlize Heading PSE shift');
        ylabel('normlize Rotation PSE shift');
        
        % in 4-target fine task
        FourT_H_PSE_shift_temp = nor_diffBias_task{1}(:,1);
        FourT_R_PSE_shift_temp = nor_diffBias_task{1}(:,2);
        FourT_H_PSE_shift = FourT_H_PSE_shift_temp(select_4T_task);
        FourT_R_PSE_shift = FourT_R_PSE_shift_temp(select_4T_task);
        
        
        figure(206);clf
        set(206,'Name', 'Compare Heading and Rotation PSE shift in 4-target task', 'color','w');
        max_color_value = 20;
        jet_color = colormap(jet(max_color_value));
        color_index4 = ceil(FourT_spiral_index*max_color_value/2 + max_color_value/2); % 将spiral index（-1~1）对应到max_color_value（1-20）个color上去
        % 越Rotation，spiral index（-1）越小，color_index越小（1），颜色越蓝；
        % 越Heading，spiral index（1）越大，color index越大（max_color_value20），颜色越红
        for i = 1:length(select_4T_task)
            try
                selected_color = jet_color(color_index4(i),:);
                plot(FourT_H_PSE_shift(i),FourT_R_PSE_shift(i),'o','color',selected_color,'markerfacecolor',selected_color);
            catch
                plot(FourT_H_PSE_shift(i),FourT_R_PSE_shift(i),'ko');
            end
            hold on
        end
        colorbar('Ticks',[0,1],'TickLabels',{'Rotation','Heading'});
        xlim([-4 10]);
        ylim([-4 10]);
        plot([-4 10],[-4 10],'k:');
        xlabel('normlize Heading PSE shift');
        ylabel('normlize Rotation PSE shift');
        
        
        
        %% 电刺激后 0 degree附近猴子四种choice选择的变化
        % 如何定义“改变”了？
        for i = 1:length(group_choice_array_0{1})
            %             change_thr = sum(group_choice_array_0{1}(i,:))/10; % change threshold   % 更合理？
            change_thr = std([group_choice_array_0{1}(i,:), group_choice_array_0{2}(i,:)])/2;
            
            % 增加定义为1，减少定义为-1，不变定义为0
            change_choice_num = group_choice_array_0{2}(i,:) - group_choice_array_0{1}(i,:);
            change_choice_sign = sign(change_choice_num);
            change_choice_sign(abs(change_choice_num)<change_thr) = 0; % not change above the threshold
            change_choice_code(i,:) = change_choice_sign; % save it
        end
        
        % 认为定义几种电刺激可能引起的变化情况
        for n = 1:length(change_choice_code)
            A = sum(change_choice_code(n,[1:2])); % L + R
            B = sum(change_choice_code(n,[3:4])); % CW + CCW
            unique_code = length(unique(change_choice_code(n,:)));
            if ~any([A B] - [1 -1]) || ~any([A B] - [-1 1])
                situList(n,1) = 1; % 一个motion type内1种direction增加，另一个motion type的一个direction减少
            elseif (A==0 && B==0 && unique_code == 3)
                situList(n,1)=2; % 一个motion type内，一个direction增加，另一个减少
            elseif ~any([A B] - [0 -1]) || ~any([A B] - [-1 0])
                situList(n,1)=3; % 一个motion type内，一个direction增加，另一个减少，另一个motiontype内其中一个direction也减少
            elseif ~any([A B] - [1 -2]) || ~any([A B] - [-2 1])
                situList(n,1)=4; % 一个motion type内，一个direction增加，另一个motion type两个direction都减少
            elseif ~any([A B] - [0 -2]) || ~any([A B] - [-2 0])
                situList(n,1)=5; % 一个direction增加，其余3个direction减少
            elseif ~any([A B] - [0 1]) || ~any([A B] - [1 0])
                situList(n,1)=6; % 两个不同motion type的direction增加，两外两个一个不变一个减少
            elseif (A==0 && B==0 && unique_code == 2)
                situList(n,1)=7; % 两个不同motion type的direction增加，另外两个都减少
            elseif ~any([A B] - [-2 2]) || ~any([A B] - [2 -2])
                situList(n,1)=8; % 两个同样motion type的direction增加，另外零个减少
            elseif (A==0 && B==0 && unique_code == 1)
                situList(n,1)=9; % 都没有变化
            else
                situList(n,1)=nan;
            end
        end
        xbins = [1:9];
        figure(201);
        set(201,'Name', 'Stim type in 2 axes in all 4-target task', 'color','w');
        hist(situList,xbins);
        
        figure(202)
        set(202,'Name', 'Stim type in 2 axes in all 4-target task (seperate cell type)', 'color','w');
        % 区分cell type看电刺激引起的类型： heading cell/rotation cell/spiral cell
        for ct = 1:3 % heading cell/rotation cell/spiral cell
            nn1 = 1;
            for n = 1:length(change_choice_code)
                if FourT_cell_type(n) == ct
                    situList_cell_type{ct}(nn1,1) = situList(n,1);
                    nn1 = nn1 + 1;
                end
            end
            subplot(1,3,ct)
            hist(situList_cell_type{ct},xbins);
        end
        
        %% 电刺激引起不同cell type在0度附近选Heading target 与 rotatino target数量的变化
        figure(203)
        set(203,'Name', 'Choice change around 0 across cell type (normalized)', 'color','w');
        for ct = 1:3 % heading cell/rotation cell/spiral cell
            control = group_choice_array_HR_0{1}(FourT_cell_type == ct,:);
            stim = group_choice_array_HR_0{2}(FourT_cell_type == ct,:);
            
            choose_h = [control(:,1) stim(:,1)];
            choose_r = [control(:,2) stim(:,2)];
            
            choose_ctrl = control(:,2) ./ control(:,1) ; % control choose R/H
            choose_stim =  stim(:,2) ./  stim(:,1); % stim: choose R/H
            
            % normalize control and stim
            choose_h_nor = normr(choose_h);
            choose_r_nor = normr(choose_r);
            
            temp = [];
            temp = [choose_ctrl choose_stim];
            temp_nor = normr(temp);
            
            edges = 10:10:70;
            
            choose_h_mean = mean(choose_h_nor);
            choose_r_mean = mean(choose_r_nor);
            choose_h_sem = std(choose_h_nor) / sqrt(length(choose_h_nor));
            choose_r_sem = std(choose_r_nor) / sqrt(length(choose_r_nor));
            
            subplot(3,3,ct)
            plot([1,2],choose_h_nor,'r.-');
            hold on
            plot([1,2],choose_r_nor,'b.-');
            if ct ==1
                ylabel('Choose H or Choose R');
            end
            
            try
                subplot(3,3,ct+3)
                shadedErrorBar([1 2],choose_h_mean,choose_h_sem,'r',0.5);
                hold on
                shadedErrorBar([1 2],choose_r_mean,choose_r_sem,'b',0.5);
            end
            
            subplot(3,3,ct+6)
            plot([1,2],temp_nor,'.-');
            if ct == 1
                ylabel('Choose H/R');
                xlabel('Translation Cell');
            elseif ct == 2
                xlabel('Rotation Cell');
            elseif ct == 3
                xlabel('Spiral Cell');
            end
            
            % pair ttest 检验 cell type内的显著性
            [~,p(ct)] = ttest(choose_ctrl,choose_stim);
        end
        
        
        
        
        % 电刺激引起不同spiral index在0度附近/所有condition 选Heading target 与 rotatino target数量的变化
        % 横坐标：spiral index  纵坐标：Hstim/Hcon  Rstim/Rcon
        group_HR{1} = group_choice_array_HR;
        group_HR{2} = group_choice_array_HR_0;
        
        %         FourT_spiral_index = [];
        %         FourT_spiral_index = cell2mat({group_result(select_4T_task).neu_sen_index_coarse})';
        
        %                 modulation_index{1} = [group_result.ahead_slope_index]'; % 斜率绝对值
        %         modulation_index{2} = [group_result.d_prime_index_fine]';
        %         modulation_index{3} = [group_result.d_prime_index_coarse]';
        %         modulation_index{4} = [group_result.d_prime_maxDiff_index]'; % 差值最大的两个点,在Group_SpiT中去掉>6std的异常值
        %         modulation_index{5} = [group_result.d_prime_maxAxis_index]'; % 差值最大的一条轴（180相对两点）,在Group_SpiT中去掉>6std的异常值
        %         modulation_index{6} = [group_result.d_prime_maxFR_index]'; % FR最大点及其对侧点,在Group_SpiT中去掉>6std的异常值
        %         modulation_index{7} = [group_result.DDI_inc_index]';
        %         modulation_index{8} = [group_result.DDI_exc_index]';
        %         modulation_index{9} = [group_result.HTI_index]';
        %         modulation_index{10} = [group_result.neu_sen_index_fine]';
        %         modulation_index{11} = [group_result.neu_sen_index_coarse]';
        %         modulation_index{12} = [group_result.spiral_index]';
        %
        figure(204);clf
        set(204,'Name', 'Choice change across spiral index (all condition or around 0)', 'color','w');
        col = [];
        col{1} = 'ro';
        col{2} = 'bo';
        for j = 1:2 % all condition/ condition around 0
            for i = 1:2 % h or r
                subplot(2,2,i+(j-1)*2)
                choice_change = [group_HR{j}{2}(:,i) ./ group_HR{j}{1}(:,i)];
                plot(FourT_spiral_index,choice_change,col{i}); % stim/control
                xlim([-1 1]);
                xlabel('Spiral Index');
                if i == 1
                    ylabel('Hstim/Hcon');
                else
                    ylabel('Rstim/Rcon');
                end
                hold on
                plot([-1 1],[1 1],'k:');
                plot_corr_line(FourT_spiral_index,[group_HR{j}{2}(:,1) ./ group_HR{j}{1}(:,1)],'MethodOfCorr','Pearson');
                
                % 对spiral index分bin处理，bin内取平均值
                bins = 0.2;
                xlimit = [-1 1];
                % xlimit = [-0.5 0.5];
                bins_num = ceil((xlimit(2)-xlimit(1))/bins);
                for k = 1:bins_num
                    bins_left(k) = xlimit(1)+bins*(k-1);
                    bins_mid(k) = bins_left(k)+bins/2;
                    select = logical(FourT_spiral_index>bins_left(k) & FourT_spiral_index<=bins_left(k)+bins);
                    if sum(select)==0
                        mean_choice_change(k) = nan;
                        sem_choice_change(k) = nan;
                    else
                        mean_choice_change(k) = mean(choice_change(select));
                        sem_choice_change(k) = std(choice_change(select))/sqrt(sum(select));
                    end
                end
                hold on
                errorbar(bins_mid,mean_choice_change,sem_choice_change,'k','linewidth',1.5);
            end
        end
        
        % 卡方检验 电刺激效果 (对每个细胞进行)
        % 条件 at least 80% of the expected frequencies exceed 5 and all the expected frequencies exceed 1
        for n = 1:length(loadM_data)
            p_chi_HR(n,1) = group_result(n).p_chi_HR{1};
            p_chi_all(n,1) = group_result(n).p_chi_all{1};
            p_chi_0(n,1) = group_result(n).p_chi_0{1};
        end
        
        %         p_chi_HR(isnan(p_chi_HR)) = [];
        %         p_chi_all(isnan(p_chi_all)) = [];
        %         p_chi_0(isnan(p_chi_0)) = [];
        
        nansum(p_chi_HR<0.05)
        
        nansum(p_chi_all<0.05)
        nansum(p_chi_0<0.05)
    end

    function f4p5(debug) % 'Four choice changed in stimulation'
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        stimtask_type_this{1} = 1; % 4 targets fine task
        stimtask_type_this{2} = 4; % 4 targets coarse task
        
        % ---------------------- get data ----------------------------
        % choice_array: eveyy cell, R-L-CCW-CW
        for fc = 1:2 % fine and coarse
            for n = 1:length(loadM_data)
                for cs = 1:2 % control or stim
                    % select near zero degree choice
                    heading_part_num = size(group_result(n).choice_array{stimtask_type_this{fc}}{1,cs},2);
                    rotation_part_num = size(group_result(n).choice_array{stimtask_type_this{fc}}{2,cs},2);
                    
                    % have zero degree?
                    if rem(heading_part_num,2) == 0 % 偶数，no zero degree
                        selectH = [heading_part_num/2,heading_part_num/2+1]; % [1 2 3 4 5 6]中的[3 4]
                    else
                        selectH = [floor(heading_part_num/2):floor(heading_part_num/2)+2]; % [1 2 3 4 5 6 7]中的[3 4 5]
                    end
                    if rem(rotation_part_num,2) == 0 % 偶数，no zero degree
                        selectR = [rotation_part_num/2,rotation_part_num/2+1]; % [1 2 3 4 5 6]中的[3 4]
                    else
                        selectR = [floor(rotation_part_num/2),floor(rotation_part_num/2)+2]; % [1 2 3 4 5 6 7]中的[3 5]，中间0度【4】与heading part重复了
                    end
                    
                    for i = 1:4 % R L CCW CW
                        % choice around 0 degree
                        group_choice_array_around0{fc,cs}(n,i) = sum(group_result(n).choice_array{stimtask_type_this{fc}}{1,cs}(i,selectH)) + ...
                            sum(group_result(n).choice_array{stimtask_type_this{fc}}{2,cs}(i,selectR)); % i:R-L-CCW-CW
                        
                        % choice across all condition
                        group_choice_array{fc,cs}(n,i) = sum(group_result(n).choice_array{stimtask_type_this{fc}}{1,cs}(i,:)) + ... % Heading part choice
                            sum(group_result(n).choice_array{stimtask_type_this{fc}}{2,cs}(i,:)); % Rotation part choice  % R-L-CCW-CW
                    end
                end
                % pref right/CCW(code = 2): after - before < 0
                % pref left/CW(code = 1): after - before > 0
                pref_direction = global_pref_code(n,1:2); % h,r,mt preference
                
                % find the prefer and null direction, R-L-CCW-CW sort as [TransPref(1) TransNull(2) RotPref(3) RotNull(4)]
                prefSort(n,:) = [3-pref_direction(1) pref_direction(1) 5-pref_direction(2) pref_direction(2)+2];
                
                choice_change_prop{fc}(n,:) = (group_choice_array_around0{fc,2}(n,:) - group_choice_array_around0{fc,1}(n,:)) ./ sum(group_choice_array_around0{fc,1}(n,:));  % (stim-contrl)./sum, R-CW-L-CCW
                
                % rearrange as [TransPref TransNull RotPref RotNull]
                if ~isnan(prefSort(n,1))
                    choice_change_PrefNull{fc}(n,:) = choice_change_prop{fc}(n,prefSort(n,:));
                else
                    choice_change_PrefNull{fc}(n,:) = choice_change_prop{fc}(n,:);
                end
            end
            
            % ***** cluster order 1: according tuning preference: [TransPref(1) TransNull(2) RotPref(3) RotNull(4)]
            % 统计每种情况的个数
            mask{fc} = logical(~isnan(choice_change_PrefNull{fc}(:,1)) & methods_of_select{1});
            data = [];
            data = choice_change_PrefNull{fc}(mask{fc},:);
            
            % 人为分类
            thresh = 0.11;
            %         thresh = std(data(:));
            %         datain = ctrs;
            cluster{fc}{1} = find(data(:,1)>thresh & data(:,2)<-thresh & abs(data(:,3))<thresh & abs(data(:,4))<thresh); % TransPref+↑，TransNull-↓,RotPref, RotNull (Only Between Translation)
            cluster{fc}{2} = find(data(:,3)>thresh & data(:,4)<-thresh & abs(data(:,1))<thresh & abs(data(:,2))<thresh); % TransPref，TransNull,RotPref+↑, RotNull-↓ (Only Between Rotation)
            cluster{fc}{3} = find(data(:,1)>thresh & data(:,2)<-thresh & abs(data(:,3))<thresh & data(:,4)<-thresh); % TransPref+↑，TransNull-↓,RotPref, RotNull-↓ (grab from Null)
            cluster{fc}{4} = find(abs(data(:,1))<thresh & data(:,2)<-thresh & data(:,3)>thresh & data(:,4)<-thresh); % TransPref, TransNull-↓,RotPref+↑,RotNull-↓ (grab from Null)

            cluster{fc}{5} = find(data(:,1)>thresh & abs(data(:,2))<thresh & abs(data(:,3))<thresh & abs(data(:,4))<thresh); % TransPref+↑, TransNull ,RotPref,RotNull (Only TransPref)
            cluster{fc}{5} = [];
            temp = find(data(:,1)>thresh & data(:,2)<-thresh & data(:,3)<-thresh & data(:,4)<-thresh); % TransPref+↑, TransNull-,RotPref-,RotNull-↓ (Only TransPref2) % 并入上一种cluster
            cluster{fc}{5} = [cluster{fc}{5};temp];
            
            cluster{fc}{6} = find(abs(data(:,1))<thresh & abs(data(:,2))<thresh & data(:,3)>thresh & abs(data(:,4))<thresh); % TransPref, TransNull,RotPref+↑,RotNull (Only RotPref)
            cluster{fc}{6} = [];
            temp = find(data(:,1)<-thresh & data(:,2)<-thresh & data(:,3)>thresh & data(:,4)<-thresh); % TransPref-↑, TransNull-,RotPref+,RotNull-↓ (Only RotPref) % 并入上一种cluster
            cluster{fc}{6} = [cluster{fc}{6};temp];
            
            cluster{fc}{7} = find(data(:,1)>thresh & data(:,2)<-thresh & data(:,3)>thresh & data(:,4)<-thresh); % TransPref+↑，TransNull-↓,RotPref+↑，RotNull-↓ (prefer-null)
            cluster{fc}{8} = find(data(:,1)>thresh & data(:,2)>thresh & data(:,3)<-thresh & data(:,4)<-thresh); % TransPref+↑，TransNull+↑,RotPref-↓，RotNull-↓ (between Motion type)
            cluster{fc}{9} = find(data(:,1)<-thresh & data(:,2)<-thresh & data(:,3)>thresh & data(:,4)>thresh); % TransPref-↓，TransNull-↓,RotPref+↑，RotNull+↑ (between Motion type)

            cluster{fc}{10} = find(abs(data(:,1))<thresh & abs(data(:,2))<thresh & abs(data(:,3))<thresh & abs(data(:,4))<thresh); % TransPref，TransNull,RotPref,RotNull (totally no effect)
            
            unit_num{fc} = cellfun(@length,cluster{fc});
            
            cluster_indx{fc} = nan(length(data),1);
            for i = 1:length(cluster{fc})
                cluster_indx{fc}(cluster{fc}{i})=i;
            end
            
            % ***** cluster order 2: according final choice change (use cross style，roll和translation两条正交轴)
            %  定义右target为增加的target（和下为增加target）
            group{1,fc} = logical(cluster_indx{fc} == 1 | cluster_indx{fc} == 2); % 1. microstimulation only affects choices within one type of flow-pattern
            group{2,fc} = logical(cluster_indx{fc} == 3 | cluster_indx{fc} == 4 | cluster_indx{fc} == 5 | cluster_indx{fc} == 6); % 2. microstimulation increases the choices of preferred direction in one of the flow-pattern and decreases/not change the choices of the other three targets
            group{3,fc} = logical(cluster_indx{fc} == 7); % 3. microstimulation increase the preferred direction and decrease the non-preffered direction choices for each flow-pattern.
            group{4,fc} = logical(cluster_indx{fc} == 8 | cluster_indx{fc} == 9); % 4. microstimulation increases the choices of one type of flow-pattern and decreases the choice of the other type of flow-pattern.  
        end
        group_unit_num = cellfun(@sum,group); % 第一列fine，第二列coarse
        
        % 合并fine and coarse task
        choice_change_PrefNull{3} = [choice_change_PrefNull{1}(mask{1},:);choice_change_PrefNull{2}(mask{2},:)];
        mask{3} = true(length(choice_change_PrefNull{3}),1);
        
        %                 % 卡方检验 电刺激效果 (对每个细胞进行)
        %         % 条件 at least 80% of the expected frequencies exceed 5 and all the expected frequencies exceed 1
        %         for n = 1:length(loadM_data)
        %             p_chi_HR(n,1) = group_result(n).p_chi_HR{1};
        %             p_chi_all(n,1) = group_result(n).p_chi_all{1};
        %             p_chi_0(n,1) = group_result(n).p_chi_0{1};
        %         end
        %
        %         %         p_chi_HR(isnan(p_chi_HR)) = [];
        %         %         p_chi_all(isnan(p_chi_all)) = [];
        %         %         p_chi_0(isnan(p_chi_0)) = [];
        %
        %         nansum(p_chi_HR<0.05)
        %         nansum(p_chi_all<0.05)
        %         nansum(p_chi_0<0.05)
        %
        
        
        % ------------------------ Plot ------------------------------
        tl{1} = 'Fine task';
        tl{2} = 'Coarse task';
        tl{3} = 'Fine + Coarse';
        % 先看总体情况: choice change across translation preference or rotation
        % reference
        set(figure(figN),'name','Choice change across T/R preference','pos',[50,50,1000,400]); clf
        for fc = 1:3
            subplot(1,3,fc)
            plot([1:4]',choice_change_PrefNull{fc}(mask{fc},:)','.-','color',c{99});
            
            hold on
            choice_change_mean = nanmean(choice_change_PrefNull{fc}(mask{fc},:));
            choice_change_sem = nanstd(choice_change_PrefNull{fc}(mask{fc},:))/sum(mask{fc});
            errorbar([1:4],choice_change_mean,choice_change_sem,'r-o','linewidth',3);
            
            set(gca,'xtick',[1:4],'xticklabel',{'TransPref','TransNull','RotPref','RotNull'});
            xlim([0.9 4.1]);
            ylimit = ylim;
            text(2,ylimit(2) - ylimit(2)/10,sprintf('N=%g',sum(mask{fc})));
            title(tl{fc});
            
            % ttest
            % pref and null
            [~,p_prefnull(1)] = ttest(choice_change_PrefNull{fc}(mask{fc},1),choice_change_PrefNull{fc}(mask{fc},2));
            [~,p_prefnull(2)] = ttest(choice_change_PrefNull{fc}(mask{fc},3),choice_change_PrefNull{fc}(mask{fc},4));
            
            for i = 1:2 % T or R
                if p_prefnull(i)<0.05
                    hold on
                    line([(i-1)*2+1,(i-1)*2+1],[choice_change_mean((i-1)*2+1)+0.2 ylimit(2)-ylimit(2)/3],'color','k','linestyle','-','linewidth',1.5);
                    line([(i-1)*2+2,(i-1)*2+2],[choice_change_mean((i-1)*2+2)+0.2 ylimit(2)-ylimit(2)/3],'color','k','linestyle','-','linewidth',1.5);
                    line([(i-1)*2+1,(i-1)*2+2],[ylimit(2)-ylimit(2)/3 ylimit(2)-ylimit(2)/3],'color','k','linestyle','-','linewidth',1.5);
                    if p_prefnull(i)<0.001
                        text((i-1)*2+1.2,ylimit(2)-ylimit(2)/4,'***','FontSize',25,'color','k');
                    elseif p_prefnull(i)<0.005
                        text((i-1)*2+1.3,ylimit(2)-ylimit(2)/4,'**','FontSize',25,'color','k');
                    else
                        text((i-1)*2+1.4,ylimit(2)-ylimit(2)/4,'*','FontSize',25,'color','k');
                    end
                end
            end
            
            % comparison to 0 (have stimulation effect?)
            for i = 1:4
                [~,p_Zero(i)] = ttest(choice_change_PrefNull{fc}(mask{fc},i));
                
                % p value marker
                if p_Zero(i)<0.001
                    text(i-0.33,choice_change_mean(i)+0.1,'***','FontSize',30,'color','r');
                elseif p_Zero(i)<0.01
                    text(i-0.2,choice_change_mean(i)+0.1,'**','FontSize',30,'color','r');
                elseif p_Zero(i)<0.05
                    text(i-0.1,choice_change_mean(i)+0.1,'*','FontSize',30,'color','r');
                end
            end
        end
        figN = figN + 1;
        
        % manual cluster 1: according tuning preference: [TransPref(1) TransNull(2) RotPref(3) RotNull(4)]
        for fc = 1:2
            set(figure(figN+fc-1),'name','Cluster for Four choice changed','pos',[900+fc*50,50+fc*50,900,500]); clf
            mask{fc} = logical(~isnan(choice_change_PrefNull{fc}(:,1)) & methods_of_select{1});
            data = [];
            data = choice_change_PrefNull{fc}(mask{fc},:);
            
            for i = 1:length(cluster{fc})
                subplot(3,ceil(length(cluster{fc})/3),i)
                hold on
                if unit_num{fc}(i)>0
                    plot([1:4]',data(cluster{fc}{i},:)','k-');
                    if unit_num{fc}(i)>1
                        plot([1:4]',mean(data(cluster{fc}{i},:))','ro','markerfacecolor','r');
                    end
                end
                plot([1 4],[0 0],'k--');
                plot([1 4],[-thresh -thresh],'k:');
                plot([1 4],[thresh thresh],'k:');
                ylim([-0.7 0.7]);
                title(num2str(i));
                text(1.2,0.6,sprintf('N=%g',unit_num{fc}(i)));
            end
            
            % other
            subplot(3,ceil(length(cluster{fc})/3),3*ceil(length(cluster{fc})/3))
            plot([1:4]',data(isnan(cluster_indx{fc}),:)');hold on
            plot([1 4],[0 0],'k:');
            ylim([-0.7 0.7]);
            title('Other');
            text(1.2,0.6,sprintf('N=%g',length(data)-sum(unit_num{fc})));
            suptitle(tl{fc});
            SetFigure(10);
        end
        figN = figN + 2;
        
        % ***** cluster order 2: according final choice change (use cross style，roll和translation两条正交轴)
        for fc = 1:2
            set(figure(figN+fc-1),'name','Cluster2 for Four choice changed','pos',[900+fc*50,50+fc*50,500,500]); clf
            mask{fc} = logical(~isnan(choice_change_PrefNull{fc}(:,1)) & methods_of_select{1});
            data = [];
            data = choice_change_PrefNull{fc}(mask{fc},:); % 以0为基准，0等于选择没有变化
            data = data + 1; % 以1为基准
            
            for i = 1:4
                subplot(2,2,i)
                hold on
                if i == 1 || i == 2 || i == 3% 第1.2.3种group，只有一个最大值为prefer direction
                    if group_unit_num(i,fc)>0
                        data_temp = data(group{i,fc},:);
                        % 调整顺序为 pref-non-null-null
                        [~, Imax] = max(data_temp,[],2);
                        adjust_data = nan(size(data_temp));
                        adjust_data(Imax == 1,:) = data_temp(Imax == 1,:);
                        adjust_data(Imax == 3,:) = data_temp(Imax == 3,[3 4 1 2]);
                        
                        hold on
                        plot([-1.5 1.5],[0 0],'k--');
                        plot([0 0],[-1.5 1.5],'k--');
                        for j = 1:length(adjust_data)
                            plot([adjust_data(j,1) 0 -adjust_data(j,2) 0 adjust_data(j,1)], [0 adjust_data(j,3) 0 -adjust_data(j,4) 0],'k.-'); % microstimulation
                        end
                        plot([1 0 -1 0 1], [0 1 0 -1 0],'g.-'); % normal
                        
                        text(-1.5,1.5,sprintf('n = %g',group_unit_num(i,fc)));
                        axis square
                    end
                elseif i == 4 % 没有 group{4,fc}
                    hold on
                    plot([-1.5 1.5],[0 0],'k--');
                    plot([0 0],[-1.5 1.5],'k--');
                    axis square
                    text(-1.5,1.5,sprintf('n = %g',group_unit_num(i,fc)));
                end
            end
            SetFigure(10);
            suptitle(tl{fc});
        end
        figN = figN + 2;
        
        
        
        
        
        
        %         data_zscore = zscore(data,1,2); % zscore会保持相对0的大小，但是与0之间的差别会放大
        %         data_zscore = data;
        
        % kmean等cluster只关注曲线的形状，而不关注大于或者小于0，感觉不太适合
        % clunum = 20;
        %         % 聚类分析
        %         % 1. hierarchical clustering
        %         % 1.1 clusterdata 一个函数计算
        %         T = clusterdata(temp_data,'linkage','ward','savememory','off','maxclust',clunum);
        %         figure
        %         for i = 1:clunum
        %             subplot(3,ceil(clunum/3),i)
        %             plot([1:4]',temp_data(T==i,:)','k.-');
        %             hold on
        %             plot([1 4],[0 0],'k:');
        %             ylim([-0.5 0.5]);
        %         end
        %
        %         % 1.2 分步计算
        %         corrDist = pdist(temp_data,'corr');
        %         clusterTree = linkage(corrDist,'average');
        %         clusters = cluster(clusterTree,'maxclust',clunum);
        %         % plot cluser together
        %         figure
        %         for i = 1:clunum
        %             subplot(3,ceil(clunum/3),i)
        %             plot([1:4]',temp_data((clusters == i),:)');
        %             hold on
        %             plot([1 4],[0 0],'k:');
        %             ylim([-0.5 0.5]);
        % %             axis tight
        %         end
        %         suptitle('Hierarchical Clustering');
        
        %         % 2. kmean (will difference to upper algorithm)
        %         [idx,ctrs] = kmeans(data_zscore,clunum,'Distance','cityblock','Replicates',10,...
        %             'disp','final','MaxIter',1000,'EmptyAction','singleton');
        %         figure
        %         for i = 1:clunum
        %             subplot(3,ceil(clunum/3),i)
        %             plot([1:4]',data_zscore(idx==i,:)');
        %             hold on
        %             plot([1 4],[0 0],'k:');
        % %             ylim([-0.5 0.5]);
        %             title(num2str(i));
        % %             axis tight
        %         end
        %         suptitle('K-Means Clustering');
        %
        %         % just plot the controids
        %         figure
        %         for i = 1:clunum
        %             subplot(3,ceil(clunum/3),i)
        %             plot([1:4]',ctrs(i,:)');
        %             hold on
        %             plot([1 4],[0 0],'k:');
        % %             ylim([-0.5 0.5]);
        %             xlim([1 4]);
        %             title(num2str(i));
        % %             axis tight
        % %             axis off
        %         end
        %         suptitle('K-Means Clustering Centroids');
        %
        %         % return to non-zscore data
        %         figure
        %         for i = 1:clunum
        %             subplot(3,ceil(clunum/3),i)
        %             plot([1:4]',data(idx==i,:)');
        %             hold on
        %             plot([1 4],[0 0],'k:');
        %             %             ylim([-0.5 0.5]);
        %             xlim([1 4]);
        %             title(num2str(i));
        %         end
        
        
        %         clu_ind_all = [];
        %         for i = 1:length(cluster)
        %             clu_ind_all = [clu_ind_all cluster{i}];
        %         end
        %         clu_ind_all = sort(clu_ind_all);
        
        %         % select the other and kmeans
        %         clunum = 6;
        %         datain2 = data(isnan(cluster_index),:);
        %         [idx,ctrs] = kmeans(datain2,clunum,'Distance','cityblock','Replicates',10,...
        %             'disp','final','MaxIter',1000,'EmptyAction','singleton');
        %         figure
        %         for i = 1:clunum
        %             subplot(3,ceil(clunum/3),i)
        %             plot([1:4]',datain2(idx==i,:)');
        %             hold on
        %             plot([1 4],[0 0],'k:');
        %             %             ylim([-0.5 0.5]);
        %             title(num2str(i));
        %         end
        %         suptitle('K-Means Clustering');
        
        
        
        
        
        
    end

    function f4p6(debug) % 'Translation and Rotatin choice changed in stimulation'
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        is_separate_stimAmp = 0;
        select_stars_type = [0 1];
        
        % stimtask_type, always choose 4-target fine and coarse task
        stimtask_type_this{1} = 1; % 4 targets fine task
        stimtask_type_this{2} = 4; % 4 targets coarse task
        stimtask_type_this{3} = [1 4]; % combine fine and coarse
        
        % ---------------------- get data ----------------------------
        % 1. 总体看Translation and Rotation choice在电刺激前后选择的变化情况 （不区分细胞）
        % 其实只需要看其中一个的选择变化就可以，因为Translation choice change总是等于 -Rotationchoice change
        % choice_array: eveyy cell, R-L-CCW-CW
        for fc = 1:2 % fine and coarse
            for n = 1:length(loadM_data)
                for cs = 1:2 % control or stim
                    % select near zero degree choice
                    heading_part_num = size(group_result(n).choice_array{stimtask_type_this{fc}}{1,cs},2);
                    rotation_part_num = size(group_result(n).choice_array{stimtask_type_this{fc}}{2,cs},2);
                    
                    % have zero degree?
                    if rem(heading_part_num,2) == 0 % 偶数，no zero degree
                        selectH = [heading_part_num/2,heading_part_num/2+1]; % [1 2 3 4 5 6]中的[3 4]
                    else
                        selectH = [floor(heading_part_num/2):floor(heading_part_num/2)+2]; % [1 2 3 4 5 6 7]中的[3 4 5]
                    end
                    if rem(rotation_part_num,2) == 0 % 偶数，no zero degree
                        selectR = [rotation_part_num/2,rotation_part_num/2+1]; % [1 2 3 4 5 6]中的[3 4]
                    else
                        selectR = [floor(rotation_part_num/2),floor(rotation_part_num/2)+2]; % [1 2 3 4 5 6 7]中的[3 5]，中间0度【4】与heading part重复了
                    end
                    
                    % separate foure choice
                    for i = 1:4 % R L CCW CW
                        % choice around 0 degree
                        group_choice_array_around0{fc,cs}(n,i) = sum(group_result(n).choice_array{stimtask_type_this{fc}}{1,cs}(i,selectH)) + ...
                            sum(group_result(n).choice_array{stimtask_type_this{fc}}{2,cs}(i,selectR)); % i:R-L-CCW-CW
                        
                        % choice across all condition
                        group_choice_array{fc,cs}(n,i) = sum(group_result(n).choice_array{stimtask_type_this{fc}}{1,cs}(i,:)) + ... % Heading part choice
                            sum(group_result(n).choice_array{stimtask_type_this{fc}}{2,cs}(i,:)); % Rotation part choice  % R-L-CCW-CW
                    end
                    
                    % group Translation(L+R) and Rotation(CW+CCW) choice
                    group_choice_array_TR_around0{fc,cs}(n,1) = group_choice_array_around0{fc,cs}(n,1) + group_choice_array_around0{fc,cs}(n,2);
                    group_choice_array_TR_around0{fc,cs}(n,2)= group_choice_array_around0{fc,cs}(n,3) + group_choice_array_around0{fc,cs}(n,4);
                end
                TR_choice_change_prop{fc}(n,:) = (group_choice_array_TR_around0{fc,2}(n,:) - group_choice_array_TR_around0{fc,1}(n,:)) ./ sum(group_choice_array_TR_around0{fc,1}(n,:)); % (stim-contrl)./sum, Translation/Rotation
            end
            mask{fc} = logical(~isnan(TR_choice_change_prop{fc}(:,1)) & methods_of_select{1,1});
        end
        
        % 合并fine and coarse task
        %         TR_choice_change_prop{3} = [TR_choice_change_prop{1}(mask{1},:);TR_choice_change_prop{2}(mask{2},:)];
        %         mask{3} = true(length(TR_choice_change_prop{3}),1);
        
        TR_choice_change_prop{3} = [TR_choice_change_prop{1};TR_choice_change_prop{2}];
        mask{3} = [mask{1};mask{2}];
        
        % ------------------------ Plot ------------------------------
        % 1. 总体看Translation and Rotation choice在电刺激前后选择的变化情况 （不区分细胞）
        % 其实只需要看其中一个的选择变化就可以，因为Translation choice change总是等于 -Rotationchoice change
        tl{1} = 'Fine task';
        tl{2} = 'Coarse task';
        tl{3} = 'Fine + Coarse';
        set(figure(figN),'name','Translation and Rotation choice changed in stimulation','pos',[50,50,1000,300]); clf
        for fc = 1:3
            subplot(1,3,fc)
            plot([1:2]',TR_choice_change_prop{fc}(mask{fc},:)','.-','color',c{99});
            hold on
            choice_change_mean = nanmean(TR_choice_change_prop{fc}(mask{fc},:));
            choice_change_std = nanstd(TR_choice_change_prop{fc}(mask{fc},:));
            choice_change_sem = nanstd(TR_choice_change_prop{fc}(mask{fc},:))/sum(mask{fc});
            errorbar([1:2],choice_change_mean,choice_change_sem,'r-o','linewidth',3);
            set(gca,'xtick',[1:2],'xticklabel',{'Tran Choice changed','Rot Choice changed'});
            xlim([0.9 2.1]);
            ylimit = ylim;
            text(1.5,ylimit(2)-ylimit(2)/10,sprintf('N=%g',sum(mask{fc})));
            title(tl{fc});
            
            % ttest
            [~,p] = ttest(TR_choice_change_prop{fc}(mask{fc},1),TR_choice_change_prop{fc}(mask{fc},2));
            text(1.3,ylimit(2)-ylimit(2)/4,sprintf('ttest p=%0.3g',p));
        end
        figN = figN + 1;
    end


    function f4p7(debug) % 'Choice changed across spiral index'
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        
        % stimtask_type, always choose 4-target fine and coarse task
        stimtask_type_this{1} = 1; % 4 targets fine task
        stimtask_type_this{2} = 4; % 4 targets coarse task
        
        % ---------------------- get data ----------------------------
        select = logical(methods_of_select{1,1});
        
        modulation_index{1} = [group_result.ahead_slope_index]'; % 斜率绝对值
        modulation_index{2} = [group_result.d_prime_index_fine]';
        modulation_index{3} = [group_result.d_prime_index_coarse]';
        modulation_index{4} = [group_result.d_prime_maxDiff_index]'; % 差值最大的两个点,在Group_SpiT中去掉>6std的异常值
        modulation_index{5} = [group_result.d_prime_maxAxis_index]'; % 差值最大的一条轴（180相对两点）,在Group_SpiT中去掉>6std的异常值
        modulation_index{6} = [group_result.d_prime_maxFR_index]'; % FR最大点及其对侧点,在Group_SpiT中去掉>6std的异常值
        modulation_index{7} = [group_result.DDI_inc_index]';
        modulation_index{8} = [group_result.DDI_exc_index]';
        modulation_index{9} = [group_result.HTI_index]';
        modulation_index{10} = [group_result.neu_sen_index_fine]';
        modulation_index{11} = [group_result.neu_sen_index_coarse]';
        
        d_prime_coarse = cell2mat({group_result.d_prime_coarse}');
        modulation_index{13} = abs(d_prime_coarse(:,1)); % heading d prime coarse
        modulation_index{14} = abs(d_prime_coarse(:,2)); % rotation d prime coarse
        
        neu_sen_coarse = cell2mat({group_result.neu_sen_coarse}');
        modulation_index{15} = neu_sen_coarse(:,1); % heading neural sensitivity
        modulation_index{16} = neu_sen_coarse(:,2); % rotation neural sensitivity
        
        modulation_index{17} = [group_result.pair_TRDI]';
        modulation_index{18} = [group_result.nonpair_TRDI]';
        modulation_index{19} = [group_result.nonpair_TRDI2]';
        modulation_index{20} = [group_result.SVM_index]';
        modulation_index{21} = [group_result.AUC_index]';

        
        % tuning index
        spiral_index{1} = cell2mat({group_result.spiral_index}');
        
        spiral_index{2} = cell2mat({group_result.spiral_index_angle}'); % 相关性较差
        spiral_index{3} = cell2mat({group_result.SVM_index}');
        spiral_index = cellfun(@(x) x(select),spiral_index,'uniformoutput',0);
        
        % choice change
        % choice_change_pro(CCI)为正表示该motion type choice在电刺激后比电刺激前选择增增加了，两列分别是T和R，proportion = (stim-ctrl)/ctrl   详细见CueSwitch_performance.m
        % choice_change_index为正表示在电刺激后比电刺激前选择translation target增加了,index = sign * sqrt(chi2), 详细见CueSwitch_performance.m
        for fc = 1:2
            ind = 3;
            choice_change_pro_diff_4T{fc} = choice_change_pro_diff{stimtask_type_this{fc},ind}; % CCI
            choice_change_pro_in0_4T{fc} = choice_change_pro_in0{stimtask_type_this{fc},ind}; 
            choice_change_pro_around0_4T{fc} = choice_change_pro_around0{stimtask_type_this{fc},ind}; 
            
            chi_square_p_diff_4T{fc} = chi_square_p_diff{stimtask_type_this{fc}}(:,ind) ; % p value for chi-square test
            chi_square_p_in0_4T{fc} = chi_square_p_in0{stimtask_type_this{fc}}(:,ind);
            chi_square_p_around0_4T{fc} = chi_square_p_around0{stimtask_type_this{fc}}(:,ind);
        end
        
        
        %         choice_change_pro_diff_4T = cellfun(@(x) x(select,:),choice_change_pro_diff_4T,'uniformoutput',0);
        %         choice_change_pro_in0_4T = cellfun(@(x) x(select,:),choice_change_pro_in0_4T,'uniformoutput',0);
        %         choice_change_pro_around0_4T = cellfun(@(x) x(select,:),choice_change_pro_around0_4T,'uniformoutput',0);
        %
        %         choice_change_index_diff_4T = cellfun(@(x) x(select),choice_change_index_diff_4T,'uniformoutput',0);
        %         choice_change_index_in0_4T = cellfun(@(x) x(select),choice_change_index_in0_4T,'uniformoutput',0);
        %         choice_change_index_around0_4T = cellfun(@(x) x(select),choice_change_index_around0_4T,'uniformoutput',0);
        %
        %         chi_square_p_diff_4T = cellfun(@(x) x(select),chi_square_p_diff_4T,'uniformoutput',0);
        %         chi_square_p_in0_4T = cellfun(@(x) x(select),chi_square_p_in0_4T,'uniformoutput',0);
        %         chi_square_p_around0_4T = cellfun(@(x) x(select),chi_square_p_around0_4T,'uniformoutput',0);
        
        % select and group together: CCI = (stimR-ctrlR)/ctrlR, (stimT-ctrlT)/ctrlT   
        %  (ind=3时，第一列是R chocie变化，第二列是T choice变化）
        CCI{1} = cellfun(@(x) x(select,:),choice_change_pro_diff_4T,'uniformoutput',0);
        CCI{2} = cellfun(@(x) x(select,:),choice_change_pro_in0_4T,'uniformoutput',0);
        CCI{3} = cellfun(@(x) x(select,:),choice_change_pro_around0_4T,'uniformoutput',0);
        
        % do not use this index
        %         choice_change_index{1} = cellfun(@(x) x(select),choice_change_index_diff_4T,'uniformoutput',0);
        %         choice_change_index{2} = cellfun(@(x) x(select),choice_change_index_in0_4T,'uniformoutput',0);
        %         choice_change_index{3} = cellfun(@(x) x(select),choice_change_index_around0_4T,'uniformoutput',0);
        
        %         % use a new choice_change_index: choice_change_index = Tchoice_change_proportion - Rchoice_change_proportion;
        %         % index = 0: T and R chocie change similar
        %         for j = 1:3 % difficult trial, 0 degree, around 0 degree
        %             for fc = 1:2 % 1:only fine, 2:only coarse
        %                 % i = 1:2 % R or T (第一列是R chocie变化，第二列是T choice变化）
        %                 T-R_CCI{j}{fc} = CCI{j}{fc}(:,2) - CCI{j}{fc}(:,1);  % T-R
        %             end
        %         end
                 
        chi_square_p{1} = cellfun(@(x) x(select),chi_square_p_diff_4T,'uniformoutput',0);
        chi_square_p{2} = cellfun(@(x) x(select),chi_square_p_in0_4T,'uniformoutput',0);
        chi_square_p{3} = cellfun(@(x) x(select),chi_square_p_around0_4T,'uniformoutput',0);
        
        % behavioral performance
        mt_cr = [];
        for fc = 1:length(stimtask_type_this)
            mt_cr{fc} = ctrl_CR{stimtask_type_this{fc}}(:,3); % control task CR
        end
        mt_cr = cellfun(@(x) x(select,:),mt_cr,'uniformoutput',0);
        
        % choice-change index, CCI: stim/ctrl % 第一列T，第二列R
        for fc = 1:2
            TRchoice_change{fc} = choice_change_around0{stimtask_type_this{fc}}(select,:);
        end
        
        
        % ------------------------ Plot ------------------------------
        % % Choice-change index, CCI: stim-ctrl / ctrl
        % select with behavioral performance
        mt_performance_cri = [65 100]; % only select performance > 65% and < 100% correct rate 
        
%         for j = 1:3  % difficult trial, 0 degree, around 0 degree
        for j = 3
            set(figure(figN),'name','Translation and Rotation choice changed proportion across spiral index','pos',[600*j-600,50,630,950]); clf
            ha = tight_subplot(3,2,[.1 .1],[.07 0.07],[.2 .05],1);
            for i = 1:2 % R or T (第一列是R chocie变化，第二列是T choice变化）
                for fc = 1:3 % 1:only fine, 2:only coatse, 3: combine
                    axes(ha(i+(fc-1)*2))
                    if fc~=3
                        mt_perform = mt_cr{fc};
                        perform_select = logical(mt_perform>=mt_performance_cri(1) & mt_perform<=mt_performance_cri(2));
                        tuning_index = spiral_index{1}(perform_select);
                        choice_change = CCI{j}{fc}(perform_select,i);
                        sig_site = logical(chi_square_p{j}{fc}<0.05);
                        sig_site = sig_site(perform_select);
                    else % combine fine and coarse
                        mt_perform = [mt_cr{1};mt_cr{2}];
                        perform_select = logical(mt_perform>=mt_performance_cri(1) & mt_perform<=mt_performance_cri(2));
                        tuning_index = [spiral_index{1};spiral_index{1}];
                        tuning_index = tuning_index(perform_select);
                        choice_change = [CCI{j}{1}(:,i);CCI{j}{2}(:,i)];
                        choice_change = choice_change(perform_select);
                        sig_site = logical([chi_square_p{j}{1}; chi_square_p{j}{2}]<0.05);
                        sig_site = sig_site(perform_select);
                    end
                    hold on
                    plot(tuning_index,choice_change,'o','color',c{3-i});
                    plot(tuning_index(sig_site),choice_change(sig_site),'o','color',c{3-i},'markerfacecolor',c{3-i});
                    plot([-1 1],[0 0],'k:');
                    ylabel('Choice change');
                    xlabel('Spiral index');
                    ylim([-1.5 1.5]);
                    ylimit = ylim;
                    ylimit = max(abs(ylimit));
                    ylim([-ylimit ylimit]);
                    
                    % pearson correlation
%                     [~,p_pear] = plot_corr_line(tuning_index,choice_change,'MethodOfCorr','Pearson','FittingMethod',2); % Type 2: Perpendicular error (solved by PCA，垂直距离最小)
                    [~,p_spear] = plot_corr_line(tuning_index,choice_change,'MethodOfCorr','Spearman','FittingMethod',2);
%                     text(-0.9,ylimit,sprintf('p. corr. p = %0.3g',p_pear));
                    text(-0.9,ylimit-ylimit/6,sprintf('s. corr. p = %0.3g',p_spear));
                    
                    [r_spear] = corr(tuning_index,choice_change,'type','spearman');
                    text(-0.9,ylimit,sprintf('s. corr. R = %0.3g',r_spear));
                    
%                     if i == 1 && fc ==1
                        text(-0.9,ylimit-ylimit/3,sprintf('CR = [%g-%g]%%, N = %g',mt_performance_cri(1),mt_performance_cri(2), sum(perform_select)));
%                     end
                        
                    
                    if i==1 && fc == 1
                        text(-1.8,0,'Fine task','rotation',90,'HorizontalAlignment','center');
                    elseif i == 1 && fc == 2
                        text(-1.8,0,'Coarse task','rotation',90,'HorizontalAlignment','center');
                    elseif i == 1 && fc == 3
                        text(-1.8,0,'Fine+Coarse task','rotation',90,'HorizontalAlignment','center');
                    end
                end
            end
           
            if j == 1 % difficult trial, 0 degree, around 0 degree
                suptitle('Difficult trial');
            elseif j == 2
                suptitle('0 degree');
            else
                suptitle('Around 0 degree');
            end
            SetFigure(13);
            figN = figN + 1;
        end
        %         [~,~,~,til3] = SetTitle(stimtask_type_this,monkey_included_for_analysis);
        %         til_this = 'Around 0° choice changed across spiral index (4-Target)';
        %         til_this = strcat(til_this,32,til3);
        %         suptitle(til_this);
        
        
%         % choice change index
%         for j = 1:3 % difficult trial, 0 degree, around 0 degree
%             set(figure(figN),'name','Translation and Rotation choice changed index across spiral index','pos',[300*j+800,50,315,950]); clf
%             ha = tight_subplot(3,1,[.1 .1],[.07 0.07],[.2 .05],1);
%             for fc = 1:3 % 1:only fine, 2:only coatse, 3: combine
%                 axes(ha(fc))
%                 if fc~=3
%                     tuning_index = spiral_index{1};
%                     choice_change = choice_change_index{j}{fc};
%                     sig_site = logical(chi_square_p{j}{fc}<0.05);
%                 else % combine fine and coarse
%                     tuning_index = [spiral_index{1};spiral_index{1}];
%                     choice_change = [choice_change_index{j}{1};choice_change_index{j}{2}];
%                     sig_site = logical([chi_square_p{j}{1}; chi_square_p{j}{2}]<0.05);
%                 end
%                 hold on
%                 plot(tuning_index,choice_change,'o','color',c{99});
%                 plot(tuning_index(sig_site),choice_change(sig_site),'o','color',c{99},'markerfacecolor',c{99});
%                 plot([-1 1],[0 0],'k:');
%                 ylabel('Choice changed index');
%                 xlabel('Spiral index');
%                 ylimit = ylim;
%                 ylimit = max(abs(ylimit));
%                 ylim([-ylimit ylimit]);
%                 
%                 % pearson correlation
%                 [~,p_pear] = plot_corr_line(tuning_index,choice_change,'MethodOfCorr','Pearson');
%                 [~,p_spear] = plot_corr_line(tuning_index,choice_change,'MethodOfCorr','Spearman');
%                 text(-0.9,ylimit,sprintf('p. corr. p = %0.3g',p_pear));
%                 text(-0.9,ylimit-ylimit/6,sprintf('s. corr. p = %0.3g',p_spear));
%                 
%                 if fc == 1
%                     text(-1.8,0,'Fine task','rotation',90,'HorizontalAlignment','center');
%                 elseif fc == 2
%                     text(-1.8,0,'Coarse task','rotation',90,'HorizontalAlignment','center');
%                 elseif fc == 3
%                     text(-1.8,0,'Fine+Coarse task','rotation',90,'HorizontalAlignment','center');
%                 end
%             end
%             if j == 1 % difficult trial, 0 degree, around 0 degree
%                 suptitle('Difficult trial');
%             elseif j == 2
%                 suptitle('0 degree');
%             else
%                 suptitle('Around 0 degree');
%             end
%             SetFigure(13);
%             figN = figN + 1;
%         end
        
        % choice change: stim / ctrl
        set(figure(figN),'name','Translation and Rotation choice changed (stim/ctrl) across spiral index','pos',[300*4+800,50,315,950]); clf
        ha = tight_subplot(3,1,[.1 .1],[.07 0.07],[.2 .05],1);
        for fc = 1:3 % 1:only fine, 2:only coatse, 3: combine
            axes(ha(fc))
            if fc~=3
                mt_perform = mt_cr{fc};
                perform_select = logical(mt_perform>=mt_performance_cri(1) & mt_perform<=mt_performance_cri(2));
                tuning_index = spiral_index{1}(perform_select);
                choice_change = TRchoice_change{fc}(perform_select,:); 
                sig_site = logical(chi_square_p{j}{fc}<0.05);
                sig_site = sig_site(perform_select);
            else % combine fine and coarse
                mt_perform = [mt_cr{1};mt_cr{2}];
                perform_select = logical(mt_perform>=mt_performance_cri(1) & mt_perform<=mt_performance_cri(2));
                tuning_index = [spiral_index{1};spiral_index{1}];
                tuning_index = tuning_index(perform_select);
                choice_change = [TRchoice_change{1};TRchoice_change{2}]; 
                choice_change = choice_change(perform_select,:);
                sig_site = logical([chi_square_p{j}{1}; chi_square_p{j}{2}]<0.05);
                sig_site = sig_site(perform_select);
            end
            hold on
            for hr = 1:2
                if hr == 1
                    plot(tuning_index(~sig_site),choice_change(~sig_site,hr),'o','color','m');
                else
                    plot(tuning_index(~sig_site),choice_change(~sig_site,hr),'o','color','g');
                end
                plot(tuning_index(sig_site),choice_change(sig_site,hr),'o','color',c{hr},'markerfacecolor',c{hr});
            end
            axis square
            ylabel('Choice changed');
            xlabel('Spiral index');
            ylimit = ylim;
            ylimit = max(abs(ylimit));
            xlim([-1 1]);
            ylim([0 ylimit]);
            
            % pearson correlation
            for hr = 1:2
                [r_pear(fc,hr),p_pear(fc,hr)] = plot_corr_line(tuning_index,choice_change(:,hr),'MethodOfCorr','Pearson','LineStyles',{'-','color',c{hr}},'LineWidth',2);
                [r_spear(fc,hr),p_spear(fc,hr)] = plot_corr_line(tuning_index,choice_change(:,hr),'MethodOfCorr','Spearman','LineStyles',{'-','color',c{hr}},'LineWidth',2);
            end
            
            %             text(-0.9,ylimit,sprintf('p. corr. p = %0.3g',p_pear));
            %             text(-0.9,ylimit-ylimit/6,sprintf('s. corr. p = %0.3g',p_spear));
            
            if fc == 1
                text(-1.8,0,'Fine task','rotation',90,'HorizontalAlignment','center');
            elseif fc == 2
                text(-1.8,0,'Coarse task','rotation',90,'HorizontalAlignment','center');
            elseif fc == 3
                text(-1.8,0,'Fine+Coarse task','rotation',90,'HorizontalAlignment','center');
            end
        end
        suptitle('Around 0 degree');
        SetFigure(13);
        figN = figN + 1;
    end

    function f4p16(debug) % 'Choice changed across cell type'
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        % stimtask_type, always choose 4-target fine and coarse task
        stimtask_type_this{1} = 1; % 4 targets fine task
        stimtask_type_this{2} = 4; % 4 targets coarse task
        
        % ---------------------- get data ----------------------------
        select = logical(methods_of_select{1,1});
        cell_type = [group_result.cell_type];
        cell_type = cell_type(select)';
    
        spiral_index = cell2mat({group_result.spiral_index}');
        spiral_index = spiral_index(select);
        
        % behavioral performance
        mt_cr = [];
        for fc = 1:length(stimtask_type_this)
            mt_cr{fc} = ctrl_CR{stimtask_type_this{fc}}(:,3); % control task CR
        end
        mt_cr = cellfun(@(x) x(select,:),mt_cr,'uniformoutput',0);

        % choice change
        % choice_change_pro为正表示该motion type choice在电刺激后比电刺激前选择增增加了，两列分别是T和R，proportion = (stim-ctrl)/ctrl   详细见CueSwitch_performance.m
        % choice_change_index为正表示在电刺激后比电刺激前选择translation target增加了,index = sign * sqrt(chi2), 详细见CueSwitch_performance.m
        for fc = 1:2
            ind = 3;
            choice_change_pro_diff_4T{fc} = choice_change_pro_diff{stimtask_type_this{fc},ind}; % 两列分别是T和R
            choice_change_pro_in0_4T{fc} = choice_change_pro_in0{stimtask_type_this{fc},ind}; % 两列分别是T和R
            choice_change_pro_around0_4T{fc} = choice_change_pro_around0{stimtask_type_this{fc},ind}; % 两列分别是T和R
            
            choice_change_index_diff_4T{fc} = choice_change_index_diff{stimtask_type_this{fc}}(:,ind);
            choice_change_index_in0_4T{fc} = choice_change_index_in0{stimtask_type_this{fc}}(:,ind);
            choice_change_index_around0_4T{fc} = choice_change_index_around0{stimtask_type_this{fc}}(:,ind);
            
            chi_square_p_diff_4T{fc} = chi_square_p_diff{stimtask_type_this{fc}}(:,ind) ; % p value for chi-square test
            chi_square_p_in0_4T{fc} = chi_square_p_in0{stimtask_type_this{fc}}(:,ind);
            chi_square_p_around0_4T{fc} = chi_square_p_around0{stimtask_type_this{fc}}(:,ind);
        end
        
        % select and group together: (stimR-ctrlR)/ctrlR, (stimT-ctrlT)/ctrlT   
        choice_change_proportion{1} = cellfun(@(x) x(select,:),choice_change_pro_diff_4T,'uniformoutput',0);
        choice_change_proportion{2} = cellfun(@(x) x(select,:),choice_change_pro_in0_4T,'uniformoutput',0);
        choice_change_proportion{3} = cellfun(@(x) x(select,:),choice_change_pro_around0_4T,'uniformoutput',0);
        chi_square_p{1} = cellfun(@(x) x(select),chi_square_p_diff_4T,'uniformoutput',0);
        chi_square_p{2} = cellfun(@(x) x(select),chi_square_p_in0_4T,'uniformoutput',0);
        chi_square_p{3} = cellfun(@(x) x(select),chi_square_p_around0_4T,'uniformoutput',0);
        
        
        
        % ------------------------ Plot ------------------------------
        % chocie change proportion
        % select with behavioral performance
        mt_performance_cri = [0 100]; % only select performance > 65% and < 100% correct rate
        
%         for j = 1:3  % difficult trial, 0 degree, around 0 degree
                    for j = 3
            set(figure(figN),'name','Translation and Rotation choice changed across cell type','pos',[600*j-600,50,630,950]); clf
            ha = tight_subplot(3,2,[.1 .1],[.07 0.07],[.2 .05],1);
            for i = 1:2 % R or T (第一列是R chocie变化，第二列是T choice变化）
                for fc = 1:3 % 1:only fine, 2:only coatse, 3: combine
                    axes(ha(i+(fc-1)*2))
                    if fc~=3
                        mt_perform = mt_cr{fc};
                        perform_select = logical(mt_perform>=mt_performance_cri(1) & mt_perform<=mt_performance_cri(2));
                        cell_type_this = cell_type(perform_select);
                        choice_change = choice_change_proportion{j}{fc}(perform_select,i);
                        sig_site = logical(chi_square_p{j}{fc}<0.05);
                        sig_site = sig_site(perform_select);
                    else % combine fine and coarse
                        mt_perform = [mt_cr{1};mt_cr{2}];
                        perform_select = logical(mt_perform>=mt_performance_cri(1) & mt_perform<=mt_performance_cri(2));
                        cell_type_this = [cell_type;cell_type];
                        cell_type_this = cell_type_this(perform_select);
                        choice_change = [choice_change_proportion{j}{1}(:,i);choice_change_proportion{j}{2}(:,i)];
                        choice_change = choice_change(perform_select);
                        sig_site = logical([chi_square_p{j}{1}; chi_square_p{j}{2}]<0.05);
                        sig_site = sig_site(perform_select);
                    end
                    
                    hold on
                    choice_change_celltype = nan(300,3);
                    ct = [2 3 1]; % r cell, spi cell, t cell
                    for cti = 1:3 % t cell, r cell, spi cell
                        choice_change_celltype(1:sum(cell_type_this == ct(cti)),cti) = choice_change(cell_type_this == ct(cti));
                        choice_change_mean_celltype(cti) = nanmean(choice_change_celltype(:,cti));
                        choice_change_median_celltype(cti) = nanmedian(choice_change_celltype(:,cti));
                        choice_change_sem_celltype(cti) = nanstd(choice_change_celltype(:,cti)) / sqrt(sum(cell_type_this == ct(cti)));
                    end
                    plot([1 2 3],choice_change_celltype,'k.')
                    %                     boxplot(choice_change_celltype); % r, spi, t
                    errorbar([1 2 3],choice_change_mean_celltype,choice_change_sem_celltype,'o-','color',c{3-i}); % mean
                    %                     errorbar([1 2 3],choice_change_median_celltype,choice_change_sem_celltype,'ko--'); %v meadian
                    
                    % 与0比较
                    for cti = 1:3 % t cell, r cell, spi cell
                        [~, p(cti)] = ttest(choice_change_celltype(~isnan(choice_change_celltype(:,cti)),cti));
                    end
                    text(0.6,-0.8,sprintf('p = %.3g  %.3g  %.3g',p(1),p(2),p(3)));
                    
                    xlim([0.5 3.5]);
                    ylim([-1.5 1.5]);
                    xticks([1 2 3]);
                    xticklabels({'R cell','Spi cell','T cell'});
                    
                    if i==1 && fc == 1
                        text(-0.3,0,'Fine task','rotation',90,'HorizontalAlignment','center');
                    elseif i == 1 && fc == 2
                        text(0-0.3,0,'Coarse task','rotation',90,'HorizontalAlignment','center');
                    elseif i == 1 && fc == 3
                        text(-0.3,0,'Fine+Coarse task','rotation',90,'HorizontalAlignment','center');
                    end
                    
                    if fc == 1
                        if i == 1
                            title('Rotation choice change');
                        elseif i == 2
                            title('Translation choice change');
                        end
                    end
                end
            end
            
            if j == 1 % difficult trial, 0 degree, around 0 degree
                suptitle('Difficult trial');
            elseif j == 2
                suptitle('0 degree');
            else
                suptitle('Around 0 degree');
            end
            SetFigure(13);
            figN = figN + 1;
        end 
    end

    function f4p9(debug) % 'near 0° choice changed'
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        select_stars_type = [0 1];
        
        % stimtask_type, always choose 4-target fine and coarse task
        stimtask_type_this{1} = 1; % 4 targets fine task
        stimtask_type_this{2} = 4; % 4 targets coarse task
        
        % ---------------------- get data ----------------------------
        select = logical(methods_of_select{1});
        
        % tuning index
        %         spiral_index = cell2mat({group_result.spiral_index}');
        %         spiral_index = cell2mat({group_result.spiral_index_angle}'); % 相关性较差
        %         spiral_index = cell2mat({group_result.SVM_index}');
        for fc = 1:2 % fine and coarse
            chi_square_p{fc} = chi_square_p_around0{stimtask_type_this{fc}}(:,[1 2]);
            for n = 1:length(loadM_data)
                if select(n)
                    for cs = 1:2 % control or stim
                        % select near zero degree choice
                        heading_part_num = size(group_result(n).choice_array{stimtask_type_this{fc}}{1,cs},2);
                        rotation_part_num = size(group_result(n).choice_array{stimtask_type_this{fc}}{2,cs},2);
                        
                        % have zero degree?
                        if rem(heading_part_num,2) == 0 % 偶数，no zero degree
                            selectH = [heading_part_num/2,heading_part_num/2+1]; % [1 2 3 4 5 6]中的[3 4]
                        else
                            selectH = [floor(heading_part_num/2):floor(heading_part_num/2)+2]; % [1 2 3 4 5 6 7]中的[3 4 5]
                        end
                        if rem(rotation_part_num,2) == 0 % 偶数，no zero degree
                            selectR = [rotation_part_num/2,rotation_part_num/2+1]; % [1 2 3 4 5 6]中的[3 4]
                        else
                            selectR = [floor(rotation_part_num/2),floor(rotation_part_num/2)+2]; % [1 2 3 4 5 6 7]中的[3 5]，中间0度【4】与heading part重复了
                        end
                        
                        % separate foure choice，选每个target的个数
                        for i = 1:4 % R L CCW CW
                            % choice around 0 degree
                            group_choice_array_around0{cs}(i) = sum(group_result(n).choice_array{stimtask_type_this{fc}}{1,cs}(i,selectH)) + ...
                                sum(group_result(n).choice_array{stimtask_type_this{fc}}{2,cs}(i,selectR)); % i:R-L-CCW-CW
                        end
                        
                        % group Translation(L+R) and Rotation(CW+CCW) choice
                        group_choice_array_TR_around0{cs}(1) = group_choice_array_around0{cs}(1) + group_choice_array_around0{cs}(2); % 选Translation target个数
                        group_choice_array_TR_around0{cs}(2)= group_choice_array_around0{cs}(3) + group_choice_array_around0{cs}(4); % 选Rotation target个数
                    end
                    
                    % ************************************************************************************************
                    % 相比于control trial，stinm trial选择Translation或选择Rotation target的变化：(stim-contrl)./contrl
                    % ************************************************************************************************
                    %                     TR_choice_change_index_around0{fc}(n,:) = (group_choice_array_TR_around0{2} - group_choice_array_TR_around0{1}) ./ sum(group_choice_array_TR_around0{1});
                    TR_choice_change_index_around0{fc}(n,:) = (group_choice_array_TR_around0{2} - group_choice_array_TR_around0{1}) ./ (group_choice_array_TR_around0{1});
                else
                    TR_choice_change_index_around0{fc}(n,:) = nan(1,2);
                end
            end
        end
        
        %         % 人为设定 |TR_choice_change_index_around0|>0.3 or chi_square_p<0.05
        %         for fc = 1:2 % fine and coarse
        %             select = any(chi_square_p{fc}<0.05,2);
        %             a = find(select)
        %             a(17)
        %             a(3)
        %             choice_change_this = TR_choice_change_index_around0{fc}(select,:);
        %
        %         end
        
    end

    function f4p8(debug) % 'Choice changed between control and stim and chi-square test'
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        is_separate_stimAmp = 0;
        select_stars_type = [0 1];
        
        % ---------------------- get data ----------------------------
        % select for stars_type: 0, 1, [0 1](combine)
        stars_type = cell2mat({group_result.stars_type})';
        select1 = false(size(stars_type)); % select for stars_type
        for st = 1:length(select_stars_type)
            select1_temp = logical(stars_type == select_stars_type(st));
            select1 = logical(select1_temp | select1);
        end
        select = logical(methods_of_select{1} & select1);
        select_this = repmat(select,length(cell2mat(stimtask_type)),1);
        
        choice_change_index_diff_all = [];
        choice_change_index_in0_all = [];
        choice_change_index_around0_all = [];
        
        chi_square_p_diff_all = [];
        chi_square_p_in0_all = [];
        chi_square_p_around0_all = [];
        
        choice_change_pro_diff_all = cell(1,3);
        choice_change_pro_in0_all = cell(1,3);
        choice_change_pro_around0_all = cell(1,3);
        
        for fc = 1:length(stimtask_type) % fine or coarse，需要combine一起，所以最后的index没有fr
            if ~isempty(stimtask_type{fc})
                for i = 1:length(stimtask_type{fc}) % 具体task, 需要combine一起，所以最后的index没有i
                    choice_change_index_diff_all = [choice_change_index_diff_all;choice_change_index_diff{stimtask_type{fc}(i)}];
                    choice_change_index_in0_all = [choice_change_index_in0_all;choice_change_index_diff{stimtask_type{fc}(i)}];
                    choice_change_index_around0_all = [choice_change_index_around0_all;choice_change_index_diff{stimtask_type{fc}(i)}];
                    
                    chi_square_p_diff_all = [chi_square_p_diff_all;chi_square_p_diff{stimtask_type{fc}(i)}];
                    chi_square_p_in0_all = [chi_square_p_in0_all;chi_square_p_in0{stimtask_type{fc}(i)}];
                    chi_square_p_around0_all = [chi_square_p_around0_all;chi_square_p_around0{stimtask_type{fc}(i)}];
                    
                    for ind = 1:3
                        choice_change_pro_diff_all{ind} = [choice_change_pro_diff_all{ind};choice_change_pro_diff{stimtask_type{fc}(i),ind}];
                        choice_change_pro_in0_all{ind} = [choice_change_pro_in0_all{ind};choice_change_pro_in0{stimtask_type{fc}(i),ind}];
                        choice_change_pro_around0_all{ind} = [choice_change_pro_around0_all{ind};choice_change_pro_around0{stimtask_type{fc}(i),ind}];
                    end
                end
            end
        end
        
        choice_change_index{1} = choice_change_index_diff_all(select_this,:); % difficult condition
        choice_change_index{2} = choice_change_index_in0_all(select_this,:); % zero condition
        choice_change_index{3} = choice_change_index_around0_all(select_this,:); % around zero condition
        
        chi_square_p{1} = chi_square_p_diff_all(select_this,:); % difficult condition
        chi_square_p{2} = chi_square_p_in0_all(select_this,:); % zero condition
        chi_square_p{3} = chi_square_p_around0_all(select_this,:); % around zero condition
        
        for ind = 1:3
            choice_change_pro{1,ind} = choice_change_pro_diff_all{ind}(select_this,:); % difficult condition
            choice_change_pro{2,ind} = choice_change_pro_in0_all{ind}(select_this,:); % zero condition
            choice_change_pro{3,ind} = choice_change_pro_around0_all{ind}(select_this,:); % around zero condition
        end
        
        % ------------------------ Plot ------------------------------
        bin_num = 8; % 偶数时，0是边界处
        for i = 1:3 % diff, 0, around 0
            for ind = 1:3 % h r mt
                xlimit_all{i}(ind) = max(abs(choice_change_index{i}(:,ind)));
                xlimit_all{i}(ind) = ceil(xlimit_all{i}(ind)*10)/10;
                xbins{i,ind} = linspace(-xlimit_all{i}(ind),xlimit_all{i}(ind),bin_num);
                
                S_unit{i,ind} = chi_square_p{i}(:,ind)<0.05;
                NS_unit{i,ind} = chi_square_p{i}(:,ind)>=0.05;
            end
        end
        
        for i = 1:3
            if i == 1
                set(figure(figN),'Position',[100,650 900,200], 'Name', 'Choice change in difficult trial');
            elseif i == 2
                set(figure(figN+1),'Position',[100,350 900,200], 'Name', 'Choice change in 0 degree');
            elseif i == 3
                set(figure(figN+2),'Position',[100,50 900,200], 'Name', 'Choice change around 0 degree');
            end
            
            for ind = 1:3
                if ~isnan(xlimit_all{i}(ind))
                    nbins = hist(choice_change_index{i}(:,ind),xbins{i,ind});
                    histS = hist(choice_change_index{i}(S_unit{i,ind},ind),xbins{i,ind});
                    histNS = hist(choice_change_index{i}(NS_unit{i,ind},ind),xbins{i,ind});
                    
                    subplot(1,3,ind)
                    hbars = bar(xbins{i,ind},[histS' histNS'],1,'stacked','LineWidth',1.5);
                    set(hbars,'EdgeColor','k','FaceColor',c{ind});
                    set(hbars(2),'FaceColor','none');
                    
                    
                    xlim([-xlimit_all{i}(ind)-xlimit_all{i}(ind)/3 xlimit_all{i}(ind)+xlimit_all{i}(ind)/3]);
                    ylimit_all = max(nbins);
                    ylim([0 ylimit_all+ylimit_all/5]);
                    xlabel('choice change index');
                    ylabel('Cases');
                    
                    xtick_temp = get(gca,'xtick');
                    xtick_this = max(abs(xtick_temp));
                    set(gca,'xtick',[-xtick_this 0 xtick_this]);
                    if ind == 1
                        set(gca,'xticklabel',{'More Left','0','More Right'});
                    elseif ind == 2
                        set(gca,'xticklabel',{'More CW','0','More CCW'});
                    elseif ind == 3
                        set(gca,'xticklabel',{'More Rot','0','More Trans'});
                    end
                    
                    SetFigure(10);
                end
            end
        end
        figN = figN + 3;
        
        % chocie change across Spiral index only for MT task
        spiral_index = cell2mat({group_result.spiral_index}');
        spiral_index = repmat(spiral_index,length(cell2mat(stimtask_type)),1);
        spiral_index = spiral_index(select_this);
        % ------------------------ Plot ------------------------------
        set(figure(figN),'Position',[1050,50 600,900], 'Name', 'Choice change proportion');
        for i = 1:3
            for ind = 3 % only mt task
                for lr = 1:2 % 1 for Rotation, 2 for Translation
                    choice_change_index_this = choice_change_pro{i,ind}(:,lr);
                    spiral_index_this = spiral_index;
                    % choice_change_index_this = choice_change_pro{i,ind}(chi_square_p{i}(:,ind)<0.05,lr);
                    % spiral_index_this = spiral_index(chi_square_p{i}(:,ind)<0.05);
                    
                    subplot(3,2,lr+(i-1)*2)
                    plot(spiral_index_this,choice_change_index_this,'o','color',c{lr});
                    hold on
                    plot([-1 1],[0 0],'k:');
                    ylabel('Choice changed');
                    xlabel('Spiral index');
                    ylimit = ylim;
                    ylimit = max(abs(ylimit));
                    ylim([-ylimit ylimit]);
                    
                    % pearson correlation
                    [~,p_pear] = plot_corr_line(spiral_index_this,choice_change_index_this,'MethodOfCorr','Pearson');
                    [~,p_spear] = plot_corr_line(spiral_index_this,choice_change_index_this,'MethodOfCorr','Spearman');
                    text(-0.9,ylimit,sprintf('p. corr. p = %0.3g',p_pear));
                    text(-0.9,ylimit-ylimit/6,sprintf('s. corr. p = %0.3g',p_spear));
                    
                    if i == 1
                        title('In difficult degree');
                    elseif i == 2
                        title('In 0 degree');
                    elseif i == 3
                        title('Around 0 degree');
                    end
                end
            end
        end
        figN = figN + 3; 
    end

    function f4p15(debug) % 'Choice changed with some T-R ratio'
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        
        % stimtask_type, always choose 4-target fine and coarse task
        stimtask_type_this{1} = 1; % 4 targets fine task
        stimtask_type_this{2} = 4; % 4 targets coarse task
        
        % ---------------------- get data ----------------------------
        select = logical(methods_of_select{1,1});
        
        modulation_index{1} = cell2mat({group_result.fr_axis_ratio}');
        modulation_index{1}(modulation_index{1}>50) = nan;
        modulation_index{2} = cell2mat({group_result.vec_axis_ratio}');
        modulation_index{3} = cell2mat({group_result.d_prime_ratio}');
        modulation_index = cellfun(@(x) x(select),modulation_index,'uniformoutput',0);
        
        % choice change
        % choice_change_pro为正表示该motion type choice在电刺激后比电刺激前选择增增加了，两列分别是T和R，proportion = (stim-ctrl)/ctrl   详细见CueSwitch_performance.m
        % choice_change_index为正表示在电刺激后比电刺激前选择translation target增加了,index = sign * sqrt(chi2), 详细见CueSwitch_performance.m
        for fc = 1:2
            ind = 3;
            choice_change_pro_diff_4T{fc} = choice_change_pro_diff{stimtask_type_this{fc},ind}; % 两列分别是T和R
            choice_change_pro_in0_4T{fc} = choice_change_pro_in0{stimtask_type_this{fc},ind}; % 两列分别是T和R
            choice_change_pro_around0_4T{fc} = choice_change_pro_around0{stimtask_type_this{fc},ind}; % 两列分别是T和R
            
            choice_change_index_diff_4T{fc} = choice_change_index_diff{stimtask_type_this{fc}}(:,ind);
            choice_change_index_in0_4T{fc} = choice_change_index_in0{stimtask_type_this{fc}}(:,ind);
            choice_change_index_around0_4T{fc} = choice_change_index_around0{stimtask_type_this{fc}}(:,ind);
            
            chi_square_p_diff_4T{fc} = chi_square_p_diff{stimtask_type_this{fc}}(:,ind) ;
            chi_square_p_in0_4T{fc} = chi_square_p_in0{stimtask_type_this{fc}}(:,ind);
            chi_square_p_around0_4T{fc} = chi_square_p_around0{stimtask_type_this{fc}}(:,ind);
        end
        % select and group together
        choice_change_proportion{1} = cellfun(@(x) x(select,:),choice_change_pro_diff_4T,'uniformoutput',0);
        choice_change_proportion{2} = cellfun(@(x) x(select,:),choice_change_pro_in0_4T,'uniformoutput',0);
        choice_change_proportion{3} = cellfun(@(x) x(select,:),choice_change_pro_around0_4T,'uniformoutput',0);
        chi_square_p{1} = cellfun(@(x) x(select),chi_square_p_diff_4T,'uniformoutput',0);
        chi_square_p{2} = cellfun(@(x) x(select),chi_square_p_in0_4T,'uniformoutput',0);
        chi_square_p{3} = cellfun(@(x) x(select),chi_square_p_around0_4T,'uniformoutput',0);
        
        % behavioral performance
        mt_cr = [];
        for fc = 1:length(stimtask_type_this)
            mt_cr{fc} = ctrl_CR{stimtask_type_this{fc}}(:,3); % control task CR
        end
        mt_cr = cellfun(@(x) x(select,:),mt_cr,'uniformoutput',0);

        % select with behavioral performance
        mt_performance_cri = [65 100]; % only select performance > 65% and < 100% correct rate
        
        select_index = 3;
        
        % ------------------------ Plot ------------------------------
        for j = 1:3 % difficult trial, 0 degree, around 0 degree
            set(figure(figN),'name','Translation and Rotation choice changed across T-R ratio','pos',[600*j-600,50,630,950]); clf
            ha = tight_subplot(3,2,[.1 .1],[.07 0.07],[.2 .05],1);
            for i = 1:2 % R or T (第一列是R chocie变化，第二列是T choice变化）
                for fc = 1:3 % 1:only fine, 2:only coatse, 3: combine
                    axes(ha(i+(fc-1)*2))
                    if fc~=3
                        mt_perform = mt_cr{fc};
                        perform_select = logical(mt_perform>=mt_performance_cri(1) & mt_perform<=mt_performance_cri(2));
                        tuning_index = modulation_index{select_index}(perform_select);
                        choice_change = choice_change_proportion{j}{fc}(perform_select,i);
                        sig_site = logical(chi_square_p{j}{fc}<0.05);
                        sig_site = sig_site(perform_select);
                    else % combine fine and coarse
                        mt_perform = [mt_cr{1};mt_cr{2}];
                        perform_select = logical(mt_perform>=mt_performance_cri(1) & mt_perform<=mt_performance_cri(2));
                        tuning_index = [modulation_index{select_index};modulation_index{select_index}];
                        tuning_index = tuning_index(perform_select);
                        choice_change = [choice_change_proportion{j}{1}(:,i);choice_change_proportion{j}{2}(:,i)];
                        choice_change = choice_change(perform_select);
                        sig_site = logical([chi_square_p{j}{1}; chi_square_p{j}{2}]<0.05);
                        sig_site = sig_site(perform_select);
                    end
                    % log10 scale
                    tuning_index = log10(tuning_index);
                    hold on
                    plot(tuning_index,choice_change,'o','color',c{3-i});
                    plot(tuning_index(sig_site),choice_change(sig_site),'o','color',c{3-i},'markerfacecolor',c{3-i});
                    plot([-1 1],[0 0],'k:');
                    ylabel('Choice change');
                    if select_index == 1
                        xlabel('log(FR axis ratio)');
                    elseif select_index == 2
                        xlabel('log(Vec-sum proj. axis ratio)');
                    else
                        xlabel('log(d prime ratio)');
                    end
                    xlimit = xlim;
                    ylimit = ylim;
                    ylimit = max(abs(ylimit));
                    ylim([-ylimit ylimit]);
                    
                    % pearson correlation
                    [~,p_pear] = plot_corr_line(tuning_index,choice_change,'MethodOfCorr','Pearson');
                    [~,p_spear] = plot_corr_line(tuning_index,choice_change,'MethodOfCorr','Spearman');
                    text(-0.9,ylimit,sprintf('p. corr. p = %0.3g',p_pear));
                    text(-0.9,ylimit-ylimit/6,sprintf('s. corr. p = %0.3g',p_spear));

                    text(-0.9,ylimit-ylimit/3,sprintf('CR = [%g-%g]%%, N = %g',mt_performance_cri(1),mt_performance_cri(2), sum(perform_select)));
                                        
                    if i==1 && fc == 1
                        text((xlimit(1)-(xlimit(2)-xlimit(1))/10),0,'Fine task','rotation',90,'HorizontalAlignment','center');
                    elseif i == 1 && fc == 2
                        text((xlimit(1)-(xlimit(2)-xlimit(1))/10),0,'Coarse task','rotation',90,'HorizontalAlignment','center');
                    elseif i == 1 && fc == 3
                        text((xlimit(1)-(xlimit(2)-xlimit(1))/10),0,'Fine+Coarse task','rotation',90,'HorizontalAlignment','center');
                    end
                end
            end
            if j == 1 % difficult trial, 0 degree, around 0 degree
                suptitle('Difficult trial');
            elseif j == 2
                suptitle('0 degree');
            else
                suptitle('Around 0 degree');
            end
            SetFigure(13);
            figN = figN + 1;
        end
    end

    function f4p10(debug) % 'Motion type PSE shift'
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        effect_type = 1; % 1: PSE, 2: threshold
        
        is_separate_stimAmp = 0; % 默认为0，combine 20 and 40 stim amplitude
        select_stars_type = [0 1]; % 0: dot stimulus; 1: triangle stimulus
        methods_of_select = {select_typical, 'Typical cells'};
        
        % ---------------------- get data ----------------------------
        % must select 4-target task
        if sum(stimtask_type{1}==1)+sum(stimtask_type{2}==4)==0 % fine 4T和coarse 4T都没选
            disp('********** Please choose 4-Target task! **********');
            return
        end
        stimtask_type_this = stimtask_type;
        stimtask_type_this{1}(stimtask_type_this{1}~=1) = []; % only left 4-Target task
        stimtask_type_this{2}(stimtask_type_this{2}~=4) = [];
        
        % get the MT PSE shift, direction according to spiral index("MT preference"): positive mean prefer motion type
        % is_original_bias = 1; % 0: bias positive = prefer direction; 1: only original bias in MT;  2: original bias(after-before); 3： normalize PSE； 4：normlize PSE shift and original bias: after-before
        if length(cell2mat(stimtask_type_this))==1 % only select fine or coarse, 不使用normalize PSE shift
            is_original_bias = 0; % adjusted PSE: bias positive = prefer direction
            [diff_thistask_temp1, ~, ~] = get_stim_data_uA(stimtask_type_this,methods_of_select{1,1},is_separate_stimAmp,is_original_bias,select_stars_type);
        else
            is_original_bias = 3; % normalize adjusted PSE
            [diff_thistask_temp1, ~, ~] = get_stim_data_uA(stimtask_type_this,methods_of_select{1,1},is_separate_stimAmp,is_original_bias,select_stars_type);
        end
        mt_bias_adj = diff_thistask_temp1{effect_type}(:,3); % >0: bias to prefer; <0 bias to non-prefer

        
        % get the MT PSE shift, direction according to bias to T or R: positive mean prefer Translation
        if length(cell2mat(stimtask_type_this))==1 % only select fine or coarse, 不使用normalize PSE shift
            is_original_bias = 1; % original PSE
            [diff_thistask_temp2, S_unit, NS_unit] = get_stim_data_uA(stimtask_type_this,methods_of_select{1,1},is_separate_stimAmp,is_original_bias,select_stars_type);
        else
            is_original_bias = 4; % normalize original PSE
            [diff_thistask_temp2, S_unit, NS_unit] = get_stim_data_uA(stimtask_type_this,methods_of_select{1,1},is_separate_stimAmp,is_original_bias,select_stars_type);
        end
        mt_bias_ori = diff_thistask_temp2{effect_type}(:,3); % >0: bias to T; <0 bias to R
        
        S_unit = S_unit{effect_type}(:,3);
        NS_unit = NS_unit{effect_type}(:,3);
        %         S_unit_ori = S_unit{effect_type}(:,3); % 和上面的是一样的
        %         NS_unit_ori = NS_unit{effect_type}(:,3); % 和上面的是一样的
        
        % group mt bias
        mt_bias{1} = mt_bias_adj;
        mt_bias{2} = mt_bias_ori;
        til{1} = 'Adjusted motion type PSE shift';
        til{2} = 'Original motion type PSE shift';
        
        % ------------------------ Plot ------------------------------
        xlimit_all = ceil(max(abs([mt_bias_adj;mt_bias_ori]))*10)/10;
        xbins = linspace(-xlimit_all,xlimit_all,9);
        
        set(figure(figN),'Position',[50,50 1400,700], 'Name', 'Motion type PSE shift');
        for i = 1:2
            subplot(1,2,i)
            histogram(mt_bias{i},xbins,'facecolor','w');
            hold on
            histogram(mt_bias{i}(S_unit),xbins,'facecolor',c{3});
            
            % mean+sem
            mean_mt = nanmean(mt_bias{i});
            sem_mt = nanstd(mt_bias{i}) / sqrt(sum(S_unit|NS_unit));
            % ttest
            [~,p_mt] = ttest(mt_bias{i},0);
            
            axis tight
            xlim([-xlimit_all xlimit_all]);
            ylimit = ylim;
            plot(mean_mt,ylimit(2)*1.1,'kv');
            text(-xlimit_all/1.1,ylimit(2),sprintf('%.3g ± %.3g (p = %.3g)',mean_mt,sem_mt,p_mt));
            ylim([0 ylimit(2)*1.2]);

            title(til{i});
            if i == 1
                xla = strcat('Bias to non-prefer',repmat(32,1,20),'Bias to prefer');
                xlabel(xla);
            else
                xla = strcat('Bias to Rot.',repmat(32,1,20),'Bias to Trans.');
                xlabel(xla);
            end
        end
        
        SetFigure(13)
        [til] = SetTitle(stimtask_type,monkey_included_for_analysis);
        str = ['Motion type PSE shift in ',til];
        suptitle(str);
        figN = figN + 1;
    end

    
%%
    function f5p1(debug) % 'Preferred direction'
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        stimtask_type_this = 1; % 4 targets fine task
        
        % stimCueS stimHonly stimRonly stimCo stimCoH stimCoR
        %    1         2         3       4        5      6
        select_task = false(size(loadM_data,2),1);
        stimtask_type_this = cell2mat(stimtask_type);
        for n = 1:length(loadM_data)
            for i = 1:length(stimtask_type_this)
                if ischar(group_result(n).stim_FILE{stimtask_type_this(i)}) % have done one of the selected task
                    select_task(n) = true;
                end
            end
        end
        select_task(~methods_of_select{1,1}) = false;
        
        cell_type = [group_result.cell_type]';
        
        % VectorSum result
        pref_direc = cell2mat({group_result.pref_direc}');
        
        spiral_index = cell2mat({group_result.spiral_index}');
        
        % ------------------------ Plot ------------------------------
        set(figure(figN),'Position',[50,50 1200,700], 'Name', 'Prefer Direction Distribution');
        xbins{1} = linspace(-180,180,11);
        xbins{2} = linspace(-180,180,11);
        xbins{3} = linspace(-1,1,11);
        xtl{1} = {[-180:90:180]};
        xtl{2} = {'Con' 'CW' 'Exp' 'CCW' 'Con'};
        xtl{3} = {[-1:0.5:1]};
        yt(1) = {'Translation preference'};
        yt(2) = {'Rotation preference'};
        yt(3) = {'Spiral index'};
        
        % 预处理决定distribution ylim 范围
        for d = 1:2 % heading pref and rotation pref
            nbins(d,:) = hist(pref_direc(select_task,d),xbins{d});
        end
        ylimit = ceil(max(nbins(:))/10)*10;
        
        % 预处理决定polarplot rlim 范围
        h_temp = figure(9999);
        for d = 1:2 % heading pref and rotation pref
            h = polarhistogram(HeadingToazi(pref_direc(select_task,d))*pi/180,16);
            rbins(d,:) = h.Values;
        end
        rlimit = max(rbins(:));
        close(h_temp);
        
        for d = 1:3 % heading pref, rotation pref, spiral index
            subplot(2,3,d)
            if d ~=3
                % bar plot
                nbins = hist(pref_direc(select_task,d),xbins{d});
                hold on
                h = bar(xbins{d},nbins,'facecolor',c{d});
                xlim([-205 205]);
                ylim([0 ylimit+ylimit/3]);
                set(gca,'XTick',[-180:90:180]);
                set(gca,'xticklabel',xtl{d});
                ylabel('Cases');
                xlabel(yt(d));
                
                % 检验tuning的分布情况
                az = HeadingToazi(pref_direc(methods_of_select{1,1},d));
                az = az';
                % 1. circular data uniformity test: Hodges-Ajne test
                p_HA(d) = circ_otest(deg2rad(az));  % h=1,p<0.05 不是均匀分布
                
                % 2. modality test   p>0.05为符合该mode
                [mode(d), p_mod{d}] = modalityTestForYong(az,[1 2],1000);
                
                text(-170,ylimit+ylimit/3.2,sprintf('Hodges-Ajne test,p=%0.3g',p_HA(d)),'fontsize',10); % circular data uniformity test
                text(-170,ylimit+ylimit/5,sprintf('modality test, p_u_n_i=%0.2g',p_mod{d}(1)),'fontsize',10);
                text(-170,ylimit+ylimit/10,sprintf('modality test, p_b_i=%0.2g',p_mod{d}(2)),'fontsize',10);
                
                % rose plot
                subplot(2,3,3+d)
                h = polarhistogram(HeadingToazi(pref_direc(select_task,d))*pi/180,16,'facecolor',c{d},'facealpha',0.5);
                if d == 1
                    set(gca,'ThetaTick',[0:45:315],'ThetaTickLabel',{'90','45','0','-45','-90','-135','180','135'});
                else
                    set(gca,'ThetaTick',[0:45:315],'ThetaTickLabel',{'CCW','45','0','-45','CW','-135','180','135'});
                end
                rlim([0 rlimit]);
            else
                % spiral index
                histogram(spiral_index(select_task),xbins{d},'facecolor',c{d});
                set(gca,'xticklabel',xtl{d});
                xlabel(yt(d));
                
                mean_spiral_index = mean(spiral_index(select_task));
                [~,p] = ttest(spiral_index(select_task),0);
                ylimit = ylim;
                text(-0.9,ylimit(2)-ylimit(2)/18,sprintf('mean = %2.2g',mean_spiral_index));
                text(-0.9,ylimit(2)-ylimit(2)/10,sprintf('ttest,p = %2.2g',p));
                text(-0.9,ylimit(2)-ylimit(2)/7,sprintf('N = %g',sum(select_task)));
                xlim([-1 1]);
                ylabel('Cases');
            end
        end
        suptitle(SetTitle(stimtask_type,monkey_included_for_analysis));
        SetFigure(13);
        figN = figN + 1;
        
        
        set(figure(figN),'Position',[50,50 700,700], 'Name', 'Prefer Direction with RF location');
        for i = 1:3
            hold on
            plot(pref_direc(select_task & cell_type == i,1),pref_direc(select_task & cell_type == i,2),'ko','markersize',15,'markerfacecolor',c{i}); % all unit
            cell_num(i) = sum(select_task & cell_type == i);
        end
        
        cell_num_all = sum(cell_num);
        %         plot(pref_direc(select_task,1),pref_direc(select_task,2),'ko','markersize',11); % all unit
        hold on
        xlim([-180 180]);ylim([-180 180]);
        set(gca,'xtick',-180:90:180);
        set(gca,'ytick',-180:90:180);
        set(gca,'yticklabel',{'Con' 'CW' 'Exp' 'CCW' 'Con'});
        
        text(-150,150,sprintf('N = %g',cell_num_all));
        
        xlabel('Translation preference');
        ylabel('Rotation preference');
        
        % 4 * 4等分
        dir_range = [-180:90:180];
        for i = 2:4
            plot([dir_range(i) dir_range(i)],[-180 180],'k:');
            plot([-180 180],[dir_range(i) dir_range(i)],'k:');
        end
        SetFigure(13);
        figN = figN + 1;
    end

    function f5p2(debug) % RF in stim
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        
        % 做了stim task的unit的 RF
        select_stim_task = false(length(loadM_data),1);
        for n = 1:length(loadM_data)
            for i = 1:6 % 6 stim tasks
                if ischar(group_result(n).stim_FILE{i}) % have done the fine 4 target task
                    select_stim_task(n) = true;
                end
            end
        end
        
        rf_temp =  {group_result.RF}';
        rf_temp_stim = rf_temp(select_stim_task & methods_of_select{1,1});
        have_rf = ~cellfun(@isempty,rf_temp_stim);
        rf_temp_stim(~have_rf) = {[nan nan nan nan]};
        rf = cell2mat(rf_temp_stim);
        rf_size = rf(:,3) .* rf(:,4);
        rf_dia = 2.* sqrt(rf_size / pi); % 等效圆直径
        rf_dia2 = sqrt(rf_size); % 等效正方形边长
        
        % ------------------------ Plot ------------------------------
        set(figure(figN),'Position',[500,250, 1000,500], 'Name', 'MST RF (stim task)');
        subplot(1,2,1)
        axis equal; box on; axis([-60 60 -60 60]);
        line([-60 60],[0 0],'color','k','LineStyle',':'); hold on;
        line([0 0],[-60 60],'color','k','LineStyle',':');
        set(gca,{'xtick','ytick'},{-60:20:60,-60:20:60});
        xlabel('degree');
        ylabel('degree');
        title('MST RF (stim task)');
        
        for i = 1:length(rf)
            if ~isnan(rf(i,1))
                rectangle('position',[rf(i,1)-rf(i,3)/2 rf(i,2)-rf(i,4)/2, rf(i,3) rf(i,4)],...
                    'Curvature',[0.3 0.3],'EdgeColor',[0.6 0.6 0.6],'LineWidth',1);
                hold on
            end
        end
        
        text(-55,55,sprintf('N=%g',sum(have_rf)));
        
        % RF distribution
        subplot(1,2,2)
        xbins = [0:10:110];
        h = hist(rf_dia,xbins);
        hold on;
        h_h = bar(xbins,h,1,'facecolor','w');
        xlabel('Equivalent diameter of receptive field (degree)');
        ylabel('Number of cases');
        set(gca,'xlim',[0 120]);
        % mean
        mean_d = nanmean(rf_dia);
        top = max(h);
        plot([mean_d mean_d],[0 top+1],'k:');
        plot(mean_d,top+2,'kv','markerfacecolor','k');
        text(mean_d+4,top+2,num2str(mean_d),'fontsize',15);
        
        SetFigure(13);
        figN = figN+1;
    end

    function f5p3(debug) %'RF type and stimulation effect'
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        stimtask_type_this_fine = [1 2 3];
        stimtask_type_this_coarse = [4 5 6];
        effect_type = 1; % 1: PSE, 2: threshold
        
        % ---------------------- get data ----------------------------
        %         is_separate_stimAmp = 0;
        %         is_original_bias = 0; %  正值：bias to prefer direction
        %
        %         [diff_temp_fine, S_unit_temp_fine, NS_unit_temp_fine, ~] = get_stim_data_uA(stimtask_type_this_fine,methods_of_select,is_separate_stimAmp,is_original_bias);
        %         [diff_temp_coarse, S_unit_temp_coarse, NS_unit_temp_coarse, ~] = get_stim_data_uA(stimtask_type_this_coarse,methods_of_select,is_separate_stimAmp,is_original_bias);
        %         PSE_shift{1} = diff_temp_fine{effect_type};
        %         PSE_shift{2} = diff_temp_coarse{effect_type};
        %         S_unit{1} = S_unit_temp_fine{effect_type};
        %         S_unit{2} = S_unit_temp_coarse{effect_type};
        %         NS_unit{1} = NS_unit_temp_fine{effect_type};
        %         NS_unit{2} = NS_unit_temp_coarse{effect_type};
        
        TVectorSumALL = []; RVectorSumALL = [];
        TVectorSumMidd = []; RVectorSumMidd = [];
        TVectorSumQuad = []; RVectorSumQuad = [];
        
        TRAngleALL = nan(length(loadM_data),1);
        TRAngleMidd = nan(length(loadM_data),1);
        TRAngleQuad = nan(length(loadM_data),1);
        
        cd('D:\user\desktop\test202010');
        namelist = dir('D:\user\desktop\test202010');
        len = length(namelist);
        for ii = 1:len
            %             try
            file_name11{ii}=namelist(ii).name;
            load(file_name11{ii});
            %             catch
            %                 disp('1111111111');
            %             end
        end
        TRAngleALL(TRAngleALL== 0.00000) = nan;
        
        %         TRAngleALL_thistask = [];
        %         for i = 1:length(stimtask_type_this_fine) % fine or coarse都一样，重复3就行了
        %             TRAngleALL_thistask = [TRAngleALL_thistask;TRAngleALL(methods_of_select{1,1},:)];
        %         end
        
        % 将RF内optic flow vector分为2大类：
        RFtype = nan(size(TRAngleALL));
        %         RFtype(TRAngleALL<=45) = 1; % similar
        %         RFtype(TRAngleALL>45 & TRAngleALL<=135) = 2; % vertical
        %         RFtype(TRAngleALL>135) = 3; % opposite
        RFtype(TRAngleALL<=45) = 1; % similar
        RFtype(TRAngleALL>45) = 2; % not similar
        
        
        % ---------------------- get data ----------------------------
        % ===============  Condition base  ====================
        % stimCueS stimHonly stimRonly stimCo stimCoH stimCoR
        %    1         2         3       4        5      6
        % =====================================================
        % 1. 2-target task
        f2T_temp(:,1) = diffBias_task{2}(methods_of_select{1,1},1); % fine two-target heading task
        f2T_temp(:,2) = diffBias_task{3}(methods_of_select{1,1},2); % fine two-target rotation task
        c2T_temp(:,1) = diffBias_task{5}(methods_of_select{1,1},1); % coarse two-target heading task
        c2T_temp(:,2) = diffBias_task{6}(methods_of_select{1,1},2); % coarse two-target rotation task
        
        % 2. 4-target task
        f4T_temp = diffBias_task{1}(methods_of_select{1,1},1:2); % fine two-target heading task
        c4T_temp = diffBias_task{4}(methods_of_select{1,1},1:2); % coarse two-target heading task
        
        % find which do heading and rotation simultaneously
        select_f2T = sum(isnan(f2T_temp),2)==0;
        select_c2T = sum(isnan(c2T_temp),2)==0;
        select_f4T = sum(isnan(f4T_temp),2)==0;
        select_c4T = sum(isnan(c4T_temp),2)==0;
        
        hr_comp{1} = f2T_temp(select_f2T,:);
        hr_comp{2} = c2T_temp(select_c2T,:);
        hr_comp{3} = f4T_temp(select_f4T,:);
        hr_comp{4} = c4T_temp(select_c4T,:);
        
        set(figure(figN),'Position',[100,100 1000,500], 'Name', '');
        for rft = 1:2
            fine2target{rft} = f2T_temp(select_f2T & RFtype(methods_of_select{1,1})==rft,:);
            fine4target{rft} = f4T_temp(select_f4T & RFtype(methods_of_select{1,1})==rft,:);
            coarse2target{rft} = c2T_temp(select_c2T & RFtype(methods_of_select{1,1})==rft,:);
            coarse4target{rft} = c4T_temp(select_c4T & RFtype(methods_of_select{1,1})==rft,:);
            
            meanf2t{rft} = mean(fine2target{rft});
            meanf4t{rft} = mean(fine4target{rft});
            meanc2t{rft} = mean(coarse2target{rft});
            meanc4t{rft} = mean(coarse4target{rft});
            
            errf2t{rft} = std(fine2target{rft}) / sqrt(size(fine2target{rft},1));
            errf4t{rft} = std(fine4target{rft}) / sqrt(size(fine4target{rft},1));
            errc2t{rft} = std(coarse2target{rft}) / sqrt(size(coarse2target{rft},1));
            errc4t{rft} = std(coarse4target{rft}) / sqrt(size(coarse4target{rft},1));
            
            [~,pvalue(1,rft)] = ttest(fine2target{rft}(:,1),fine2target{rft}(:,2));
            [~,pvalue(2,rft)] = ttest(fine4target{rft}(:,1),fine4target{rft}(:,2));
            [~,pvalue(3,rft)] = ttest(coarse2target{rft}(:,1),coarse2target{rft}(:,2));
            [~,pvalue(4,rft)] = ttest(coarse4target{rft}(:,1),coarse4target{rft}(:,2));
            
            subplot(4,2,rft)
            errorbar(meanf2t{rft},errf2t{rft});
            xlim([0.5 2.5]);
            ylabel('Fine 2-T');
            if rft == 1
                title('similar RF vectorsum')
            else
                title('different RF vectorsum')
            end
            
            subplot(4,2,rft+2)
            errorbar(meanf4t{rft},errf4t{rft});
            xlim([0.5 2.5]);
            ylabel('Fine 4-T');
            
            subplot(4,2,rft+4)
            errorbar(meanc2t{rft},errc2t{rft});
            xlim([0.5 2.5]);
            ylabel('Coarse 2-T');
            
            subplot(4,2,rft+6)
            errorbar(meanc4t{rft},errc4t{rft});
            xlim([0.5 2.5]);
            ylabel('Coarse 4-T');
        end
        
        
        
        
        %              subplot(2,2,rft)
        %             xlimit = xlim; ylimit = ylim;
        %             Limit = max(max(abs(ylimit)),max(abs(xlimit)));
        %             xlim([-Limit Limit])
        %             ylim([-Limit Limit]);
        %             xlimit = xlim; ylimit = ylim;
        %
        %             hold on
        %             plot([0 0],ylimit,'k:');
        %             plot(xlimit,[0 0],'k:');
        %         end
        %
        %         % ------------------------ Plot ------------------------------
        %         set(figure(figN),'Position',[100,100 700,900], 'Name', '');
        %         for rft = 1:3
        %            for fc = 1:2 % fine or coarse
        %                subplot(2,3,rft+(fc-1)*3)
        %
        %            end
        %         end
    end

    function f5p4(debug) %'Some tuning index in four-target task'
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        stimtask_type_this = [1 4];% 4 targets fine task and 4 targets coarse task
        % ---------------------- get data ----------------------------
        select_4T_task = false(length(loadM_data),1);
        for n = 1:length(loadM_data)
            for fc = 1:2
                if ischar(group_result(n).stim_FILE{stimtask_type_this(fc)}) % have done the fine 4 target task
                    select_4T_task(n) = true;
                end
            end
        end
        
        select_4T_task = logical(select_4T_task & methods_of_select{1,1});
        d_prime_index = cell2mat({group_result.d_prime_index_coarse}');
        spiral_index = cell2mat({group_result.spiral_index}');
        
        % ------------------------ Plot ------------------------------
        set(figure(figN),'Position',[100,100 800,400], 'Name', 'Spiral Index and D prime index in 4-T task');
        edges = linspace(-1,1,11);
        subplot(1,2,1)
        histogram(spiral_index(select_4T_task),edges);
        mean_spiral_index = mean(spiral_index(select_4T_task));
        [~,p] = ttest(spiral_index(select_4T_task),0);
        ylimit = ylim;
        text(-0.9,ylimit(2)-ylimit(2)/18,sprintf('mean = %2.2g',mean_spiral_index));
        text(-0.9,ylimit(2)-ylimit(2)/10,sprintf('p = %2.2g',p));
        text(-0.9,ylimit(2)-ylimit(2)/7,sprintf('N = %2.2g',sum(select_4T_task)));
        xlim([-1 1]);
        title('Spiral Index in 4-T task');
        ylabel('Cases');
        xlabel('Spiral Index');
        SetFigure(13);
        
        subplot(1,2,2)
        histogram(d_prime_index(select_4T_task),edges);
        mean_d_index = mean(d_prime_index(select_4T_task));
        [~,p] = ttest(d_prime_index(select_4T_task),0);
        ylimit = ylim;
        text(-0.9,ylimit(2)-ylimit(2)/18,sprintf('mean = %2.2g',mean_d_index));
        text(-0.9,ylimit(2)-ylimit(2)/10,sprintf('p = %2.2g',p));
        text(-0.9,ylimit(2)-ylimit(2)/7,sprintf('N = %2.2g',sum(select_4T_task)));
        xlim([-1 1]);
        title('D prime Index in 4-T task');
        ylabel('Cases');
        xlabel('D prime Index');
        figN = figN+1;
        
        set(figure(figN),'Position',[1000,100 500,500], 'Name', 'Comparison Spiral Index and D prime index in 4-T task');
        plot(spiral_index(select_4T_task),d_prime_index(select_4T_task),'ko');
        xlabel('Spiral Index');
        ylabel('D prime Index');
        xlim([-1 1]);ylim([-1 1]);
        hold on
        plot([-1 1],[-1 1],'k:');
        SetFigure(13)
        figN = figN+1;
    end

    function f5p5(debug) % 'Cell type in stim task'
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        
        % ---------------------- get data ----------------------------
        cell_type = [group_result.cell_type];
        % 做了stim task的unit的 cell type
        select_stim_task = false(length(loadM_data),1);
        stimtask_type_this = cell2mat(stimtask_type);
        for n = 1:length(loadM_data)
            for i = 1:length(stimtask_type_this)
                if ischar(group_result(n).stim_FILE{stimtask_type_this(i)}) % have done the fine 4 target task
                    select_stim_task(n) = true;
                end
            end
        end
        
        select = logical(methods_of_select{1,1} & select_stim_task);
        cell_type = cell_type(select)';
        
        heading_cell_num = sum(cell_type == 1);
        rotation_cell_num = sum(cell_type == 2);
        spiral_cell_num = sum(cell_type == 3);
        no_tuning_cell_num = sum(cell_type == 4) + sum(cell_type == 5); % combine cell type 4 and 5: no tuning cell
        
        cell_all_num = heading_cell_num+rotation_cell_num+spiral_cell_num+no_tuning_cell_num;
        mt_tuning_num = sum([group_result(select_all).mt_type] == 1);
        all_active_num = heading_cell_num + rotation_cell_num + spiral_cell_num + no_tuning_cell_num;
        
        % ------------------------ Plot ------------------------------
        set(figure(figN),'name','Cell type in stim task','pos',[200,300 800,500]); clf
        pp = pie([heading_cell_num rotation_cell_num spiral_cell_num no_tuning_cell_num],[0 0 0 0]);
        legend('Translation cell','Roll cell','Spiral cell','Untuned cell','location','bestoutside');
        
        axes('position',[0.67 0.05 0.3 0.5]);
        text(0,1,sprintf('Translation Cell = %g (%2.4g%%)',heading_cell_num, heading_cell_num/cell_all_num*100));
        text(0,0.8,sprintf('Roll Cell = %g (%2.4g%%)',rotation_cell_num, rotation_cell_num/cell_all_num*100));
        text(0,0.6,sprintf('Spiral Cell = %g (%2.4g%%)',spiral_cell_num, spiral_cell_num/cell_all_num*100));
        text(0,0.4,sprintf('Untuned = %g (%2.4g%%)',no_tuning_cell_num, no_tuning_cell_num/cell_all_num*100));
        text(0,0.2,sprintf('N = %g',cell_all_num));
        axis off;
        
        title(SetTitle(stimtask_type,monkey_included_for_analysis));
        SetFigure(15);figN = figN+1;
        
    end

    function f6p1(debug) % 'Psychometric function' in Two-target fine/coarse task (non-stim trials)
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        % t = 2 stimHonly
        % t = 3 stimRonly
        % t = 5 stimCoH
        % t = 6 stimCoR
        stimtask_type_this = [2 3 5 6];
        
        % ---------------------- get data ----------------------------
        % 提取condition list, may be more than one condition list in one specific task
        % 1. 2 target fine task, h and r
        condition_list_temp{1} = {unique_condition_everysite_h{:,2}}'; % 2T,f-h
        condition_list_temp{2} = {unique_condition_everysite_r{:,3}}'; % 2T,f-r
        
        % 2. 2 target coarse task
        condition_list_temp{3} = {unique_condition_everysite_co{:,5}}'; % 2T,co-h
        condition_list_temp{4} = {unique_condition_everysite_co{:,6}}'; % 2T,co-r
        
        stimtask_type_this = [2 3 5 6]; % stimCueS, stimCueS, stimCo, stimCo, stimCueS, stimCo
        ind_this = [1 2 1 2]; % h r h r mt mt
        
        condition_raw = arrayfun(@(x) nan(length(loadM_data),13), [1:6],'UniformOutput',0); % 足够大的矩阵放下多有condition,每个site最大13个condition
        psy_correct_raw = arrayfun(@(x) nan(length(loadM_data),13), [1:6],'UniformOutput',0); % 足够大的矩阵放下多有condition,每个site最大13个condition
        
        at_least_rep = 15; % condition at least have been done for 10 site will be considered
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        use_finally_cond = 1; % use_finally_cond=0: use all condition to plot and fitting; use_finally_cond=1: use finally condition to plot and fitting
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        plot_this = false(1,4);
        for t = 1:4 % f-h, f-r, co-h, co-r
            cond_num = 0;
            for n = 1:length(loadM_data)
                if ~isnan(condition_list_temp{t}{n}(1)) && sum(condition_list_temp{t}{n} == 0)>0 && methods_of_select{1}(n)==1 % remove nan and remove not include zero degree session
                    cond_num = cond_num + 1; % next condition list
                    condition_list_str{cond_num,1} = num2str(condition_list_temp{t}{n}); % 转为string,求unique
                    
                    
                    % 原始的condition及psy_correct数据，按照各自的condition排列
                    condition_raw{t}(n,1:length(condition_list_temp{t}{n})) = condition_list_temp{t}{n}; % 保留数字类型，填充到大矩阵里
                    psy_correct_raw{t}(n,1:length(condition_list_temp{t}{n})) = group_result(n).Stim_Psy_correct{stimtask_type_this(t),ind_this(t)}{1}; % 每一行中第一小行是电刺激前的，第二小行是电刺激后的 Stim_Psy_correct{}{1}: ctrl; Stim_Psy_correct{}{2} stim
                end
            end
            
            % 提取uniqe_condition: 不同的uniqe_condition合并成一个完整的uniqe_condition（每个condition至少重复了10个site）
            % 例如两个condition list: -6 -3 0 3 6;  -4 -2 0 2 4
            % 合并成： -6 -4 -3 -2 0 2 3 4 6, 原来的condition list重排为：
            % -6 [] -3 [] 0 [] 3 [] 6;
            % [] -4 [] -2 0 2 [] 4 []
            if cond_num~=0
                plot_this(t) = true;
                tbl = tabulate(condition_raw{t}(:)); % fine Frequency
                unique_condition_new{t} = tbl(tbl(:,2)>=at_least_rep,1); % 1st column: the unique values of x ; % 2st: the number of instances of each value
                
                % 按照unique_condition排列，重建原始的condition_raw和psy_correct_raw （这两个矩阵是一一对应的）
                psy_correct_new{t} = nan(length(loadM_data),length(unique_condition_new{t}));
                
                for n = 1:length(loadM_data)
                    if ~isnan(psy_correct_raw{t}(n,1))
                        for i = 1:length(unique_condition_new{t})
                            select = condition_raw{t}(n,:) == unique_condition_new{t}(i);
                            if sum(select)~=0
                                psy_correct_new{t}(n,i) = psy_correct_raw{t}(n,select);
                            end
                        end
                    end
                end
                
                % save site selection
                psy_correct_plot{t} = find(any(~isnan(psy_correct_new{t}),2));
                psy_correct_new{t}(~any(~isnan(psy_correct_new{t}),2),:) = []; % remove all nan site
                
                % find the finally used condition
                first_cond = condition_raw{t}(~isnan(condition_raw{t}(:,1)),1);
                unique_first_cond = unique(first_cond)';
                % find the location of this first_cond
                this_cond_find = []; this_cond_loc = []; this_cond_num = [];
                for i = 1:length(unique_first_cond)
                    this_cond_find{i} = find(first_cond == unique_first_cond(i));
                    this_cond_loc(i) = mean(this_cond_find{i}); % 该condition平均出现的位置
                    this_cond_num(i) = length(this_cond_find{i});
                end
                select = this_cond_num>=10; % only select this_cond_num>=10
                % find the finally used condition
                finally_used_first_cond = unique_first_cond(this_cond_loc==max(this_cond_loc(select))); % location越大，也就是越靠后的文件
                finally_used_cond_range{t} = [finally_used_first_cond -finally_used_first_cond]; % 总是对称的
            else
                % do not do this task, 0 cond_num
                plot_this(t) = false;
            end
        end
        
        % ------------------------ Plot ------------------------------
        set(figure(figN),'Position',[100,50 800,800], 'Name', 'Psychometric function in Two-target control task');clf
        for t = 1:4
            subplot(2,2,t)
            if plot_this(t)
                xh = min(unique_condition_new{t}):0.01:max(unique_condition_new{t});
                hold on
%                 for i = 1:length(psy_correct_plot{t})
%                     plot(xh, cum_gaussfit(psy_perf_ctrl{stimtask_type_this(t),ind_this(t)}(psy_correct_plot{t}(i),:),xh),'-','color',c{99},'linewidth',0.7); % control
%                     %                 plot(xh, cum_gaussfit(psy_perf_stim{stimtask_type_this(t),ind_this(t)}(psy_correct_plot{t}(i),:),xh),'-','color',c{99},'linewidth',0.7); % stim
%                 end
                
                rep_num = sum(~isnan(psy_correct_new{t}),1); % 每个condition重复site个数
                % fitting
                if use_finally_cond
                    % fitting of finally used condtion
                    select_cond = logical(unique_condition_new{t}>=finally_used_cond_range{t}(1) & unique_condition_new{t}<=finally_used_cond_range{t}(2));
                    fit_data_psy_cum = [];
                    fit_data_psy_cum(:,1) = unique_condition_new{t}(select_cond);
                    fit_data_psy_cum(:,2) = nanmean(psy_correct_new{t}(:,select_cond))';
                    fit_data_psy_cum(:,3) = rep_num(select_cond)'; % simply this condtion repetition count
                else
                    fit_data_psy_cum = [];
                    fit_data_psy_cum(:,1) = unique_condition_new{t};
                    fit_data_psy_cum(:,2) = nanmean(psy_correct_new{t})';
                    fit_data_psy_cum(:,3) = rep_num'; % simply this condtion repetition count
                end
                
                
                method = 0; % 0=Maximum likelihood (default), 1=Square error
                tolerance = 0; % 0=No tolerance (default), other = error tolerance of the max abs (heading)
                [bb,tt] = cum_gaussfit_max1(fit_data_psy_cum,method,tolerance);
                psy_perf_all{t} = [bb, tt];
                plot(xh, cum_gaussfit(psy_perf_all{t},xh), '-','color',c{ind_this(t)},'linewidth',3);
                % mean+sem
                errorbar(unique_condition_new{t},nanmean(psy_correct_new{t}),nanstd(psy_correct_new{t})./sqrt(rep_num),'k.','markersize',10);
                
                % finally_used_cond_range
                plot([finally_used_cond_range{t}(1) finally_used_cond_range{t}(1)],[0 1],'k:');
                plot([finally_used_cond_range{t}(2) finally_used_cond_range{t}(2)],[0 1],'k:');
                
                if use_finally_cond
                    xlim([finally_used_cond_range{t}(1)*1.1 finally_used_cond_range{t}(2)*1.1]); % range of finally used condtion
                else
                    xlim([min(unique_condition_new{t})*1.1 max(unique_condition_new{t})*1.1]); % range of all condition
                end
                
                ylim([0,1]);
                if min(unique_condition_new{t})~=finally_used_cond_range{t}(1)
                    set(gca,'xtick',[min(unique_condition_new{t}) finally_used_cond_range{t}(1) finally_used_cond_range{t}(1)/2 0 finally_used_cond_range{t}(2)/2 finally_used_cond_range{t}(2) max(unique_condition_new{t})]);
                else
                    set(gca,'xtick',[finally_used_cond_range{t}(1) finally_used_cond_range{t}(1)/2 0 finally_used_cond_range{t}(2)/2 finally_used_cond_range{t}(2)]);
                end
                set(gca,'xTickLabel',num2str(get(gca,'xTick')','%.1f'));
                %                 set(gca,'xtick',[finally_used_cond_range{t}(1) 0 finally_used_cond_range{t}(2)]);
                set(gca,'ytick',[0:0.25:1]);
                %                 set(gca,'ytick',[0:0.5:1]);
                
                xlimit = xlim;
                text(min(xlimit/1.1),0.9,sprintf('\\it\\mu_{psy} = \\rm%0.3g\\circ',psy_perf_all{t}(1)),'fontsize',10);
                text(min(xlimit/1.1),0.75,sprintf('\\it\\sigma_{psy} = \\rm%0.3g\\circ',psy_perf_all{t}(2)),'fontsize',10);
                text(min(xlimit/1.1),0.5,sprintf('n = %g',length(psy_correct_new{t})),'fontsize',10);
                
                if t==1
                    title('Fine trans.');
                elseif t==2
                    title('Fine rot.');
                elseif t==3
                    title('Coarse trans.');
                elseif t==4
                    title('Coarse rot.');
                end
            else
                text(0.1,0.1,'Did not do this task');
            end
        end
        SetFigure(13)
        suptitle('Psychometric function in Two-target control task');
        figN = figN+1;
    end

    function f6p2(debug) % 'Psychometric function' in Four-target fine/coarse task (non-stim trials)
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        % ===============  Condition base  ====================
        % stimCueS stimHonly stimRonly stimCo stimCoH stimCoR
        %    1         2         3       4        5      6
        % ---------------------- get data ----------------------------
        % 提取condition list, may be more than one condition list in one specific task
        % 1. 4 target fine task, h r mt
        condition_list_temp{1} = {unique_condition_everysite_h{:,1}}'; % 4T,f-h
        condition_list_temp{2} = {unique_condition_everysite_r{:,1}}'; % 4T,f-r
        condition_list_temp{3} = {unique_condition_everysite_mt{:,1}}'; % 4T,f-mt
        
        % 2. 4 target coarse task
        condition_list_temp{4} = {unique_condition_everysite_co{:,4}}'; % 4T,co-h
        condition_list_temp{5} = {unique_condition_everysite_co{:,4}}'; % 4T,co-r
        condition_list_temp{6} = {unique_condition_everysite_mt{:,4}}'; % 4T,co-mt

        stimtask_type_this = [1 1 1 4 4 4]; % stimCueS, stimCueS, stimCueS, stimCo, stimCo, stimCo
        ind_this = [1 2 3 1 2 3]; % h r h r mt mt
        
        condition_raw = arrayfun(@(x) nan(length(loadM_data),13), [1:6],'UniformOutput',0); % 足够大的矩阵放下多有condition,每个site最大13个condition
        psy_correct_raw = arrayfun(@(x) nan(length(loadM_data),13), [1:6],'UniformOutput',0); % 足够大的矩阵放下多有condition,每个site最大13个condition
        
        at_least_rep = 10; % condition at least have been done for 10 site will be considered
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        use_finally_cond = 1; % use_finally_cond=0: use all condition to plot and fitting; use_finally_cond=1: use finally condition to plot and fitting
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for t = 1:6 % f-h, f-r,f-mt, co-h, co-r, , co-mt
            cond_num = 1;
            for n = 1:length(loadM_data)
                if ~isnan(condition_list_temp{t}{n}(1)) && sum(condition_list_temp{t}{n} == 0)>0 && methods_of_select{1}(n)==1 % remove nan and remove not include zero degree session
                    condition_list_str{cond_num,1} = num2str(condition_list_temp{t}{n}); % 转为string,求unique
                    cond_num = cond_num + 1; % next condition list
                    
                    % 原始的condition及psy_correct数据，按照各自的condition排列
                    condition_raw{t}(n,1:length(condition_list_temp{t}{n})) = condition_list_temp{t}{n}; % 保留数字类型，填充到大矩阵里
                    psy_correct_raw{t}(n,1:length(condition_list_temp{t}{n})) = group_result(n).Stim_Psy_correct{stimtask_type_this(t),ind_this(t)}{1}; % 每一行中第一小行是电刺激前的，第二小行是电刺激后的 Stim_Psy_correct{}{1}: ctrl; Stim_Psy_correct{}{2} stim
                    psy_correct_true_raw{t}(n,1:length(condition_list_temp{t}{n})) = group_result(n).Stim_Psy_correct_true{stimtask_type_this(t),ind_this(t)}{1}; % true error
                end
            end
            
            % 提取uniqe_condition: 不同的uniqe_condition合并成一个完整的uniqe_condition（每个condition至少重复了10个site）
            % 例如两个condition list: -6 -3 0 3 6;  -4 -2 0 2 4
            % 合并成： -6 -4 -3 -2 0 2 3 4 6, 原来的condition list重排为：
            % -6 [] -3 [] 0 [] 3 [] 6;
            % [] -4 [] -2 0 2 [] 4 []
            tbl = tabulate(condition_raw{t}(:)); % find Frequency
            unique_condition_new{t} = tbl(tbl(:,2)>=at_least_rep,1); % 1st column: the unique values of x ; % 2st: the number of instances of each value
            
            % 按照unique_condition排列，重建原始的condition_raw和psy_correct_raw （这两个矩阵是一一对应的）
            psy_correct_new{t} = nan(length(loadM_data),length(unique_condition_new{t}));
            psy_correct_true_new{t} = nan(length(loadM_data),length(unique_condition_new{t}));
            
            for n = 1:length(loadM_data)
                if ~isnan(psy_correct_raw{t}(n,1)) % raw中第一列不为nan
                    for i = 1:length(unique_condition_new{t})
                        select = condition_raw{t}(n,:) == unique_condition_new{t}(i); % 只选取unique_condition_new对应角度的psy_correct （可能只有0度对应，所以只有1个值）
                        if sum(select)==1 % select this degree
                            psy_correct_new{t}(n,i) = psy_correct_raw{t}(n,select);
                            psy_correct_true_new{t}(n,i) = psy_correct_true_raw{t}(n,select); % true error
                        end
                    end
                end
            end
            
            % save site selection
            plot_this = any(~isnan(psy_correct_new{t}),2);
            psy_correct_new{t}(~plot_this,:) = []; % remove all nan site
            psy_correct_true_new{t}(~plot_this,:) = []; % remove all nan site
            
            % find the finally used condition
            first_cond = condition_raw{t}(~isnan(condition_raw{t}(:,1)),1);
            unique_first_cond = unique(first_cond)';
            % find the location of this first_cond
            this_cond_find = []; this_cond_loc = []; this_cond_num = [];
            for i = 1:length(unique_first_cond)
                this_cond_find{i} = find(first_cond == unique_first_cond(i));
                this_cond_loc(i) = mean(this_cond_find{i}); % 该condition平均出现的位置
                this_cond_num(i) = length(this_cond_find{i});
            end
            select = this_cond_num>=10; % only select this_cond_num>=10
            % find the finally used condition
            finally_used_first_cond = unique_first_cond(this_cond_loc==max(this_cond_loc(select))); % location越大，也就是越靠后的文件
            finally_used_cond_range{t} = [finally_used_first_cond -finally_used_first_cond]; % 总是对称的
        end
        
        % ------------------------ Plot ------------------------------
        % for paper plot, only show curve
        Curve_only = 0;
        
        set(figure(figN),'Position',[100,50 1400,800], 'Name', 'Psychometric function in Four-target control task');clf
        for t = 1:6
            subplot(2,3,t)
            xh = min(unique_condition_new{t}):0.01:max(unique_condition_new{t});
%             % plot every site curve
%             for i = 1:length(psy_correct_plot{t})
%                 hold on
%                 plot(xh, cum_gaussfit(psy_perf_ctrl{stimtask_type_this(t),ind_this(t)}(psy_correct_plot{t}(i),:),xh),'-','color',c{99},'linewidth',0.7); % control
%                 %                 plot(xh, cum_gaussfit(psy_perf_stim{stimtask_type_this(t),ind_this(t)}(psy_correct_plot{t}(i),:),xh),'-','color',c{99},'linewidth',0.7); % stim
%             end
            
            rep_num = sum(~isnan(psy_correct_new{t}),1); % 每个condition重复site个数
            rep_true_num = sum(~isnan(psy_correct_true_new{t}),1); % 每个condition重复site个数
            
            % fitting
            if use_finally_cond
                % fitting of finally used condtion
                select_cond = logical(unique_condition_new{t}>=finally_used_cond_range{t}(1) & unique_condition_new{t}<=finally_used_cond_range{t}(2));
                fit_data_psy_cum = [];
                fit_data_psy_cum(:,1) = unique_condition_new{t}(select_cond);
                fit_data_psy_cum(:,2) = nanmean(psy_correct_new{t}(:,select_cond))';
                fit_data_psy_cum(:,3) = rep_num(select_cond)'; % simply this condtion repetition count
                
                fit_data_psy_cum_true = [];
                fit_data_psy_cum_true(:,1) = unique_condition_new{t}(select_cond);
                fit_data_psy_cum_true(:,2) = nanmean(psy_correct_true_new{t}(:,select_cond))';
                fit_data_psy_cum_true(:,3) = rep_true_num(select_cond)'; % simply this condtion repetition count
            else
                fit_data_psy_cum = [];
                fit_data_psy_cum(:,1) = unique_condition_new{t};
                fit_data_psy_cum(:,2) = nanmean(psy_correct_new{t})';
                fit_data_psy_cum(:,3) = rep_num'; % simply this condtion repetition count
                
                fit_data_psy_cum_true = [];
                fit_data_psy_cum_true(:,1) = unique_condition_new{t};
                fit_data_psy_cum_true(:,2) = nanmean(psy_correct_true_new{t})';
                fit_data_psy_cum_true(:,3) = rep_true_num'; % simply this condtion repetition count
            end
            
            method = 0; % 0=Maximum likelihood (default), 1=Square error
            tolerance = 0; % 0=No tolerance (default), other = error tolerance of the max abs (heading)
            [bb,tt] = cum_gaussfit_max1(fit_data_psy_cum,method,tolerance);
            psy_perf_all{t} = [bb, tt];
            
            [bb,tt] = cum_gaussfit_max1(fit_data_psy_cum_true,method,tolerance);
            psy_perf_true_all{t} = [bb, tt];

            plot(xh, cum_gaussfit(psy_perf_all{t},xh), '-','color',c{ind_this(t)},'linewidth',3);
            hold on
            plot(xh, cum_gaussfit(psy_perf_true_all{t},xh), '-','color',[0.5 0.5 0.5],'linewidth',2);
            
            if ~Curve_only
                % mean+sem
                hold on
                errorbar(unique_condition_new{t},nanmean(psy_correct_new{t}),nanstd(psy_correct_new{t})./sqrt(rep_num),'k.','markersize',10);
                errorbar(unique_condition_new{t},nanmean(psy_correct_true_new{t}),nanstd(psy_correct_true_new{t})./sqrt(rep_true_num),'.','color',[.5 .5 .5],'markersize',10);
            end

            % finally_used_cond_range
            hold on
            plot([finally_used_cond_range{t}(1) finally_used_cond_range{t}(1)],[0 1],'k:');
            plot([finally_used_cond_range{t}(2) finally_used_cond_range{t}(2)],[0 1],'k:');
            plot([0 0],[0 1],'k:');
            plot([finally_used_cond_range{t}(1) finally_used_cond_range{t}(2)],[0.5 0.5],'k:');
            plot([finally_used_cond_range{t}(1) finally_used_cond_range{t}(2)],[0.25 0.25],'k:');
            
            if use_finally_cond
                xlim([finally_used_cond_range{t}(1)*1.1 finally_used_cond_range{t}(2)*1.1]); % range of finally used condtion
            else
                xlim([min(unique_condition_new{t})*1.1 max(unique_condition_new{t})*1.1]); % range of all condition
            end
            
            ylim([0,1]);
            if min(unique_condition_new{t})~=finally_used_cond_range{t}(1)
                set(gca,'xtick',[min(unique_condition_new{t}) finally_used_cond_range{t}(1) finally_used_cond_range{t}(1)/2 0 finally_used_cond_range{t}(2)/2 finally_used_cond_range{t}(2) max(unique_condition_new{t})]);
            else
                set(gca,'xtick',[finally_used_cond_range{t}(1) finally_used_cond_range{t}(1)/2 0 finally_used_cond_range{t}(2)/2 finally_used_cond_range{t}(2)]);
            end
            set(gca,'ytick',[0:0.25:1]);
            
            if Curve_only
                set(gca,'xtick',[finally_used_cond_range{t}(1) 0 finally_used_cond_range{t}(2)]);
                set(gca,'xTickLabel',num2str(get(gca,'xTick')','%.1f'));
                set(gca,'ytick',[0:0.5:1]);
            end

            xlimit = xlim;
            
            % use 25% PSE instead of 50% PSE for 4 target true error
            PSE_true25{t} = norminv(0.25,bb,tt);
            
            if ~Curve_only
                str = [sprintf('\\it\\mu_{psy} = \\rm%2.3g\\circ',psy_perf_all{t}(1)) newline ...
                    sprintf('\\it\\sigma_{psy} = \\rm%2.3g\\circ',psy_perf_all{t}(2)) newline ...
                    sprintf('True\\it\\mu_{25%%} = \\rm%2.3g\\circ',PSE_true25{t}) newline ...  % use 25% PSE instead of 50% PSE for 4 target true error
                    sprintf('True\\it\\sigma_{50%%} = \\rm%2.3g\\circ',psy_perf_true_all{t}(2)) ... % use fitting threshold for 4 target true error
                    ];
                text(min(xlimit/1.1),0.7,str,'fontsize',10);

                if t==1
                    title('Fine trans.');
                elseif t==2
                    title('Fine rot.');
                elseif t==3
                    title('Fine motion type');
                elseif t==4
                    title('Coarse trans.');
                elseif t==5
                    title('Coarse rot.');
                elseif t==6
                    title('Coarse motion type');
                end
            end
            
            axis square
        end
        SetFigure(13)
        if ~Curve_only
            [~,~,~,monkey] = SetTitle(stimtask_type,monkey_included_for_analysis);
            str = ['Psychometric function in Four-target control task' monkey];
            suptitle(str);
        end
        figN = figN+1;
    end

    function f6p6(debug) % 'Comparison of 2T and 4T task performance'
        if debug  ; dbstack;   keyboard;      end
        
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        effect_type = 1; % 1: PSE, 2: threshold
        select_stars_type = [0 1]; % 0: dot stimulus; 1: triangle stimulus
        is_separate_stimAmp = 0; % 默认为0，combine 20 and 40 stim amplitude
        is_original_bias = 0;  % 1: only original bias in MT;  2: original bias(after-before);   3: bias positive = prefer direction
        
        % ---------------------- get data ----------------------------
        stimtask_type_this{1,1} = [1]; % fine 4T
        stimtask_type_this{2,1} = [4]; % coarse 4T
        stimtask_type_this{1,2} = [2 3]; % fine 2T
        stimtask_type_this{2,2} = [5 6]; % coarse 2T
        
        % compare: correct rate, bias, threshold, neuron-by-neuron
        for fc = 1:2
            mt_cr{fc} = ctrl_CR{stimtask_type_this{fc,1}}(:,3);
            T_cr{fc} = ctrl_CR{stimtask_type_this{fc,2}(1)}(:,1);
            R_cr{fc} = ctrl_CR{stimtask_type_this{fc,2}(2)}(:,2);
        end
        
        mt_cr = cellfun(@(x) x(methods_of_select{1},:),mt_cr,'uniformoutput',0);
        T_cr = cellfun(@(x) x(methods_of_select{1},:),T_cr,'uniformoutput',0);
        R_cr = cellfun(@(x) x(methods_of_select{1},:),R_cr,'uniformoutput',0);
        
        % ------------------------ Plot ------------------------------
        % 1. CR
        set(figure(figN),'name','CR for 2T vs 4T (neuron-by-neuron)','pos',[100,50,600,700]); clf
        % find neurons with both 2T and 4T，不一定需要2T-T task和2T-R task同时做
        for fc = 1:2
            for hr = 1:2 
                select = [mt_cr{fc} T_cr{fc} R_cr{fc}];
                select_this = find(~any(isnan(select(:,(1:1+hr))),2)); % 区分T and R
                
                subplot(2,2,hr+(fc-1)*2);
                plot(select(select_this,1),select(select_this,1+hr),'ko');
                xlim([50 100]);
                ylim([50 100]);
                xlabel('CR for 4T');
                ylabel('CR for 2T');
                axis square;
                SameScale(1);
                
                if fc == 1
                    if hr == 1
                        title('fine T task');
                    else
                        title('fine R task');
                    end
                else
                    if hr == 1
                        title('coarse T task');
                    else
                        title('coarse R task');
                    end
                end
            end
        end
        SetFigure(13)
        suptitle('CR for 2T vs 4T (neuron-by-neuron)');
        figN = figN+1; 
    end


    function f6p3(debug) % 'Microstimulation site location'
        if debug  ; dbstack;   keyboard;      end
        
        % ------------------ default parameter -----------------------
        % always separate monkey
        if length(monkey_included_for_analysis)>1
            disp('********** Please choose ONE monkey! **********');
            return
        end
        
        methods_of_select = {select_typical, 'Typical cells'};
        effect_type = 1;
        select_stars_type = [0 1]; % 0: dot stimulus; 1: triangle stimulus
        is_separate_stimAmp = 0; % 默认为0，combine 20 and 40 stim amplitude
        is_original_bias = 0;  % 1: only original bias in MT;  2: original bias(after-before);   3: bias positive = prefer direction
        
        % ---------------------- get data ----------------------------
        % at least one of six stim task do microstimulation: include 2 and
        % 4 target task, fine and coarse task
        for tt = 1:6 % six stim task
            do_stim_temp(:,tt) = any(~isnan(diffBias_task{tt}),2);
        end
        do_stim = any(do_stim_temp,2);
        
        % stim site information
        GridX = cell2mat({group_result.GridX}');
        GridY = cell2mat({group_result.GridY}');
        main_depth = [group_result.main_depth]';
        guidetube = [group_result.guidetube]';
        offset = [group_result.offset]';
        % file name
        cellID = [group_result.cellID]';
        
        % overlap all stimulated site in one example MRI Frame
        if monkey_included_for_analysis == 7
            hemisphere = 1; % left
            rowNo = 17; % use this MRI frame
            MRI_path = 'Z:\Data\MOOG\forDrawMapping\Ringbell\';
            MRI_offset_new = {[-1.5, 95, 74.9862],[0, 0]}; % first adjust by mapping
            gridRange = [10 23; 1 25];  % Grid plotting ragne = [xLow xHigh; yLow yHigh]
            maxZ = 300;   % Drawing range (in 100 um)
            AP0 = 18; % Ringbell_left_AP0 = 18
            
        elseif monkey_included_for_analysis == 16
            hemisphere = 1; % left
            rowNo = 13; % use this MRI frame
            MRI_path = 'Z:\Data\MOOG\forDrawMapping\Arthas\';
            MRI_offset_new = {[-0.5, 230, 216.22],[0, -2]}; % first adjust by mapping
            gridRange = [10 23; 1 25];  % Grid plotting ragne = [xLow xHigh; yLow yHigh]
            maxZ = 300;   % Drawing range (in 100 um)
            AP0 = 14; % Arthas_left_AP0 = 14
        end
                
        % ------------------------ Plot ------------------------------
        % plot MRI frame
        set(figure(figN),'name','Microstimulation site on MRI'); clf
        try
            MRI = imread([MRI_path 'row' num2str(rowNo) '.bmp']);
        catch
            MRI = imread([MRI_path 'row' num2str(rowNo) '.jpg']);
        end
        ratioyx = size(MRI,1)/size(MRI,2) * 10 * 0.8; % Auto keep ratio !!!   HH20180606
        x_lims = [MRI_offset_new{1}(1) - MRI_offset_new{1}(3)/2, MRI_offset_new{1}(1) + MRI_offset_new{1}(3)/2];
        y_lims = [MRI_offset_new{1}(2) - ratioyx * MRI_offset_new{1}(3)/2, MRI_offset_new{1}(2) + ratioyx * MRI_offset_new{1}(3)/2];
        MRI = flip(MRI,2); % 预先左右翻转MRI图片，后面和axes一起再反转过来
        h_MRI = image(x_lims + rowNo * MRI_offset_new{2}(1), y_lims + rowNo * MRI_offset_new{2}(2), MRI);
        set(h_MRI,'AlphaData',1);
        
        % Frame
        xlim([gridRange(2,1) gridRange(2,2)+1]);
        ylim([-30 maxZ]);
        
        figurePosition = [ 700  0  459  1000 ];
        aspectRatio = (range(ylim) * 100) / (range(xlim) * 1200);  % grid interval = 0.8 mm grid上两个孔圆心间距离
        set(gcf,'Position',[figurePosition(1:2) figurePosition(4)/aspectRatio figurePosition(4)]);
        daspect([1/0.8 10 1]); % This is betther man
        set(gca,'xdir','rev');
        
        % Set y scale to mm
        ytick_temp = 0:50:200;
        set(gca,'ytick',ytick_temp);
        set(gca,'yticklabel',ytick_temp/10);
        ylabel('Depth (mm)');
        
        xlabel('Grid Y No. (x 0.8 mm)');
        hemisphere_text = {'Left','Right'};
        AP_this = AP0 - rowNo;
        title(sprintf('Monkey %g, %s[%g], AP \\approx %g',monkey_included_for_analysis, hemisphere_text{hemisphere},rowNo,AP_this)); 
        
        % plot stimulated site
        select_do_stim = logical(methods_of_select{1,1} & do_stim);
        unique_cell = unique(cellID(select_do_stim));

        for j = 1:length(unique_cell)
            select_cell = logical(cellID == unique_cell(j) & select_do_stim);
            cellID_this = cellID(select_cell);
            GridX_this = GridX(select_cell); % row
            GridY_this = GridY(select_cell);
            
            depth_this = main_depth(select_cell);
            gt_this = guidetube(select_cell);
            offset_this = offset(select_cell);
            offSet = round((gt_this  + offset_this - 2.0) * 100);  % Related to Guide Tube 2.0 cm!!
            xx = GridY_this;
            yy = offSet + round(depth_this/100);
            hold on
            h = scatter(xx,yy,30,'markeredgecolor','w','markerfacecolor','k','markerfacealpha',0.2,'markeredgealpha',(7-abs(GridX_this-rowNo))/7); % near row17, more clearer
            h.Tag = 'stimsite';
        end
        box off
        
        % show whole MRI frame
        if monkey_included_for_analysis == 7
            xlim([-40 40]);
            ylim([-200 400]);
        elseif monkey_included_for_analysis == 16
            xlim([-40 50]);
            ylim([-200 400]);
        end
        
        
        
        

        
%         % ------------------------ Plot ------------------------------
%         unique_GridX = unique(GridX);
%   
%         for i = 1:length(unique_GridX) % 每一行一张图
%             %           for i = 5 % for testing
%             select_x = logical(GridX == unique_GridX(i) & methods_of_select{1,1} & do_stim);
%             unique_cell = unique(cellID(select_x));
%             if ~isempty(unique_cell)
%                 % save all marker location
%                 xx_all = [];
%                 yy_all = [];
%                 select_actual_plot = false(size(select_x));
%                 
%                 % wait for plot
%                 for j = 1:length(unique_cell)
%                     select_cell = logical(cellID == unique_cell(j) & select_x);
%                     cellID_this = cellID(select_cell);
%                     GridX_this = unique_GridX(i);
%                     GridY_this = GridY(select_cell);
%                     depth_this = main_depth(select_cell);
%                     gt_this = guidetube(select_cell);
%                     offset_this = offset(select_cell);
%                     offSet = round((gt_this  + offset_this - 2.0) * 100);  % Related to Guide Tube 2.0 cm!!
%                     xx = GridY_this;
%                     yy = offSet + round(depth_this/100);
%                     xx_all = [xx_all;xx];
%                     yy_all = [yy_all;yy];
%                     select_actual_plot(select_cell) = true;
%                 end
%                 
%                 % plot marker together
%                 if monkey_included_for_analysis==7
%                     DrawMapping_Ringbell(1, [unique_GridX(i) 1], figN); % for Ringbell
%                 elseif monkey_included_for_analysis==16
%                     DrawMapping_Arthas(1, [unique_GridX(i) 1], figN); % for Ringbell
%                 else
%                     keyboard
%                 end
%                 figure(figN);
%                 %                 h_all(hn) = scatter(xx_all,yy_all,50*ones(size(xx_all)),'g','filled','visible','on');hold on;
%                 h_all = plot(xx_all,yy_all,'o','markersize',6,'markerfacecolor','g','color','w','visible','on');
%                 h_all.Tag = 'stimsite';
%                 % Show individual cell selected from the figure. HH20150424
%                 %                 set([gca h_all(hn)],'ButtonDownFcn',{@Show_individual_cell, h_all(hn), select_actual_plot});
%                 set(figN,'WindowButtonDownFcn',{@Show_individual_cell, select_actual_plot});
%                 figN = figN + 1;
%             end
%         end
%                         
%         % Show individual cell selected from the figure. @HH20150424
%         % changed by Lwh 20210831, 调用较慢。
%         function Show_individual_cell(~,~, select_for_this, couples)
%             %         function Show_individual_cell(hObject,~,h_line, select_for_this, couples)
%             if nargin < 5
%                 couples = 1; % HH20150515 Show coupled cells in the figure
%             end
%             
%             % get handle from figure
%             handles_this = guihandles(gca);
%             allX = handles_this.stimsite.XData;
%             allY = handles_this.stimsite.YData;
%             
%             try
%                 h_marker = handles_this.showselect;
%                 t_showcell = handles_this.showcellname;
%             catch
%                 h_marker = [];
%                 t_showcell = [];
%             end
% 
%             n_single_group = length(allX)/couples;
%             
%             % Select cell from figure
%             pos = get(gca,'currentPoint'); 
% %             pos = get(hObject,'currentPoint'); 
%             posX = pos(1,1); posY = pos(1,2);
%             [min_dis,id] = min(abs(((posX-allX)/range(xlim)).^2+((posY-allY)/range(ylim)).^2));
%             if min_dis > (range(xlim)^2+range(ylim)^2)/100 +inf ; return; end
%             
%             id = mod(id-1,n_single_group) + 1; % Deal with coupled cells
%             ori_cell_no = find(cumsum(select_for_this)==id,1);
%             
%             % Plotting
%             if ~isempty(h_marker) ; try delete(h_marker); catch ; end ;end
%             all_inds = mod(1:length(allX),n_single_group) == mod(id,n_single_group); % Deal with coupled cells
%             h_marker = plot(allX(all_inds),allY(all_inds),'x','color','m','markersize',8,'linew',3);
%             h_marker.Tag = 'showselect';
%             
%             cellID1 = [group_result.cellID]';
%             if ~isempty(t_showcell) ; try delete(t_showcell); catch ; end ;end
%             t_showcell = text(10,250,num2str(cellID1(ori_cell_no)),'fontsize',13,'color','w','FontWeight','bold');
%             t_showcell.Tag = 'showcellname';
%             disp(cellID1(ori_cell_no));
%         end
    end

    function f6p4(debug) % 'Motion type behavioral performance'
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        stimtask_type_this{1} = [1];
        stimtask_type_this{2} = [4];
        
        % behavioral performance
        mt_cr = [];
        for fc = 1:length(stimtask_type_this)
            mt_cr{fc} = ctrl_CR{stimtask_type_this{fc}}(:,3); % control task CR
        end
        mt_cr = cellfun(@(x) x(methods_of_select{1,1},:),mt_cr,'uniformoutput',0);

        % remove nan
        for fc = 1:2
            mt_cr{fc}(isnan(mt_cr{fc})) = [];
        end
        
        % ------------------------ Plot ------------------------------
        set(figure(figN),'Position',[50,50 1400,500], 'Name', 'Motion type CR');
        binsize = 5;
        xedges = 40:binsize:100;
        xedges_mid = xedges+binsize/2;
        xedges_mid(end) = [];
        
        for i = 1:length(stimtask_type_this)
            subplot(1,2,i)
            histogram(mt_cr{i},xedges,'normalization','count');
            binN{i} = histcounts(mt_cr{i},xedges,'normalization','probability');
            axis tight
            xlim([min(xedges) max(xedges)]);
            ylimit = ylim;
            
            % mean
            mean_cr(i) = nanmean(mt_cr{i});
            sem_cr(i) = nanstd(mt_cr{i}) / sqrt (length(mt_cr{i}(~isnan(mt_cr{i}))));
            
            % medain
            median_cr(i) = nanmedian(mt_cr{i});
            
            hold on
            plot([mean_cr(i) mean_cr(i)],[0 ylimit(2)*1.1],'k--');
            text(mean_cr(i),ylimit(2)*1.1,num2str(mean_cr(i)),'horizontalalign','center');
            ylim([0 ylimit(2)*1.2]);
            
            % t test with chance level 50%
            [~,p(i)] = ttest(mt_cr{i},50);
            
            text(45, ylimit(2),sprintf('n=%g',length(mt_cr{i})));
            
            if i == 1
                title('Fine motion type');
            else
                title('Coarse motion type');
            end
        end
        
        SetFigure(13)
        [til] = SetTitle(stimtask_type_this,monkey_included_for_analysis);
        str = ['Motion type behavioral performance: correct rate' newline til];
        suptitle(str);
        figN = figN + 1;
        
        %         % boxplot
        %         set(figure(figN),'Position',[200,50 400,400], 'Name', 'Motion type CR');clf
        %         plot([ones(size(mt_cr{1})) 2*ones(size(mt_cr{2}))],[mt_cr{1} mt_cr{2}],'k.');
        %         hold on
        %         boxplot([mt_cr{1} mt_cr{2}]);
        %         xticklabels({'Fine task','Coarse task'});
        %         ylim([45 100]);
        %         SetFigure(15);
        %         figN = figN + 1;
        
        
        % 插值 for smooth histogram
        set(figure(figN),'Position',[200,100 750,700], 'Name', 'Motion type CR'); clf
        binN_interp = [];
        for i = 1:length(stimtask_type_this)
            subplot(1,2,i)
            % remove fist and last several zero
            select = [find(binN{i}~=0,1,'first'): find(binN{i}~=0,1,'last')];
            binN_adj = binN{i}(select);
            xedges_mid_adj = xedges_mid(select);
            
            % add first and last 0 for interpolation
            binN_adj = [0 binN_adj 0];
            xedges_mid_adj = [min(xedges_mid_adj)-binsize/2 xedges_mid_adj max(xedges_mid_adj)+binsize/2];
            xedges_mid_interp = min(xedges_mid_adj):0.1:max(xedges_mid_adj);
            
            if sum(binN_adj)~=0
                binN_interp = interp1(xedges_mid_adj,binN_adj,xedges_mid_interp,'pchip');
                maxN = max(binN_interp);
                
                hold on
                plot(binN_interp,xedges_mid_interp,'-','color',c{3});
                plot(-binN_interp,xedges_mid_interp,'-','color',c{3}); % copy for other side
                
                % add cross = meadian
                plot([-0.2 0.2],[median_cr(i) median_cr(i)],'-','color',c{3}); % plot meadian
                plot([0 0],[median_cr(i)-2.5 median_cr(i)+2.5],'-','color',c{3});
                text(0.1,median_cr(i)+2,num2str(median_cr(i)));
                
                ylim([40 100]);
                xlim([-0.5 0.5]);
                xlabel('Cases/N');
                ylabel('Correct rate (%)');
                
                % t test
                text(-0.4,95,sprintf('t test, p=%2.2g',p(i)));
                text(-0.4, 90,sprintf('n=%g',length(mt_cr{i})));
                
                [~,Two_or_Four,fine_or_coarse,monkey_this] = SetTitle(stimtask_type_this{i},monkey_included_for_analysis);
                title([Two_or_Four fine_or_coarse]);
            end
        end
        str = [];
        str = ['Motion type behavioral performance: Correct Rate'  monkey_this];
        suptitle(str);
        
        SetFigure(15);
        figN = figN + 1;
    end

    function f6p5(debug) % Translation and Rotation task behavioral performance (control task)
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};

        % task                stimCueS stimHonly stimRonly stimCo stimCoH stimCoR
        % stimtask_type_code     1         2         3       4       5      6
        
        % behavioral performance
        % fine and coarse
        % 2-T task: T and R
        % 4-T task: T and R
        % 2+4-T task: T and R
        stimtask_type_this{1,1} = [2 3]; % fine 2-T
        stimtask_type_this{2,1} = [5 6]; % coarse 2-T
        stimtask_type_this{1,2} = [1]; % fine 4-T
        stimtask_type_this{2,2} = [4]; % coarse 4-T
        stimtask_type_this{1,3} = [1 2 3]; % fine 2+4-T
        stimtask_type_this{2,3} = [4 5 6]; % coarse 2+4-T
        
        for fc = 1:2 % fine or coarse
            for tf = 1:2 % two- or four-target
                task_num = length(stimtask_type_this{fc,tf});
                T_cr{fc,tf} = nan(length(loadM_data),2); % prepare 3 site for each cell deal to [2 3] or [1 2 3]
                R_cr{fc,tf} = nan(length(loadM_data),2); % prepare 3 site for each cell deal to [2 3] or [1 2 3]
                T_cr_true{fc,tf} = nan(length(loadM_data),2); % prepare 3 site for each cell deal to [2 3] or [1 2 3]
                R_cr_true{fc,tf} = nan(length(loadM_data),2); % prepare 3 site for each cell deal to [2 3] or [1 2 3]
                
                for i = 1:task_num
                    T_cr{fc,tf}(:,i) = ctrl_CR{stimtask_type_this{fc,tf}(i)}(:,1);
                    R_cr{fc,tf}(:,i) = ctrl_CR{stimtask_type_this{fc,tf}(i)}(:,2);
                    
                    T_cr_true{fc,tf}(:,i) = ctrl_CR_true{stimtask_type_this{fc,tf}(i)}(:,1);
                    R_cr_true{fc,tf}(:,i) = ctrl_CR_true{stimtask_type_this{fc,tf}(i)}(:,2);
                end
            end
        end
        
        % select
        T_cr = cellfun(@(x) x(methods_of_select{1,1},:),T_cr,'uniformoutput',0);
        R_cr = cellfun(@(x) x(methods_of_select{1,1},:),R_cr,'uniformoutput',0);
        T_cr_true = cellfun(@(x) x(methods_of_select{1,1},:),T_cr_true,'uniformoutput',0);
        R_cr_true = cellfun(@(x) x(methods_of_select{1,1},:),R_cr_true,'uniformoutput',0);
        
        % combine task type
        T_cr = cellfun(@(x) x(:),T_cr,'uniformoutput',0);
        R_cr = cellfun(@(x) x(:),R_cr,'uniformoutput',0);
        T_cr_true = cellfun(@(x) x(:),T_cr_true,'uniformoutput',0);
        R_cr_true = cellfun(@(x) x(:),R_cr_true,'uniformoutput',0);
        
        % remove nan
        for fc = 1:2 % fine or coarse
            for tf = 1:2 % two- or four-target
                T_cr{fc,tf}(isnan(T_cr{fc,tf})) = [];
                R_cr{fc,tf}(isnan(R_cr{fc,tf})) = [];
                
                T_cr_true{fc,tf}(isnan(T_cr_true{fc,tf})) = [];
                R_cr_true{fc,tf}(isnan(R_cr_true{fc,tf})) = [];
            end
            
            % combine 2T and 4T
            T_cr{fc,3} = [T_cr{fc,1};T_cr{fc,2}];
            R_cr{fc,3} = [R_cr{fc,1};R_cr{fc,2}];
            T_cr_true{fc,3} = [T_cr_true{fc,1};T_cr_true{fc,2}];
            R_cr_true{fc,3} = [R_cr_true{fc,1};R_cr_true{fc,2}];
        end

        % ------------------------ Plot ------------------------------
        set(figure(figN),'Position',[50,50 1400,500], 'Name', 'Translation and Rotation CR'); clf
        % bar plot
        binsize = 5;
        xedges = 0:binsize:100;
        xedges_mid = xedges+binsize/2;
        xedges_mid(end) = [];
        for fc = 1:2 % fine or coarse
            for tf = 1:3 % two- or four-target or 2+4 T
                subplot(2,3,(fc-1)*3+tf)
                hold on
                
                for hr = 1:2
                    if hr == 1
                        hold on
                        % normal
                        histogram(T_cr{fc,tf},xedges,'normalization','probability','facecolor',c{hr});
                        mean_cr{fc,tf}(hr) = mean(T_cr{fc,tf});
                        median_cr{fc,tf}(hr) = median(T_cr{fc,tf});
                        sem_cr{fc,tf}(hr) = std(T_cr{fc,tf}) / sqrt(length(T_cr{fc,tf}));
                        binN{fc,tf}(hr,:) = histcounts(T_cr{fc,tf},xedges,'normalization','probability');
                        N_this{fc,tf}(hr,:) = length(T_cr{fc,tf});
                        
                        % true error
                        histogram(T_cr_true{fc,tf},xedges,'normalization','probability','facecolor',c{99});
                        mean_cr_true{fc,tf}(hr) = mean(T_cr_true{fc,tf});
                        median_cr_true{fc,tf}(hr) = median(T_cr_true{fc,tf});
                        sem_cr_true{fc,tf}(hr) = std(T_cr_true{fc,tf}) / sqrt(length(T_cr_true{fc,tf}));
                        binN_true{fc,tf}(hr,:) = histcounts(T_cr_true{fc,tf},xedges,'normalization','probability');
                        N_this_true{fc,tf}(hr,:) = length(T_cr_true{fc,tf});
                        
                        xlabel('Correct rate');
                        ylabel('Cases/N');
                        
                    else
                        hold on
                        % normal
                        histogram(R_cr{fc,tf},xedges,'normalization','probability','facecolor',c{hr});
                        mean_cr{fc,tf}(hr) = mean(R_cr{fc,tf});
                        median_cr{fc,tf}(hr) = median(R_cr{fc,tf});
                        sem_cr{fc,tf}(hr) = std(R_cr{fc,tf}) / sqrt(length(R_cr{fc,tf}));
                        binN{fc,tf}(hr,:) = histcounts(R_cr{fc,tf},xedges,'normalization','probability');
                        N_this{fc,tf}(hr,:) = length(R_cr{fc,tf});
                        
                         % true error
                        histogram(R_cr_true{fc,tf},xedges,'normalization','probability','facecolor',c{99});
                        mean_cr_true{fc,tf}(hr) = mean(R_cr_true{fc,tf});
                        median_cr_true{fc,tf}(hr) = median(R_cr_true{fc,tf});
                        sem_cr_true{fc,tf}(hr) = std(R_cr_true{fc,tf}) / sqrt(length(R_cr_true{fc,tf}));
                        binN_true{fc,tf}(hr,:) = histcounts(R_cr_true{fc,tf},xedges,'normalization','probability');
                        N_this_true{fc,tf}(hr,:) = length(R_cr_true{fc,tf});
                    end
                    
                    axis tight
                    xlim([min(xedges) max(xedges)]);
                    ylimit = ylim;
                    
                    plot([mean_cr{fc,tf}(hr) mean_cr{fc,tf}(hr)],[0 ylimit(2)*1.1],'-','color',c{hr});
                    text(mean_cr{fc,tf}(hr),ylimit(2)*1.1,sprintf('mean = %.2g (± %.2g), median = %.2g', mean_cr{fc,tf}(hr), sem_cr{fc,tf}(hr), median_cr{fc,tf}(hr)),'horizontalalign','center','color',c{hr});
                    
                    ylim([0 0.5]);
                end
                
                if ~isempty(T_cr{fc,tf}) && ~isempty(R_cr{fc,tf})
                    % ttest for h and r
                    [~,p_ttest{fc,tf}] = ttest2(T_cr{fc,tf},R_cr{fc,tf});
                    
                    % ranksum test
                    p_ranksum{fc,tf} = ranksum(T_cr{fc,tf},R_cr{fc,tf});
                    p_ranksum_true{fc,tf} = ranksum(T_cr_true{fc,tf},R_cr_true{fc,tf});
                else
                    p_ttest{fc,tf} = nan;
                    p_ranksum{fc,tf} = nan;
                    p_ranksum_true{fc,tf} = nan;
                    binN{fc,tf} = zeros(size(binN{fc,tf}));
                    binN_true{fc,tf} = zeros(size(binN_true{fc,tf}));
                end
            end
        end
        [~,~,~,monkey_this] = SetTitle(stimtask_type_this{fc,tf},monkey_included_for_analysis);
        str = ['T and R behavioral performance: Correct Rate ('  monkey_this ,')'];
        suptitle(str);
        SetFigure(15);
        figN = figN + 1;
        
        % 插值
        set(figure(figN),'Position',[200,100 750,700], 'Name', 'Translation and Rotation CR'); clf
        binN_interp = [];
        for fc = 1:2 % fine or coarse
            for tf = 1:3 % two- or four-target or 2+4 T
                subplot(2,3,(fc-1)*3+tf)
                hold on
                for hr = 1:2 % h or r
                    % remove fist and last several zero
                    select = [find(binN{fc,tf}(hr,:)~=0,1,'first'): find(binN{fc,tf}(hr,:)~=0,1,'last')];
                    binN_adj = binN{fc,tf}(hr,select);
                    xedges_mid_adj = xedges_mid(select);
                    
                    select2 = [find(binN_true{fc,tf}(hr,:)~=0,1,'first'): find(binN_true{fc,tf}(hr,:)~=0,1,'last')];
                    binN_adj_true = binN_true{fc,tf}(hr,select2);
                    xedges_mid_adj_true = xedges_mid(select2);

                    % add first and last 0 for interpolation
                    binN_adj = [0 binN_adj 0];
                    xedges_mid_adj = [min(xedges_mid_adj)-binsize/2 xedges_mid_adj max(xedges_mid_adj)+binsize/2];
                    xedges_mid_interp = min(xedges_mid_adj):0.1:max(xedges_mid_adj);
                    
                    binN_adj_true = [0 binN_adj_true 0];
                    xedges_mid_adj_true = [min(xedges_mid_adj_true)-binsize/2 xedges_mid_adj_true max(xedges_mid_adj_true)+binsize/2];
                    xedges_mid_interp_true = min(xedges_mid_adj_true):0.1:max(xedges_mid_adj_true);
                    
                    % plot true error
                    if sum(binN_adj_true)~=0
                        binN_interp_true = interp1(xedges_mid_adj_true,binN_adj_true,xedges_mid_interp_true,'spline'); % spline % pchip
                        maxN = max(binN_interp_true);
                        
                        plot(sign(1.5-hr).*binN_interp_true,xedges_mid_interp_true,'-','color',c{99}); % sign(1.5-hr): t plot on positive, r plot on negative
                        plot(sign(1.5-hr).*[maxN/2-0.1 maxN/2+0.1],[median_cr_true{fc,tf}(hr) median_cr_true{fc,tf}(hr)],'-','color',c{99}); % plot median
                        plot(sign(1.5-hr).*[maxN/2 maxN/2],[median_cr_true{fc,tf}(hr)-2.5 median_cr_true{fc,tf}(hr)+2.5],'-','color',c{99});
                        
                        text(sign(1.5-hr).*0.25,median_cr_true{fc,tf}(hr),sprintf('%.2g',median_cr_true{fc,tf}(hr)),'horizontalalignment','center');
                    end
                    
                    % plot normal
                    if sum(binN_adj)~=0
                        binN_interp = interp1(xedges_mid_adj,binN_adj,xedges_mid_interp,'spline'); % spline % pchip
                        maxN = max(binN_interp);
                        
                        plot(sign(1.5-hr).*binN_interp,xedges_mid_interp,'-','color',c{hr}); % sign(1.5-hr): t plot on positive, r plot on negative
                        plot(sign(1.5-hr).*[maxN/2-0.1 maxN/2+0.1],[median_cr{fc,tf}(hr) median_cr{fc,tf}(hr)],'k-'); % plot median
                        plot(sign(1.5-hr).*[maxN/2 maxN/2],[median_cr{fc,tf}(hr)-2.5 median_cr{fc,tf}(hr)+2.5],'k-');
                        
                        text(sign(1.5-hr).*0.25,median_cr{fc,tf}(hr),sprintf('%.2g',median_cr{fc,tf}(hr)),'horizontalalignment','center','color',c{hr});
                        text(sign(1.5-hr).*0.3,95,sprintf('N = %g',N_this{fc,tf}(hr,:)),'horizontalalignment','center');
                    end
                end
                
                plot([0 0],[0 100],'k-');
                if tf == 1
                    ylabel('Correct rate');
                end
                xlabel('Cases/N');
                ylim([0 100]);
                yticks([0:25:100]);
                
                % show t-test
                str = [];
                str = ['ranksum' newline sprintf('p=%.2g',p_ranksum{fc,tf}) newline sprintf('p=%.2g',p_ranksum_true{fc,tf})];
                text(0.1,30,str);
                [~,Two_or_Four,fine_or_coarse,monkey_this] = SetTitle(stimtask_type_this{fc,tf},monkey_included_for_analysis);
                title([fine_or_coarse Two_or_Four])
            end
        end
        str = [];
        str = ['T and R behavioral performance: Correct Rate'  monkey_this];
        suptitle(str);
        
        
        SetFigure(15);
        figN = figN + 1;
    end

    function f6p7(debug) %  'Eyetrace during VisualOn period'
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        
        % separate for monkey
        % always separate monkey
        if length(monkey_included_for_analysis)>1
            disp('********** Please choose ONE monkey! **********');
            return
        end
        % combine all task type and all trial
        stimtask_type_this{1} = [1 2 3];
        stimtask_type_this{2} = [4 5 6];
        
        % ---------------------- get data ----------------------------
        % combine all task type and all trial
        eyetraceX = [];
        eyetraceY = [];
        eye_dia = [];
        
        progressbar('Finding Eyetrace...');
        num_this = 1;
        for n = 1:length(loadM_data)
            if methods_of_select{1,1}(n)
                for t = 1:6
                    eyeX_temp = group_result(n).eye_data_VisualOn_X{t}; % % 每一行一个trial，每一列一个bin
                    eyeY_temp = group_result(n).eye_data_VisualOn_Y{t};
                    
                    if ~isnan(eyeX_temp)
                        eyetraceX = [eyetraceX;eyeX_temp];
                        eyetraceY = [eyetraceY;eyeY_temp];
                    end
                end
                
                % 计数
                progressbar(num_this/sum(methods_of_select{1,1}));
                num_this = num_this + 1;
            end
        end
        
        set(figure(figN),'Position',[200,100 1200,600], 'Name', 'Eyetrace distribution'); clf
        subplot(1,2,1)
        histogram(eyetraceX);
        xlabel('Eye position X (°)');
        xlim([-5 5]);
        set(gca,'xtick',[-5:5]);
        subplot(1,2,2)
        histogram(eyetraceY);
        xlabel('Eye position Y (°)');
        xlim([-5 5]);
        set(gca,'xtick',[-5:5]);
        SetFigure(15)
        [til] = SetTitle(stimtask_type_this,monkey_included_for_analysis);
        suptitle(['Eye position distribution' newline til]);
        figN = figN + 1;
        
        set(figure(figN),'Position',[1000,100 700,700], 'Name', 'Eyetrace plot'); clf
        plot(eyetraceX,eyetraceY,'k.');
        xlabel('Eye position X (°)');
        ylabel('Eye position Y (°)');
        SameScale
        set(gca,'xtick',[-6:6],'ytick',[-6:6]);
        
        % proportion of eye position in diameter 2 ([-1 1]);
        eye_dia_temp = sqrt(eyetraceX.^2 + eyetraceX.^2); % 距离中心点距离
        eye_dia = eye_dia_temp(:);
        pro_2 = sum(eye_dia<1)/length(eye_dia)*100;
        text(-4.5,4.5,sprintf('Position<[-1 1] = %.4g%%',pro_2));
        
        pro_1 = sum(eye_dia<0.5)/length(eye_dia)*100;
        text(-4.5,4,sprintf('Position<[-0.5 0.5] = %.4g%%',pro_1));
        
        SetFigure(15)
        box on
        title(['Eye position' newline til]);
        figN = figN + 1;
    end

    function f6p8(debug) % 'μ distribution in 4AFC'
        if debug  ; dbstack;   keyboard;      end
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        stimtask_type_temp{1} = [1 2 3];
        stimtask_type_temp{2} = [4 5 6];
%         stimtask_type_this = cell2mat(stimtask_type);
        stimtask_type_this = cell2mat(stimtask_type_temp);
        
        % normalize PSE for control trial
        % combine fine/coarse task
        PSE_ctrl{1} = [nor_PSE{1}(methods_of_select{1,1},:); nor_PSE{2}(methods_of_select{1,1},:); nor_PSE{3}(methods_of_select{1,1},:)]; % fine: 1, 2, 3
        PSE_ctrl{2} = [nor_PSE{4}(methods_of_select{1,1},:); nor_PSE{5}(methods_of_select{1,1},:); nor_PSE{6}(methods_of_select{1,1},:)]; % coarse: 4, 5, 6
        
        % ---------------------- get data ----------------------------
        set(figure(figN),'Position',[200,100 600,1200], 'Name', 'Flow-pattern PSE distribution in 4AFC'); clf
        for fc = 1:2
            for hr = 1:3
                subplot(3,2,fc+(hr-1)*2)
                PSE_this = PSE_ctrl{fc}(:,hr);
                PSE_this(isnan(PSE_this)) = [];
                
                if hr~=3
                    histogram(PSE_this);
                else
                    xbin = linspace(-2,2,41);
                    histogram(PSE_this,xbin);
                end
                
                
                PSE_mean = nanmean(PSE_this);
                PSE_sem = nanstd(PSE_this)/sqrt(length(PSE_this));
                % t test with 0
                [~,p] = ttest(PSE_this,0);
                
                str = [sprintf('n = %g',length(PSE_this)) newline ...
                    sprintf('mean±SEM = %2.2g ± %2.2g', PSE_mean, PSE_sem) newline ...
                    sprintf('t test, p = %2.2g', p)];
                
                axis tight
                ylimit = ylim;
                ylim([0 ylimit(2)+ylimit(2)/10]);
                text(-1.9,ylimit(2),str);
                if fc == 1 && hr == 1
                    title('Fine version');
                    xlabel('PSE (°)');
                elseif  fc == 2 && hr == 1
                    title('Coarse version');
                    xlabel('PSE (% coherence)');
                end
                if hr == 3
                    xlabel('PSE'); % motion axes is normalized
                end
                
                ylabel('Caese');
                
                hold on
                plot(PSE_mean,ylimit(2)+ylimit(2)/12,'kv');
                plot([0 0],[0 ylimit(2)+ylimit(2)/10],'k-');
                
                axis square
            end
        end
        
        SetFigure(13)
        [til] = SetTitle(stimtask_type_temp,monkey_included_for_analysis);
        suptitle(['Flow-pattern non-stim. PSE distriubution in 4AFC' newline til]);
        figN = figN + 1;
    end


%%
    function [diff_thistask, S_unit, NS_unit, unique_stimAmp_this, select_save] = get_stim_data_uA(stimtask_type_this,methods_of_select_this,is_separate_stimAmp,is_original_bias,select_stars_type)
        
        % 区分旧视觉刺激（三角形，stars_type=1）还是新视觉刺激(圆点，stars_type=0), 因为Translation task中不受影响，只影响Rotation task，所以分开讨论。
        % ====================  Condition base  ========================
        % task                stimCueS stimHonly stimRonly stimCo stimCoH stimCoR
        % stimtask_type_code     1         2         3       4       5      6
        % stars_type            [0]      [0 1]      [0]     [0]    [0 1]   [0]
        % ===============================================================
        %
        % 根据stimtask_type获得需要combine的task
        % stimtask_type{1}: fine task部分，stimtask_type{2}: coarse task部分
        % stimtask_type{1}=[1 2 3], stimtask_type{2}=[4 5 6]：2&4-target fine&coarse task
        % stimtask_type{1}=[1 2 3], stimtask_type{2}=[]     ：2&4-target fine task
        % stimtask_type{1}=[2 3],   stimtask_type{2}=[]     ：2-target fine task
        % stimtask_type{1}=[1],     stimtask_type{2}=[]     ：4-target fine task
        % stimtask_type{1}=[],      stimtask_type{2}=[4 5 6]：2&4-target coarse task
        % stimtask_type{1}=[],      stimtask_type{2}=[5 6]  ：2-target coarse task
        % stimtask_type{1}=[],      stimtask_type{2}=[4]    ：4-target coarse task
        % stimtask_type{1}=[2 3],   stimtask_type{2}=[5 6]  ：2-target fine&coarse task
        % stimtask_type{1}=[1],     stimtask_type{2}=[4]    ：4-target fine&coarse task
        %
        % parameters:
        % is_separate_stimAmp  0: combine difference stim amplitude;  1:separate
        % is_original_bias  0: bias positive = prefer direction; 1: only
        % original bias in MT;  2: original bias(after-before);  3:normalize bias (divided by the threshold under non-stim trials   (Xuefei Yu, 2018))
        % select_stars_type   0: dot stimulus(correct in rotation)    1: triangle stimulus (wrong in rotation)
        %
        % data structure, example for diff_thistask:
        % 1. diffBias_thistask{1,fr},diffBias_thistask{2,fr}: seperated by 'is_separate_stimAmp' (20 or 40 uA)
        % 2. in each cell like diffBias_thistask{1,1}: combine difference selected stimtask_type_this, and in fine
        %       2.1 only selected by 'is_separate_stimAmp' and 'select_stars_type', length = length(dataset) * stimtask_type_this num
        %       2.2 for example: 300 dataset, stimtask_type_this = [4 5 6], so that length(diffBias_thistask{1}) = 300 * 3 = 900;
        %       2.3 后续处理如果要进行筛选，只需要重构数组: reshape(diffBias_thistask{1},length(dataset),length(stimtask_type_this))
        %       2.4 后续处理如果不需要筛选，只需要去掉NAN即可
        if nargin == 1
            methods_of_select_this = true(length(group_result),1);
        elseif  nargin == 2
            is_separate_stimAmp = 0;
            is_original_bias = 0;
            select_stars_type = [0 1]; % select all star type
        elseif nargin == 3
            is_original_bias = 0;
            select_stars_type = [0 1]; % select all star type
        elseif nargin == 4
            select_stars_type = [0 1]; % select all star type
        end
        
        % ---------------------- get data ----------------------------
        % all type stim amplitude in 6 stim protocol
        unique_stimAmp_this = unique(stimAmp(:));
        unique_stimAmp_this = unique_stimAmp_this(~isnan(unique_stimAmp_this));
        if ~is_separate_stimAmp % 不分开不同的电刺激强度
            unique_stimAmp_this = 9999;
        end
        
        % 如果輸入的是uA的准確值，則只画该uA的unit
        if is_separate_stimAmp~=0 && is_separate_stimAmp~=1
            unique_stimAmp_this = is_separate_stimAmp;
        end
        
        % 区分旧视觉刺激（三角形，stars_type=1）还是新视觉刺激(圆点，stars_type=0)
        stars_type = cell2mat({group_result.stars_type})';
        
        % preparation for output
        diffBias_thistask = cell(length(unique_stimAmp_this),1);
        diffThresh_thistask = cell(length(unique_stimAmp_this),1);
        diffBias_p_thistask = cell(length(unique_stimAmp_this),1);
        diffThresh_p_thistask = cell(length(unique_stimAmp_this),1);
        select_save = cell(length(unique_stimAmp_this),1);
        
        % select for stars_type: 0, 1, [0 1](combine)
        select1 = false(size(stars_type)); % select for stars_type
        for st = 1:length(select_stars_type)
            select1_temp = logical(stars_type == select_stars_type(st));
            select1 = logical(select1_temp | select1);
        end
        
        % select for stimtask_type_this and stimAmp
        for amp = 1:length(unique_stimAmp_this)
            for fc = 1:length(stimtask_type_this) % fine or coarse，需要combine一起，所以最后的index没有fr
                if ~isempty(stimtask_type_this{fc})
                    for i = 1:length(stimtask_type_this{fc}) % 具体task, 需要combine一起，所以最后的index没有j
                        if is_separate_stimAmp % 不分开不同的电刺激强度
                            select2 = stimAmp(:,stimtask_type_this{fc}(i)) == unique_stimAmp_this(amp); % select for stimAmp and stimtask_type_this
                        else
                            select2 = true(size(stimAmp,1),1);
                        end
                        
                        select = logical(select2 & select1 & methods_of_select_this); % select for stimAmp, stimtask_type_this and stars_type
                        temp = []; temp2 = []; temp3 = []; temp4 = [];
                        
                        % combine selected stimtask_type, seperate or combine difference stimAmp
                        % Bais
                        if is_original_bias == 0 % 默认，正值为bias to prefer direction
                            temp = diffBias_task{stimtask_type_this{fc}(i)};
                        elseif is_original_bias==1 % only original bias in MT
                            temp = diffBias_task_originalmt{stimtask_type_this{fc}(i)};
                        elseif is_original_bias==2 % original bias: after-before
                            temp = diffBias_task_origin{stimtask_type_this{fc}(i)};
                        elseif is_original_bias == 3 % normlize PSE shift (divided by the threshold under non-stim trials
                            temp = nor_diffBias_task{stimtask_type_this{fc}(i)};
                        elseif is_original_bias == 4 % normlize PSE shift and original bias: after-before
                            temp = nor_diffBias_task_origin{stimtask_type_this{fc}(i)};
                        elseif is_original_bias == 5 % general error level PSE
                            temp = diffBias_true_task{stimtask_type_this{fc}(i)};
                        elseif is_original_bias == 6 % normlize general error level PSE
                            temp = nor_diffBias_true_task{stimtask_type_this{fc}(i)};
                        end
                        
                        temp(~select,:) = nan; % remove undesired unit
                        diffBias_thistask{amp} = [diffBias_thistask{amp};temp]; % combine selected stimtask_type_this
                        
                        
                        if sum(is_original_bias == [0 1 2 3 4])>0 % normal psy
                            % threshold
                            temp2 = diffThresh_task{stimtask_type_this{fc}(i)};
                            temp2(~select,:) = nan;
                            diffThresh_thistask{amp} = [diffThresh_thistask{amp};temp2];
                            
                            % p value
                            temp3 = diffBias_p_task{stimtask_type_this{fc}(i)};
                            temp3(~select,:) = nan;
                            diffBias_p_thistask{amp} = [diffBias_p_thistask{amp};temp3];
                            
                            temp4 = diffThresh_p_task{stimtask_type_this{fc}(i)};
                            temp4(~select,:) = nan;
                            diffThresh_p_thistask{amp} = [diffThresh_p_thistask{amp};temp4];
                            
                        elseif sum(is_original_bias == [5 6])>0 % general error level
                            % threshold
                            temp2 = diffThresh_true_task{stimtask_type_this{fc}(i)};
                            temp2(~select,:) = nan;
                            diffThresh_thistask{amp} = [diffThresh_thistask{amp};temp2];
                            
                            % p value
                            temp3 = diffBias_p_true_task{stimtask_type_this{fc}(i)};
                            temp3(~select,:) = nan;
                            diffBias_p_thistask{amp} = [diffBias_p_thistask{amp};temp3];
                            
                            temp4 = diffThresh_p_true_task{stimtask_type_this{fc}(i)};
                            temp4(~select,:) = nan;
                            diffThresh_p_thistask{amp} = [diffThresh_p_thistask{amp};temp4];
                        end

                        select_save{amp} = [select_save{amp}; select]; % select for stimAmp, stimtask_type_this and stars_type
                    end
                end
                
                S_PSE{amp} = logical(diffBias_p_thistask{amp}<0.05);
                NS_PSE{amp} = logical(diffBias_p_thistask{amp}>=0.05);
                S_Thr{amp} = logical(diffThresh_p_thistask{amp}<0.05);
                NS_Thr{amp} = logical(diffThresh_p_thistask{amp}>=0.05);
                
                % pool bias and threshold together
                diff_thistask{amp,1} = diffBias_thistask{amp};
                diff_thistask{amp,2} = diffThresh_thistask{amp};
                S_unit{amp,1} = S_PSE{amp};
                S_unit{amp,2} = S_Thr{amp};
                NS_unit{amp,1} = NS_PSE{amp};
                NS_unit{amp,2} = NS_Thr{amp};
            end
        end
    end

    function plot_stim_effect_distribution_uA(stimtask_type_this,monkey_this,effect_type,diff_thistask,S_unit,NS_unit,unique_stimAmp_this,is_separate_stimAmp,is_original_bias)
%         plot_stim_effect_distribution_uA(stimtask_type,monkey_included_for_analysis,effect_type,diff_thistask,S_unit,NS_unit,unique_stimAmp_this,is_separate_stimAmp,is_original_bias);
        
        if  nargin == 6
            unique_stimAmp_this = 9999; % combine 20 and 40 uA
            is_separate_stimAmp = 0;
            is_original_bias = 0;  % 0: bias positive = prefer direction; 1: only original bias in MT;  2: original bias(after-before); 3： normalize PSE; 4 normalize original bias
        end
        
        % ------------------ default parameter -----------------------
        pic_name_this = SetTitle(stimtask_type_this,monkey_this);
        pic_name_this = strcat("u-Stim effect in ",pic_name_this);
        
        
        % prepare for plotting
        %         diff_thistask{amp,1} = diffBias_thistask{amp};
        %         diff_thistask{amp,2} = diffThresh_thistask{amp};
        
        bin_num = 8; % 偶数时，0是边界处
        
        % 统一不同stimAmp的 xlimit, xbins, ylimit
        for i = 1:3
            diff_all_stimAmp = [];
            for amp = 1:length(unique_stimAmp_this) % 合并20/40uA的数据
                diff_all_stimAmp = [diff_all_stimAmp;diff_thistask{amp,effect_type}(:,i)];
            end
            diff_all_stimAmp(isnan(diff_all_stimAmp)) = [];
            if ~isempty(diff_all_stimAmp)
                % 四分位距(IQR) 判断outlier
                IQR = iqr(diff_all_stimAmp); % = q3-q1
                q1 = prctile(diff_all_stimAmp,25);
                q3 = prctile(diff_all_stimAmp,75);
                
                outlier_mult = 4;
                upbound = q3+outlier_mult*IQR; % q3+4*IQR 且取整
                if upbound>1
                    upbound = ceil(upbound);
                elseif upbound>0 && upbound<=1
                    upbound = ceil(upbound*10)/10;
                else
                    keyboard
                end
                
                downbound = q1-outlier_mult*IQR; % q1-4*IQR 且取整
                if downbound<-1
                    downbound = floor(downbound);
                elseif downbound<0 && downbound>=-1
                    downbound = ceil(downbound*10)/10;
                else
                    keyboard
                end

                reansonable_range{i} = [downbound upbound]; % 全体数据q1-3*IQR  q3+3*IQR 且取整
                
                outlier_num = sum(diff_all_stimAmp>reansonable_range{i}(2)) + sum(diff_all_stimAmp<reansonable_range{i}(1));

                % have outlier
                if outlier_num>0
                    normal_value = diff_all_stimAmp(diff_all_stimAmp>=reansonable_range{i}(1) & diff_all_stimAmp<=reansonable_range{i}(2)); % 去掉outlier
                else
                    normal_value = diff_all_stimAmp;
                end
                
                % use normal_value to find the bin size and bin edge
                % adjust the normal binedge to "整数”， exclude outlier
                % 为了画图好看，规则：
                % 1. >24: 取30，40，50等
                % 2. >2: 取最接近的3/4的倍数, 例如2.1=3, 2.5=3, 5=6, 7=8
                % 3. >=1 & <=2: 取[1, 1.5, 2]，例如1=1, 1.1=1.5, 1.6=2
                % 4. <1: 取0.9,0.8,0.7...等
                normal_value_max = max(abs(normal_value));
                
                if normal_value_max > 24 % 较大的x轴
                    normal_edge{i} = ceil(normal_value_max/10)*10; % 取30，40，50等
                elseif  normal_value_max >2 % 取最接近的3/4的倍数；目标：使得xtick包括0一共3或4个数字，且每个bin的宽度比较易读
                    normal_edge{i} = near_multiple(normal_value_max,[3 4]);
                elseif  normal_value_max >=1 && normal_value_max <=2 % use 1, 1.5
                    set_max = [1 1.5 2];
                    temp = set_max - normal_value_max;
                    normal_edge{i} = set_max(temp == min(temp(temp>=0))); % 最接近set_max且小于等于set_max的值
                else %<1, for motion type task, typical PSE is very small
                    normal_edge{i} = ceil(normal_value_max*10)/10; % 取0.9,0.8,0.7...等
                end
                
                % 使每个bin的宽度比较易读
                if normal_edge{i}==6 || normal_edge{i}==9 || normal_edge{i}==15 % 12 暂时用10个bin
                    bin_size{i} = normal_edge{i}/9;
                else
                    bin_size{i} = normal_edge{i}/10;
                end

                %                 % original range, include outlier,取整
                %                 ii = 1;
                %                 xlimit_max{i} = 0;
                %                 while xlimit_max{i} < max(abs(diff_all_stimAmp))
                %                     xlimit_max{i} = bin_size{i}*ii;
                %                     ii = ii + 1;
                %                 end
                %                 binedge{i} = -xlimit_max{i}: bin_size{i}: xlimit_max{i}; % original edge

                binedge{i} = -normal_edge{i}: bin_size{i}: normal_edge{i}; % normal data edge
                binmid{i} = binedge{i}(1:end-1) + bin_size{i}/2; % for bar plot
            else
                binedge{i} = [];
            end
            
            nbins = [];
            if ~isempty(diff_all_stimAmp)
                for amp1 = 1:length(unique_stimAmp_this) % 合并分开20/40uA的数据
                    try
                        nbins = [nbins,histcounts(diff_thistask{amp1,effect_type}(:,i),binedge{i})];
                    catch
                        keyboard
                    end
                end
                ylimit_all{i} = max(nbins);
            else
                ylimit_all{i} = [];
            end
        end
        
        % ------------------------ Plot ------------------------------
        set(figure(figN),'name',pic_name_this,'pos',[100,200 1500,400]); clf
        % plotting
        % 1. 不区分电刺激强度 （原来的画法）
        % 2. 区分电刺激强度，40uA颜色，透明度,bar宽度
        Only_plotBar = 0;
        for i = 1:3
            for amp = 1:length(unique_stimAmp_this)
                subplot(1,3,i)
                hold on
                if ~isempty(binedge{i})      
                    % deal with outlier
                    diffthis = diff_thistask{amp,effect_type}(:,i);
                    
                    % first plot the normal data
%                     select_normal = logical(diffthis<=reansonable_range{i}(2) & diffthis>=reansonable_range{i}(1));
                    select_normal = logical(diffthis<=normal_edge{i} & diffthis>=-normal_edge{i});
                    
                    histS = histcounts(diffthis(S_unit{amp,effect_type}(:,i) & select_normal),binedge{i});
                    histNS = histcounts(diffthis(NS_unit{amp,effect_type}(:,i) & select_normal),binedge{i});
                    hbars = bar(binmid{i},[histS' histNS'],1,'stacked','LineWidth',1.5);
                    xlimit_max = max(binedge{i});
                    
                    % then plot the outlier, separate positive outlier and negative outlier
                    select_pos_out = logical(diffthis>normal_edge{i});
                    select_neg_out = logical(diffthis<-normal_edge{i});
                    
                    % cout the num and plot together
                    out_dis_bin = 2.5; % n binsize of outlier away from ori max bar
                    out_dis = (out_dis_bin+1) * bin_size{i}; % distance between outlier bar mid position and ori max bar mid position
                    
                    % positive outlier
                    if sum(select_pos_out)>0
                        histS_out = sum(S_unit{amp,effect_type}(:,i) & select_pos_out);
                        histNS_out = sum(NS_unit{amp,effect_type}(:,i) & select_pos_out);
                        
                        % location for plot outlier
                        out_loc = [max(binmid{i})+out_dis max(binmid{i})+out_dis+bin_size{i}]; % add a fake zero bar for easy to plot stack bar
                        hold on
                        hbars_pos_out = bar(out_loc,[histS_out histNS_out;0 0],1,'stacked','LineWidth',1.5);
                        set(hbars_pos_out,'EdgeColor','k','FaceColor',c{i});
                        set(hbars_pos_out(2),'FaceColor','y');
                    end
                    
                    % negative outlier
                    if sum(select_neg_out)>0
                        histS_out = sum(S_unit{amp,effect_type}(:,i) & select_neg_out);
                        histNS_out = sum(NS_unit{amp,effect_type}(:,i) & select_neg_out);
                        
                        % location for plot outlier
                        out_loc = [min(binmid{i})-out_dis-bin_size{i} min(binmid{i})-out_dis]; % add a fake zero bar for easy to plot stack bar
                        hold on
                        hbars_neg_out = bar(out_loc,[0 0;histS_out histNS_out;],1,'stacked','LineWidth',1.5);
                        set(hbars_neg_out,'EdgeColor','k','FaceColor',c{i});
                        set(hbars_neg_out(2),'FaceColor','y');
                    end
                    
                    if amp == 1
                        set(hbars,'EdgeColor','k','FaceColor',c{i});
%                         set(hbars(2),'FaceColor','none');
                        set(hbars(2),'FaceColor','y'); % for segement
                    else
                        set(hbars,'barwidth',0.6,'EdgeColor','k','FaceColor',1-c{i}); % 显著
                        set(hbars(2),'barwidth',0.6,'FaceColor','k'); % 不显著
                    end
                    
                    % set the xticklabel
                    % auto set xtick, 使得xtick包括0一共3或4个数字
                    if mod(normal_edge{i},3)==0 % 优先分3个tick
                        xtick_this = [-normal_edge{i}:normal_edge{i}/3:normal_edge{i}];
                    elseif mod(normal_edge{i},4)==0 % 其次分4个tick
                        xtick_this = [-normal_edge{i}:normal_edge{i}/4:normal_edge{i}];
                    elseif normal_edge{i} == 2 || normal_edge{i} == 1.5
                        xtick_this = [-normal_edge{i}:0.5:normal_edge{i}];
                    else % <1
                        xtick_this = [-normal_edge{i}:0.1:normal_edge{i}];
                    end
                    
                    if normal_edge{i} == 12 % 分4个tick
                        xtick_this = [-normal_edge{i}:normal_edge{i}/4:normal_edge{i}];
                    end
                    
                    xticks(xtick_this);
                    
%                     % 调整增加outlier bar后的xlim，会缩小有outlier的每个bar的宽度
%                     if sum(select_pos_out)+ sum(select_neg_out) > 0
%                         xlimit_max = max(binedge{i})+out_dis;
%                     end
%                     xlim([-xlimit_max*1.1 xlimit_max*1.1])
                    
                    % 没有outlier时两边也预留位置，但每个bar宽度一样
                    xlimit_max = max(binedge{i})+out_dis;
                    xlim([-xlimit_max xlimit_max]);
                    
                    % set xticklabel
                    if sum(select_pos_out)>0
                        xtick_ori = get(gca,'xTick')';
                        xticklabel_ori = get(gca,'xTicklabel');
                        set(gca,'xtick',[get(gca,'xtick'), max(binmid{i})+out_dis]); % set xtick
                        xticklabel_this = [xticklabel_ori;{['>' num2str(normal_edge{i})]}]; % set xticklabel
                        set(gca,'xticklabel',xticklabel_this);
                    end
                    if sum(select_neg_out)>0
                        xtick_ori = get(gca,'xTick');
                        xticklabel_ori = get(gca,'xTicklabel');
                        set(gca,'xtick',[min(binmid{i})-out_dis, get(gca,'xtick')]);
                        xticklabel_this = [{['<' num2str(-normal_edge{i})]};xticklabel_ori];
                        set(gca,'xticklabel',xticklabel_this);
                    end
                    
                    % auto set yticklabel
                    if ylimit_all{i}>25 % 较大的y轴,十位取3,4,5的倍数
                        ylimit_max = near_multiple(ylimit_all{i}/10,[3 4 5])*10;
                        if mod(ylimit_max/10,4)==0 % 优先分4个tick
                            ytick_this = [0:ylimit_max/4:ylimit_max];
                        elseif mod(ylimit_max/10,3)==0 % 其次分3个tick
                            ytick_this = [0:ylimit_max/3:ylimit_max];
                        elseif mod(ylimit_max/10,5)==0 % 其次分5个tick
                            ytick_this = [0:ylimit_max/5:ylimit_max];
                        end
                    else
                        % set the yticklabel: 4/5倍数
                        ylimit_max = near_multiple(ylimit_all{i},[4 5]);
                        if mod(ylimit_max,4)==0 % 优先分4个tick
                            ytick_this = [0:ylimit_max/4:ylimit_max];
                        elseif mod(ylimit_max,5 )==0 % 其次分5个tick
                            ytick_this = [0:ylimit_max/5:ylimit_max];
                        end
                    end
                    yticks(ytick_this);
                    ylim([0 ylimit_max*1.1]);
                    
                    % for segment
                    set(gca,'XColor','c');
                    set(gca,'YColor','m');
                    
                    hold on
                    plot([0,0],[0,ylimit_all{i}*1.1],'k-','linew',1.5);
                     
                    % mean+sem
                    temp = diff_thistask{amp,effect_type}(S_unit{amp,effect_type}(:,i)|NS_unit{amp,effect_type}(:,i),i); % same stimAmp
                    psy_mean = nanmean(temp);
                    psy_std = nanstd(temp);
                    psy_sem = nanstd(temp) / sqrt(length(temp));
                    
                    psy_median = nanmedian(temp);
                    
                    % Shapiro-Wilk test for normality 正态性检验 N：3~5000
                    [H_norm, p_norm] = SWtest(temp, 0.05); % H0: is composite norality
                    
                    %ttest for significant
                    [~,p_ttest] = ttest(temp);
                    %                     t = stats.tstat;
                    %                     diff_saze{i} = temp;
                    
                    %                     % 正态性不通过，用Wilcoxon signed-rank test 符号秩检验
                    %                     p_signrank = signrank(temp);
                    
                    
                    % 分布不一定对称，不能使用Wilcoxon signed-rank test，改为使用sign test
                    p_signtest = signtest(temp);
                    
                    
                    
                    % only plot the bar
                    if Only_plotBar
                        ylabel([]);
                        xticks([]);
                        yticks([]);
                    end
                    
                    % add Asterisks marker
                    if p_signtest<0.001
                        asterisk = '***';
                    elseif p_signtest<0.01
                        asterisk = '**';
                    elseif p_signtest<0.05
                        asterisk = '*';
                    else
                        asterisk = 'n.s.';
                    end
                    
                    %annotation
                    if amp == 1
                        plot(psy_median,ylimit_max+ylimit_max/30,'vk','linewidth',1,'markerfacecolor',c{i});  %三角标记
                        % plot mean±sem
                        %                         if ~Only_plotBar
                        %                             text(psy_mean,ylimit_max+ylimit_max/7,sprintf('%2.2g ± %2.2g (p(t) = %2.2g, p(signrank) = %2.2g)',psy_mean, psy_sem,p_psy_mean, p_signrank),'fontsize',10,'color',c{i},'HorizontalAlign','center'); % ttest 结果
                        %                         end
                        % plot median
                        text(psy_median+xlimit_max/20,ylimit_max+ylimit_max/15,sprintf('%2.2g (p(signtest) = %2.2g)',psy_median, p_signtest),'fontsize',10,'color',c{i});
                        
                        % t test
                        text(psy_median+xlimit_max/10,ylimit_max-ylimit_max/15,sprintf('p(ttest) = %2.2g)',p_ttest),'fontsize',10,'color',c{i})
                        % add Asterisks marker
                        text(psy_median,ylimit_max+ylimit_max/15,asterisk,'fontsize',15,'color',c{i},'HorizontalAlign','center');
                            
                        
                        
                    else
                        plot(psy_median,ylimit_max+ylimit_max/30,'vk','linewidth',1,'markerfacecolor',1-c{i});  %三角标记
                        %                         if ~Only_plotBar
                        %                             text(psy_mean,ylimit_max+ylimit_max/7,sprintf('%2.2g ± %2.2g (p = %2.2g), p(signrank) = %2.2g)',psy_mean, psy_sem,p_psy_mean,p_signrank),'fontsize',10,'color','k','HorizontalAlign','center'); % ttest 结果
                        %                         end
                        text(psy_median,ylimit_max+ylimit_max/20,sprintf('%2.2g (p(signtest) = %2.2g)',psy_median, p_signtest),'fontsize',10,'color',c{i},'HorizontalAlign','center');
                    end
                    
                    % unit number
                    if ~Only_plotBar
                        if is_separate_stimAmp
                            text(-xlimit_max,ylimit_max-ylimit_max*amp/15,sprintf('%g uA, N = %g',unique_stimAmp_this(amp),sum(S_unit{amp,effect_type}(:,i)|NS_unit{amp,effect_type}(:,i))));
                        else
                            text(-xlimit_max,ylimit_max-ylimit_max*amp/15,sprintf('20 and 40 uA, N = %g',sum(S_unit{amp,effect_type}(:,i)|NS_unit{amp,effect_type}(:,i))));
                        end
                        
%                         ylabel('Cases');
                    end
                end
            end
            if i == 1
                title('Translation');
            elseif i == 2
                title('Roll');
            else
                title('Flow pattern');
            end
        end
        
        SetFigure(10);
        if effect_type == 1
            if is_original_bias == 0 % 0: bias positive = prefer direction; 1: only original bias in MT;  2: original bias(after-before); 3： normalize PSE; 4 normalize original bias
                bigtitle = char(strcat(pic_name_this,"(PSE shift)"));
            elseif is_original_bias == 1
                bigtitle = char(strcat(pic_name_this,"(PSE shift and original bias in MT)"));
            elseif is_original_bias == 2
                bigtitle = char(strcat(pic_name_this,"(original PSE shift)"));
            elseif is_original_bias == 3
                bigtitle = char(strcat(pic_name_this,"(normalized PSE shift)"));
            elseif is_original_bias == 4
                bigtitle = char(strcat(pic_name_this,"(normalized original PSE shift)"));
            elseif is_original_bias == 5
                bigtitle = char(strcat(pic_name_this,"(general error PSE shift)"));
            elseif is_original_bias == 6
                bigtitle = char(strcat(pic_name_this,"(normalized general error PSE shift)"));
            end
        else
            bigtitle = char(strcat(pic_name_this,"(Threshold changde)"));
        end
        if ~Only_plotBar
            suptitle(bigtitle);
        end
        figN = figN+1;
    end
end
