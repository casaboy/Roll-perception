function function_handles = Group_CueS(loadM_data)
%% Extract some data to global value for easier access
% 将常用的参数提取出来
global group_result;

% ================================================
%    only select eye center cell and spike channel 1
%
% ===============================================


% Add flexible monkey mask here (but I've still decided to choose monkey for analysis below). HH20150723
monkey_included_for_loading = [7 16];

for n = 1:length(loadM_data)
    % CueS data
    main_depth_loc = [];
    main_depth_loc = [(loadM_data(n).Data.Depth_code)]==0;
    select_cues_data = loadM_data(n).Data(main_depth_loc).SpikeChan(1).centerEye_Protocol{2};  % CueS4T
    if isempty(select_cues_data) % main depth 没有做 cues task
        % 找出是否在别的深度也做了 cues task
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
            for m = 1:length(depth_temp)
                select_cues_data = loadM_data(n).Data(depth_ori == depth_temp(m)).SpikeChan(1).centerEye_Protocol{2};
                if ~isempty(select_cues_data)
                    break % 找到有
                end
            end
        end
    end
    
    % window number in partial correlation over time
    timewin = 200; % in ms
    stepwin = 10; % in ms
    binStartOffset = 100; % in ms
    binStopOffset = 100; % in ms
    
    try
        group_result(n).actul_stim_time = select_cues_data.actul_stim_time;
    catch
        group_result(n).actul_stim_time = 9999;
    end
    time_start = 0-binStartOffset;
    time_stop = group_result(n).actul_stim_time+binStopOffset;
    t_centers = time_start + timewin/2 : stepwin : time_stop - timewin/2;
    winnum = length(t_centers);
    
    group_result(n).cellID = loadM_data(n).cellID;
    if isempty(select_cues_data) % 仍然没找到该unit做了spiral tuning
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        winnum = 102;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        group_result(n).CueS_FILE = nan;
        group_result(n).CueS_repN = nan;
        group_result(n).unique_coherence = nan;
        group_result(n).SUMU = nan;
        group_result(n).Psy_correct = {nan(1,7),nan(1,7),nan(1,7)};
        group_result(n).Bias_psy = nan(1,3);
        group_result(n).Thresh_psy = nan(1,3);
        
        group_result(n).Neuro_correct_with0 = {nan(1,7),nan(1,7),nan(1,7)};
        group_result(n).Bias_neu_with0 = nan(1,3);
        group_result(n).Thresh_neu_with0 = nan(1,3);
        
        group_result(n).Neuro_correct_anti = {nan(1,7),nan(1,7),nan(1,7)};
        group_result(n).Bias_neu_anti = nan(1,3);
        group_result(n).Thresh_neu_anti = nan(1,3);
        
        group_result(n).p_local_tun_all_choice = nan(1,3);
        group_result(n).p_local_tun_tasktype_choice = nan(1,3);
        group_result(n).line_re = nan(1,3);
        group_result(n).line_p = nan(1,3);
        group_result(n).local_pref_code = nan(1,3);
        group_result(n).global_pref_code = loadM_data(n).global_pref_code; % 1: 左，CW，rotation；  2: 右，CCW，heading
        
        group_result(n).partial_corr_coef_sen = nan(1,3);
        group_result(n).partial_corr_coef_choice = nan(1,3);
        group_result(n).partial_corr_p_sen = nan(1,3);
        group_result(n).partial_corr_p_choice = nan(1,3);
        group_result(n).p_effect_sen = nan(1,3);
        group_result(n).p_effect_choice = nan(1,3);
        
        group_result(n).partial_corr_coef_time_sen = {nan(1,winnum) nan(1,winnum) nan(1,winnum)};
        group_result(n).partial_corr_coef_time_choice = {nan(1,winnum) nan(1,winnum) nan(1,winnum)};
        group_result(n).partial_corr_p_time_sen = {nan(1,winnum) nan(1,winnum) nan(1,winnum)};
        group_result(n).partial_corr_p_time_choice = {nan(1,winnum) nan(1,winnum) nan(1,winnum)};
        group_result(n).partial_corr_coef_all = nan(1,4);
        group_result(n).partial_corr_p_all = nan(1,4);
        
        group_result(n).partial_corr_coef_all_time = {nan(winnum,4)};
        group_result(n).partial_corr_p_all_time = {nan(winnum,4)};
        
        group_result(n).dead_ahead = nan(1,3);
        group_result(n).CP_allcond = {nan(1,7),nan(1,7),nan(1,7)};
        group_result(n).CP_center = nan(1,3);
        group_result(n).CP_grand = nan(1,3);
        group_result(n).CP_center_p_perm = nan(1,3);
        group_result(n).CP_grand_p_perm = nan(1,3);
        group_result(n).CP_center_p_ttest = nan(1,3);
        group_result(n).CP_grand_p_ttest = nan(1,3);
        
        group_result(n).selection_0 = nan(1,4);
        group_result(n).selection_around_0 = nan(1,4);
        group_result(n).task_selection = {nan(7,1),nan(7,1),nan(7,1)};
    else
        group_result(n).CueS_FILE = select_cues_data.FILE;
        group_result(n).CueS_repN = select_cues_data.repetition;
        group_result(n).repN_zero = select_cues_data.zero_repetition;
        group_result(n).unique_cond = select_cues_data.unique_cond;
        group_result(n).unique_coherence = select_cues_data.unique_coherence_coarse;
        group_result(n).SUMU = select_cues_data.SUMU;
        
        group_result(n).Psy_correct = squeeze(select_cues_data.Psy_correct)';
        group_result(n).Bias_psy = cell2mat(squeeze(select_cues_data.Bias_psy)');
        group_result(n).Thresh_psy = cell2mat(squeeze(select_cues_data.Thresh_psy)');
        
        group_result(n).Neuro_correct_with0 = squeeze(select_cues_data.Neuro_correct_with0)';
        group_result(n).Bias_neu_with0 = cell2mat(squeeze(select_cues_data.Bias_neu_with0)');
        group_result(n).Thresh_neu_with0 = cell2mat(squeeze(select_cues_data.Thresh_neu_with0)');
        
        group_result(n).Neuro_correct_anti = squeeze(select_cues_data.Neuro_correct_anti)';
        group_result(n).Bias_neu_anti = cell2mat(squeeze(select_cues_data.Bias_neu_anti)');
        group_result(n).Thresh_neu_anti = cell2mat(squeeze(select_cues_data.Thresh_neu_anti)');
        
        group_result(n).p_local_tun_all_choice = cell2mat(squeeze(select_cues_data.p_local_tun_all_choice)');
        group_result(n).p_local_tun_tasktype_choice = cell2mat(squeeze(select_cues_data.p_local_tun_tasktype_choice)');
        group_result(n).line_re = cell2mat(squeeze(select_cues_data.line_re)');
        group_result(n).line_p = cell2mat(squeeze(select_cues_data.line_p)');
        group_result(n).local_pref_code = cell2mat(squeeze(select_cues_data.local_pref_code)'); % 1: 左，CW，rotation；  2: 右，CCW，heading
        group_result(n).global_pref_code = loadM_data(n).global_pref_code; % 1: 左，CW，rotation；  2: 右，CCW，heading
        
        partial_corr_coef = cell2mat(squeeze(select_cues_data.partial_corr_coef));
        partial_corr_p = cell2mat(squeeze(select_cues_data.partial_corr_p));
        group_result(n).partial_corr_coef_sen = partial_corr_coef(:,1)';
        group_result(n).partial_corr_coef_choice = partial_corr_coef(:,2)';
        group_result(n).partial_corr_p_sen = partial_corr_p(:,1)';
        group_result(n).partial_corr_p_choice = partial_corr_p(:,2)';
        
        p_effect = cell2mat(squeeze(select_cues_data.p_effect));
        group_result(n).p_effect_sen = p_effect(:,1)';
        group_result(n).p_effect_choice = p_effect(:,2)';
        
        partial_corr_coef_time_temp = squeeze(select_cues_data.partial_corr_coef_time);
        partial_corr_p_time_temp = squeeze(select_cues_data.partial_corr_p_time);
        for m = 1:3 % h r mt
            group_result(n).partial_corr_coef_time_sen{m} = partial_corr_coef_time_temp{m}(:,1)';
            group_result(n).partial_corr_coef_time_choice{m} = partial_corr_coef_time_temp{m}(:,2)';
            group_result(n).partial_corr_p_time_sen{m} = partial_corr_p_time_temp{m}(:,1)';
            group_result(n).partial_corr_p_time_choice{m} = partial_corr_p_time_temp{m}(:,2)';
        end
        
        % heading, rotation, switch_index, choice
        group_result(n).partial_corr_coef_all = cell2mat(squeeze(select_cues_data.partial_corr_coef_all));
        group_result(n).partial_corr_p_all = cell2mat(squeeze(select_cues_data.partial_corr_p_all));
        group_result(n).partial_corr_coef_all_time = select_cues_data.partial_corr_coef_all_time;
        group_result(n).partial_corr_p_all_time = select_cues_data.partial_corr_p_all_time;
        group_result(n).p_effect_all = cell2mat(squeeze(select_cues_data.p_effect_all));
        
        group_result(n).dead_ahead = select_cues_data.i_0;
        group_result(n).CP_allcond = squeeze(select_cues_data.CP_allcond)';
        group_result(n).CP_center = cell2mat(squeeze(select_cues_data.CP_center)');
        group_result(n).CP_grand = cell2mat(squeeze(select_cues_data.CP_grand)');
        group_result(n).CP_center_p_perm = cell2mat(squeeze(select_cues_data.CP_center_p_perm)');
        group_result(n).CP_grand_p_perm = cell2mat(squeeze(select_cues_data.CP_grand_p_perm)');
        group_result(n).CP_center_p_ttest = cell2mat(squeeze(select_cues_data.CP_center_p_ttest)');
        group_result(n).CP_grand_p_ttest = cell2mat(squeeze(select_cues_data.CP_grand_p_ttest)');
        
        % L R CW CCW
        group_result(n).selection_0 = cell2mat(squeeze(select_cues_data.selection_0)');
        group_result(n).selection_around_0 = cell2mat(squeeze(select_cues_data.selection_around_0)');
        group_result(n).task_selection = squeeze(select_cues_data.task_selection)';
    end
end

%% Cell Selection and Cell Counter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
select_all = [];
select_typical = [];
select_no_typical = [];
loop_num = [];
selection_num = [];
t_criterion_txt = [];
cell_selection();

    function cell_selection(typical_cell_selection_num)  % Cell Selection and Cell Counter
        if nargin < 1
            typical_cell_selection_num = 4; % Default,select all cell type
        end
        
        %         if typical_cell_selection_num == 2
        %             loop_num = [1:3];
        %         else
        %             loop_num = typical_cell_selection_num - 2;
        %         end
        
        selection_num = typical_cell_selection_num;
        
        area = {loadM_data.Area}';
        select_all_all_monkey =  [group_result.CueS_repN]' >= 10 & cellfun(@(x) strcmp(x,'MST'),area); % 最少10rep
        
        typical_cell_selection_criteria = ...
            {%  Logic                                   Notes
            % Bottom-line
            'Bottom-line (all)', select_all_all_monkey; % Just all bottom-line cells
            % Condition base
            'Low Coherence (<25 Coh)',select_all_all_monkey & [group_result.unique_coherence]'<25;  % early experience
            'High Coherence (≥25 Coh)',select_all_all_monkey & [group_result.unique_coherence]'>=25; % later experience
            
            % cell type base
            'Tuning Cell',select_all_all_monkey & (([loadM_data.cell_type] == 1) | ([loadM_data.cell_type] == 2) | ([loadM_data.cell_type] == 3))';
            'Translation Cell', select_all_all_monkey & ([loadM_data.cell_type] == 1)';
            'Rotation Cell', select_all_all_monkey & ([loadM_data.cell_type] == 2)';
            'Spiral Cell', select_all_all_monkey & ([loadM_data.cell_type] == 3)';
            'Ohter Cell', select_all_all_monkey & ([loadM_data.cell_type] ~=3)' & ([loadM_data.cell_type] ~=2)' & ([loadM_data.cell_type] ~=1)';
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
    end

% -------- Update/refresh some related datasets that are influenced by cell_selection ---------


%% Final Preparation

% =========== Data for common use ============
function_handles = [];

% ================ Miscellaneous ===================
set(0,'defaultAxesColorOrder',[41 89 204; 248 28 83; 14 153 46]/255);

colors = [41 89 204; 248 28 83; 14 153 46]/255;
c{1} = [248 28 83]/255;
c{2} = [41 89 204]/255;
c{3} = [14 153 46]/255;

% c{1} = 'r';c{2} = 'b';c{3} = 'g';
c{4}='c'; c{5}='w';
c{11} = [0.8 0.5 0.5];
c{12} = [0.5 0.5 0.8];
c{13} = [0.5 0.8 0.5];


transparent = get(findall(gcbf,'tag','transparent'),'value');
p_critical = 0.05;
figN = 200;

%% ====================================== Function Handles =============================================%
% Today I reconstruct the group analysis codes in a more logical and effcient way:
% 1. I divided each figure into small nested functions in Group_HD
% 2. Group_HD returns a cell called "function_handles" back to GROUP_GUI (with all the data available for each nested function)
% 3. Now I can choose which figure(s) to plot or which function to debug in GROUP_GUI
% 4. I feel so happy.
% @HH20150425;

function_handles = {
    'Partial correlation',{
    'Partial correlation',@f1p1;
    'Partial correlation across time',@f1p2;
    'Partial correlation for all condition',@f1p3;
    'Partial correlation across cell',@f1p4;
    }
    'Choice Probability',{
    'CP across task type',@f2p1;
    'Grand CP and neuronal threshold',@f2p2;
    'CP and Partial correlation',@f2p3;
    'Grand CP and Center CP',@f2p4;
    }
    'Other',{
    'Local and Global preferred direction',@f3p1;
    }
    'NoShow',{@cell_selection};
    };

% ---------------------- get data ----------------------------
% ------------------ default parameter -----------------------
% ------------------------ Plot ------------------------------

%% ====================================== Function Definitions =============================================%
    function f1p1(debug) % partial correlation (effect cell)
        if debug  ; dbstack;   keyboard;      end
        % ---------------------- get data ----------------------------
        % two-way anova to see whether heading or choice has a significant effect on FR 后面分析只用至少有一个p<0.05的细胞
        p_effect_sen = cell2mat({group_result.p_effect_sen}');
        p_effect_choice = cell2mat({group_result.p_effect_choice}');
        effect_cell = logical(p_effect_sen < 0.05 | p_effect_choice < 0.05); % Use ANOVA2 p value (Dora)
        effect_cell_num = sum(effect_cell);
        
        partial_corr_coef_sen = cell2mat({group_result.partial_corr_coef_sen}');
        partial_corr_coef_choice = cell2mat({group_result.partial_corr_coef_choice}');
        partial_corr_p_sen = cell2mat({group_result.partial_corr_p_sen}');
        partial_corr_p_choice = cell2mat({group_result.partial_corr_p_choice}');
        
        %         % 取log
        %         type_neg = logical(R_parcor_type<0);
        %         choice_neg = logical(R_parcor_choice<0);
        %         R_parcor_type_log = nan(size(R_parcor_type));
        %         R_parcor_choice_log = nan(size(R_parcor_choice));
        %
        %         R_parcor_type_log(~type_neg) = log10(R_parcor_type(~type_neg));
        %         R_parcor_type_log(type_neg) = -log10(-R_parcor_type(type_neg));
        %
        %         R_parcor_choice_log(~choice_neg) = log10(R_parcor_choice(~choice_neg));
        %         R_parcor_choice_log(choice_neg) = -log10(-R_parcor_choice(choice_neg));
        
        % ------------------ default parameter -----------------------
        methods_of_select = {select_typical, 'Typical cells'};
        set(figure(figN),'name','Partial correlation','pos',[100,200 1700,500]); clf
        
        % ------------------------ Plot ------------------------------
        % Task condition and choice partial correlation
        for i = 1:3 % h r mt
            subplot(1,4,i)
            plot([-1,1],[0,0],'k-');
            hold on
            plot([0,0],[-1,1],'k-');
            
            find_for_partial = [];sen_sig = [];choice_sig = [];
            find_for_partial = methods_of_select{1,1} & effect_cell(:,i);
            sen_sig = partial_corr_p_sen(:,i)<0.05;
            choice_sig = partial_corr_p_choice(:,i)<0.05;
            
            temp_x = partial_corr_coef_sen(find_for_partial,i);
            temp_y = partial_corr_coef_choice(find_for_partial,i);
            
            temp_x_all_sig = partial_corr_coef_sen(find_for_partial & sen_sig & choice_sig,i);
            temp_y_all_sig = partial_corr_coef_choice(find_for_partial & sen_sig & choice_sig,i);
            
            temp_x_one_sig = partial_corr_coef_sen(find_for_partial & (sen_sig | choice_sig) & ~(sen_sig & choice_sig),i);
            temp_y_one_sig = partial_corr_coef_choice(find_for_partial & (sen_sig | choice_sig) & ~(sen_sig & choice_sig),i);
            
            %             plot(temp_x,temp_y,'o','color',c{i},'markersize',9);
            plot(temp_x_one_sig,temp_y_one_sig,'o','color',c{i},'markerfacecolor',c{i+10},'markersize',9); % "One sig" cells
            plot(temp_x_all_sig,temp_y_all_sig,'o','color',c{i},'markerfacecolor',c{i},'markersize',9); % "All sig" cells
            
            axis equal;axis([-1 1 -1 1]);
            set(gca,'xtick',[-1:0.5:1],'xticklabel',[-1:0.5:1]);
            set(gca,'ytick',[-1:0.5:1],'yticklabel',[-1:0.5:1]);
            ylabel('Choice partial corr. (R)');
            if i==1
                xlabel('Translation partial corr. (R)');
                title('Translation and choice partial correlation');
            elseif i==2
                xlabel('Rotation partial corr. (R)');
                title('Rotation and choice partial correlation');
            else
                xlabel('Motion type partial corr. (R)');
                title('Motion type partical correlation');
            end
            
            % Show individual cell selected from the figure. HH20150424
            select_actual_plot = logical(find_for_partial & (sen_sig | choice_sig));
            xx = partial_corr_coef_sen(select_actual_plot,i);
            yy = partial_corr_coef_choice(select_actual_plot,i);
            h_all = plot(xx,yy,'visible','off');hold on;
            
            text(0.4, 0.8, sprintf('N = %0.5g',sum(select_actual_plot))); % show actual plot number
            
            % Draw line if significant
            [~,ppp] = corr(xx,yy,'type','Pearson');
            if ppp < 0.05
                coeff = pca([xx yy]);
                linPara(1) = coeff(2) / coeff(1);
                linPara(2) = mean(yy)- linPara(1) *mean(xx);
                
                % -- Plotting
                xxx = linspace(min(xx),max(xx),150);
                yyy = linPara(1) * xxx + linPara(2);
                plot(xxx,yyy,'k--','linew',2);
            end
            
            set([gca h_all],'ButtonDownFcn',{@Show_individual_cell, h_all, select_actual_plot});
        end
        
        subplot(1,4,4)
        for i = 1:3
            num_corr = sum(find_for_partial);
            both_corr_sig(i) = sum(partial_corr_p_sen(find_for_partial,i)<0.05 & partial_corr_p_choice(find_for_partial,i)<0.05) / num_corr ;
            sen_corr_sig(i) = sum(partial_corr_p_sen(find_for_partial,i)<0.05 & partial_corr_p_choice(find_for_partial,i)>0.05) / num_corr ;
            choice_corr_sig(i) = sum(partial_corr_p_sen(find_for_partial,i)>0.05 & partial_corr_p_choice(find_for_partial,i)<0.05) / num_corr ;
        end
        h3 = bar([both_corr_sig;sen_corr_sig;choice_corr_sig]','stacked');
        set(gca,'xticklabel',{'Translation','Rotation','MT'});
        ylabel('Significant correlation proportion');
        xlim([0.5 3.5]);
        ylim([0 1]);
        legend('Both','Sensory','Choice','location','Northeast');
        box off
        
        suptitle(sprintf('Partial correlation, any signigicant out of  %s (N = %g),',...
            t_criterion_txt,sum(find_for_partial)));
        
        SetFigure(15);figN = figN+1;
    end

    function f1p2(debug) % Partial correlation across time (all cell)
        %         (Zaidel 2017, Figure 6)  HH20180821
        if debug  ; dbstack;   keyboard;      end
        % ---------------------- get data ----------------------------
        % 组合 type 和 choice 方便后面处理, 取平方
        for nn = 1:length(group_result)
            for i = 1:3 % h r mt
                partial_corr_coef_time_flipped{1,i}(nn,:) = group_result(nn).partial_corr_coef_time_sen{i}.^2; % sen
                partial_corr_coef_time_flipped{2,i}(nn,:) = group_result(nn).partial_corr_coef_time_choice{i}.^2; % choice
                
                 partial_corr_coef_time{1,i}(nn,:) = group_result(nn).partial_corr_coef_time_sen{i}; % sen
                partial_corr_coef_time{2,i}(nn,:) = group_result(nn).partial_corr_coef_time_choice{i}; % choice
            end
        end
        
        actul_stim_time = cell2mat({group_result.actul_stim_time})';
        actul_stim_time(actul_stim_time==9999) = []; % 去掉没有cue task的cell
        unique_actul_stim_time = unique(actul_stim_time);
        if length(unique_actul_stim_time)>1
            keyboard
        end
        
%         heading_500ms_sen = mean(partial_corr_coef_time{1,1}(:,1:51),2);
%         heading_500ms_choice = mean(partial_corr_coef_time{2,1}(:,1:51),2);
%         
%         figure
%         plot(heading_500ms_sen,heading_500ms_choice,'ko');
%         [r_square,ppp,h] = plot_corr_line(heading_500ms_sen,heading_500ms_choice)

        % ------------------ default parameter -----------------------
        set(figure(figN),'name','Partial correlations as a function of time','pos',[100,200, 1500,750]); clf
        methods_of_select = { select_typical, 'Typical cells'};
        
        % window number in partial correlation over time
        timewin = 200; % in ms
        stepwin = 10; % in ms
        binStartOffset = 100; % in ms
        binStopOffset = 100; % in ms
        
        time_start = 0 - binStartOffset; % visual on: 0 ms
        time_stop = unique_actul_stim_time + binStopOffset; % visual off: (actul_stim_time) ms
        t_centers = time_start + timewin/2 : stepwin : time_stop - timewin/2;
        winnum = length(t_centers);
        %         smoothFactor = 50;
        
        % ------------------------ Plot ------------------------------
        for j = 1:2 %type correlation, choice correlation
            for i = 1:3 % heading, rotation, moting type
                % 每个time window的 R mean, R err
                select = ~isnan(partial_corr_coef_time_flipped{j,i}(:,1)) & methods_of_select{1,1};
                aver_this = nanmean(partial_corr_coef_time_flipped{j,i}(select,:));
                sem_this = nanstd(partial_corr_coef_time_flipped{j,i}(select,:)) / sqrt(sum(select));
                
                subplot(2,3,i)
                hold on
                if j == 1 % sen
                    shadedErrorBar(t_centers,aver_this,sem_this,{'color',c{i},'markerfacecolor','r','linew',2});
                    %                     ylimit(i) = max(aver_this+sem_this)/6+max(aver_this+sem_this);
                    ylimit(i) = 0.07;
                    ylim([0 ylimit(i)]);
                else
                    shadedErrorBar(t_centers,aver_this,sem_this,{'k','markerfacecolor','k','linew',2});
                end
                plot([0 0],[0 ylimit(i)],'k--','linew',1.5); % visual on and off marker
                plot([unique_actul_stim_time unique_actul_stim_time],[0 ylimit(i)],'k--','linew',1.5);
                xlim([-50 unique_actul_stim_time+50]);
                xlabel('Time[ms]');
                if i==1
                    ylabel('Avg. Partial Corr. [R^2]');
                    title('Translation task');
                elseif i==2
                    title('Rotation task');
                else
                    title('Motion type');
                end
                text(50, ylimit(i)-ylimit(i)/10, sprintf('N = %g',sum(select))); % show actual plot number
                
                % Pool all sensory and choice partial correlation
                subplot(2,3,j+3)
                hold on
                shadedErrorBar(t_centers,aver_this,sem_this,{'color',c{i},'markerfacecolor','r','linew',2});
                plot([0 0],[0 ylimit(i)],'k--','linew',1.5); % visual on and off marker
                plot([unique_actul_stim_time unique_actul_stim_time],[0 ylimit(i)],'k--','linew',1.5);
                xlim([-50 unique_actul_stim_time+50]);
                xlabel('Time[ms]');
                
                if j==1
                    ylabel('Avg. Partial Corr. [R^2]');
                    ylim([0 0.07]);
                    title('Sensory partial R^2');
                elseif j == 2
                    ylim([0 0.025]);
                    title('Choice partial R^2');
                end
            end
        end
        
        
        suptitle(sprintf('Partial correlation over time, %s',t_criterion_txt));
        
        SetFigure(15);figN = figN+1;
    end

    function f1p3(debug) % Partial correlation for all condition
        if debug  ; dbstack;   keyboard;      end
        % ---------------------- get data ----------------------------
        % two-way anova to see whether heading or choice has a significant effect on FR 后面分析只用至少有一个p<0.05的细胞
        p_effect_all = cell2mat({group_result.p_effect_all}');
        effect_cell = logical(p_effect_all(:,1) < 0.05 | p_effect_all(:,2) < 0.05 | p_effect_all(:,3) < 0.05 | p_effect_all(:,4) < 0.05);
        
        % heading, rotation, switch_index, choice
        partial_corr_coef_all = cell2mat({group_result.partial_corr_coef_all}');
        partial_corr_p_all = cell2mat({group_result.partial_corr_p_all}');
        
        for nn = 1:length(group_result)
            for i = 1:4 % heading, rotation, switch index, choice
                partial_corr_coef_all_time{i}(nn,:) = group_result(nn).partial_corr_coef_all_time{1}(:,i).^2;
                partial_corr_p_all_time{i}(nn,:) = group_result(nn).partial_corr_p_all_time{1}(:,i).^2;
            end
        end
        
        actul_stim_time = cell2mat({group_result.actul_stim_time})';
        actul_stim_time(actul_stim_time==9999) = []; % 去掉没有cue task的cell
        unique_actul_stim_time = unique(actul_stim_time);
        if length(unique_actul_stim_time)>1
            keyboard
        end
        
        % ------------------ default parameter -----------------------
        set(figure(figN),'name','Partial correlations for all condition','pos',[100,200, 1500,750]); clf
        methods_of_select = { select_typical, 'Typical cells'};
        
        % window number in partial correlation over time
        timewin = 200; % in ms
        stepwin = 10; % in ms
        binStartOffset = 100; % in ms
        binStopOffset = 100; % in ms
        
        time_start = 0 - binStartOffset; % visual on: 0 ms
        time_stop = unique_actul_stim_time + binStopOffset; % visual off: (actul_stim_time) ms
        t_centers = time_start + timewin/2 : stepwin : time_stop - timewin/2;
        winnum = length(t_centers);
        
        % ------------------------ Plot ------------------------------
        % across time
        for i = 1:4
            select = ~isnan(partial_corr_coef_all_time{i}(:,1)) & methods_of_select{1,1};
            aver_this = nanmean(partial_corr_coef_all_time{i}(select,:));
            sem_this = nanstd(partial_corr_coef_all_time{i}(select,:)) / sqrt(sum(select));
            
            hold on
            h(i) = shadedErrorBar(t_centers,aver_this,sem_this,{'color',c{i},'markerfacecolor','r','linew',2});
        end
        legend([h(1).mainLine h(2).mainLine h(3).mainLine h(4).mainLine],'Translation degree','Rotation degree','Switch Index','Choice');
        
        SetFigure(15);figN = figN+1;
    end

    function f1p4(debug) % Partial correlation across cell
        if debug  ; dbstack;   keyboard;      end
        % ---------------------- get data ----------------------------
        % 对于每一个细胞都有：
        % 控制 choice，heading 与FR
        % 控制 choice，rotation 与FR
        % 控制 choice，motion type 与FR
        % 控制 heading，choice_heading 与FR
        % 控制 rotation，choie_rotation 与FR
        % 控制 motion tpye，choice_mt 与FR
        % cell-by-cell暂时只看两种sensory之间的和两种choice之间的
        
        % two-way anova to see whether heading or choice has a significant effect on FR 后面分析只用至少有一个p<0.05的细胞
        p_effect_sen = cell2mat({group_result.p_effect_sen}');
        p_effect_choice = cell2mat({group_result.p_effect_choice}');
        effect_cell = logical(p_effect_sen < 0.05 | p_effect_choice < 0.05); % Use ANOVA2 p value (Dora)
        effect_cell_num = sum(effect_cell);
        
        partial_corr_coef_sen = cell2mat({group_result.partial_corr_coef_sen}');
        partial_corr_coef_choice = cell2mat({group_result.partial_corr_coef_choice}');
        
        R_parcor{1,1} = partial_corr_coef_sen(:,1); % heading sensory
        R_parcor{1,2} = partial_corr_coef_sen(:,2); % rotation sensory
        R_parcor{2,1} = partial_corr_coef_choice(:,1); % heading choice
        R_parcor{2,2} = partial_corr_coef_choice(:,2); % heading choice
        
        % ------------------ default parameter -----------------------
        set(figure(figN),'pos',[100,100, 850,800], 'Name', 'Comparision of Translation task and Rotation task in partial correlation'); clf
        methods_of_select = { select_typical, 'Typical cells'};
        
        % ------------------------ Plot ------------------------------
        for i = 1:2 % h or r
            subplot(2,2,i)
            plot(R_parcor{i,1}(methods_of_select{1,1}),R_parcor{i,2}(methods_of_select{1,1}),'ko'); % 两种sensory或两种choice比较{i,1},{i,2}
            hold on
            xlim([-1 1]);
            ylim([-1 1]);
            if i==1
                xlabel('r-translation');
                ylabel('r-rotation');
                title('r-sensory cell-by-cell');
            else
                xlabel('r-choice_{t}');
                ylabel('r-choice_{r}');
                title('r-choice cell-by cell');
            end
        end
        
        % 将坐标轴取绝对值
        for i = 1:2
            subplot(2,2,i+2)
            plot(abs(R_parcor{i,1}(methods_of_select{1,1})),abs(R_parcor{i,2}(methods_of_select{1,1})),'ko'); % 两种sensory或两种choice比较{i,1},{i,2}
            
            if i==1
                xlim([0 1]);
                ylim([0 1]);
                xlabel('r-translation');
                ylabel('r-rotation');
                title('|r-sensory|');
            else
                xlim([0 0.5]);
                ylim([0 0.5]);
                xlabel('r-choice_{t}');
                ylabel('r-choice_{r}');
                title('|r-choice|');
            end
        end
        
        suptitle(sprintf('Partial correlation, %s',t_criterion_txt));
        SetFigure(15);figN = figN+1;
    end

    function f2p1(debug) % CP across task type
        if debug  ; dbstack;   keyboard;      end
        % ---------------------- get data ----------------------------
        Grand_CP = cell2mat({group_result.CP_grand}');
        Center_CP = cell2mat({group_result.CP_center}');
        Grand_CP_p = cell2mat({group_result.CP_grand_p_perm}');
        Center_CP_p = cell2mat({group_result.CP_center_p_perm}');
        
        % prefer direction
        global_pref_code = cell2mat({group_result.global_pref_code}');
        local_pref_code = cell2mat({group_result.local_pref_code}');
        
        % whether tuning
        p_stim = cell2mat({group_result.p_stim}');
        global_pref_p = cell2mat({group_result.global_pref_p}');
        local_pref_p = cell2mat({group_result.p_local_tun_tasktype_choice}');
        
        % Translation and rotation CP sign according to local tuning
        % motion type CP sign according to global tuning (有些细胞并没有global tuning，暂时使用local tuning)
        % re-sign mt CP
        Grand_CP(local_pref_code(:,3) == 2,3) = 1 -  Grand_CP(local_pref_code(:,3) == 2,3);
        Center_CP(local_pref_code(:,3) == 2,3) = 1 -  Center_CP(local_pref_code(:,3) == 2,3);
        
        mt_pref_code = global_pref_code(:,3);
        mt_pref_code(isnan(mt_pref_code)) = local_pref_code(isnan(mt_pref_code));
        Grand_CP(mt_pref_code==2,3) = 1 - Grand_CP(mt_pref_code==2,3);
        Center_CP(mt_pref_code==2,3) = 1 - Center_CP(mt_pref_code==2,3);
        
        % group together
        CP{1} = Grand_CP;
        CP{2} = Center_CP;
        CP_p{1} = Grand_CP_p;
        CP_p{2} = Center_CP_p;
        
        % SUMU
        SUMU = cell2mat({group_result.SUMU}'); % SU=1; MU=2
%         SUMU_mask = logical(SUMU==1);
        SUMU_mask = logical(SUMU==1 | SUMU==2);
        
        % ------------------ default parameter -----------------------
        set(figure(figN),'pos',[100,100,950,700], 'Name', 'CP across task type'); clf
        methods_of_select = { select_typical, 'Typical cells'};
        xbins = [];
        xbins = linspace(0,1,30);
        
        % ------------------------ Plot ------------------------------
        for j = 1:2 % grand, center CP
            for i = 1:3 % task type: h r mt
                subplot(2,3,i+(j-1)*3)
                % 只包含有当前motion type的global tuning的细胞，并且以local tuning判断正负
                temp = CP{j}(p_stim<0.05 & global_pref_p(:,i)<0.05 & methods_of_select{1,1} & SUMU_mask,i);
                temp_sig = CP{j}(p_stim<0.05 & global_pref_p(:,i)<0.05 & CP_p{j}(:,i)<0.05 & methods_of_select{1,1} & SUMU_mask,i);
                cpN = hist(temp, xbins);   % for all cell type
                cp_sigN = hist(temp_sig,xbins);  % for p<0.05
                
                h1 = bar(xbins,cpN',1,'facecolor','w','edgecolor',c{i});
                hold on;
                h2 = bar(xbins,cp_sigN',1,'facecolor',c{i},'edgecolor',c{i});
                alpha(.5);
                set(gca,'xlim',[0 1]);
                set(gca,'XTickLabelMode','auto');
                set(gca,'YTickLabelMode','auto');
                ylabel('Number of neurons');
                
                if i==2
                    if j==1
                        xlabel('Grand CP');
                    else
                        xlabel('Center CP');
                    end
                end
                
                %                 ylim([0 30]);
                ymax = max(cpN);
                ylimit = ymax+ymax/10;
                ylim([0 ylimit])
                
                [~,p_CP] = ttest(CP{j}(p_stim<0.05 & global_pref_p(:,i)<0.05 & methods_of_select{1,1} & SUMU_mask,i),0.5); % 与0.5做ttest
                mean_cp = nanmean(CP{j}(p_stim<0.05 & global_pref_p(:,i)<0.05 & methods_of_select{1,1} & SUMU_mask,i)); % 均值
                num_sig_larger = sum(CP{j}(p_stim<0.05 & global_pref_p(:,i)<0.05 & CP_p{j}(:,i)<0.05 & methods_of_select{1,1} & SUMU_mask,i) > 0.5); % CP显著且大于0.5的数量
                
                %annotation
                plot(mean_cp,ylimit,'vk','linewidth',1,'markerfacecolor',c{i});  %三角标记
                text(0.6,ymax,sprintf('%2.2g (p = %2.2g)',mean_cp,p_CP),'fontsize',12,'color',c{i});  % ttest 结果
                text(0.05,ymax,sprintf('N = %g',sum(cpN)),'fontsize',5,'color',c{i});  % ttest 结果
            end
        end
        SetFigure(13);figN = figN+1;
    end

    function f2p2(debug) % 'Grand CP and neuronal threshold'
        if debug  ; dbstack;   keyboard;      end
        % ---------------------- get data ----------------------------
        % prefer direction
        global_pref_p = cell2mat({group_result.global_pref_p}');
        p_stim = cell2mat({group_result.p_stim}');
        global_pref_code = cell2mat({group_result.global_pref_code}');
        
        % Grand CP
        local_pref_code = cell2mat({group_result.local_pref_code}');
        mt_pref_code = global_pref_code(:,3);
        mt_pref_code(isnan(mt_pref_code)) = local_pref_code(isnan(mt_pref_code));
        
        Grand_CP = cell2mat({group_result.CP_grand}');
        Grand_CP_p = cell2mat({group_result.CP_grand_p_perm}');
        Grand_CP(local_pref_code(:,3) == 2,3) = 1 -  Grand_CP(local_pref_code(:,3) == 2,3);
        Grand_CP(mt_pref_code==2,3) = 1 - Grand_CP(mt_pref_code==2,3);
        
        % neuronal threshold
        Thresh_neu_anti = cell2mat({group_result.Thresh_neu_anti}');
        
        % ------------------ default parameter -----------------------
        set(figure(figN),'pos',[100,100,1500,450], 'Name', 'Grand CP and neuronal threshold'); clf
        methods_of_select = { select_typical, 'Typical cells'};
        
        % ------------------------ Plot ------------------------------
        for i = 1:3
            temp = logical(p_stim<0.05 & global_pref_p(:,i)<0.05 & methods_of_select{1,1});
            temp_sig = logical(p_stim<0.05 & global_pref_p(:,i)<0.05 & methods_of_select{1,1} & Grand_CP_p(:,i)<0.05);
            
            subplot(1,3,i)
            plot(Thresh_neu_anti(temp,i),Grand_CP(temp,i),'o','color',c{i},'markersize',10);
            hold on
            plot(Thresh_neu_anti(temp_sig,i),Grand_CP(temp_sig,i),'o','color',c{i},'markerfacecolor',c{i},'markersize',10);
            xlimit = xlim;
            xlimit = xlimit(2);
            plot([0 xlimit],[0.5 0.5],'k--');
            ylim([0 1]);
            xlabel('Neuronal Threshold');
            ylabel('Grand CP');
            
            % Draw line if significant
            xx = Thresh_neu_anti(temp,i);
            yy = Grand_CP(temp,i);
            nanCP = logical(isnan(xx)|isnan(yy)); % 去掉没有CP的site
            xx = xx(~nanCP);
            yy = yy(~nanCP);
            
            [~,ppp] = corr(xx,yy,'type','Pearson');
            if ppp < 0.05
                coeff = pca([xx yy]);
                linPara(1) = coeff(2) / coeff(1);
                linPara(2) = mean(yy)- linPara(1) *mean(xx);
                
                % -- Plotting
                xxx = linspace(min(xx),max(xx),150);
                yyy = linPara(1) * xxx + linPara(2);
                plot(xxx,yyy,'k--','linew',2);
            end
            
        end
        
        SetFigure(15);figN = figN+1;
    end

    function f2p3(debug) % Grand CP and Partial correlation
        if debug  ; dbstack;   keyboard;      end
        % ---------------------- get data ----------------------------
        % CP
        local_pref_code = cell2mat({group_result.local_pref_code}');
        global_pref_code = cell2mat({group_result.global_pref_code}');
        mt_pref_code = global_pref_code(:,3);
        mt_pref_code(isnan(mt_pref_code)) = local_pref_code(isnan(mt_pref_code));
        Grand_CP = cell2mat({group_result.CP_grand}');
        Grand_CP_p = cell2mat({group_result.CP_grand_p_perm}');
        Grand_CP(local_pref_code(:,3) == 2,3) = 1 -  Grand_CP(local_pref_code(:,3) == 2,3);
        Grand_CP(mt_pref_code==2,3) = 1 - Grand_CP(mt_pref_code==2,3);
        
        % Partical correlation
        p_effect_sen = cell2mat({group_result.p_effect_sen}');
        p_effect_choice = cell2mat({group_result.p_effect_choice}');
        effect_cell = logical(p_effect_sen < 0.05 | p_effect_choice < 0.05); % Use ANOVA2 p value (Dora)
        effect_cell_num = sum(effect_cell);
        
        partial_corr_coef_sen = cell2mat({group_result.partial_corr_coef_sen}');
        partial_corr_coef_choice = cell2mat({group_result.partial_corr_coef_choice}');
        partial_corr_p_sen = cell2mat({group_result.partial_corr_p_sen}');
        partial_corr_p_choice = cell2mat({group_result.partial_corr_p_choice}');
        
%         
%         find_for_partial = [];sen_sig = [];choice_sig = [];
%         find_for_partial = methods_of_select{1,1} & effect_cell(:,i);
%         sen_sig = partial_corr_p_sen(:,i)<0.05;
%         choice_sig = partial_corr_p_choice(:,i)<0.05;
%         
%         temp_x = partial_corr_coef_sen(find_for_partial,i);
%         temp_y = partial_corr_coef_choice(find_for_partial,i);
%         
%         temp_x_all_sig = partial_corr_coef_sen(find_for_partial & sen_sig & choice_sig,i);
%         temp_y_all_sig = partial_corr_coef_choice(find_for_partial & sen_sig & choice_sig,i);
%         
%         temp_x_one_sig = partial_corr_coef_sen(find_for_partial & (sen_sig | choice_sig) & ~(sen_sig & choice_sig),i);
%         temp_y_one_sig = partial_corr_coef_choice(find_for_partial & (sen_sig | choice_sig) & ~(sen_sig & choice_sig),i);
        
        keyboard;
        
        
        % ------------------ default parameter -----------------------
        set(figure(figN),'pos',[100,100,1500,450], 'Name', 'Grand CP and neuronal threshold'); clf
        methods_of_select = { select_typical, 'Typical cells'};
        
        % ------------------------ Plot ------------------------------
    end

    function f2p4(debug) %'Grand CP and Center CP'
        if debug  ; dbstack;   keyboard;      end
        % ---------------------- get data ----------------------------
        % ------------------ default parameter -----------------------
        % ------------------------ Plot ------------------------------
    end

    function f3p1(debug) %'Local and Global preferred direction'
        if debug  ; dbstack;   keyboard;      end
        % 比较45度内的global prefer direction
        % ---------------------- get data ----------------------------
        
        % ------------------ default parameter -----------------------
        set(figure(figN),'pos',[100,100,1500,450], 'Name', 'Local and Global preferred direction'); clf
        methods_of_select = { select_typical, 'Typical cells'};
        % ------------------------ Plot ------------------------------
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
        path_temp = 'Z:\Data\Tempo\Batch\Arthas_CueS';
        
        open_name_temp = group_result(ori_cell_no).CueS_FILE;
        dot_loc = find(open_name_temp=='.'); % 找到.的位置
        open_name_temp = open_name_temp(1:dot_loc-1);
        
        open_name = strcat(path_temp,'\',open_name_temp,'_5_CueS_fig_10.bmp');
        figure(9999);clf;
        i = imread(open_name);
        i = imrotate(i,0);
        imshow(i,'InitialMagnification','fit');
        guidata(gcbo,h_marker);
    end

end

