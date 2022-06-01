% population analysis
% Protocol: 1->spiral tuning, 2-> cue switching, 3-> cue preformance, 4-> microstimulation
% protcol = 11 -> forward motion and back motion 201905 ask for Gu
% Lwh 20170919
% Lwh 20210815 等待优化： 不同result的匹对不再需要顺序一样，只要内容一样就行

function data_all = popul_load_data_LWH(FnameCode, Protocol,monkey_choose)

% progress bar, 进度条
if Protocol == 1
    progressbar('Load SpiT .mat files');
elseif Protocol == 2
    progressbar('Load CueS .mat files');
elseif Protocol == 4
    progressbar('Load Stim .mat files');
elseif Protocol == 5
    progressbar('Load AfterStim SpiT.mat files');
end

% deal with two monkey
monkey_this = [];
if length(monkey_choose)==1 % 1 for Ringbell, 2 fro Arthas
    if monkey_choose==1
        monkey_this{1} = 'Ringbell';
    elseif monkey_choose==2
        monkey_this{1} = 'Arthas';
    end
else % both two monkey
    monkey_this{1} = 'Ringbell';
    monkey_this{2} = 'Arthas';   
end

if exist('Z:\BaiduNetdiskWorkspace\data\Batch')
    path_temp = 'Z:\BaiduNetdiskWorkspace\data\Batch\';
elseif exist('Z:\Data\TEMPO\BATCH\')
    path_temp = 'Z:\Data\TEMPO\BATCH\';
end

switch FnameCode
    case 1 % choose files manually
        
    case 2 % load all mat files from interested folder
        switch Protocol
            case 1
                for i = 1:length(monkey_this)
                    pathname{i} = strcat(path_temp,monkey_this{i},'_SpiT');
                    show_temp = 'Spiral Tuning data...';
                end
            case 2
                for i = 1:length(monkey_this)
                    pathname{i} = strcat(path_temp,monkey_this{i},'_CueS');
                    show_temp = 'CueSwitching CP data...';
                end
            case 3
                for i = 1:length(monkey_this)
                    pathname{i} = strcat(path_temp,monkey_this{i},'_CueS_performance');
                    show_temp = 'Performance data...';
                end
            case 4
                for i = 1:length(monkey_this)
                    pathname{i} = strcat(path_temp,monkey_this{i},'_microstimulation');
                    show_temp = 'Microstimulation data...';
                end
            case 5
                for i = 1:length(monkey_this)
                    pathname{i} = strcat(path_temp,monkey_this{i},'_AfterStimSpiT');
                    show_temp = 'AfterStim Tuning data...';
                end
            case 11
                pathname{1} = 'Z:\Data\TEMPO\BATCH\BatchFile_2coh_tuning_MST';
                disp('Load data...');
            case 12
                pathname{1} = 'Z:\Data\Tempo\Batch\motion_parallax_in_tuning (Data from Lihua Yang)\BatchFile_withoutMotionParallaxvis(2D)-sort6';
                disp('Load data...');
            case 13
                pathname{1} = 'Z:\Data\Tempo\Batch\motion_parallax_in_tuning (Data from Lihua Yang)\BatchFile_withMotionParallax_vis(3D)-sort9';
                disp('Load data...');
            case 14
                pathname{1} = 'Z:\Data\Tempo\Batch\data_from_Gu\FisherInformation992allvissig_new_adapted_by_Lwh';
                disp('Load data...');
            case 15
                pathname{1} = 'Z:\Data\Tempo\Batch\FEF_tuning (Data from Gu)\Batch_FEF_visual_Gu';
                disp('Load data...');
        end
        
        % file num in select monkey
        filenum_all = 0;
        for n = 1:length(monkey_choose)
            cd(pathname{n});
            matfile = dir([pathname{n},'\*m*_*.mat']); %所要读取文件的类型
            filenum_all = length(matfile) + filenum_all;
        end
        
        % load .mat file
        file_index = 0;
        if filenum_all > 0
            for n = 1:length(monkey_choose) % 1:sheet1→Ringbell, 2:sheet2→Arthas
                show_this = strcat('Load',32,monkey_this{n},32,show_temp);
                disp(show_this);
                
                cd(pathname{n});
                matfile = dir([pathname{n},'\*m*_*.mat']); %所要读取文件的类型
                
                dircell = struct2cell(matfile);
                ori_name = dircell(1,:);
                temp_name = sort_nat(ori_name);  % 文件名按照自然顺序排序
                
                filenum_this = length(temp_name); % this folder file number
                if filenum_this>0
                    match_index = []; match_index_temp= [];
                    
                    for i = 1:filenum_this % loop for every folder
                        if i > 1
                            fname_data = fieldnames(data_all(1));
                            temp = load(char(temp_name(i)));
                            fname_this = fieldnames(temp.result);
                            data_field_num = max([length(fname_data),length(fname_this)]);
                            
                            % compare data file in previous and this file,
                            % 保持每一个data的结构都是一致的，如果只选择部分文件跑了新的batch，可能会导致与原来的data结构不同，此处报错
                            compare_data_field = cell(data_field_num,3);
                            for f = 1:data_field_num
                                try
                                    compare_data_field{f,1} = fname_data{f};
                                    compare_data_field{f,2} = fname_this{f};
                                    compare_data_field{f,3} = strcmp(fname_data{f}, fname_this{f});
                                catch
                                    disp('************** Data_all Struct do not match ! *****************');
                                    keyboard
                                end
                            end
                            if sum([compare_data_field{:,3}])~=data_field_num
                                disp('************** Data_all Struct do not match ! *****************');
                                disp(compare_data_field);
                                disp('************** Data_all Struct do not match ! *****************');
                                keyboard
                            end
                        end
                        
                        try
                            temp = load(char(temp_name(i)));
                            file_index = file_index + 1;
                            eval(['data_all(file_index) = ','temp.result;']);
                        catch
                            temp_name(i)
                            keyboard
                        end
                        % 计数
                        progressbar(file_index/filenum_all);
                    end
                    %             else
                    %                 data_all = data_all; % if one of the monkeys or all monkeys do not have this task data
                    %                 % 如果第一只猴子没有data，第二只有data，如何处理？？ Lwh 20211231
                end
            end
            
        else % 所有猴子都没有data
            data_all = struct;
        end
        
        
        
        disp('---------------------------------------------------------------------------------------------');
        disp('                                   Pack data ALL SUCCESS!                                    ');
        disp('---------------------------------------------------------------------------------------------');
end
