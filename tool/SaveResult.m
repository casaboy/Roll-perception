function SaveResult(config, result, IF_FIGURE, IF_EXCEL)
% Output text for xls; save figures, .mat and .dat file (BATCH only)
% HH20141124
% Lwh201709
% Lwh20210401 excel control

if nargin==3
    IF_EXCEL = true;
end
    
persistent XlsData; % Hold some xls data to prevent read xls file repeatly during batch processing. HH20150724

%%%%%%%%%%%%%%%%%%%%%  Output   HH20140510 / HH20140621 / HH20141003 /HH20141124 %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  remove test option Lwh 202007    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEMPO_GUI: config.batch_flag = []; % original is 'test'
% BATCH_GUI: config.batch_flag = 'Z:\Data\Tempo\Batch\201709014_Ringbell_SpiT.m';  % SpiT task

if ~isempty(config.batch_flag)  % Figures and raw data (always in "result" structure)
    outpath = [config.batch_flag(1:end-2) '\'];
    
    % Check directory
    if ~exist(outpath,'dir')
        mkdir(outpath);
    end
    
    l = length(result.FILE);
    if (result.FILE(l-3:l) == '.htb')	% .htb extension already there
        if ~isempty(result.SpikeChan)
            savefilename = [outpath [result.FILE(1:end-4) '_' num2str(result.SpikeChan)] '_' config.suffix];
        else
            savefilename = [outpath [result.FILE(1:end-4)] '_' config.suffix];
        end
        
    else
        if ~isempty(result.SpikeChan)
            savefilename = [outpath [result.FILE(1:end) '_' num2str(result.SpikeChan)] '_' config.suffix];
        else
            savefilename = [outpath [result.FILE(1:end)] '_' config.suffix];
        end
    end
    
    % Overwrite or Append
    if exist([savefilename '.mat'],'file')
        if isfield(config,'append') && config.append  % Load already existing fields before deletion
            temp = load([savefilename '.mat']);
            result_old = temp.result;
            result_new = result;
            
            % Merge two structures
            M = [fieldnames(result_new)' fieldnames(result_old)' ; struct2cell(result_new)' struct2cell(result_old)' ];
            [~, rows] = unique(M(1,:), 'stable'); % Handle Conflicting ('stable' selects the first occurrence, i.e., result_new)
            M = M(:, rows);
            %              result = struct(M{:}); % This is wrong because we have cells in some fields...
            for i = 1:size(M,2)
                result.(M{i*2-1})= M{i*2};
            end
            disp('Appending fields in .mat file...');
        else
            delete([savefilename '*.*']);
            disp('Overwrite the .mat file...');
        end
    end
    
    % Save raw data
    save([savefilename '.mat'],'result');
    disp('Saving to .mat finished...');
   
    % Save figures
    if IF_FIGURE
    for ff = 1:length(config.save_figures)
        config.save_figures(ff).PaperPositionMode  = 'auto';
        orient(config.save_figures(ff),'portrait'); % 
        set(config.save_figures(ff),'Visible','on');
        %         print(config.save_figures(ff),'-dbitmap',[savefilename '_fig_' num2str(config.save_figures(ff)) '.bmp']);
        saveas(config.save_figures(ff),[savefilename '_fig_' num2str(config.save_figures(ff).Number),'.bmp']);
%         if ~strcmp(config.batch_flag,'test.m')  % HH20160415
%             close(config.save_figures(ff));
%         end
    end
    end
end


% Print part of data to texts (clipboard / .dat file / Data hub "Result.xlsm")
if IF_EXCEL
    if exist('Z:\Data\MOOG\Results\Result_LJY.xlsm')
        excel_path = 'Z:\Data\MOOG\Results\Result_LJY.xlsm';
    elseif exist('Z:\Data\MOOG\Results\Result_Lwh.xlsm')
        excel_path = 'Z:\Data\MOOG\Results\Result_Lwh.xlsm';
    end
    
% try
if ~isempty(config.batch_flag)
    % Print to file
    outfile = [outpath config.suffix '.dat'];
    printHead = 0;
    if (exist(outfile, 'file') == 0)   % file does not yet exist
        printHead = 1;
    end
    
    fid = fopen(outfile, 'a');
    % This line controls the output format
    
    if (printHead)
        fprintf(fid, ['FILE\t ' config.sprint_once_contents '|\t']);
        
        for ll = 1:length(config.sprint_loop_contents)
            fprintf(fid,[config.sprint_loop_contents{ll} '|\t']);
        end
        fprintf(fid, '\r\n');
    end
    
    fprintf(fid,'%s\t',[result.FILE '_' num2str(result.SpikeChan)]);
else  % Print to screen
    fid = 1;
end

toClip = [];
%     toXls = {};

% Print once
sprint_once_marker = [];
for i = 1:length(config.sprint_once_marker)
    sprint_once_marker = [sprint_once_marker '%' config.sprint_once_marker(i) '\t'];
end

% Print once
if ~isempty(config.sprint_once_marker)
    eval(['buff = sprintf(sprint_once_marker,' config.sprint_once_contents ');']);
    %         eval(['buff_toXls = {' config.sprint_once_contents '};' ]); % HH20150724
    
    fprintf(fid, '%s', buff);
    toClip = [toClip sprintf('%s', buff)];
end

% Print loops
if strcmpi(config.suffix,'CueS')
    for n = 1:length(result.unique_condition) % normal �� microstimulation �ֿ�
        if n == 1 % normal
            loop = config.sprint_loop_contents;
            loop_marker = config.sprint_loop_marker;
        elseif n == 2 % microsti only output psy function parameters
            loop = config.sprint_loop_contents(1:2);
            loop_marker = config.sprint_loop_marker(1:2);
        end
        
        for ll = 1:length(loop)
            sprint_loop_marker = [];
            for i = 1:length(loop_marker{ll})
                sprint_loop_marker = [sprint_loop_marker '%' loop_marker{ll}(i) '\t '];
            end
            
            for m = 1:length(result.unique_coherence)
                for k=1:length(result.unique_stim_type)
                    for ind=1:3 % heading, rotation, motion type
                        try
                            eval(['buff = sprintf(sprint_loop_marker,' loop{ll} ');']);
                        catch
                            buff = sprintf(sprint_loop_marker,nan(1,sum(sprint_loop_marker=='%')));
                        end
                        fprintf(fid, '%s', buff);
                        toClip = [toClip sprintf('%s', buff)];
                    end
                end
            end
        end
    end
elseif strcmpi(config.suffix,'SpiT')
    loop = config.sprint_loop_contents;
    loop_marker = config.sprint_loop_marker;
    for ll = 1:length(loop)
        sprint_loop_marker = [];
        for i = 1:length(loop_marker{ll})
            sprint_loop_marker = [sprint_loop_marker '%' loop_marker{ll}(i) '\t '];
        end
        
        for ind = 1:3
            if ind~=3
                r=2; % 2.�exclude exp/con
            else
                r=1;% 1.�includ exp/con
            end
            try
                eval(['buff = sprintf(sprint_loop_marker,' loop{ll} ');']);
            catch
                buff = sprintf(sprint_loop_marker,nan(1,sum(sprint_loop_marker=='%')));
            end
            fprintf(fid, '%s', buff);
            toClip = [toClip sprintf('%s', buff)];
        end 
    end
elseif strcmpi(config.suffix,'CueS_performance')
    for ll = 1:length(config.sprint_loop_contents)
        sprint_loop_marker = [];
        for i = 1:length(config.sprint_loop_marker{ll})
            sprint_loop_marker = [sprint_loop_marker '%' config.sprint_loop_marker{ll}(i) '\t '];
        end
        last_step = result.last_step;
        for ind = 1:3 % h r mt
            try
                eval(['buff = sprintf(sprint_loop_marker,' config.sprint_loop_contents{ll} ');']);
            catch
                buff = sprintf(sprint_loop_marker,nan(1,sum(sprint_loop_marker=='%')));
            end
%              fprintf(fid, '%s', buff);
            toClip = [toClip sprintf('%s', buff)];
        end
    end
end

% if strcmpi(config.suffix,'CueS')
%     %toClip(isnan(toClip))=[];  %remove NaN     Here NaN is string, not the logical NaN
%     toClip = toClip(1:(end-21));  %remove the last four 'NaN'
% end

fprintf(fid, '\r\n');

if ~isempty(config.batch_flag)  % Close .dat file
    fclose(fid);
end

toClip = [toClip(1:end-1) sprintf('\r\n')]; % Delete the last '\t' for clipboard
%     clipboard('copy',toClip);
%     Turn [NaN]s (number) into 'NaN's (string)
%     for ll = 1:length(toXls)
%         if isnan(toXls{ll})
%             toXls{ll} = 'NaN';
%         end
%     end

% --- Save back to xls ---
% (finally I decide to save somethings back to xls for easier and better visualization in xls). HH20150724
sheet_num = 1;
xlRange = 'A1:DG3000';
if ~isempty(config.batch_flag) && isfield(config,'xls_column_begin')
    % Turn toClip into toXls (Separate strings by TAB and reoranize into cells)
    toXls = textscan(toClip,'%s','Delimiter','\t');
    toXls = toXls{1}';
    
    % Read xls if needed. (only for the first file in BATCH mode)
    if isempty(XlsData) || strcmp(config.batch_flag,'test.m')  % If we are in test mode, we reload xls each time. HH20160415
        XlsData = ReadXls(excel_path,sheet_num,3,xlRange);
        if isempty(XlsData)
            disp('Check Excel path, sheet, header!!!');
            keyboard
        end
    end
    
    % Locate where to paste "toClip"
    row = find(strcmp(XlsData.txt(:,XlsData.header.FileNo),result.FILE)) + 3;
    if isempty(row)
        row = find(strcmp(XlsData.txt(:,XlsData.header.FileNo),result.FILE(1:end-4))) + 3;
    end

    if  isempty(config.xls_column_begin)
        disp('No saving back to .xls...');
    elseif ~isempty(row) 
        if numel(row) > 1 % More than one SpikeChans
            row = intersect(row,find(XlsData.num(:,XlsData.header.Chan1) == result.SpikeChan)+3);
        else % Only one SpikeChan or no SpikeChan (Training sessions)
        end
        
        column_begin = XlsData.header.(config.xls_column_begin);
        column_end = XlsData.header.(config.xls_column_end);
        
        % Write "toClip" into excel file
        if length(toXls) == column_end - column_begin + 1
            column_begin_name = num2ExcelName(column_begin);
            range_name = [column_begin_name num2str(row)];
            try
                xlswrite(excel_path, toXls,sheet_num,range_name);  % Speed-up of xlswrite
                disp('Writing to .xls finished...');
            catch
                try
                    options.Interpreter = 'tex';
                    options.Default = 'Yes';
                    button=questdlg('\fontsize{15}\color{red} Waiting for manual save and close the Excel ! Then you can press "YES..."','Warning!','Yes',options);
                    if strcmp(button,'Yes')
                        % Shows a list of active excel processes
                        !tasklist /FI "IMAGENAME eq excel.exe"
                        xlswrite(excel_path, toXls,sheet_num,range_name);  % Speed-up of xlswrite
                        % system('taskkill /F /IM EXCEL.EXE'); % STOPS ALL EXCEL PROCESSES [USE WITH CAUTION] 自动杀掉EXCEL程序
                        try
                            excelApp = actxserver('Excel.Application'); % Starts an Excel Application
                            disp('Writing to .xls finished...');
                        catch
                            disp('You do not close the EXCEL!');
                            keyboard
                        end
                    end
                catch
                    disp('Writing to .xls failed :<');
                    %keyboard
                end
            end
            
            % write the "Coherence" in specific column, Lwh 20200803
            coherence_column = num2ExcelName(XlsData.header.coherence);
            coherence_range_name = [coherence_column num2str(row)];
            coherence = num2str(result.unique_coherence);
            toXls_co = textscan(coherence,'%s','Delimiter','\t');
            xlswrite(excel_path, toXls_co{1},sheet_num,coherence_range_name);
        else
            disp('Size not match when write back to xls...');
            keyboard
        end
    else
        disp('No file entry found in .xls...');
    end
end
end

if strcmp(config.batch_flag,'test.m')
    assignin('base','result',result); 
end

end

function [col_str] = num2ExcelName(num_loc) % Convert xls column number to column name (such as "A", "AB", ...)
test = 2;
old = 0;
x = 0;
while test >= 1
    old = 26^x + old;
    test = num_loc/old;
    x = x + 1;
end
num_letters = x - 1;
str_array = zeros(1,num_letters);
for i = 1:num_letters
    loc = floor(num_loc/(26^(num_letters-i)));
    num_loc = num_loc - (loc*26^(num_letters-i));
    str_array(i) = char(65 + (loc - 1));
end
col_str = strcat(str_array(1:length(str_array)));
end