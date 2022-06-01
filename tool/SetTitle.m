% title setting, Lwh 202102
% example
% title_select = {[2] [1] monkey_included_for_analysis} (monkey_included_for_analysis=7)
% "2-Target Fine task (Ringbell)"
%
% title_select = {[2] [1] monkey_included_for_analysis} (monkey_included_for_analysis=[7 16])
% "2-Target task (Two Monkey)"
%
% [til] = titleset(stimtask_type,monkey_included_for_analysis);
% str = ['Correlation between motion type PSE shift and spiral index' newline til];
% title(str);

function [til,Two_or_Four,fine_or_coarse,monkey] = SetTitle(stimtask_type_temp,monkey_included_for_analysis)

% stimtask_type should be two part for fine task selection and coarse task
% selection
% e.g.
% simtask_type{1} = [2 3];
% simtask_type{2} = [4 5];

% for title
tartitle{1} = ["2-Target "];
tartitle{2} = ["4-Target "];
tartitle{3} = ["2&4-Target "];

fctitle{1} = ["Fine "];
fctitle{2} = ["Coarse "];
fctitle{3} = ["Fine & Coarse "];

mktitle{1} = ["(Ringbell)"];
mktitle{2} = ["(Arthas)"];
mktitle{3} = ["(Two Monkey)"];


if ~isempty(stimtask_type_temp) % fine or caorse, 2 or 4-AFC
    % turn to cell
    if iscell(stimtask_type_temp)
        stimtask_type_temp = cell2mat(stimtask_type_temp);
    end
    
    select_fine_part = ismember(stimtask_type_temp,[1 2 3]);
    stimtask_type{1} = stimtask_type_temp(select_fine_part);
    
    select_coarse_part = ismember(stimtask_type_temp,[4 5 6]);
    stimtask_type{2} = stimtask_type_temp(select_coarse_part);
    
    % fine or coarse
    if ~isempty(stimtask_type{1}) && ~isempty(stimtask_type{2})
        fine_or_coarse = fctitle{3};
    elseif ~isempty(stimtask_type{1}) && isempty(stimtask_type{2})
        fine_or_coarse = fctitle{1};
    elseif isempty(stimtask_type{1}) && ~isempty(stimtask_type{2})
        fine_or_coarse = fctitle{2};
    end
    
    % 2 or 4 target
    num = max(cellfun(@(x) length(x),stimtask_type));
    if num == 3
        Two_or_Four = tartitle{3}; % 例如[1 2 3],[4 5 6]分别是4-T fine and coarse task
    elseif num == 2
        Two_or_Four = tartitle{1}; % 例如[2 3],[5 6]分别是2-Target fine and coarse task
    elseif num == 1
        Two_or_Four = tartitle{2}; % 例如[1],[4]分别是4-Target fine and coarse task
    end
    
else % for tuning task
    Two_or_Four = 'Tuning ';
    fine_or_coarse = [];
end


% which monkey, 一定会有的
if length(monkey_included_for_analysis)==2 % two monkey1
    monkey = mktitle{3};
elseif length(monkey_included_for_analysis)==1 % one monkey
    if monkey_included_for_analysis==7
        monkey = mktitle{1};
    elseif monkey_included_for_analysis==16
        monkey = mktitle{2};
    else
        disp('Please add monkey');
        keyboard
    end
else
    disp('Please add monkey');
    keyboard
end

til = strcat(Two_or_Four,fine_or_coarse,"task ",monkey);
til = char(til);

Two_or_Four = char(Two_or_Four);
fine_or_coarse = char(fine_or_coarse);
monkey = char(monkey);

end
