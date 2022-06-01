 function [MI,normalMI] = nmi(A,B) %NMI Normalized mutual information
% http://en.wikipedia.org/wiki/Mutual_information
% http://nlp.stanford.edu/IR-book/html/htmledition/evaluation-of-clustering-1.html
% Author: http://www.cnblogs.com/ziqiao/   [2011/12/13] if length( A ) ~= length( B)

if length( A ) ~= length( B)
    error('length( A ) must == length( B)');
end

total = length(A);
A_ids = unique(A);
B_ids = unique(B);

% Mutual information
MI = 0;
for idA = 1:length(A_ids)
    for idB = 1:length(B_ids)
         idAOccur = find( A == A_ids(idA) ); 
         idBOccur = find( B == B_ids(idB) );
         idABOccur = intersect(idAOccur,idBOccur); 
         
         px = length(idAOccur)/total;
         py = length(idBOccur)/total;
         pxy = length(idABOccur)/total;
         
         % ��
         MI = MI + pxy*log2(pxy/(px*py)+eps); % eps : the smallest positive number
    end
end

% Normalized Mutual information
Hx = 0; % Entropies
for idA = 1:length(A_ids)
    idAOccurCount = length( find( A == A_ids(idA) ) );
    Hx = Hx - (idAOccurCount/total) * log2(idAOccurCount/total + eps);
end
Hy = 0; % Entropies
for idB = 1:length(B_ids)
    idBOccurCount = length( find( B == B_ids(idB) ) );
    Hy = Hy - (idBOccurCount/total) * log2(idBOccurCount/total + eps);
end

normalMI = 2 * MI / (Hx+Hy);
% normalMI = MI / sqrt(Hx*Hy); another version of NMI
end


% Example :  
% (http://nlp.stanford.edu/IR-book/html/htmledition/evaluation-of-clustering-1.html)
% A = [1 1 1 1 1 1   2 2 2 2 2 2    3 3 3 3 3];
% B = [1 2 1 1 1 1   1 2 2 2 2 3    1 1 3 3 3];
% nmi(A,B) 
% ans =  0.3646