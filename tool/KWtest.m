function [P1, P2, X2, F] = KWtest(X,alpha)
% add output by Lwh 20220314
% P1：Chi-squared distribution (same as original function "kruskalwallis")
% P2: F approximations: less conservative and this means that it refuses the H0 more than you want
% input:
% [group1_data 1
%  group1_data 1
%  group1_data 1
%  group2_data 2
%  group2_data 2
%  group3_data 3]

% Kruskal-Wallis' nonparametric analysis of variance.
%(This file is applicable for equal or unequal sample sizes and if there are tied ranks.)
%
%   Syntax: function [KWtest] = KWtest(X,alpha) 
%      
%     Inputs:
%          X - data matrix (Size of matrix must be n-by-2; data=column 1, sample=column 2). 
%       alpha - significance level (default = 0.05).
%     Outputs:
%          - Sample medians vector.
%          - Whether or not the equality among samples was met.
%
%    Example: From the example 10.11 of Zar (1999, p. 199-200), to test the 
%             equality among samples by the Kruskal-Wallis' nonparametric 
%             analysis of variance with a significance level = 0.05.
%
%                                 Sample
%                   ---------------------------------
%                       1       2       3       4
%                   ---------------------------------
%                      7.68    7.71    7.74    7.71
%                      7.69    7.73    7.75    7.71
%                      7.70    7.74    7.77    7.74
%                      7.70    7.74    7.78    7.79
%                      7.72    7.78    7.80    7.81
%                      7.73    7.78    7.81    7.85
%                      7.73    7.80    7.84    7.87
%                      7.76    7.81            7.91
%                   ---------------------------------
%                                       
%       Data matrix must be:
%    X=[7.68 1;7.69 1;7.70 1;7.70 1;7.72 1;7.73 1;7.73 1;7.76 1;7.71 2;7.73 2;7.74 2;7.74 2;
%    7.78 2;7.78 2;7.80 2;7.81 2;7.74 3;7.75 3;7.77 3;7.78 3;7.80 3;7.81 3;7.84 3;7.71 4;7.71 4;
%    7.74 4;7.79 4;7.81 4;7.85 4;7.87 4;7.91 4];
%
%     Calling on Matlab the function: 
%             KWtest(X)
%
%       Answer is:
%
% The number of samples are: 4
%
%  Sample    Size        Median          Sum of Ranks
%    1        8          7.7100             55.00
%    2        8          7.7600            132.50
%    3        7          7.7800            145.00
%    4        8          7.8000            163.50
%   
% Kruskal-Wallis' test for difference among samples X2=11.9435, F= 5.9531
% With X2, the assumption of equality among samples was not met.
% With F, the assumption of equality among samples was not met.
%

%  Created by A. Trujillo-Ortiz and R. Hernandez-Walls
%             Facultad de Ciencias Marinas
%             Universidad Autonoma de Baja California
%             Apdo. Postal 453
%             Ensenada, Baja California
%             Mexico.
%             atrujo@uabc.mx
%
%  April 26, 2003.
%
%  To cite this file, this would be an appropriate format:
%  Trujillo-Ortiz, A. and R. Hernandez-Walls. (2003). KWtest: Kruskal-Wallis' 
%    nonparametric analysis of variance. A MATLAB file. [WWW document]. URL http://
%    www.mathworks.com/fileexchange/loadFile.do?objectId=3361&objectType=FILE
%
%  References:
% 
%  Zar, J. H. (1999), Biostatistical Analysis (2nd ed.).
%           NJ: Prentice-Hall, Englewood Cliffs. p. 195-200. 
%

if nargin < 2,
   alpha = 0.05;
end 

k=max(X(:,2));
% fprintf('The number of samples are:%2i\n\n', k); % remopve by Lwh 20220314

sample=[];
indice=X(:,2);
for i=1:k
   Xe=find(indice==i);
   eval(['x' num2str(i) '=X(Xe,1);']);
   eval(['x= x' num2str(i) ';']);
   sample=[sample;x];
end

N=[];
for i=1:k
   eval(['n' num2str(i) '=length(x' num2str(i) ');']);
   eval(['x= n' num2str(i) ';']);
   N=[N,x];
end

Md=[];
for i=1:k
   eval(['Md' num2str(i) '=median(x' num2str(i) ');']);
   eval(['xMd= Md' num2str(i) ';']);
   Md=[Md;xMd];
end

%Ranks of the samples procedure.
Y=X(:,1);
sz=sort(Y);
tiescor=0;
for i=1:k
   eval(['ka=n' num2str(i) ';']);
   for idx=1:ka
      eval(['tmp=find(sz==x' num2str(i) '(idx));']);
      nties=length(tmp);
      if nties > 1
        eval(['rank' num2str(i) '(idx)=mean(tmp);']);
        tiescor=tiescor+(nties.^2-1)./(sum(N)*(sum(N)-1));
      else
       eval(['rank' num2str(i) '(idx)=tmp;']);
      end
   end
end

R=[];  %sum of the ranks.
for i=1:k
  eval(['R' num2str(i) '=sum(rank' num2str(i) ');']); 
  eval(['x= R' num2str(i) ';']);
  R=[R,x];
end

C=[];
for i=1:k
   eval(['c' num2str(i) '=((R' num2str(i) ').^2)/n' num2str(i) ';']);
   eval(['x= c' num2str(i) ';']);
   C=[C,x];
end

% remopve by Lwh 20220314
% disp(' Sample    Size            Median            Sum of Ranks')
% for i=1:k
%    fprintf('   %d       %2i         %11.4f         %11.2f\n',i,N(i),Md(i),R(i))
% end
% disp(' ')

H=((12/(sum(N)*(sum(N)+1)))*sum(C))-(3*(sum(N)+1));  %Kruskal-Wallis' statistic.

n=sum(N);
y=sz;
c=[];
m=1;
m2=1;
while m<=n
   b=y(m);
   cont=1;
   c(m2,1)=b;
   c(m2,2)=cont;
   if m<n;
      while y(m)==y(m+1);
         cont=cont+1;
         c(m2,2)=cont;
         m=m+1;
         if m==n;
            break;
         end
      end
   end
   m=m+1;
   m2=m2+1;
end
x=c(:,1);
f=c(:,2);
t3=f.^3;
t=t3-f;
T=sum(t);  %Total number of ties.

CT=1-(T/((sum(N)^3)-sum(N)));  %correction factor.
Hc=H/CT;  %Kruskal-Wallis' corrected statistic.

X2=Hc; %Chi-squared-statistic approximation.
v=k-1; %degrees of freedom.
P1=1-chi2cdf(X2,v);  %probability associated to the Chi-squared-statistic.

F=((n-k)*Hc)/((k-1)*(n-1-Hc));  %F-statistic approximation.
P2=1-fcdf(F,v,n-k-1);  %probability associated to the F-statistic.

% remopve by Lwh 20220314
% fprintf('Kruskal-Wallis'' test for difference among samples X2=%3.4f, F=%7.4f\n', X2,F);
% 
% if P1 >= alpha;
%    fprintf('With X2, the assumption of equality among samples was met.\n');
% else
%    fprintf('With X2, the assumption of equality among samples was not met.\n');
% end
% 
% if P2 >= alpha;
%    fprintf('With F, the assumption of equality among samples was met.\n');
% else
%    fprintf('With F, the assumption of equality among samples was not met.\n');
% end

