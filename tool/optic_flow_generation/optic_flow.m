% by LZQ

function [varargout] = optic_flow(varargin)
%% default input parameters:
options = struct('plot',1,'density',30,'duration',400,'angular_velocity',0,'velocity',[0,0.5,0],'ang_mode',1);

%% read input parameters
optionNames = fieldnames(options);
if mod(length(varargin),2) == 1
    error('Please provide propertyName/propertyValue pairs')
end
for pair = reshape(varargin,2,[])    % pair is {propName; propValue}
    if any(strcmp(pair{1}, optionNames))
        options.(pair{1}) = pair{2};
    else
        error('%s is not a recognized parameter name', pair{1})
    end
end
%%
t = options.duration; % ms
density = options.density;
angular_v = options.angular_velocity;
velocity = options.velocity; % 像素/ms ?

rang = velocity*t;
vol = (rang(1)+1)*(rang(2)+1)*(rang(3)+1);% 运动经过的空间
n = vol*density; % 运动空间内的点个数
x = (rand(1,n)-0.5)*(rang(1)+10)+rang(1)/2;
y = (rand(1,n)-0.5)*(rang(2)+10)+rang(2)/2;
z = (rand(1,n)-0.5)*(rang(3)+10)+rang(3)/2;
if options.plot==1
    figure
end
switch options.ang_mode
    case 0
        axs = -0.5:0.001:0.5;
        s = length(axs);
        stim = zeros(s,s,t);
        
    case 1
        axs = -pi/6:0.001:pi/6;
        s = length(axs);
        stim = zeros(s,s,t);
end
for i = 1:t
    [x_,y_,z_] = movecam(x,y,z,velocity,norm(velocity)*i);
    [x_,y_,z_] = rotatecam(x_,y_,z_,angular_v*i);
    [t,p] = camproject(x_,y_,z_);
    
    switch options.ang_mode
        case 0
            [a,b] = retina2ISO(t,p);
            [~,i1] = arrayfun(@(x)min(abs(axs-x)),a);
            [~,i2] = arrayfun(@(x)min(abs(axs-x)),b);
            tmp = stim(:,:,i);
            idx = sub2ind(size(tmp),i1,i2);
            tmp(idx)=1;
            stim(:,:,i) = tmp;
            if options.plot==1
                plot(a,b,'.');
                xlim([-0.5,0.5]);
                ylim([-0.5,0.5]);
            end
        case 1
            [~,i1] = arrayfun(@(x)min(abs(axs-x)),t);
            [~,i2] = arrayfun(@(x)min(abs(axs-x)),p);
            tmp = stim(:,:,i);
            idx = sub2ind(size(tmp),i1,i2);
            tmp(idx)=1;
            stim(:,:,i) = tmp;
            if options.plot==1
                
                plot(t,p,'.');
                xlim([pi/3,pi*2/3]);
                ylim([-pi/6,pi/6]);
            end
    end
    pause(0.01)
end
%%
if nargout>=1
    varargout{1} = stim;
    if nargout>=2
    varargout{2} = axs;
    end
end
if options.plot==1             
close
end