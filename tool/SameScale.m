function h_diag = SameScale(diag,fixLim)

% Lwh 202101, for x,y same scale
% diag = 1; plot the diagonal line 对角线
% fixLim: 1: keep original lim; 0: a bit larger than original limit

if nargin<1
    diag = 0;
    fixLim = 0;
elseif nargin<2
    fixLim = 0;
end

axis square
xlims = xlim;
ylims = ylim;
xlims_new = [min(xlims(1),ylims(1)) max(xlims(2),ylims(2))];

if ~fixLim
    xlims_new = [xlims_new(1)-range(xlims_new)/20 xlims_new(2)+range(xlims_new)/20];
end

axis([xlims_new xlims_new]);
xtick_set = get(gca,'xtick');
set(gca,'ytick',xtick_set);

if diag
    hold on
    h_diag = plot([xlims_new(1) xlims_new(2)],[xlims_new(1) xlims_new(2)],'k-');
end

% keep this tick not change, Lwh 202102
xticks('manual');
yticks('manual');

% axis equal

end