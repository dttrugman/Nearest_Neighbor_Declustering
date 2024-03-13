function allx(X)

% Sets new x-axis to all current subplots.
%
% Usage
%     allx(X)
%
% Input:
%     X is a two-element vector: [Xmin Xmax]
%


ch=get(gcf,'Children');
for i=1:length(ch)
   set(ch(i),'XLim',X);
end
