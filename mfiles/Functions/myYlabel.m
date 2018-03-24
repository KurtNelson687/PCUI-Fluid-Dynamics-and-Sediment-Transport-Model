function myYlabel(label, fontSize)
% Function: myYlbel
%
% Author: Kurt Nelson
%
% Purpose: creates ylabel with Latexfont
%
% Inputs:
% 1) label - x label
% 2) fontSize - fontsize of label
%%
ylab = ylabel(label);
set(ylab,'interpreter','Latex','FontSize',fontSize)
end

