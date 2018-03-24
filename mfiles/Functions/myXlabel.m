function myXlabel(label, fontSize)
% Function: myXlbel
%
% Author: Kurt Nelson
%
% Purpose: creates xlabel with Latexfont
%
% Inputs:
% 1) label - x label
% 2) fontSize - fontsize of label
%%
xlab = xlabel(label);
set(xlab,'interpreter','Latex','FontSize',fontSize)
end

