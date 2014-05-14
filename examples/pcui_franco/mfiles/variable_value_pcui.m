function value = variable_value_pcui(variable,ftext)
%
% Filename    : variable_value_pcui.m
% Author      : Goncalo Gil
% Description : Searches a file containing variable assignments with equal
%               signs and extracts just the value given the variable name. 
%               Designed for Fortran type assignments in decimal or 
%               exponential form.
%                                       
% Author      : Goncalo Gil, Stanford University
% email       : gilg@stanford.edu
%

% REGEX expression that finds line containing variable
expr = [variable '\s*=\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)'];
temp = regexp(ftext, expr, 'match');
% REGEX expression that extracts the value from the extracted string 'temp'
expr = ['(?<=\s)[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?'];
temp = regexp(temp, expr, 'match');
value = str2double(temp{1});