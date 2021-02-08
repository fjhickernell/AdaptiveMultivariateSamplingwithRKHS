function outString = cleanStringFJH(inString)
% This function takes strings of MATLAB formatted numbers and replaces
% certain characters by others that look better for a table
outString = strrep(inString,'E','\\text{E}');
outString = strrep(outString,'-0','{-}');
outString = strrep(outString,'+0','');
outString = strrep(outString,'e','\\text{E}');
