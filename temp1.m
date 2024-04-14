% temporary - used to check code to read the default directory from the
% config.json file - once it works I'll incorporate it into bpp_Res2.m
clear
text = fileread('config.json');
data = jsondecode(text);
%file_lists=dir(fullfile(string(data.folder_name), '*.*'));
folder_name = string(data.folder_name);
disp(folder_name)