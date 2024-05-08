% temporary - used to check code to read the default directory from the
% config.json file - once it works I'll incorporate it into bpp_Res2.m
clear
jtext = fileread('config.json');
jdata = jsondecode(jtext);
%file_lists=dir(fullfile(string(data.folder_name), '*.*'));
folder_name = jdata.folder_name;
disp(folder_name)