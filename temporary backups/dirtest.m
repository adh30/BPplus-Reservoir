% testing json etc
clear
% get user directory
if ispc; userdir=getenv('USERPROFILE');
else; userdir= getenv('HOME');
end
% read json
S=readstruct("config.json", FileType="json");
A=convertStringsToChars(S.folder_name);
datadir=[userdir A];
% datadir=[userdir '\BPPdata'];
if not(isfolder(datadir))
   % set data folder
   dirnew = uigetdir(userdir);
   S.folder_name=dirnew;
   writestruct(S,"config.json",PreserveInfAndNaN=false);
   % type config.json
else
    disp(['Data folder is ' datadir])
end
