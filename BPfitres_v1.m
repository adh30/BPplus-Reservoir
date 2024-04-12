%% BPfitres_v1 - batch analysis of BPplus data using kreservoir_v14 by KHP
%  Copyright 2019 Alun Hughes
%  This software is distributed under under the terms of the GNU General Public License
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%  http://www.gnu.org/licenses/gpl.html

%% Versions
% v1 First stable version(03/01/20) based on fitres_v5 and uses
% kreservoir_v14
%%
function  [Tn_av, Pinf_av, P_av,Pr_av,Pn_av,...
    fita_av, fitb_av,rsq_av]=BPfitres_v1(p_av,samplerate)
%%
% analyse the average beat
% *****NB  peripheral pulse starts at the foot!!!

% all data
T_av=length(p_av)/samplerate;

% prevent upturn in pressure affecting fit (different methods compared
% but settled on method a (other methods retained but commented out

% (a) exclude diff (p)>0 in diastole
diffp=diff(p_av);
cut=find(diffp<0,1,'last');
P_av=p_av(1:cut);

% fit average beat using kreservoir
[Pr_av,fita_av,fitb_av,Pinf_av, Tn_av,Pn_av]=kreservoir_v14(P_av,T_av,samplerate);

% Calculate R^2 (cofficient determination) for Pr fit in diastole
xn=round(Tn_av*samplerate);         % parameters for length diastole for R2
C = corrcoef(P_av(xn:end),Pr_av(xn:end)); % calculate correlation for diastolic fit
rsq_av = C(1,2)^2;
end