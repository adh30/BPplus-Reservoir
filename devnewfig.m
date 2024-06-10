close all
PP=aoPr_av-min(aoP_av);
subplot(2,2,3)
plot(Time_av, PP); hold on;
if (lsys <= dicroticNotchIdx)
    area(Time_av(1:dicroticNotchIdx),PP(1:dicroticNotchIdx),'FaceColor', [1 0.8 1]);
    area(Time_av(1:lsys),PP(1:lsys),'FaceColor', [0.8 0.8 1]);
    %area1=area(PP(1:lsys));
else
    area(Time_av(1:lsys),PP(1:lsys),'FaceColor', [0.8 0.8 1]);
    area(Time_av(1:dicroticNotchIdx),PP(1:dicroticNotchIdx),'FaceColor', [1 0.8 1]);
end
% x = [lsys/samplerate,lsys/samplerate];
% y = [0,PP(lsys)];
% plot(x,y);
% x = [dicroticNotchIdx/samplerate,dicroticNotchIdx/samplerate];
% y = [0,PP(dicroticNotchIdx)];
% plot(x,y);
%sevr_title = sprintf('%s SEVR=%3.0f %% cAIx=%3.0f %% cAP=%3.0f %% cPP=%3.0f %%',filename,ao_sevr*100,ao_ai,ao.sbp-aoPi,ao.pp);
title('SEVR')
%cAIx_txt = sprintf('\\leftarrow T1',aoTmax);
cAIx_txt = '\leftarrow T1';
offset = 0.005;
text((ao_Ti)+offset,PP(ao_Ti*samplerate),cAIx_txt,'FontSize', 6);
% if to catch when lsys==dicroticNotchIdx so they dont overwrite.
if lsys~=dicroticNotchIdx
    maxnegdPdt_txt = '\leftarrow max(-dP/dt)';
    text(((lsys+1)/samplerate),PP(lsys),maxnegdPdt_txt,'FontSize', 6);
else
    maxnegdPdt_txt = 'max(-dP/dt) =';
    text(((lsys+1)/samplerate),PP(lsys)+.075*(max(PP)),maxnegdPdt_txt,'FontSize', 6);
end
dicroticNotchIdx_txt = '\leftarrow notch';
text(dicroticNotchIdx+offset,PP(dicroticNotchIdx),dicroticNotchIdx_txt,'FontSize', 6);
text((lsys/(2*samplerate)),10,'SPTI','FontSize', 6, HorizontalAlignment='center');
text((lsys+lsys/2)/samplerate,10,'DPTI','FontSize', 6);