close all
PP=aoPr_av-min(aoP_av);
plot(Time_av, PP); hold on;
if (lsys <= dicroticNotchIdx)
    area(Time_av(1:dicroticNotchIdx),PP(1:dicroticNotchIdx),'FaceColor', [1 0.8 1]);
    area(Time_av(1:lsys),PP(1:lsys),'FaceColor', [0.8 0.8 1]);
    %area1=area(PP(1:lsys));
else
    area(Time_av(1:lsys),PP(1:lsys),'FaceColor', [0.8 0.8 1]);
    area(Time_av(1:dicroticNotchIdx),PP(1:dicroticNotchIdx),'FaceColor', [1 0.8 1]);
end