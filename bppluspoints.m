% replicate index points in BP+
sr=0.200;
t=(1:length(ba_p_av))/sr;
plot(t, ba_p_av,'Color',[0.3,0.75,0.93]); hold on; 
plot(t, ao_p_av,'Color',[0.5,0.5,0.5])
td1=(1:length(ba_p_av)-1)/sr;
plot(td1, diff(ba_p_av)*10,'Color',"k")
td2=(1:length(ba_p_av)-2)/sr;
plot(td2,diff(diff(ba_p_av)*10)*10,"Color",'m')
