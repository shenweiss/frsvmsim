close all
x = simoutput.max_r;
y1 = simoutput.vSOZRR;
y2 = simoutput.vRonSRR;
y3 = simoutput.vFRRR;
y4 = simoutput.mRDRRD; 
y5 = simoutput.mgRR;
y6 = simoutput.urmLE;

plot(x,y1,'k');
hold on
addaxis(x,y2,'r');
addaxis(x,y3,'g');
addaxis(x,y4,'b');
addaxis(x,y5,'c');
addaxis(x,y6,'m');

addaxislabel(1,'vSOZRR');
addaxislabel(2,'vRonSRR');
addaxislabel(3,'vFRRR');
addaxislabel(4,'mRDRRD');
addaxislabel(5,'mgRR');
addaxislabel(6,'urmLE');