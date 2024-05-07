close all
x = simoutput.max_r;
y1 = simoutput.percent_r;
y2 = simoutput.novel_r;


plot(x,y1,'k');
hold on
addaxis(x,y2,'r');

addaxislabel(1,'% real resected nodes');
addaxislabel(2,'% unresected nodes');
