function PlotDs(time,X,Y,Z)
%% Documentation
% Function to plot out selected part of the synthetic seismograms
% Added by Xiaoran Chen on 08/16/2020
% Matlab R2016a

MAX=max([X, Y, Z]);
MIN=min([X, Y, Z]);
shift_max = abs(MAX)*0.2;
shift_min = abs(MIN)*0.2;

figure
h1=subplot(3,1,1);
plot(time,X,'b','LineWidth',2.5);
% axis([time(1) time(end) MIN-shift_min MAX+shift_max]);
set(gca,'XTicklabel',[]);
set(gca,'linewidth',2);
ylabel('R','FontSize',16,'linewidth',3);

h2=subplot(3,1,2);
plot(time,Y,'b','LineWidth',2.5);
% axis([time(1) time(end) MIN-shift_min MAX+shift_max]);
set(gca,'XTicklabel',[]);
set(gca,'linewidth',2);
ylabel('T','FontSize',16','linewidth',3);

h3=subplot(3,1,3);
plot(time,Z,'b','LineWidth',2.5);
% axis([time(1) time(end) MIN-shift_min MAX+shift_max]);
set(gca,'linewidth',2);
xlabel('Time(s)','FontSize',16','linewidth',3);
ylabel('Z','FontSize',16','linewidth',3);

p1 = get(h1, 'Position');
p2 = get(h2, 'Position');
p3 = get(h3, 'Position');

p2(2) = p3(2)+p3(4);
p1(2) = p2(2)+p2(4);
set(h1, 'pos', p1);
set(h2, 'pos', p2);
savefig(sprintf('%d.pdf',baz));






