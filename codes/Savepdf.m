function Savepdf(time,X,Y,Z,baz,K,name,Nmodel,targetdir)
%% Documentation
% Function to plot and save all the files in .pdf format
% Added by Xiaoran Chen on 08/16/2020
% Matlab R2016a

%% Edition starts from here
X=X(1:K);
Y=Y(1:K);
Z=Z(1:K);
time=time(1:K);

figure
h1=subplot(3,1,1);
plot(time,X,'b','LineWidth',2.5);
set(gca,'XTicklabel',[]);
set(gca,'linewidth',2);
ylabel('R','FontSize',16,'linewidth',3);

h2=subplot(3,1,2);
plot(time,Y,'b','LineWidth',2.5);
set(gca,'XTicklabel',[]);
set(gca,'linewidth',2);
ylabel('T','FontSize',16','linewidth',3);

h3=subplot(3,1,3);
plot(time,Z,'b','LineWidth',2.5);
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
saveas(gcf,[targetdir '/Output/' sprintf('%s-%d-%0.f.pdf',name, Nmodel, baz)]);






