%Exemplo: Filtro de Kalman

close all; clear all;

% controlled signal
randn("seed", 10);
var_ruido_proc = 0.05;
var_ruido_obs =0.25;
T1 = [1:500];
T2 = [501:2501];
y1 = sin(2*pi*1/50*T1);
y2 = 3*sin(2*pi*1/100*T2);
y = [y1, y2];
yk = y + sqrt(var_ruido_obs) * randn(1,length(y));

%%%%%%%%%%%%%%%%%%
%DEKF

shift = 2500;
window = shift;
current = 1;
x_estim = zeros(size(y));
while (current + window <= length(yk))
  chunk = yk(current:current+window-1);
  if (current == 1)
    [xk_estim_pos, Pxk_pos, wk_estim_pos, Pwk_pos]= dekf(chunk, var_ruido_proc, var_ruido_obs);
  else
    [xk_estim_pos, Pxk_pos, wk_estim_pos, Pwk_pos]= dekf(chunk, var_ruido_proc, var_ruido_obs, xk_estim_pos(:,:,end), Pxk_pos(:,:,end), wk_estim_pos(:,:,end), Pwk_pos(:,:,end));
  end
  x_estim(current:current+window-1) = xk_estim_pos(1,1,:);
  current = current + shift;
end

%%%%%%%%%%%%%%%%%%
% Results: states and theirs estimates

final_pos = length(x_estim);
squared_error = (y(1:final_pos) - x_estim).^2;
MSE = mean(squared_error)

norm_diff_wk = zeros(size(wk_estim_pos,3) - 1,1);
norm_wk = zeros(size(wk_estim_pos,3),1);
for i=2:size(wk_estim_pos,3)
  diff_wk = wk_estim_pos(:,1,i) - wk_estim_pos(:,1,i-1);
  norm_diff_wk(i-1) = sqrt(diff_wk' * diff_wk);
  norm_wk(i) = sqrt(wk_estim_pos(:,1,i)' * wk_estim_pos(:,1,i));
end

h = figure(1);
%title(["Response of an Underdamped Single Degree of\n"...
%      "Freedom System Subjected to an Initial Excitation"],...
%      'FontName','/usr/share/fonts/dejavu/DejaVuSerif-Italic.ttf',...
%     'FontSize',8);

FN = findall(h,'-property','FontName');
set(FN,'FontName','/usr/share/fonts/dejavu/DejaVuSerifCondensed.ttf');
FS = findall(h,'-property','FontSize');
set(FS,'FontSize',16);

subplot(2,1,1)
hold on
plot(yk(1:final_pos),'ko');
plot(y(1:final_pos), 'r','LineWidth',4);
plot(x_estim,'b-.','LineWidth',4);
set (gca, "fontsize", 16)
xlim([0,2500])
L=legend('Observations', 'True States', 'Estimated', 'location',  'north', 'orientation', 'horizontal');
FL1= findall(L,'-property','FontName');
set(FL1,'FontName','/usr/share/fonts/msttcore/cour.ttf');
set(L,'FontSize',16);
ylabel('Value','FontSize',16);

subplot(2,1,2)
plot(norm_diff_wk);
set (gca, "fontsize", 16)
xlim([0,2500])
xlabel('Sample','FontSize',16);
ylabel('Value','FontSize',16);
L=legend('Neural net delta weight');
FL1= findall(L,'-property','FontName');
set(FL1,'FontName','/usr/share/fonts/msttcore/cour.ttf');
set(L,'FontSize',16);


H = 9; W = 12;
set(h,'PaperUnits','inches')
set(h,'PaperOrientation','portrait');
set(h,'PaperSize',[H,W])
set(h,'PaperPosition',[0,0,W,H])
print(h,'-dpng','-color','impacto_snr.png');



h = figure(2);
FN = findall(h,'-property','FontName');
set(FN,'FontName','/usr/share/fonts/dejavu/DejaVuSerifCondensed.ttf');
FS = findall(h,'-property','FontSize');
set(FS,'FontSize',16);

hold on
plot(yk(1:500),'ko','LineWidth',2);
plot(y(1:500), 'r','LineWidth',4);
plot(x_estim(1:500),'b--','LineWidth',4);
set (gca, "fontsize", 16)
xlim([0,500])
ylim([-4,6])
L=legend('Observations', 'True States', 'Estimated', 'location',  'north', 'orientation', 'horizontal');
FL1= findall(L,'-property','FontName');
set(FL1,'FontName','/usr/share/fonts/msttcore/cour.ttf');
set(L,'FontSize',16);
ylabel('Value','FontSize',16);
xlabel('Sample','FontSize',16);

H = 9; W = 12;
set(h,'PaperUnits','inches')
set(h,'PaperOrientation','portrait');
set(h,'PaperSize',[H,W])
set(h,'PaperPosition',[0,0,W,H])
print(h,'-dpng','-color','impacto_snr_trecho1.png');



h = figure(3);
FN = findall(h,'-property','FontName');
set(FN,'FontName','/usr/share/fonts/dejavu/DejaVuSerifCondensed.ttf');
FS = findall(h,'-property','FontSize');
set(FS,'FontSize',16);

hold on
plot([501:1000], yk(501:1000),'ko','LineWidth',2);
plot([501:1000],y(501:1000), 'r','LineWidth',4);
plot([501:1000],x_estim(501:1000),'b--','LineWidth',4);
set (gca, "fontsize", 16)
L=legend('Observations', 'True States', 'Estimated', 'location',  'north', 'orientation', 'horizontal');
FL1= findall(L,'-property','FontName');
set(FL1,'FontName','/usr/share/fonts/msttcore/cour.ttf');
set(L,'FontSize',16);
ylabel('Value','FontSize',16);
xlabel('Sample','FontSize',16);

H = 9; W = 12;
set(h,'PaperUnits','inches')
set(h,'PaperOrientation','portrait');
set(h,'PaperSize',[H,W])
set(h,'PaperPosition',[0,0,W,H])
print(h,'-dpng','-color','impacto_snr_trecho2.png');


h = figure(4);
FN = findall(h,'-property','FontName');
set(FN,'FontName','/usr/share/fonts/dejavu/DejaVuSerifCondensed.ttf');
FS = findall(h,'-property','FontSize');
set(FS,'FontSize',16);

hold on
plot([2001:2500], yk(2001:2500),'ko','LineWidth',2);
plot([2001:2500], y(2001:2500), 'r','LineWidth',4);
plot([2001:2500], x_estim(2001:2500),'b--','LineWidth',4);
set (gca, "fontsize", 16)
L=legend('Observations', 'True States', 'Estimated', 'location',  'north', 'orientation', 'horizontal');
FL1= findall(L,'-property','FontName');
set(FL1,'FontName','/usr/share/fonts/msttcore/cour.ttf');
set(L,'FontSize',16);
ylabel('Value','FontSize',16);
xlabel('Sample','FontSize',16);

H = 9; W = 12;
set(h,'PaperUnits','inches')
set(h,'PaperOrientation','portrait');
set(h,'PaperSize',[H,W])
set(h,'PaperPosition',[0,0,W,H])
print(h,'-dpng','-color','impacto_snr_trecho3.png');


h = figure(5);
FN = findall(h,'-property','FontName');
set(FN,'FontName','/usr/share/fonts/dejavu/DejaVuSerifCondensed.ttf');
FS = findall(h,'-property','FontSize');
set(FS,'FontSize',16);

hold on
plot([251:750], yk(251:750),'ko','LineWidth',2);
plot([251:750], y(251:750), 'r','LineWidth',4);
plot([251:750], x_estim(251:750),'b--','LineWidth',4);
xlim([250,750])
set (gca, "fontsize", 16)
L=legend('Observations', 'True States', 'Estimated', 'location',  'north', 'orientation', 'horizontal');
FL1= findall(L,'-property','FontName');
set(FL1,'FontName','/usr/share/fonts/msttcore/cour.ttf');
set(L,'FontSize',16);
ylabel('Value','FontSize',16);
xlabel('Sample','FontSize',16);

H = 9; W = 12;
set(h,'PaperUnits','inches')
set(h,'PaperOrientation','portrait');
set(h,'PaperSize',[H,W])
set(h,'PaperPosition',[0,0,W,H])
print(h,'-dpng','-color','impacto_snr_trecho4.png');


h = figure(6);
FN = findall(h,'-property','FontName');
set(FN,'FontName','/usr/share/fonts/dejavu/DejaVuSerifCondensed.ttf');
FS = findall(h,'-property','FontSize');
set(FS,'FontSize',16);

hold on
plot(yk(1:final_pos),'ko','LineWidth',2);
plot(y(1:final_pos), 'r','LineWidth',4);
plot(x_estim,'b--','LineWidth',4);
xlim([0,2500])
set (gca, "fontsize", 16)
L=legend('Observations', 'True States', 'Estimated', 'location',  'north', 'orientation', 'horizontal');
FL1= findall(L,'-property','FontName');
set(FL1,'FontName','/usr/share/fonts/msttcore/cour.ttf');
set(L,'FontSize',16);
ylabel('Value','FontSize',16);
xlabel('Sample','FontSize',16);

H = 6; W = 12;
set(h,'PaperUnits','inches')
set(h,'PaperOrientation','portrait');
set(h,'PaperSize',[H,W])
set(h,'PaperPosition',[0,0,W,H])
print(h,'-dpng','-color','impacto_snr_tudo.png');

%subplot(3,1,3)
%plot(norm_wk, 'r')