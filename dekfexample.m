%Exemplo: Filtro de Kalman

close all; clear all;

% controlled signal
randn("seed", 10);
var_ruido_proc = 0.01;
var_ruido_obs =0.05;
T1 = [1:500];
T2 = [501:1501];
y1 = sin(2*pi*1/50*T1);
y2 = sin(2*pi*1/100*T2);
y = [y1, y2];
yk = y + sqrt(var_ruido_obs) * randn(1,length(y));

%%%%%%%%%%%%%%%%%%
%DEKF

shift = 1500;
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
% Results: states and theirs estimatives
final_pos = length(x_estim);
figure;
subplot(2,1,1)
hold on
plot(yk(1:final_pos),'ko');
plot(y(1:final_pos), 'r');
grid;
plot(x_estim,'b');
legend('Observations', 'True States', 'Estimated');

squared_error = (y(1:final_pos) - x_estim).^2;
MSE = mean(squared_error)

norm_diff_wk = zeros(size(wk_estim_pos,3) - 1,1);
norm_wk = zeros(size(wk_estim_pos,3),1);
for i=2:size(wk_estim_pos,3)
  diff_wk = wk_estim_pos(:,1,i) - wk_estim_pos(:,1,i-1);
  norm_diff_wk(i-1) = sqrt(diff_wk' * diff_wk);
  norm_wk(i) = sqrt(wk_estim_pos(:,1,i)' * wk_estim_pos(:,1,i));
end

subplot(2,1,2)
hold on
plot(norm_diff_wk);
%plot(norm_wk, 'r')