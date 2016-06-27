%Exemplo: Filtro de Kalman

close all; clear all;

% controlled signal
randn("seed", 10);
var_ruido_proc = 0.01;
var_ruido_obs =0.1;
T = [1:1026];
y = sin(2*pi*1/50*T);
yk = y + sqrt(var_ruido_obs) * randn(1,length(y));

%%%%%%%%%%%%%%%%%%
%DEKF

shift = 64;
window = 256;
current = 1;
x_estim = zeros(size(y));
while (current + window < length(yk))
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
hold on
plot(yk(1:final_pos),'ko');
plot(y(1:final_pos), 'r');
grid;
plot(x_estim,'b');
legend('Observations', 'True States', 'Estimated');

squared_error = (y(1:final_pos) - x_estim).^2;
MSE = mean(squared_error)