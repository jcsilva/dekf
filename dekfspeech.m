%Exemplo: Filtro de Kalman

close all; clear all;

% controlled signal
randn("seed", 10);

var_ruido_proc = 0;
var_ruido_obs =0.1;
T = 1026;
INPUT_NEURONS = 10;
HIDDEN_NEURONS = 4;
OUTPUT_NEURONS = 1;
total_sinapses = INPUT_NEURONS*HIDDEN_NEURONS + HIDDEN_NEURONS + HIDDEN_NEURONS*OUTPUT_NEURONS + OUTPUT_NEURONS;
ruido_proc = sqrt(var_ruido_proc) * randn(T,1);
x = zeros(INPUT_NEURONS,1);
x(1) = ruido_proc(1);
y = zeros(1,T);
y(1) = x(1);
w = randn(total_sinapses,1);
for i=2:T
  novo_x1 = sin(2*pi*1/50*i);
  x = circshift(x, 1);
  x(1) = novo_x1 + ruido_proc(i);
  y(i) = x(1);
end
yk = y + sqrt(var_ruido_obs) * randn(1,T);
var_ruido_proc = 0.01;

%%%%%%%%%%%%%%%%%%

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