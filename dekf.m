function [xk_estim_pos, Pxk_pos, wk_estim_pos, Pwk_pos] = dekf(yk, var_ruido_proc, var_ruido_obs, xk_estim_pos_ini, Pxk_pos_ini, wk_estim_pos_ini, Pwk_pos_ini)

  if (size(yk,1) ~= 1)
    yk = yk'; %yk must be a line vector
  end

  %Cenario
  INPUT_NEURONS = 10;
  HIDDEN_NEURONS = 4;
  OUTPUT_NEURONS = 1;
  total_sinapses = INPUT_NEURONS*HIDDEN_NEURONS + HIDDEN_NEURONS + HIDDEN_NEURONS*OUTPUT_NEURONS + OUTPUT_NEURONS;

  %numero de estados do sistema dinamico
  K = INPUT_NEURONS;

  %numero de observacoes
  T = length(yk);

  %Dual Extended Kalman Filter
  C = zeros(1,K);
  C(1) = 1;
  B = C';

  %initialization
  xk_estim_pos = zeros(K,1,T);
  Pxk_pos = zeros(K,K,T);
  Pxk_pri = zeros(K,K,T);
  xk_estim_pri = zeros(K,1,T);
  Pxkseq = zeros(K,K,2*T+1);

  wk_estim_pos = zeros(total_sinapses,1,T);
  Pwk_pos = zeros(total_sinapses,total_sinapses,T);
  Pwk_pri = zeros(total_sinapses,total_sinapses,T);
  wk_estim_pri = zeros(total_sinapses,1,T);
  Pwkseq = zeros(total_sinapses,total_sinapses,2*T+1);

  ind = 2;

  %initial estimates
  if (nargin == 3)
    xk_estim_pos(:,:,1) = zeros(K,1);
    Pxk_pos(:,:,1) = 0*eye(K,K);
    wk_estim_pos(:,:,1) = 0.01*randn(total_sinapses,1);
    Pwk_pos(:,:,1) = 0.01*eye(total_sinapses,total_sinapses);;
  else
    xk_estim_pos(:,:,1) = xk_estim_pos_ini;
    Pxk_pos(:,:,1) = Pxk_pos_ini;
    wk_estim_pos(:,:,1) = wk_estim_pos_ini;
    Pwk_pos(:,:,1) = Pwk_pos_ini;
  end
  Pxkseq(:,:,1) = Pxk_pos(:,:,1);
  Pwkseq(:,:,1) = Pwk_pos(:,:,1);

  %Kalman gain
  Kxk = zeros(K,1,T);
  Kwk = zeros(total_sinapses,1,T);

  %forgetting factor
  lambda = 0.9997;

  trx = zeros(T,1);
  trw = zeros(T,1);
  aux = zeros(T,1);

  for i=2:T,

    %time-update equations for the weight filter:
    wk_estim_pri(:,:,i) = wk_estim_pos(:,:,i-1);
    Pwk_pri(:,:,i) = (1.0/lambda)*Pwk_pos(:,:,i-1);

    %time-update equations for the state filter:
    xk_estim_pri(:,:,i) = F(xk_estim_pos(:,:,i-1), wk_estim_pri(:,:,i), INPUT_NEURONS, HIDDEN_NEURONS, OUTPUT_NEURONS);
    A = dFdx(xk_estim_pos(:,:,i-1), wk_estim_pri(:,:,i), INPUT_NEURONS, HIDDEN_NEURONS, OUTPUT_NEURONS);
    Pxk_pri(:,:,i) = A * Pxk_pos(:,:,i-1) * A' + B * var_ruido_proc * B';

    %measurement equations for the state filter:
    Kxk(:,:,i) = (Pxk_pri(:,:,i) * C')/(C * Pxk_pri(:,:,i) * C' + var_ruido_obs);
    xk_estim_pos(:,:,i) = xk_estim_pri(:,:,i) + Kxk(:,:,i) * (yk(:,i) - C*xk_estim_pri(:,:,i));
    Pxk_pos(:,:,i) = (eye(K) - Kxk(:,:,i) * C) * Pxk_pri(:,:,i);

    %measurement equations for the weight filter:
    Cw = C * dFdw(xk_estim_pri(:,:,i), wk_estim_pos(:,:,i-1), INPUT_NEURONS, HIDDEN_NEURONS, OUTPUT_NEURONS);
    Kwk(:,:,i) = (Pwk_pri(:,:,i) * Cw')/(Cw * Pwk_pri(:,:,i) * Cw' + var_ruido_proc + var_ruido_obs );
    wk_estim_pos(:,:,i) = wk_estim_pri(:,:,i) + Kwk(:,:,i) * (yk(:,i) - C*xk_estim_pri(:,:,i));
    Pwk_pos(:,:,i) = (eye(total_sinapses) - Kwk(:,:,i) * Cw) * Pwk_pri(:,:,i);
    aux(i) = trace(Pwk_pos(:,:,i));

    %grava a sequencia das matrizes de covariancia
    Pxkseq(:,:,ind) = Pxk_pri(:,:,i);
    Pxkseq(:,:,ind+1) = Pxk_pos(:,:,i);
    trx(i) = trace(Pxk_pos(:,:,i));
    Pwkseq(:,:,ind) = Pwk_pri(:,:,i);
    Pwkseq(:,:,ind+1) = Pwk_pos(:,:,i);
    trw(i) = trace(Pwk_pos(:,:,i));
    ind = ind + 2;

  end
