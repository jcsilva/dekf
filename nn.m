function [nn_out, hdn_out] = nn(x,w,n_ipt, n_hdn, n_out)
  % x = column vector (neural net input)
  % w = column vector, with w = [w1; b1; w2; b2], where:
  % w1 = 1st layer weights,
  % b1 = 1st layer bias weights,
  % w2 = 2nd layer weights
  % b2 = 2nd layer bias weights

  if (size(x,2) ~= 1)
    x = x'; % column vector
  end

  n_weights_l1 = n_ipt * n_hdn;
  n_weights_l2 = n_hdn * n_out;

  w1 = w(1:n_weights_l1);
  b1 = w(n_weights_l1+1:n_weights_l1+n_hdn);
  W1 = reshape(w1, n_hdn, n_ipt);

  w2 = w(n_weights_l1+n_hdn+1:n_weights_l1+n_hdn+n_weights_l2);
  b2 = w(n_weights_l1+n_hdn+n_weights_l2+1:end);
  W2 = reshape(w2, n_out, n_hdn);

  hdn_out = tanh(W1 * x + b1);

  nn_out = W2 * hdn_out + b2;
