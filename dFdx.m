function [J,z] = dFdx(x, w, n_ipt, n_hdn, n_out)
  %F = [f(x,w), x(k), x(k-1), ..., x(k-M+2)]'
  %f(x,w) = w2 * tanh(w1 * x + b1) + b2 (neural net with one hidden layer)

  %d=zeros(length(x0),length(x0)); % just the first line will be different from zero

  % x = column vector (neural net input)
  % w = column vector, with w = [w1; b1; w2; b2], where:
  % w1 = 1st layer weights,
  % b1 = 1st layer bias weights,
  % w2 = 2nd layer weights
  % b2 = 2nd layer bias weights
  %%%%%%%%%%%%%
  %n_weights_l1 = n_ipt * n_hdn;
  %n_weights_l2 = n_hdn * n_out;

  %w1 = w(1:n_weights_l1);
  %b1 = w(n_weights_l1+1:n_weights_l1+n_hdn);
  %W1 = reshape(w1, n_hdn, n_ipt);

  %w2 = w(n_weights_l1+n_hdn+1:n_weights_l1+n_hdn+n_weights_l2);
  %b2 = w(n_weights_l1+n_hdn+n_weights_l2+1:end);
  %W2 = reshape(w2, n_out, n_hdn);

  %K = sech(W1 * x0 + b1).^2;
  %for i=1:n_ipt
  %  d(1,i) = W2 * (W1(:,i).*K);
  %  if i < n_ipt
  %    d(i+1,i) = 1.0;
  %  end
  %end

  n_weights_l1 = n_ipt * n_hdn;
  n_weights_l2 = n_hdn * n_out;

  w1 = w(1:n_weights_l1);
  b1 = w(n_weights_l1+1:n_weights_l1+n_hdn);
  W1 = reshape(w1, n_hdn, n_ipt);

  w2 = w(n_weights_l1+n_hdn+1:n_weights_l1+n_hdn+n_weights_l2);
  b2 = w(n_weights_l1+n_hdn+n_weights_l2+1:end);
  W2 = reshape(w2, n_out, n_hdn);

  f=@(x)[W2 * tanh(W1 * x + b1) + b2; x(1);x(2);x(3);x(4);x(5);x(6);x(7);x(8);x(9)];
  [J,z]=jacobiancsd(f,x);