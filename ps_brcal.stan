data {
  int<lower=0> N;
  int<lower=0> D;
  real w_alpha;
  vector[N] A;
  matrix[N,D] X;
  vector[N] R;
}

parameters {
  real alpha0;
  vector[D] alpha;
  real<lower=0,upper=100> lambda;
}

transformed parameters{
  vector[N] psA;
  psA = inv(1+exp(-(rep_vector(alpha0, N) + X*alpha)));
}

model {
  lambda ~ gamma(0.1,0.1);
  alpha ~ double_exponential(0,1/sqrt(lambda));
  target += -w_alpha * ( A'*exp(-log(psA./(1-psA))) + (1-A)'*log(psA./(1-psA)) + (1-A)'*exp(log(psA./(1-psA))) - A'*log(psA./(1-psA)) );
}
