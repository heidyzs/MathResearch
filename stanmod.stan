
// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  vector[N] y;
  real<lower=0> sigma;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real<lower=0,upper=100> mu;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  for (n in 1:N) {
    y[n] ~ normal(mu, sigma)T[0,100];
  }
  mu ~ skew_normal(86,22,-2.2)T[0,100];
}

