# Synthetic demos (base R only)

library(mpssAIC)

set.seed(1)
n <- 150; p <- 7
X <- as.data.frame(matrix(rnorm(n*p), n, p)); names(X) <- paste0("x",1:p)
beta <- c(1.5, -1.2, 0.8, rep(0, p-3))
y <- as.numeric(as.matrix(X) %*% c(beta, rep(0, p-length(beta))) + rnorm(n, sd=1))

forest <- build_paths(X, y, family="gaussian", K=5, eps=1e-3, delta=2, L=50)
stab <- stability(X, y, B=20, resample="bootstrap",
                  build_args=list(family="gaussian", K=5, eps=1e-3, delta=2, L=50))
plaus <- plausible_models(forest, pi=stab$pi, Delta=2, tau=0.4)
print(plaus)

# Logistic
set.seed(2)
n <- 200; p <- 6
X <- as.data.frame(matrix(rnorm(n*p), n, p)); names(X) <- paste0("x",1:p)
lin <- -0.5 + 1.1*X$x1 - 0.9*X$x2 + 0.7*X$x5
prob <- 1/(1+exp(-lin))
y <- rbinom(n, 1, prob)
forest <- build_paths(X, y, family="binomial", K=4, eps=1e-3, delta=2, L=40)
stab <- stability(X, y, B=20, resample="bootstrap",
                  build_args=list(family="binomial", K=4, eps=1e-3, delta=2, L=40))
plaus <- plausible_models(forest, pi=stab$pi, Delta=2, tau=0.45, attach_fits=TRUE, X=X, y=y)
print(plaus)
if (nrow(plaus) > 0) {
  fit <- plaus$fit[[1]]
  pr <- predict(fit, type="response")
  print(confusion_metrics(y, pr))
}
