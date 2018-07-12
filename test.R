library(data.table)
library(doParallel)
library(lfe)
library(stargazer)

dt.test <- data.table(a1 = sample(c(0,0,1), 1000, replace=TRUE),
                      a2 = sample(c(0,0,1), 1000, replace=TRUE),
                      b1 = rnorm(1000),
                      b2 = rnorm(1000),
                      c1 = rnorm(1000),
                      c2 = rnorm(1000),
                      f1 = sample(1:4, 1000, replace=TRUE),
                      f2 = sample(1:4, 1000, replace=TRUE),
                      e1 = rnorm(1000),
                      e2 = rnorm(1000),
                      e3 = rnorm(1000),
                      e4 = rnorm(1000),
                      i1 = rnorm(1000),
                      i2 = rnorm(1000),
                      i3 = rnorm(1000),
                      i4 = rnorm(1000),
                      cl1 = sample(1:200, 1000, replace=TRUE),
                      cl2 = sample(1:200, 1000, replace=TRUE),
                      w1 = sample(1:4, 1000, replace=TRUE),
                      w2 = sample(1:4, 1000, replace=TRUE)
                      )

aa <- reg.combs(a1 + a2 ~ b1 + b2 | c1 + c2 ,
          fe = ~ (f1 + f2),
#          iv = ~
#            (e1 | e2 ~ i1 + i2) +
#            (e3 | e4 ~ i3 + i4),
          cl = ~ (cl1 + cl2),
          w = ~ w1 + w2,
          data = dt.test[1:500] + dt.test[501:1000] ~ .,
          reg.fn = ~ felm + logit,
          include.all = TRUE, n.cores=4)


reg.combs.models(a1 ~ b1 | c1 + c2 ,
                cl = ~ cl1,
                w = ~ w2,
                data = dt.test[1:500] + dt.test[501:1000] ~ .,
                reg.fn = ~ felm + logit, include.all = FALSE, controls.alone = FALSE, fe = ~ 0, iv = ~ 0 )

aa <- reg.combs(a1 ~ b1 + b2 | (c1 + c2) ,
                 cl = ~ (cl1 + cl2),
                 w = ~ w2,
                data =  ~ dt.test[c1 > 0] + dt.test[c1 < 0],
                reg.fn = ~ felm + logit,
                include.all = FALSE,
                controls.alone = FALSE,
                fe = ~ 0,
                iv = ~ 0)

## Test not having explicit controls
reg.combs.models(a1 + a2 + a3 + a4 ~ b1 + b2 | c1 ,
                 w = w1 + w2 ~ w3 + w4 | w5 + w6,
                 reg.fn = felm + logit ~ .,
          data = d1 + d2 ~ .)

reg.combs.models(a1 ~ b1 + b2 ,
          include.all = TRUE,
          data = ~ dt.test[1:500])

reg.combs.models(a1 ~ 1 + b1 + b2 + (b1 + b2) | (c1 + c2 + c3) + 1,
                cl = ~ (cl1 + cl2),
                w = ~ w2,
                data =  ~ dt.test[c1 > 0],
                reg.fn = ~ felm,
                fe = ~ 0,
                iv = ~ 0)

reg.combs(a1 ~ 1 + b1 + b2 + (b1 + b2) ,
                 cl = ~ (cl1 + cl2),
                 w = ~ w2,
                 data =  ~ dt.test[c1 > 0],
                 reg.fn = ~ felm,
                 fe = ~ 0,
                 iv = ~ 0, test=TRUE)
