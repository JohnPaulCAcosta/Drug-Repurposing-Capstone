###### HW 3 STAT 445 #######
### Q1
  ### part a
  set.seed(123)
  n <- 200
  U <- rnorm(n, mean = 0, sd = 1)   # confounder
  T <- rbinom(n, 1, prob = plogis(0.5 * U))  # treatment depends on U
  Y <- 2 + 1.5 * T + 2.0 * U + rnorm(n, 0, 1) # outcome depends on T and U
  mod_a <- lm(Y ~ T)
  mod_b <- lm(Y ~ T + U)

  summary(mod_a)
    # Call:
    #   lm(formula = Y ~ T)
    # 
    # Residuals:
    #   Min      1Q  Median      3Q     Max 
    # -6.5302 -1.3850 -0.1302  1.4902  5.8306 
    # 
    # Coefficients:
    #   Estimate Std. Error t value Pr(>|t|)    
    # (Intercept)   1.7646     0.2067   8.538 3.58e-15 ***
    #   T             1.9961     0.2894   6.897 6.96e-11 ***
    #   ---
    #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    # 
    # Residual standard error: 2.046 on 198 degrees of freedom
    # Multiple R-squared:  0.1937,	Adjusted R-squared:  0.1896 
    # F-statistic: 47.57 on 1 and 198 DF,  p-value: 6.964e-11
  
  summary(mod_b)
    # Call:
    #   lm(formula = Y ~ T + U)
    # 
    # Residuals:
    #   Min       1Q   Median       3Q      Max 
    # -2.73658 -0.58856  0.03362  0.61193  2.54695 
    # 
    # Coefficients:
    #   Estimate Std. Error t value Pr(>|t|)    
    # (Intercept)  2.02058    0.10304   19.61   <2e-16 ***
    #   T            1.52619    0.14481   10.54   <2e-16 ***
    #   U            1.89696    0.07695   24.65   <2e-16 ***
    #   ---
    #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    # 
    # Residual standard error: 1.015 on 197 degrees of freedom
    # Multiple R-squared:  0.8026,	Adjusted R-squared:  0.8006 
    # F-statistic: 400.6 on 2 and 197 DF,  p-value: < 2.2e-16
  
    # The estimated coefficient for T is higher in the model without the confounder, U, than the model that includes it
  
  ### part b
  # yes, omitting U biases the estimated effect of T in the direction of increasing it's value positively
  # this is because U has a positive effect on T, so when U is omitted the model attributes some of U's positive effect to T
  
  ###part c
  coefs <- matrix(NA, nrow = 1000, ncol = 2)
  colnames(coefs) <- c("bA", "bB")
  for (i in 1:1000){
    n <- 200
    U <- rnorm(n, mean = 0, sd = 1)   # confounder
    T <- rbinom(n, 1, prob = plogis(0.5 * U))  # treatment depends on U
    Y <- 2 + 1.5 * T + 2.0 * U + rnorm(n, 0, 1) # outcome depends on T and U
    bA <- coef(lm(Y ~ T))["T"]
    bB <- coef(lm(Y ~ T + U))["T"]
    coefs[i,]<- c(bA, bB)
  }
  estA <- coefs[,"bA"]
  estB <- coefs[,"bB"]
  
  summary(estA)
  summary(estB)
  
  # histograms
  par(mfrow = c(1, 2))  # side-by-side plots
  
  hist(estA, breaks = 30, col = "lightblue", border = "white",
       main = "Model A: Y ~ T",
       xlab = "Estimated Coefficient for T")
  abline(v = 1.5, col = "red", lwd = 2)  # true effect
  
  hist(estB, breaks = 30, col = "lightgreen", border = "white",
       main = "Model B: Y ~ T + U",
       xlab = "Estimated Coefficient for T")
  abline(v = 1.5, col = "red", lwd = 2)
  
  par(mfrow = c(1,1))
  
  # Confounding causes the estimate of T to be biased, omitting U makes the effect of T appear larger than it truly is
  # When U is included in the model, the bias disappears and the estimate is centered around 1.5, the true value of T
  
  
  ### part d
  
  true_tau <- 1.5
  R <- 1000
  betaU_vals <- c(0, 1, 3)    
  gamma_vals <- c(0, 0.5, 1)  
  
  one_sim <- function(n, betaU, gamma, tau=true_tau) {
    U <- rnorm(n, 0, 1)
    T <- rbinom(n, 1, prob = plogis(gamma * U))
    Y <- 2 + tau * T + betaU * U + rnorm(n, 0, 1)
    b_A <- coef(lm(Y ~ T))["T"]         
    b_B <- coef(lm(Y ~ T + U))["T"]     #
    c(b_A=b_A, b_B=b_B)
  }
  
  grid_summaries <- lapply(betaU_vals, function(bU) {
    lapply(gamma_vals, function(gm) {
      ests <- replicate(R, one_sim(n, bU, gm), simplify = "matrix")
      rownames(ests) <- c("b_A","b_B")
      est_A <- ests["b_A", ]
      est_B <- ests["b_B", ]
      data.frame(
        betaU = bU, gamma = gm,
        meanA = mean(est_A), sdA = sd(est_A), biasA = mean(est_A) - true_tau,
        meanB = mean(est_B), sdB = sd(est_B), biasB = mean(est_B) - true_tau
      )
    })
  })
  
  tab <- do.call(rbind, unlist(grid_summaries, recursive = FALSE))
  tab[order(tab$betaU, tab$gamma), ]
  print(round(tab, 4), row.names = FALSE)
   #betaU gamma meanA  sdA    biasA  meanB  sdB     biasB
    # 0   0.0   1.5029 0.1418 0.0029 1.5030 0.1416  0.0030
    # 0   0.5   1.5044 0.1386 0.0044 1.5043 0.1425  0.0043
    # 0   1.0   1.5034 0.1423 0.0034 1.5029 0.1529  0.0029
    # 1   0.0   1.5057 0.1924 0.0057 1.5008 0.1405  0.0008
    # 1   0.5   1.9700 0.1969 0.4700 1.4964 0.1462 -0.0036
    # 1   1.0   2.3173 0.1929 0.8173 1.4885 0.1566 -0.0115
    # 3   0.0   1.5026 0.4631 0.0026 1.4978 0.1456 -0.0022
    # 3   0.5   2.9288 0.4512 1.4288 1.5006 0.1501  0.0006
    # 3   1.0   3.9916 0.4202 2.4916 1.4974 0.1539 -0.0026
  # as the strength of confounding increases (larger values of beta, the effect of U on Y, and larger values of gamma, the effect of U on T), the estimated treatment effect from the model omitting U becomes increasingly biased
  # meanwhile the model including U consistently estimates T around 1.5 regardless of confouding strength
  # When either beta or gamma is 0, both models produced unbiased estimates around the true value of 1.5
  
  
### Q2
  set.seed(123)
  n <- 200
  # Confounder: age
  age <- rnorm(n, mean = 50, sd = 10)
  # Exposure: treatment (0 = control, 1 = treated), probability depends on age
  prob_trt <- plogis(0.05 * (age- 50))
  trt <- rbinom(n, size = 1, prob = prob_trt)
  # Mediator: weight after program
  weight <- 80- 5 * trt + 0.3 * age + rnorm(n, sd = 5)
  # Outcome: blood pressure
  bp <- 130- 3 * trt + 0.5 * weight + 0.2 * age + rnorm(n, sd = 8)
  # Assemble into a data frame
  meddat <- data.frame(bp = bp, trt = trt, weight = weight, age = age)
  
  ###part a
  mod_out <-lm(bp ~ trt + weight + age, data = meddat)
  mod_med <-lm(weight ~ trt + age, data = meddat)
  summary(mod_out)
  summary(mod_med)  
  
  ###part b
  coefs_out <- coef(summary(mod_out))     
  coefs_med <- coef(summary(mod_med))
  
  beta1  <- coefs_out["trt",   "Estimate"]   # direct effect (trt in outcome model)
  beta2  <- coefs_out["weight","Estimate"]   # weight effect in outcome model
  gamma1 <- coefs_med["trt",   "Estimate"]   # trt effect in mediator model
  
  direct   <- beta1
  indirect <- gamma1 * beta2
  total    <- direct + indirect
  direct; indirect; total
  # direct: -4.005527
  # indirect: -2.06614
  # total: -6.071667
  
  ###part c

  Vy <- vcov(mod_out)
  Vm <- vcov(mod_med)
 
  se_direct <- coef(summary(mod_out))["trt", "Std. Error"]  
  se_beta2  <- coef(summary(mod_out))["weight", "Std. Error"]
  cov_b1b2  <- Vy["trt", "weight"]
  var_g1    <- Vm["trt", "trt"]
  
   se_indirect <- sqrt((gamma1^2) * (se_beta2^2) + (beta2^2) * var_g1)
  
  se_total <- sqrt((se_direct^2) + 2 * gamma1 * cov_b1b2 + (gamma1^2) * (se_beta2^2) + (beta2^2) * var_g1)
  
  ci <- function(est, se) c(L = est - 1.96 * se, U = est + 1.96 * se)
  
  ci_direct   <- ci(direct, se_direct)
  ci_indirect <- ci(indirect, se_indirect)
  ci_total    <- ci(total, se_total)
  
  ci_direct
  # L         U 
  # -6.458677 -1.552378 
  ci_indirect
  # L          U 
  # -3.2852314 -0.8470479 
  ci_total 
  # L         U 
  # -8.364489 -3.778845 
  
  # The treatment lowers blood pressure by about 6.1 units on average
  # Around 2/3 of this reduction (around -2.1) is from the indirect mechanism of weight loss, i.e. the treatment causes weight loss, and weight loss leads to lower blood pressure
  # The direct effect (about -4) is the amount of blood pressure reduction that isn't explainable by weight change
  # Overall there is clear evidence of a mediation pathway through weight
  
  ### Q3
    ###part a
    library(MASS)
    data("birthwt")
    str(birthwt)
    birthwt$low <- as.integer(birthwt$bwt < 2500)
    bwt <- subset(birthwt, select = c(low, age, lwt, smoke, ht, ui))
    summary(bwt)
    
    ### part b
    # low              age             lwt       
    # Min.   :0.0000   Min.   :14.00   Min.   : 80.0  
    # 1st Qu.:0.0000   1st Qu.:19.00   1st Qu.:110.0  
    # Median :0.0000   Median :23.00   Median :121.0  
    # Mean   :0.3122   Mean   :23.24   Mean   :129.8  
    # 3rd Qu.:1.0000   3rd Qu.:26.00   3rd Qu.:140.0  
    # Max.   :1.0000   Max.   :45.00   Max.   :250.0  
    # smoke              ht                ui        
    # Min.   :0.0000   Min.   :0.00000   Min.   :0.0000  
    # 1st Qu.:0.0000   1st Qu.:0.00000   1st Qu.:0.0000  
    # Median :0.0000   Median :0.00000   Median :0.0000  
    # Mean   :0.3915   Mean   :0.06349   Mean   :0.1481  
    # 3rd Qu.:1.0000   3rd Qu.:0.00000   3rd Qu.:0.0000  
    # Max.   :1.0000   Max.   :1.00000   Max.   :1.0000  
    # 
    summary(bwt[, c("age","lwt")])
    table_smoke <- table(Smoke = bwt$smoke)
    table_ht    <- table(HT = bwt$ht)
    table_ui    <- table(UI = bwt$ui)
    
    
    by_age <- tapply(bwt$age, bwt$low, summary)
    by_lwt <- tapply(bwt$lwt, bwt$low, summary)
    tab_smoke_by_low <- with(bwt, table(low, smoke))
    tab_ht_by_low    <- with(bwt, table(low, ht))
    tab_ui_by_low    <- with(bwt, table(low, ui))
    
    by_age; by_lwt
    # $`0`
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    # 14.00   19.00   23.00   23.66   28.00   45.00 
    # 
    # $`1`
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    # 14.00   19.50   22.00   22.31   25.00   34.00 
    # 
    # $`0`
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    # 85.0   113.0   123.5   133.3   147.0   250.0 
    # 
    # $`1`
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    # 80.0   104.0   120.0   122.1   130.0   200.0 
    tab_smoke_by_low; tab_ht_by_low; tab_ui_by_low
    # smoke
 # low  0  1
    # 0 86 44
    # 1 29 30
    # ht
 # low   0   1
    # 0 125   5
    # 1  52   7
    # ui
# low   0   1
    # 0 116  14
    # 1  45  14
    par(mfrow=c(1,2))
    boxplot(age ~ low, data=bwt, names=c("Not low","Low"), main="Age by Outcome",
            ylab="Age")
    boxplot(lwt ~ low, data=bwt, names=c("Not low","Low"), main="LWT by Outcome",
            ylab="Weight at LMP")
    par(mfrow=c(1,1))
    
    op <- par(mfrow=c(1,3))
    barplot(tab_smoke_by_low, beside=TRUE, legend=TRUE,
            main="Smoking by Outcome", xlab="low (0/1)")
    barplot(tab_ht_by_low,    beside=TRUE, legend=TRUE,
            main="HT by Outcome", xlab="low (0/1)")
    barplot(tab_ui_by_low,    beside=TRUE, legend=TRUE,
            main="UI by Outcome", xlab="low (0/1)")
    par(op)
    
    # Mothers who are lighter, younger, and who smoke or have hypertension/ui seem more likely to have low weight babies
    
    ###part c
    fit_log <- glm(low ~ age + lwt + smoke + ht + ui, data=bwt, family=binomial)
    summary(fit_log)
    
    co <- coef(summary(fit_log))
    OR  <- exp(coef(fit_log))
    SE  <- co[, "Std. Error"]
    coef_log<-co[, "Estimate"]
    
    coef_log
    # (Intercept)         age         lwt       smoke 
    # 1.39979416 -0.03407314 -0.01544710  0.64753972 
    # ht          ui 
    # 1.89327417  0.88460678 
    SE
    # (Intercept)         age         lwt       smoke 
    # 1.080407199 0.033673932 0.006586788 0.336650120 
    # ht          ui 
    # 0.683392499 0.444051332 
    OR
    # (Intercept)         age         lwt       smoke 
    # 4.0543653   0.9665008   0.9846716   1.9108339 
    # ht          ui 
    # 6.6410771   2.4220318 
 
    ###part d
    ci_prof <- confint(fit_log)          # may take a moment
    OR_ci   <- exp(ci_prof)
    OR_ci[c("smoke","ht"), , drop=FALSE]
    # 2.5 %    97.5 %
    #   smoke 0.9879944  3.713646
    # ht    1.8018077 27.783502
    # Smokers have about twice the odds of a low-birth-weight baby (not as significant), while mothers with hypertension have roughly seven times the odds
   
    ###part e
    fit_q <- glm(low ~ age + lwt + smoke + ht + ui, data=bwt, family=quasibinomial)
    summary(fit_q)
    
    dispersion <- deviance(fit_q) / df.residual(fit_q)
    dispersion
    #1.157256, no evidence of overdispersion
    
    ### part f
    
    newdat <- data.frame(
      age  = c(25, 40),
      lwt  = c(120, 100),
      smoke= c(0, 1),
      ht   = c(1, 0),
      ui   = c(0, 1)
    )
    pred <- predict(fit_q, newdata=newdat, type="response")
    cbind(newdat, p_low = pred)
    # age lwt smoke ht ui     p_low
    # 1  25 120     0  1  0 0.6428115
    # 2  40 100     1  0  1 0.5060760
    #The younger mother with hypertension has a higher predicted risk (64%) of low birth weight than the older smoker with uterine irritability (51%)
    
### q4
    ###part a
    library(vcd)
    data("Arthritis")
    Arthritis$Improved2 <- ifelse(Arthritis$Improved == "Marked", "Yes",
                                  ifelse(Arthritis$Improved %in% c("None","Some"), "No", NA))
    Arthritis$Improved2 <- factor(Arthritis$Improved2, levels=c("No","Yes"))
    
    tab <- with(Arthritis, table(Treatment, Improved2))
    tab

    
    ###part b
    # Improved2
    # Treatment No Yes
    # Placebo 36   7
    # Treated 20  21
    prop.table(tab, margin=1) 
    # Improved2
    # Treatment        No       Yes
    # Placebo 0.8372093 0.1627907
    # Treated 0.4878049 0.5121951
    # barplot
    barplot(tab, beside=TRUE, legend=TRUE, main="Improvement by Treatment",
            ylab="Count")
    #Improvement occurred in 51% of treated patients versus 16% on placebo, suggesting the treatment was much more effective
    
    ###part c
    chisq <- chisq.test(tab, correct=FALSE) 
    chisq
  #     Pearson's Chi-squared test
  # 
  # data:  tab
  # X-squared = 11.53, df = 1, p-value = 0.0006847
    #The chi-square test is significant (p = 0.0006847), indicating that improvement is strongly associated with receiving treatment.
    chisq$expected
    # Improved2
    # Treatment       No      Yes
    # Placebo 28.66667 14.33333
    # Treated 27.33333 13.66667
    # expected count assumption good
    
    ### part d
    ## ---- Q4(d): Fisher's Exact test ----
    fisher <- fisher.test(tab, alternative="two.sided")
    fisher
    #     Fisher's Exact Test for Count Data
    # 
    # data:  tab
    # p-value = 0.001041
    # alternative hypothesis: true odds ratio is not equal to 1
    # 95 percent confidence interval:
    #   1.783439 17.445146
    # sample estimates:
    # odds ratio 
    #   5.284255 
    # Fisher’s test (p = 0.001) agrees with the Chi-square result, showing a strong treatment effect; Fisher’s test is probably preferred due to small sample size
    
    ###part e

    p_treated <- tab["Treated","Yes"] / sum(tab["Treated", ])
    p_placebo <- tab["Placebo","Yes"] / sum(tab["Placebo", ])
    
    odds_treated <- p_treated / (1 - p_treated)
    odds_placebo <- p_placebo / (1 - p_placebo)
    
    OR <- odds_treated / odds_placebo
    OR
    #5.4
   
    fisher$estimate        
    #odds ratio 
    #5.284255 
    fisher$conf.int        
    #1.783439, 17.445146
    #treated patients had about 5 times higher odds of improvement than those on placebo
   # The 95% CI (1.78 – 17.45) does not include 1, indicating a statistically significant treatment effect
    
    ###part f
    risk_treated <- tab["Treated","Yes"] / sum(tab["Treated", ])
    risk_placebo <- tab["Placebo","Yes"] / sum(tab["Placebo", ])
    RR <- risk_treated / risk_placebo
    RR
    #3.146341
    logRR <- log(RR)
    SE_logRR <- sqrt( (1/tab["Treated","Yes"]) - (1/sum(tab["Treated", ])) +
                        (1/tab["Placebo","Yes"]) - (1/sum(tab["Placebo", ])) )
    CI_RR <- exp(logRR + c(-1,1) * 1.96 * SE_logRR)
    CI_RR
    # 1.500052 6.599413
    #Treated patients are about 3 times as likely to improve as those on placebo.
   # Because the CI does not include 1, there is strong evidence that the treatment increases the probability of improvement.
    ###part g
    # The odds ratio compares odds (p / (1 − p)), while the risk ratio compares actual probabilities (p).
    # Both show the same direction of effect, but the risk ratio is easier to interpret — it tells us treated patients are 3 times more likely to improve rather than 5 times higher odds, which is more intuitive for most people
    # 
    