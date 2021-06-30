# A turtorial for implementing the **reverse test** to test gene-environment interaction effects in R software

The **reverse test** is a statistial approach to test the interaction effect between a genetic variant and a continuous environmental exposure on a binary disease outcome.

## Required Software

`R` software (http://www.r-project.org/) and `skedastic` package (https://cran.r-project.org/web/packages/skedastic/index.html).

## R functions

```{r}
gxe.test=function(outcome="Y",exposure="X",gene="G",confounders=NULL,data,method="reverse") {
  if (is.null(confounders)) {
    confounders.formula=""
  } else {
    confounders.formula=paste("+",paste(confounders,collapse="+"))
  }
  if (method=="logistic") {
    formula=paste(outcome,"~",exposure,"+",gene,"+I(",gene,"*",exposure,")",confounders.formula)
    m=glm(formula,data=data,family=binomial(link="logit"))
  } else if (method=="reverse") {
    formula=paste(exposure,"~",outcome,"+",gene,"+I(",gene,"*",outcome,")",confounders.formula)
    m=lm(formula,data=data)
  }
  m.summary=summary(m)
  test=m.summary$coefficients[4,c(1,2)]
  Chi2.stat=test[1]^2/test[2]^2;names(Chi2.stat)=NULL
  P.value=1-pchisq(Chi2.stat,df=1);names(P.value)=NULL
  output=list(Chi2.stat=Chi2.stat,P.value=P.value,method=method,
              formula=as.formula(formula),dataset=deparse(substitute(mydata)),model=m)
  class(output)="gxetest"
  output
}

print.gxetest=function(m) {
  if (m$method=="reverse") {
    cat("    Reverse Test for GxE Interaction\n")
    model="ordinary linear regression"
  } else if (m$method=="logistic") {
    cat("    Logistic Regression Test for GxE Interaction\n")
    model="logistic regression"
  }
  cat(paste("Dataset: ",m$dataset),"\n")
  cat("Model: ",model,"\n")
  cat(paste("Formula: ",Reduce(paste, deparse(m$formula))),"\n")
  cat("Chi-squared statistic: ",m$Chi2.stat," P-value: ",m$P.value,"\n")
}
```

The gene-environment test can be called with
```
gxe.test(outcome=,exposure=,gene=,confounders=,data=,method=)
```
The function has 6 arguments. These are
* `outcome=` (Required) Name of the outcome variable, which should be a binary datatype.
* `exposure=` (Required) Name of the exposure variable, which should be a continuous datatype.
* `gene=` (Required) Name of the genetic variant, which can be binary or ordinal.
* `confounders=` (Optional) A vector of names showing the confounding variables used in the gene-environment regression model. The default value is `NULL`, which represents no confounding variables. 
* `data=` (Required) The name of the dataset, which should be a `data.frame`.
* `method=` (Optional) `"reverse"` for the reverse test, `"logistic"` for the standard logistic regression test.

## Illustrative Example

We now use an example to illustrate the usage of the `gxe.test` function. First, we read the `gxe_simulate.csv` dataset
```
> mydata=read.csv("gxe_simulate.csv")[,-1]
> head(mydata)
          Z1 Z2          X G D
1 -0.5622591  1  0.2316188 1 1
2  0.6080418  1 -0.0249439 1 1
3 -0.1744747  1  1.7729779 2 1
4 -0.7125421  0  0.5858684 0 1
5 -1.9470882  1  2.2713262 1 1
6 -0.3435955  1  2.6412762 1 1
```
This simulated dataset can be found in this github folder, which has 2000 observations and 5 variables, including a continuous exposure (`X`), a binary outcome (`D`), an ordinal genetic variate (`G`), as well as two confounding variables (`Z1` and `Z2`). 

Now we conduct a reverse examine to whether the genetic variant `G` modifies `X`-`D` association, adjusting for both of confounders. See the code below
```
> t.reverse=gxe.test(outcome="D",
                     exposure="X",
                     gene="G",
                     confounders=c("Z1","Z2"),
                     data=mydata,
                     method="reverse")
> t.reverse
    Reverse Test for GxE Interaction
Dataset:  mydata 
Model:  ordinary linear regression 
Formula:  X ~ D + G + I(G * D) + Z1 + Z2 
Chi-squared statistic:  6.984665  P-value:  0.008221106
```
It is exhibited that the P-value by the reverse test is `0.008`, so we can reject the null hypothesis that there is no gene-environment interaction effect under a 5% significance level. If the error term in the linear regression model used in reverse test follows a constant variance normal distribution, we can calculate the Ratio of Odds Ratio (ROR), which is a metric representing the magnitude of the gene-environment interaction effect. Now we applied Shapiro-Wilk test and White test to check for the normality and homoscedasticity requirements, respectively. At first, we conduct a Shapiro-Wilk test for the residuals of the linear regression model.
```
> shapiro.test(t.reverse$model$residuals)
	Shapiro-Wilk normality test

data:  t.reverse$model$residuals
W = 0.99929, p-value = 0.6658
```
The P-value is larger than 5%, so we assume the normality requirement holds. Second, we conduct a White test to check for the homoscedasticity assumption.
```
> skedastic::white_lm(t.reverse$model,statonly=F)
# A tibble: 1 x 5
  statistic p.value parameter method       alternative
      <dbl>   <dbl>     <dbl> <chr>        <chr>      
1      16.1  0.0965        10 White's Test greater    
```
The P-value is 9.7%, which is larger than 5%, so we assume the homoscedasticity requirement hold. Now we calculate the ROR based on the reverse approach
```
> exp(t.reverse$model$coefficients[4]/var(t.reverse$model$residuals))
[1] 1.185986
```
The ROR obtained by the reverse test is 1.186. 

The `gxe.test` function can also implement the standard logistic regression test. Here is an example.
```
> t.logistic=gxe.test(outcome="D",
                      exposure="X",
                      gene="G",
                      confounders=c("Z1","Z2"),
                      data=mydata,
                      method="logistic")
> t.logistic
Logistic Regression Test for GxE Interaction
Dataset:  mydata 
Model:  logistic regression 
Formula:  D ~ X + G + I(G * X) + Z1 + Z2 
Chi-squared statistic:  7.597554  P-value:  0.005844752 
```
We observed that P-value obtained by the logistic regression test is pretty similar to it by the reverse test. The ROR obtained by the logistic regression approach is
```
> exp(t.logistic$model$coefficients[4])
[1] 1.218934
```
This is also very close to the ROR obtained by the reverse approach (1.186).

## Reference(s)

Chao Cheng, Donna Spiegelman, Zuoheng Wang, and Molin Wang. Testing Gene-Environment Interactions in the Presence of Confounders and Mismeasured Environmental
Exposures. *G3: Genes, Genomes, Genetics* (To appear).
