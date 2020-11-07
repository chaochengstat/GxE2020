# A turtorial for implementation of the **reverse test** to test gene-environment interaction effects in R software

The **reverse test** is a statistial approach to test the interaction effect between a genetic variant and a continuous environmental exposure on a binary disease outcome.

## Required Software

R software (http://www.r-project.org/). No additional R packages are required.

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
    cat("Logistic Regression Test for GxE Interaction\n")
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
mydata=read.csv("gxe_simulate.csv")[,-1]
head(mydata)
          Z1 Z2          X G D
1 -0.5622591  1  0.2316188 1 1
2  0.6080418  1 -0.0249439 1 1
3 -0.1744747  1  1.7729779 2 1
4 -0.7125421  0  0.5858684 0 1
5 -1.9470882  1  2.2713262 1 1
6 -0.3435955  1  2.6412762 1 1
```
