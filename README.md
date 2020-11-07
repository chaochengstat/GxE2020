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
