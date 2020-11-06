# A turtorial for implementation of the **reverse test** to test gene-environment interaction effects in R software

The **reverse test** is a statistial approach to test the interaction effect between a genetic variant and a continuous environmental exposure on a binary disease outcome.

## Required Software

R software (http://www.r-project.org/). No additional R packages are required.

## R functions

```{r}
gxe.test=function(outcome="Y",exposure="X",gene="G",confounders=NULL,data,method="reverse") {

if (method=="reverse") {
formula=paste()
} else if (method=="logistic") {

}
}
```
