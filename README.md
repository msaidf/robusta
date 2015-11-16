# robusta
easily replace standard error with the robust one.

In statistical analysis, it is become a standard way to report regression results using alternative standard error that is robust to heteroscedasticity and serial correlation. Sandwich package already provide a way to calculate robust standard error for various estimators. However, it is not convenient to include the produced standard error into automatically generated publication-ready table using packages such as texreg, stargazer, etc. 

RMS package has a convenient robcov function that replace the original covariance matrix with the robust one, but returning object with the same class so that it can be directly used by the automatic table generation pakcages. However this function currently handles only object created by other functions within the same package.

This package provide robcov function that will works with many object types of standard packages, such as lm, glm, plm, etc.

## Examples
```
lmo = lm(y ~ x, data)
lmo = robcov(lmo) # replace standard error with huber-white heteroscedasticity-robust standard error
texreg(lmo) # produce latex table
```

robcov provides option to alternatives standard error HC0-HC4 that is provided by sandwich package, but using easier to remember naming.
```
robcov(lmo, se='hw') # huber-white standard error, the default, equivalent to HC1 in sandwich package and is comparable to stata's robust
robcov(lmo, se='nw') # newey-west method that is robust to heteroscedasticity and serial correlation, equivalent to vcovBK
```
...

