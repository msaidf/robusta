# robusta
easily replace standard error with the robust one.
![Coffee](https://pixabay.com/static/uploads/photo/2015/03/05/14/09/coffee-660409_640.jpg)

In statistical analysis, it is become a standard way to report regression results using alternative standard error that is robust to heteroscedasticity and serial correlation. Sandwich package already provide a way to calculate robust standard error for various estimators. However, it is not convenient to include the produced standard error into automatically generated publication-ready table using packages such as texreg, stargazer, etc. 

RMS package has a convenient robcov function that replace the original covariance matrix with the robust one, but returning object with the same class so that it can be directly used by the automatic table generation pakcages. However this function currently handles only object created by other functions within the same package.

This package provide `robust` function that will works with many object types of standard packages (currently only work for lm, glm, and plm objects).

### Examples
```
lmo = lm(mpg ~ carb + cyl, mtcars)
lmo = robust(lmo) # replace standard error with huber-white heteroscedasticity-robust standard error
# or more succinctly using piping
lmo = lm(mpg ~ carb + cyl, mtcars) %>% robust
# now summarizing regression objects and generating summary table of results will use the robust standard error
summary(lmo)
texreg::screenreg(lmo) 
```

`robust` provides options to alternatives robust standard error HC0-HC4 that is provided by sandwich package, but using alternative naming to make it easier to remember.
```
robust(lmo, 'hc1') # the default, comparable to stata's robust standard error
robust(lmo) 
robust(lmo, 'NeweyWest') # newey-west method that is robust to heteroscedasticity and serial correlation
```

## Clustered
`cluster` replace standard error with cluster-proof standard error
```
lmo = lm(mpg ~ carb + cyl, mtcars) %>% cluster(am)
```

## Installation
```
# install.packages('devtools')
devtools::install_github('msaidf/robusta')
```

