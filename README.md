# BoostingExp
The experiments in Boosting Pathwise Coordinate Optimization in High Dimensions
The proposed ASP-Newton method is implemented in the lastest version of `picasso` on CRAN: `install.package`

## Dependency


**Install old picasso**
The old picasso implements the greedy selection rule and can be installed by
```R
install.packages("oldpicasso", repos = NULL, type="source")
```

**Install spams**
The old picasso implements the greedy selection rule and can be installed by
```R
install.packages("spams-R/spams", repos = NULL, type="source")
```

## Linear Regression Experiments (Table 2)

```R
source('linear_regression/x.R')
```

## Logistic Regression Experiments

```R
source('logistic_regression/x.R')
```

## Group Lasso

```R
source('group_lasso/x.R')
```

## Graph Estimation

```R
source('graph_estimation/x.R')
```
