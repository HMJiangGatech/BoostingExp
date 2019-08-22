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
source('linear_regression.R')
```

**Precision Setting**
We make every algorithm achieves the approximate KKT error around $0.001$ by setting the precision as follows for every experiments.

| Data |n  |p  |Ours | Greedy | Cyclic | APG  | Glmnet | Gcdnet | Ncvreg |
|------|---|---|-----|--------|--------|------|--------|--------|--------|
|Simulate Well-Cond.   |  100  | 1000   |   |   |   |   |   |   |   |
|Simulate Well-Cond.   |  100  | 10000  |   |   |   |   |   |   |   |
|Simulate Well-Cond.   |  1000 | 1000   |   |   |   |   |   |   |   |
|Simulate Well-Cond.   |  1000 | 10000  |   |   |   |   |   |   |   |
|Simulate Well-Cond.   |  5000 | 10000  |   |   |   |   |   |   |   |
|Simulate Well-Cond.   |  10000| 1000   |   |   |   |   |   |   |   |
|Simulate Ill-Cond.    |  100  | 1000   |   |   |   |   |   |   |   |
|Simulate Ill-Cond.    |  100  | 10000  |   |   |   |   |   |   |   |
|Simulate Ill-Cond.    |  1000 | 1000   |   |   |   |   |   |   |   |
|Simulate Ill-Cond.    |  1000 | 10000  |   |   |   |   |   |   |   |
|Simulate Ill-Cond.    |  5000 | 10000  |   |   |   |   |   |   |   |
|Simulate Ill-Cond.    |  10000| 1000   |   |   |   |   |   |   |   |

## Logistic Regression Experiments

```R
source('logistic_regression.R')
```

## Poisson Regression Experiments

```R
source('poisson_regression.R')
```

## Group Lasso

```R
source('group_lasso.R')
```

## Graph Estimation

```R
source('graph_estimation.R')
```
