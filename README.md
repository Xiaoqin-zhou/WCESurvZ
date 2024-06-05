# WCESurvZ: Weighted Composite Endpoint Survival Analysis Toolkit in R
WCESurvZ is an R package developed for weighted composite endpoint survival analysis, aiming to address the statistical analysis challenges of composite endpoints in clinical trials. The package provides a set of tools for weighted Kaplan-Meier analysis, weighted log-rank tests, and professional survival curve plotting.
## Introduction
WCESurvZ is an R package designed for weighted composite endpoint survival analysis. It provides a flexible, efficient, and professional toolkit for handling clinical trial data with composite endpoints and recurrent events. The package supports event weighting analysis, which is crucial for evaluating the clinical importance and occurrence frequency of different endpoints in a composite endpoint.
## Features
Event Weighting Analysis: WCESurvZ supports event weighting analysis, which allows assigning different weights to different events based on their clinical significance. This feature is not available in the classic survival package.
Scientific and Reasonable Testing Methods: WCESurvZ employs the weighted log-rank test with linear weight correction, which is suitable for survival analysis data. The weighted log-rank test provides more reliable and accurate results for evaluating group differences.
Professional Visualization: Built on top of the ggplot2 package, WCESurvZ generates high-quality survival curves that meet the rigorous standards of professional publications. It offers a wide range of customizable parameters, allowing users to fine-tune the appearance of the plots according to their research needs.
High-Performance Computation: WCESurvZ utilizes the Rcpp interface to leverage the faster underlying C language for key computations. This significantly improves the computational efficiency.
Support for Recurrent Events: WCESurvZ provides a module for handling recurrent event weighted survival analysis. It strictly follows the theoretical principles and mathematical formulas of recurrent event weighting and incorporates the weighted log-rank test with linear weight correction. This module is applicable to scenarios with and without recurrent events.

## Install WCESurvZ from GitHub
You can install the WCESurvZ package from GitHub using the following commands:
```R
devtools::install_github("your_username/WCESurvZ")
```

## Usage
Here is a basic example of how to use the WCESurvZ package:

library(WCESurvZ)
weights <- c(CVdeath = 1, MI = 0.55, Stroke = 0.455)
WCE_obj <- WCE_KMSurv(surv_exdata, weights)
plot(WCE_obj)

For more detailed documentation and examples, please refer to the package function reference.

## Developed by Xiaoqin Zhou and team.
