# Bayesian Network Analysis for Modeling Cervical Cancer Risk Factors

Cervical cancer is a disease that affects the cervix, the lower part of the uterus in the human female reproductive system.
Despite being easy to diagnose, it still threatens thousands of women each year, especially in areas with high poverty levels and low screening rates.
The aim of the analysis is to observe how the already known risk factors for cervical cancer affect the probability of developing it and, possibly, discover new sources of exposure.

## Dataset

The dataset took into consideration is the [Cervical Cancer Risk Factors Dataset](https://archive.ics.uci.edu/ml/datasets/Cervical%2Bcancer%2B%2528Risk%2BFactors%2529),
gathered in Caracas, Venezuela, at the *Hospital Universitario de Caracas*; it collects information about 858 female patients of the hospital, in particular their demographic data, habits and medical history.

## Technologies

* `bnlearn`: learning the graphical structure of Bayesian networks, estimate their parameters and perform some useful inference
* `gRain`: propagation in graphical independence networks
* `gRbase`: provides graphical modelling features
* `Rgraphviz`: provides plotting capabilities for R graph objects
* `dplyr`: fast and consistent tool for working with dataframes

## References

* Fernandes, K., Cardoso, J. S., & Fernandes, J. (2017). Transfer Learning with Partial Observability Applied to Cervical Cancer Screening. In Pattern Recognition and Image Analysis (pp. 243–250). Springer International Publishing. https://doi.org/10.1007/978-3-319-58838-4_27
* Heckerman, D. (2020). A Tutorial on Learning With Bayesian Networks. ArXiv. https://doi.org/10.48550/ARXIV.2002.00269
