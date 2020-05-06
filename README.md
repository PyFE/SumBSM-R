# SumBSM-R

This the set of R codes used for the numerical examples of the "Sum of all Black-Scholes-Merton models" paper ([Choi 2018])

## Abstract
Contrary to the common view that exact pricing is prohibitive owing to the curse of dimensionality, this study proposes an efficient and unified method for pricing options under multivariate Black-Scholes-Merton (BSM) models, such as the basket, spread, and Asian options. The option price is expressed as a quadrature integration of analytic multi-asset BSM prices under a single Brownian motion. Then the state space is rotated in such a way that the quadrature requires much coarser nodes than it would otherwise or low varying dimensions are reduced. The accuracy and efficiency of the method is illustrated through various numerical experiments.

## Files
* [SumBSM/blksmd.R](SumBSM/blksmd.R): Functions for pricing. Two main interface functions are 
  * `blksmd_basket`: basket or spread options
  * `blksmd_asian`: Asian options 
* [SumBSM/Table_04.R](SumBSM/Table_04.R): Table 4. Parmeter set __S1__ (spread option)
* [SumBSM/Table_05.R](SumBSM/Table_05.R): Table 5. Parmeter set __S2__ (spread option)
* [SumBSM/Table_10.R](SumBSM/Table_10.R): Table 10. Parmeter set __A1__ (discretely monitored Asian option)
* [SumBSM/Table_11.R](SumBSM/Table_11.R): Table 11. Parmeter set __A2__ (discretely monitored Asian option)
* [SumBSM/Table_12.R](SumBSM/Table_12.R): Table 12. Parmeter set __A2__ (continuously monitored Asian option)

## Reference
* __Choi, J.__ (2018). Sum of all Black-Scholes-Merton models: An efficient pricing method for spread, basket, and Asian options. *Journal of Futures Markets*, 38(6), 627–644. [ [DOI](https://doi.org/10.1002/fut.21909) | [SSRN](http://ssrn.com/abstract=2913048) | [arXiv](http://arxiv.org/abs/1805.03172) ]

[Choi 2018]: https://doi.org/10.1002/fut.21909
