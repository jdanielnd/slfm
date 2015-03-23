slfm
===
[![Build Status](https://travis-ci.org/jdanielnd/slfm.png)](https://travis-ci.org/jdanielnd/slfm)

Sparse Latent Factor Model

slfm is a set of tools to find coherent patterns in microarray data using a Bayesian sparse latent factor model. Considerable effort has been put into making slfm fast and memory efficient, turning it an interesting alternative to simpler methods in terms of execution time. It implements versions of the SLFM using both type of mixtures: using a degenerate distribution or a very concentrated normal distribution for the spike part of the mixture. It also implements additional functions to help pre-process the data and fit the model for a large number of arrays.