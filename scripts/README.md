Scripts for creating figures / performing analyses

## Linear regression of mutations against time

The model for the linear regression of number of APOBEC3-type mutations against time was as follows:
>   m<sub>i</sub> ~ Normal(μ, σ)
>   μ = α + β  (t<sub>i</sub> - tmean) 
>   α ~ Normal(11, 100)
>   β ~ Lognormal(0, 1) and
>   σ ~ Normal(0, 50),

where m<sub>i</sub> is the number of mutations for genome i,
α is the y-intercept and has a prior centred on the minimum number of mutations observed over all genomes, β is a strictly-positive evolutionary rate per year, t<sub>i</sub> is the time of collection of the sample and σ is the model error standard deviation.

### posterior estimates of parameters

Lineage A

|    |    mean   |   sd   |     2.5%    |  97.5%
|--- | --- | --- | --- | ---          
|\alpha|    25.236597| 0.7343179| 23.797360| 26.675833
|\beta|      6.176307| 0.5009553|  5.194453|  7.158162
|\sigma|     4.759052| 0.5197808|  3.740300|  5.777804
|time of origin| 2015.462 | | 2014.757 | 2016.163

Lineage B.1

|      |    mean   |   sd   |     2.5%    |  97.5%
|--- | --- | --- | --- | ---          
|\alpha |57.669440|0.09242696| 57.488287 |57.850594
|\beta  | 5.927578|1.52043863|  2.947573 | 8.907583
|\sigma | 1.768233|0.06544786|  1.639957 | 1.896508
