

Several of the "traits" were significantly correlated with population structure (in the form of PC eigenvalues). This suggests that we are better-off using the G model instead of the QG model (as the eigenvalues will be confounded with the trait values).


### GWAS

The baseline model for GWAS is as follows:

$$
\mathbf{y = 1\beta + Wm + Zu + e}
$$

where $\mathbf{y}$ is the response vector, $\mathbf{1}$ is a vector of 1s for the grand mean, $\mathbf{\beta}$ is the grand mean, $\mathbf{Z}$ is the incidence matrix modeling the random genetic background effect, $\mathbf{g}$ is the vector of random genetic background effects, $\mathbf{W}$ is the SNP marker matrix, and $\mathbf{m}$ is the vector of fixed additive SNP effects. $\mathbf{g}$ is assumed to be distributed such that $\mathbf{g} \sim N(0, \mathbf{K} \sigma^2_g )$, where $\mathbf{K}$ models the population structure among genotypes.

I use population parameters previously determined [P3D; @Kang2010] to speed up the models without reducing detection power significantly. I also the `gwas()` function in my `pbr` package. Models included K and G [@Bernardo2013b; @Chen2016].

The heritability will be calcuted for each trait across all environments by fitting the following model:

$$
y_{ij} = \mu + g_i + t_j + (gt)_{ij} + \epsilon_{ij}
$$

where $y_{ij}$ is the adjusted mean of the *i*th genotype in the *j*th environment, $\mu$ is the grand mean, $g_i$ is the random effect of the *i*th genotype, $t_j$ is the fixed effect of the $j$th environment, $(gt)_{ij}$ is the random effect of the interaction between the $i$th genotype and the $j$th environment, and $\epsilon_{ij}$ is the error associated with the observation.

We will assume the random effects are distributed such that $g \sim N(0, \sigma^2_g)$, $(gt) \sim N(0, \sigma^2_{(gt)})$, and $\epsilon \sim N(0, \mathbf{R}\sigma^2_\epsilon)$, where $\mathbf{R}$ is a diagonal matrix with elements equal to the inverse of the variances of the adjusted genotype means.



$$
H = \frac{ \sigma^2_G }{ \sigma^2_G + \frac{\sigma^2_{GE}}{r} + \frac{\sigma^2_R}{er} }
$$

for overall heritability calculations, and 

$$
H = \frac{ \sigma^2_G }{ \sigma^2_G  + \frac{\sigma^2_R}{r} }
$$


where $e$ is the number of environments and $r$ is the number of replicates. In the case of an augmented design, where the number of replicates is uneven, the following formula is used to calculate the harmonic mean of the number of replications:

$$
p_h = \frac{n}{\sum^n_{i=1} \frac{1}{p_i}}
$$

where $p_h$ is the harmonic mean, $n$ is the number of genotypes, and $p_i$ is the number of replicates of the *i*th genotype.



Use the environment means to calculate the FW stability coefficient. The model for the regression coefficients looks like:

$$
y_{ij} = \mu + G_i + (1 + b_i)E_j + \epsilon_{ij},
$$

where $y_{ij}$ is the BLUE of the *i*th genotype in the *j*th environment, $\mu$ is the overall mean, $g_i$ is the effect of the *i*th genotype, $E_j$ is the random effect of the *j*th environment, and $b_i$ is the regression coefficient of the response on the mean of the *j*th environment, andx $\epsilon_{ij}$ is the random error.

To fit FW, we first fit a regular quantitative genetic model:

$$
y_{ij} = \mu + G_i + E_j + \epsilon_{ij}
$$

to obtain the environmental effects $E_j$. We then fit the following model for every genotype

$$
y_j = \mu + (1 + b)E_j + \epsilon_{ij}
$$





