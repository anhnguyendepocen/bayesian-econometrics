# Spatial Error Model (SEM)

The SEM model is specified as,
$$
y = X\beta + u
$$
$$
u = \lambda Wu + v
$$

where $W$ is a spatial weight matrix, $u$ is a spatial error vector, and the error term $v \sim \mathcal{N}(0,\sigma^2)$.

The reduced form of the model is,
$$
y = X\beta + (I_N - \lambda W)^{-1} v
$$

So, the response $y$ is modeled as,
$$
y \sim \mathcal{N}(X\beta,\ (I_N - \lambda W)^{-1} \cdot (I_N - \lambda W^{\top})^{-1} \cdot \sigma^2)
$$

which is how I'm modeling it in Stan, but I'm taking the square root of the diagonal of the scale parameter stuff first since we do things in terms of standard deviation. 
