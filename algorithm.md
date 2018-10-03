**Non-Linear Transformation**

##Inputs:
1. A random vector $(\x,\theta)\in \R^5\times\mathbb{S}$ with GVM distribution with given parameters.
2. A diffeomorphism

##Output:
A GVM distribution that approximates the true propagated distribution.

****
#### Step 1: How to generate quadrature points with weights for a canonical GVM distribution.

1. First we fix $\kappa$. $\color{red}{\text{Figure out how this is done}}$
2. $\xi=\sqrt{3}$ and $\eta = \text{arccos}\Big(\frac{B_2(\kappa)}{2B_1(\kappa)}-1\Big)$
3. The number of $\sigma$ points needed : 2*5+3 = 13.
4. Sigma point locations for the canonical GVM: ($\boldsymbol{\z}_i,\phi_1$)
    1. $\mathcal{N}^{00} = [\text{zeros}(5,1);0]$ with $w_{00}=1-2w_{\eta0}-10*w_{\xi0}$
    2. $\mathcal{N}^{\eta0}_1 = [\text{zeros}(5,1);\eta]$ with $w_{\eta0} = \frac{B_1(\kappa)^2}{4B_1(\kappa)-B_2(\kappa)}$
    3. $\mathcal{N}^{\eta0}_2 = [\text{zeros}(5,1);-\eta]$ with $w_{\eta0} = \frac{B_1(\kappa)^2}{4B_1(\kappa)-B_2(\kappa)}$
    4. $\mathcal{N}^{\xi0}_{i1} = [[0,0,...,\xi,...,0]^T;0]$ for  $i=1(1)5$ with $w_{\xi0}=\frac{1}{6}$
    5. $\mathcal{N}^{\xi0}_{i2} = [[0,0,...,-\xi,...,0]^T;0]$ for  $i=1(1)5$ with $w_{\xi0}=\frac{1}{6}$

5. $\xi$ occurs in the $i^{th}$ location.
6. $B_p(\kappa) :=  1-\frac{I_p(\kappa)}{I_0(\kappa)}$.

#### Step 2: How to generate quadrature points with weights for a input GVM distribution.
1. Input GVM parameters:
    1. $\boldsymbol{\mu},\P,\alpha,\boldsymbol{\beta},\boldsymbol{\Gamma},\kappa$
    2. Take lower triangular Cholesky factor of $\P$ as $\A$, i.e. $\P=\A\A^T$.
    3. $\x_{\sigma_i} = \boldsymbol{\mu}+\A\boldsymbol{z}_{\sigma_i}$
    4. $\theta_{\sigma_i} = \phi_{\sigma_i}+\alpha+\boldsymbol{\beta}^T\boldsymbol{z}_{\sigma_i}+\frac{1}{2}\boldsymbol{z}_{\sigma_i}^T\boldsymbol{\Gamma}\boldsymbol{z}_{\sigma_i}$

#### Step 3: How to propagate the prior $\sigma$ points
Just use an ODE solver to propagate each sigma points to get $\{(\tilde{\x}_{\sigma_i},\tilde{\theta}_{\sigma_i})\}^{N}_{i=1}$.

#### Step 3: How to generate a GVM approximation of the propagated points
1. First calculate the weighted sample mean and covariance matrix using $\tilde{\x}_{\sigma_i}$.
2. Calculate the lower-triangular Cholesky factor of $\tilde{\P}=\tilde{\A}\tilde{\A}^T$.
3. So now we have $\tilde{\boldsymbol{\mu}}$ and $\tilde{\P}$.
4. We need to calculate $\tilde{\alpha},\tilde{\boldsymbol{\beta}},\tilde{\boldsymbol{\Gamma}}$.
5. $\color{red}{\text{How do you calculate}\ \tilde{\alpha}}$
6. Calculate $\tilde{\boldsymbol{\beta}}$ as:
$$\tilde{\A}^{-1}\A[\boldsymbol{\beta}+\A^T[-\frac{3}{2}n_0/a_0\Delta t,0,0,0,0]^T]$$
7. Calculate $\tilde{\boldsymbol{\Gamma}}$ as:
$$\tilde{\A}^{-1}\A[\boldsymbol{\Gamma}+\A^T\begin{bmatrix}\frac{15n_0}{4a_0^2\Delta t} &0 &...&0 \\ 0 & 0 & ... & 0 \\ . & . &. & .\end{bmatrix}\A]\A^T\tilde{\A}^{-T}$$
8. $a_0=\mu_1$; $n_0 = \sqrt{\frac{\mu_\oplus}{a_0^3}}$; $\Delta t = t_f-t_0$
****
#### Step 3: How to generate a OPTIMAL GVM approximation of the propagated points
1. Solve a non-linear optimization problem. The problem is:
$$(\tilde{\alpha},\tilde{\boldsymbol{\beta}},\tilde{\boldsymbol{\Gamma}}) = \arg \min_{\hat{\alpha},\hat{\boldsymbol{\beta}},\hat{\boldsymbol{\Gamma}}}\sum_{i=1}^{N}[r_i(\hat{\alpha},\hat{\boldsymbol{\beta}},\hat{\boldsymbol{\Gamma}})]^2$$

****
###Implementation using MATLAB

#### Input GVM parameters:
$$\begin{align}\boldsymbol{\mu} = ,\P,\alpha,\boldsymbol{\beta},\boldsymbol{\Gamma},\kappa\end{align}$$


****
## Appendix
Equinoctial Orbital Elements:
$$\begin{align}a &= a \\ h &= e\text{sin}(\omega+\Omega)\\ k &= e\text{cos}(\omega+\Omega)\\ p &= \text{tan}\frac{i}{2}\text{sin}(\Omega)\\q &= \text{tan}\frac{i}{2}\text{cos}(\Omega)\\l&=M+\omega+\Omega\end{align}$$

Modified Equinoctial Orbital Elements:
$$\begin{align}p &= a(1-e^2) \\ f &= e\text{cos}(\omega+\Omega)\\ g &= e\text{sin}(\omega+\Omega)\\ h &= \text{tan}\frac{i}{2}\text{cos}(\Omega)\\k &= \text{tan}\frac{i}{2}\text{sin}(\Omega)\\l&=\nu+\omega+\Omega\end{align}$$
