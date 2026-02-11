---
marp: true
theme: gaia
paginate: true
---


# Test Driven Development with pFUnit
## Using pFUnit to apply TDD to numerical Fortran development

---

<style scoped>
section { font-size: 18px; }
</style>

## KSG Mutual Information Estimator
### Algorithm

**Input:** Paired samples $(x_1,y_1)\dots(x_N,y_N)$, neighbor count $k$

Precompute $\psi(n)$ for $n = 1 \cdots N+1$

For each point $i = 1 \dots N$:

1. Compute joint distances  
   $d(j) = \max(\lvert x_i - x_j \rvert , \lvert y_i - y_j \rvert)$

2. Find the $k$ nearest neighbors (exclude itself)

3. From those neighbors define  
   $\epsilon_x = \max \lvert x_i - x_j \rvert$  
   $\epsilon_y = \max \lvert y_i - y_j \rvert$

4. Count  
   $n_x(i)$: $\lvert x_i - x_j \rvert \le \epsilon_x$ (exclude $i$)  
   $n_y(i)$: $\lvert y_i - y_j \rvert \le \epsilon_y$ (exclude $i$)

Compute the average:
$\text{avg} = \text{mean}\,[\psi(n_x(i)+1) + \psi(n_y(i)+1)]$

$MI = \psi(N+1) + \psi(k+1) - \text{avg} - \frac{1}{k}$

---

## KSG Mutual Information Estimator
### Properties

- Mutual information $\approx 0$ for independent variables  
  when $k = \frac{N}{2}$ and $N \sim 1000$

- For correlated Gaussians:
  
  $I(X;Y) = -\frac{1}{2}\log(1-\rho^2)$

  - When $\rho=0.9$, $I(X;Y) \approx 0.830366$ 
  - Estimate is most accurate when $k = (0.04 \cdot N) + 0.5$  
  - Estimate decreases as $k$ increases
- For the generalized correlation coefficient:
  
  $R(X;Y) = \sqrt{1 - e^{-2I(X;Y)}}$
- Correlation equals 1 when the correlated relationship is perfectly nonlinear (eg. $x=sin$ and $y=cos$)
- Generalized correlation equals absolute value of pearson correlation for linear correlations with negative covariance.
- Generalized correlation is zero for independent variables.


