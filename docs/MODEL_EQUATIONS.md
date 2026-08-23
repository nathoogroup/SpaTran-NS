# Full Model Equations: Stationary vs Non-stationary SPDE

This document spells out the complete stationary (M₀) and non-stationary (M₁) models with all equations substituted in. Reference: Lindgren et al. (2011), Ingebrigtsen et al. (2014).

---

## The Nugget

The **nugget** \( \epsilon(s) \) is the **measurement error** at each spot: \( \epsilon(s) \sim N(0, \sigma_\epsilon^2) \), iid across locations. It captures technical noise (e.g. sequencing, capture efficiency) and any microscale variation that cannot be resolved by the spatial model. In geostatistics it is called the "nugget effect" because it contributes variance at zero spatial lag—even at the same nominal location, repeated measurements differ. The total variance at a spot is \( \sigma_b^2 + \sigma_\epsilon^2 \) (spatial + nugget).

---

## One Big Equation

The full model in one expression (non-stationary case M₁; for M₀ set \( \theta_3 = \theta_4 = 0 \)):

\[
\boxed{
Y(s) = \mu + u(s) + \epsilon(s), \quad
\bigl(\kappa^2(s) - \Delta\bigr)\bigl(\tau(s)\,u(s)\bigr) = W(s), \quad
\tau(s) = e^{\theta_1}, \quad
\kappa(s) = \exp\bigl(\theta_2 + \theta_3\tilde x(s) + \theta_4\tilde y(s)\bigr), \quad
\epsilon(s) \sim N(0, \sigma_\epsilon^2)
}
\]

with \( \tilde x(s) = (x - \bar{x})/\text{sd}(x) \), \( \tilde y(s) = (y - \bar{y})/\text{sd}(y) \), \( W(s) \) Gaussian white noise, and \( \alpha = 2 \). The marginal variance and practical range of \( u \) are \( \sigma_b^2 = 1/(4\pi\kappa^2\tau^2) \), \( \rho = \sqrt{8}/\kappa \).

---

## 1. Observation Model

For each gene, at spatial location \( s \in \mathbb{R}^2 \):

\[
Y(s) = \mu + u(s) + \epsilon(s)
\]

where:
- \( \mu \) = intercept (mean expression)
- \( u(s) \) = Gaussian spatial random effect (Matérn GMRF)
- \( \epsilon(s) \sim N(0, \sigma_\epsilon^2) \) iid nugget (measurement noise)

---

## 2. SPDE (Stochastic Partial Differential Equation)

The spatial field \( u(s) \) is defined via the SPDE:

\[
\bigl(\kappa^2(s) - \Delta\bigr)^{\alpha/2} \bigl(\tau(s)\, u(s)\bigr) = W(s)
\]

with:
- \( \alpha = 2 \) (so \( \nu = 1 \) Matérn smoothness in 2D)
- \( \Delta \) = Laplacian
- \( W(s) \) = Gaussian white noise
- \( \tau(s), \kappa(s) \) = spatially-varying SPDE parameters

**Marginal variance** and **practical range** of \( u \) at any location where \( \tau, \kappa \) are defined:

\[
\sigma_b^2 = \frac{1}{4\pi\,\kappa^2\tau^2}, \qquad
\rho = \frac{\sqrt{8}}{\kappa}
\]

---

## 3. Stationary Model (M₀)

In M₀, \( \tau \) and \( \kappa \) are constant over space:

\[
\log\tau(s) = \theta_1 \quad \Rightarrow \quad \tau(s) = \tau = e^{\theta_1}
\]

\[
\log\kappa(s) = \theta_2 \quad \Rightarrow \quad \kappa(s) = \kappa = e^{\theta_2}
\]

**SPDE with parameters substituted:**

\[
\bigl(e^{2\theta_2} - \Delta\bigr) \bigl(e^{\theta_1}\, u(s)\bigr) = W(s)
\]

**Variance and range (constant):**

\[
\sigma_b^2 = \frac{1}{4\pi\,e^{2\theta_2}\,e^{2\theta_1}} = \frac{1}{4\pi\,e^{2(\theta_1+\theta_2)}}
\]

\[
\rho = \frac{\sqrt{8}}{e^{\theta_2}} = \sqrt{8}\,e^{-\theta_2}
\]

**INLA parameterisation (B matrices):**
- \( \mathbf{B}_\tau = [0,\; 1,\; 0] \) → \( \log\tau = 0\cdot W + 1\cdot\theta_1 + 0\cdot\theta_2 = \theta_1 \)
- \( \mathbf{B}_\kappa = [0,\; 0,\; 1] \) → \( \log\kappa = 0\cdot W + 0\cdot\theta_1 + 1\cdot\theta_2 = \theta_2 \)

---

## 4. Non-stationary Model (M₁)

In M₁, \( \tau \) is constant but \( \kappa \) varies linearly in space (Ingebrigtsen et al. 2014):

\[
\log\tau(s) = \theta_1 \quad \Rightarrow \quad \tau(s) = \tau = e^{\theta_1}
\]

\[
\log\kappa(s) = \theta_2 + \theta_3\,\tilde x(s) + \theta_4\,\tilde y(s)
\]

where \( \tilde x(s), \tilde y(s) \) are **standardised** spatial coordinates:

\[
\tilde x(s) = \frac{x(s) - \bar{x}}{\text{sd}(x)}, \qquad
\tilde y(s) = \frac{y(s) - \bar{y}}{\text{sd}(y)}
\]

So:

\[
\kappa(s) = \exp\bigl(\theta_2 + \theta_3\,\tilde x(s) + \theta_4\,\tilde y(s)\bigr)
\]

**SPDE with parameters substituted:**

\[
\bigl(\kappa^2(s) - \Delta\bigr) \bigl(e^{\theta_1}\, u(s)\bigr) = W(s)
\]

with

\[
\kappa(s) = \exp\bigl(\theta_2 + \theta_3\,\tilde x(s) + \theta_4\,\tilde y(s)\bigr)
\]

**Local variance and range** (at mesh node \( i \) with scaled coords \( \tilde x_i, \tilde y_i \)):

\[
\sigma_b^2(s_i) = \frac{1}{4\pi\,\kappa^2(s_i)\,\tau^2} = \frac{1}{4\pi\,\exp\bigl(2\theta_2 + 2\theta_3\tilde x_i + 2\theta_4\tilde y_i\bigr)\,e^{2\theta_1}}
\]

\[
\rho(s_i) = \frac{\sqrt{8}}{\kappa(s_i)} = \sqrt{8}\,\exp\bigl(-\theta_2 - \theta_3\tilde x_i - \theta_4\tilde y_i\bigr)
\]

**INLA parameterisation (B matrices, per mesh node \( i \)):**
- \( \mathbf{B}_\tau[i,] = [0,\; 1,\; 0,\; 0,\; 0] \) → \( \log\tau = \theta_1 \) (constant)
- \( \mathbf{B}_\kappa[i,] = [0,\; 0,\; 1,\; \tilde x_i,\; \tilde y_i] \) → \( \log\kappa(s_i) = \theta_2 + \theta_3\tilde x_i + \theta_4\tilde y_i \)

---

## 5. Summary Table

| Quantity | Stationary (M₀) | Non-stationary (M₁) |
|----------|------------------|---------------------|
| \( \tau(s) \) | \( e^{\theta_1} \) | \( e^{\theta_1} \) |
| \( \kappa(s) \) | \( e^{\theta_2} \) | \( \exp(\theta_2 + \theta_3\tilde x + \theta_4\tilde y) \) |
| Hyperparameters | \( \theta_1, \theta_2 \) | \( \theta_1, \theta_2, \theta_3, \theta_4 \) |
| \( \sigma_b^2 \) | \( 1/(4\pi\,e^{2\theta_1+2\theta_2}) \) | \( 1/(4\pi\,\kappa^2(s)\,\tau^2) \) (varies in space) |
| \( \rho \) | \( \sqrt{8}\,e^{-\theta_2} \) | \( \sqrt{8}/\kappa(s) \) (varies in space) |

---

## 6. Inverse Relations (for Simulation Setup)

Given target marginal variance \( \sigma_b^2 \) and practical range \( \rho \):

\[
\kappa = \frac{\sqrt{8}}{\rho}, \qquad
\tau = \frac{1}{\sqrt{4\pi\,\kappa^2\,\sigma_b^2}}
\]

So \( \theta_1 = \log\tau \), \( \theta_2 = \log\kappa \).

---

## 7. Bayes Factor

\[
\log BF_{10} = \log p(Y \mid M_1) - \log p(Y \mid M_0)
\]

INLA provides \( \log p(Y \mid M) \) via the integrated nested Laplace approximation. Both models share the same prior family on \( \theta \), so the marginal log-likelihoods are directly comparable.
