# Parallel Tempering MCMC Examples

This repository demonstrates **Parallel Tempering (Metropolis-Coupled) Markov Chain Monte Carlo (MCMC)** methods and illustrates why they are superior to standard Metropolis–Hastings (MH) sampling in the presence of **multi-modal posterior distributions**.

The examples are designed to be simple and reproducible, focusing on how coupling multiple chains at different temperatures improves **mixing**, **exploration**, and **convergence**.

---

## Motivation

Standard MCMC algorithms often fail when the target distribution is **multi-modal**. In such cases, a single chain may become trapped in one mode, resulting in poor mixing and biased inference.

**Parallel Tempering (PT-MCMC)** mitigates this issue by:
- Running multiple chains at different temperatures
- Allowing probabilistic state swaps between chains
- Enabling high-temperature chains to explore freely while the low-temperature (“cold”) chain preserves the correct target distribution

---

## General Experimental Setup

Across all examples:

- Symmetric random-walk proposal distributions are used
- Priors are weakly informative or uniform
- Each experiment runs **100,000 iterations**
- The **first 80%** of samples are discarded as burn-in
- Statistical inference is based **only on the cold chain**
---

---
## Codes

- Example1.R corresponds to the MCMC run in example 1 
- Example2.R corresponds to the MCMC run in example 2
- Example3.R corresponds to the MCMC run in example 3

Details of the examples can be found in Parallel_Tempering_MCMC_Ex.pdf

---
## Example 0: Double Well Potential

### Model

A classical double potential well distribution:

\[
U(x) = \gamma (x^2 - 1)^2, \quad \pi(x) = \exp(-U(x))
\]

The objective is to explore the distribution and estimate the parameter \( \gamma \).

### Challenge

- Increasing \( \gamma \) leads to sharper bi-modality
- Standard MH chains:
  - Become trapped in a single mode
  - Mix poorly
  - Fail to explore the full state space

### Parallel Tempering Solution

By introducing multiple chains with different temperatures and allowing **state swaps**, the sampler successfully transitions between modes, resulting in improved mixing and exploration.

---

## Example 1: Gaussian Mixture Model

### Model

A two-component Gaussian mixture:

- True parameters: \( \theta_1 = 6 \), \( \theta_2 = -6 \)
- Known variance: \( \sigma^2 = 1 \)
- Priors:
  - \( \theta_1, \theta_2 \sim \mathcal{N}(0, 100) \)

### Challenge

A single MH chain:
- Becomes stuck near one component
- Fails to recover both mixture means
- Shows slow or no convergence

### Parallel Tempering Setup

- 5 chains with temperature parameters:
  (0.1, 0.3, 0.5, 0.7, 1)

- Only **state swaps** are performed (no heat transfer)
- Inference uses the cold chain only

### Outcome

Parallel tempering achieves:
- Better mixing
- Exploration across both modes
- Accurate parameter recovery

---

## Example 2: Non-Identifiable Likelihood

### Model

Unidentifiable Gaussian with \alpha^2

### Challenge

Standard MH sampling:
- Gets trapped in one symmetric mode
- Does not converge to true parameter values
- Exhibits poor mixing

### Parallel Tempering Setup

- 5 chains with temperature parameters:
  \[
  (0.1, 0.3, 0.5, 0.7, 1)
  \]

### Outcome

Parallel tempering:
- Explores symmetric modes effectively
- Improves convergence
- Produces reliable posterior samples from the cold chain

---

## Key Takeaways

- Parallel tempering is highly effective for **multi-modal** and **non-identifiable** models
- Even low-dimensional problems can defeat standard MH algorithms
- State swapping alone can dramatically improve sampler performance
- Inference should always be based on the cold chain

---

## References

- Wilkinson, D. (2013). *Parallel tempering and Metropolis-coupled MCMC*.  
  https://darrenjw.wordpress.com/2013/09/29/parallel-tempering-and-metropolis-coupled-mcmc/

- Verity, B. (2024). *Parallel tempering*.  
  https://mrc-ide.github.io/drjacoby/articles/metropolis_coupling.html
