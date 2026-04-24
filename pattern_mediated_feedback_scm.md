# Pattern-Mediated Feedback in Simple Climate Models
## Design note for CICERO-SCM

**Date:** April 2026  
**Context:** Analysis of FaIR v2.2 posterior ensemble under SSP2-4.5.

---

## 1. Motivation

Standard simple climate models (SCMs) use a single, time-invariant climate feedback parameter $\lambda$ (W m$^{-2}$ K$^{-1}$). In a multi-box energy balance model this appears as a constant coefficient in the matrix $\mathbf{A}$ governing surface temperature evolution:

$$C_1 \frac{dT_1}{dt} = F_\text{total} - \lambda T_1 - k_2(T_1 - T_2)$$

This is a well-understood structural simplification. The problem is that GCM output, analysed via the Rugenstein & Armour (2021) differential feedback method, shows $\lambda_\text{diff}(t)$ varying by 1-2 W m$^{-2}$ K$^{-1}$ over the historical period, driven by the **pattern effect**: the spatial structure of sea surface temperature anomalies co-determines the cloud radiative feedback, so the effective feedback depends on which forcing agents are currently dominant.

Running the FaIR v2.2 841-member posterior ensemble under SSP2-4.5, $\lambda_\text{diff}(t)$ is nearly flat at $\approx -1.5$ W m$^{-2}$ K$^{-1}$ throughout 1850-2300. The reason is structural: FaIR's energy balance is linear in total forcing, and all forcing agents drive an identical SST pattern response. The model cannot represent the fact that aerosol forcing (concentrated in the Northern Hemisphere marine boundary layer) excites a different SST pattern from CO$_2$ forcing, and therefore a different cloud feedback.

Three levels of modification can address this. They are described here in increasing order of fidelity and implementation complexity.

---

## 2. Background: the forcing-composition dependence of feedback

### 2.1 Pattern effect and SST structure

Pattern effect refers to the dependence of the top-of-atmosphere (TOA) feedback on the spatial pattern of warming, not just its global mean magnitude. Clouds, particularly low marine stratocumulus, respond to local SST gradients and lower-tropospheric stability (LTS) rather than global mean temperature. This means:

- Two forcing scenarios that produce the same $\Delta T_\text{global}$ can produce different TOA imbalances if they project differently onto the SST pattern.
- Aerosol forcing (concentrated in the northern extratropics and marine boundary layer) tends to cool the tropical eastern Pacific relative to the western Pacific, opposing the El Nino-like warming pattern that suppresses low clouds. This means aerosol cooling drives SST patterns that *increase* low cloud, producing a stronger negative feedback than the same amount of CO$_2$ cooling would.
- Kawaguchi & Ceppi (2025) show that intermodel spread in the pattern effect in CMIP6 AMIP runs is dominated by spread in low cloud fraction sensitivity to LTS over the eastern subtropical oceans.

The net observable consequence is: the effective climate feedback $\lambda_\text{eff}(t)$ is less negative during periods when aerosol forcing is large (20th century industrial era), and more negative when CO$_2$ dominates (21st century onward), even if the underlying equilibrium feedback $\lambda_\text{eq}$ is constant.

### 2.2 Operational definition: differential feedback

Following Rugenstein & Armour (2021), define the differential feedback parameter at time $t$ as the slope of a regression of adjusted radiative imbalance $R = N - F$ against surface temperature $T$ in a sliding time window $[t - w, t + w]$:

$$\lambda_\text{diff}(t) = \frac{\text{Cov}(T, R)}{\text{Var}(T)}\bigg|_{[t-w,\,t+w]}$$

where $N$ is TOA radiative imbalance (positive = net energy in) and $F$ is the diagnosed ERF. In a model with a truly constant $\lambda$, this always returns $\lambda_\text{eq}$. Deviations from a constant value are a direct signal of pattern-mediated feedback. The goal of the modifications described here is for an SCM to reproduce observed/GCM-diagnosed time variation of $\lambda_\text{diff}(t)$.

---

## 3. Tier 1: Per-agent forcing efficacy (static)

### 3.1 Formulation

Replace the single total forcing $F_\text{total}$ that drives the energy balance with an efficacy-weighted sum:

$$F_\text{eff} = \sum_i E_i \, F_i$$

where $E_i$ is a static, time-invariant efficacy for forcing agent $i$. The energy balance becomes:

$$C_1 \frac{dT_1}{dt} = F_\text{eff} - \lambda T_1 - k_2(T_1 - T_2)$$

For the TOA imbalance diagnostic, continue to use the unweighted sum $F_\text{total} = \sum_i F_i$, so that the Gregory regression and $\lambda_\text{diff}$ computation reflect the pattern effect:

$$N = \lambda_\text{eq} T_1 + k_2(T_1 - T_2) + \ldots$$
$$R = N - F_\text{total} \quad \text{(observed-consistent)}$$

### 3.2 Physical interpretation

$E_i > 1$ means agent $i$ drives more warming/cooling per unit ERF than CO$_2$ does, because it excites a more sensitive SST pattern. Marvel et al. (2016) showed that aerosol forcing has $E_\text{aero} \approx 1.5$ relative to CO$_2$ in GFDL-CM3, based on single-forcing AMIP runs. Values from CMIP5/6 literature range from $\sim 1.1$ to $\sim 2.0$ depending on model and method.

The effect on $\lambda_\text{diff}(t)$: during the industrial aerosol era, $T$ is suppressed relative to $F_\text{total}$, causing the Gregory slope (regressing $R$ against $T$) to appear weaker (less negative). After aerosols peak and decline, CO$_2$ dominates and $\lambda_\text{diff}$ converges toward $\lambda_\text{eq}$. This produces exactly the time-variation seen in AMIP GCM runs.

### 3.3 Implementation requirements

For a typical multi-box SCM (upwelling-diffusion or N-box conductive):

1. **New parameters:** $E_i$ for each forcing agent (or at minimum, $E_\text{aero}$, $E_\text{CO2}$, and a common $E_\text{other}$).
2. **Forcing input:** compute $F_\text{eff} = \sum_i E_i F_i$ before passing to the energy balance step.
3. **TOA imbalance output:** continue computing $N$ from the energy balance state, not from $F_\text{eff}$, so that $N - F_\text{total}$ gives the correct Gregory intercept.
4. **No changes to the energy balance matrix:** $\mathbf{A}$ is unchanged; only the forcing input vector is modified.

This is the lowest-cost change. FaIR v2.2 already has the infrastructure for this; the calibration simply does not use it.

### 3.4 Calibration strategy

$E_i$ values can be constrained by:

- **Single-forcing AMIP pairs:** Regress $R$ against $T$ in GCM runs forced only with aerosols, and only with CO$_2$. The ratio of the two slopes (relative to a CO$_2$-only baseline) gives $E_\text{aero}$.
- **Observational constraint:** The difference between the effective climate sensitivity inferred from the historical record and the ECS inferred from abrupt-4xCO$_2$ runs is partly attributable to $E_\text{aero} \neq 1$. Lewis & Curry (2018) and Sherwood et al. (2020) discuss the implications.
- **Pattern effect diagnostics:** Andrews et al. (2022) provide CMIP6-mean pattern effect estimates that can be decomposed by forcing agent.

A minimal viable calibration target: the relationship between $\lambda_\text{diff}(t)$ in the historical period and in the 21st century, requiring $E_\text{aero}$ to reproduce the observed "step" in effective feedback across the aerosol peak.

---

## 4. Tier 2: Time-dependent per-agent efficacy

### 4.1 Motivation

In reality, efficacy is not strictly constant. It is highest when the aerosol forcing is rapidly ramping up (the SST pattern is most out of equilibrium with the forcing) and declines as the ocean adjusts. This is the same slow adjustment that $\varepsilon$ (deep-ocean efficacy) captures, but at the level of the SST spatial pattern rather than the global mean temperature departure.

### 4.2 Formulation

Replace static $E_i$ with functions of time or state:

$$F_\text{eff}(t) = \sum_i E_i(t) \, F_i(t)$$

Candidate parameterisations for $E_i(t)$:

**(a) Rate-dependent:** Efficacy scales with the rate of change of forcing:

$$E_i(t) = E_i^{\infty} + \alpha_i \left|\frac{dF_i}{dt}\right|$$

Physically: rapid forcing changes produce larger out-of-equilibrium SST patterns under thermodynamic adjustment timescales. $\alpha_i$ and $E_i^\infty$ are parameters.

**(b) State-dependent (ocean heat content):** The pattern strength decays as the deep ocean reaches quasi-equilibrium. This could be proxied by the heat content anomaly of the lower ocean layers:

$$E_i(t) = 1 + (E_i^0 - 1) \exp\!\left(-\frac{Q_\text{deep}(t)}{Q^*}\right)$$

where $Q_\text{deep}$ is the deep-ocean accumulated heat and $Q^*$ is a decay scale.

**(c) Empirical lookup:** $E_i(t)$ is prescribed as a time series derived from CMIP6 AMIP diagnostics and interpolated at runtime. This avoids the need to parameterise the functional form but limits the model to the historical period plus extrapolation.

### 4.3 Implementation requirements

1. **New parameters:** $E_i^0$, $E_i^\infty$, and $\alpha_i$ (or $Q^*$) per forcing agent.
2. **Time-loop change:** Evaluate $E_i(t)$ at each timestep before computing $F_\text{eff}$.
3. **No matrix changes:** The energy balance matrix $\mathbf{A}$ is still time-invariant.

### 4.4 Calibration strategy

This tier is significantly harder to calibrate than Tier 1 because there are more free parameters and the constraints from single-forcing AMIP runs are noisy. Bayesian calibration against the time series of $\lambda_\text{diff}(t)$ in CMIP6 models (rather than plain ECS and TCR targets) would be required.

---

## 5. Tier 3: Forcing-composition-dependent feedback

### 5.1 Formulation

This is the most physically complete formulation. Each forcing agent $i$ contributes not only to the forcing but also shifts the effective feedback parameter:

$$\lambda_\text{eff}(t) = \lambda_0 + \sum_i \Delta\lambda_i \cdot \frac{F_i(t)}{F_\text{total}(t)}$$

where:
- $\lambda_0$ is the baseline feedback (as in the current SCM)
- $\Delta\lambda_i$ is the forcing-composition feedback sensitivity for agent $i$ (units: W m$^{-2}$ K$^{-1}$)
- $F_i(t)/F_\text{total}(t)$ is the fractional contribution of agent $i$ to total forcing

When aerosols are dominant (large $|F_\text{aero}|/|F_\text{total}|$), $\lambda_\text{eff}$ is weakened by $\Delta\lambda_\text{aero}$ (a positive value, since aerosol forcing makes the feedback less negative). When CO$_2$ dominates the future, $\lambda_\text{eff} \to \lambda_0 + \Delta\lambda_\text{CO2}$, which should be close to the true equilibrium feedback.

### 5.2 Connection to Tier 1

Tiers 1 and 3 are related but not equivalent. In Tier 1, $E_i$ modifies what temperature the energy balance produces for a given forcing. In Tier 3, $\Delta\lambda_i$ modifies how efficiently that temperature radiates back to space. To first order in a single-box model:

$$F_\text{eff} - \lambda T = 0 \quad\Rightarrow\quad T = \frac{F_\text{eff}}{\lambda} = \frac{\sum_i E_i F_i}{\lambda}$$

compared to Tier 3:

$$F_\text{total} - \lambda_\text{eff} T = 0 \quad\Rightarrow\quad T = \frac{\sum_i F_i}{\lambda_0 + \sum_i \Delta\lambda_i F_i/F_\text{total}}$$

The two are equivalent when $E_i = \lambda_0/(\lambda_0 + \Delta\lambda_i)$ (approximately $1 - \Delta\lambda_i/\lambda_0$ for small $\Delta\lambda_i$). Tier 3 is more physically transparent because $\Delta\lambda_i$ directly represents a cloud feedback shift, rather than an abstract forcing amplification factor. It also makes the interaction with other pattern-effect mechanisms (e.g., deep-ocean efficacy) more interpretable.

### 5.3 The energy balance equation in matrix form

For an N-box model, the continuous-time energy balance is:

$$\mathbf{C} \frac{d\mathbf{T}}{dt} = \mathbf{A}(\lambda_\text{eff})\,\mathbf{T} + \mathbf{b}\,F_\text{total}(t)$$

where $\mathbf{C} = \text{diag}(C_1, C_2, \ldots, C_N)$ is the heat capacity matrix, $\mathbf{b} = (1, 0, \ldots, 0)^\top$ is the forcing input vector, and $\mathbf{A}$ depends on $\lambda_\text{eff}$ through its $(1,1)$ element:

$$A_{11} = -\frac{\lambda_\text{eff}(t) + k_2}{C_1}$$

In the standard (Tier 0) formulation, $\lambda_\text{eff} = \lambda_0$ is a constant pulled from calibration, and $\mathbf{A}$ is precomputed once. In Tier 3, $\lambda_\text{eff}(t)$ changes at every timestep as the forcing composition evolves, so $\mathbf{A}$ must be recomputed each step.

**Cost:** For a 3-box model, matrix construction is $O(9)$ operations and the matrix exponential (if using exact discretisation) is the main cost. Approximate schemes (forward Euler, Runge-Kutta) avoid the matrix exponential but introduce timestep sensitivity. Since SCMs already run at annual or sub-annual timesteps over centuries, even exact discretisation adds only modest overhead.

### 5.4 Implementation requirements

The following changes are needed relative to the standard SCM energy balance loop:

1. **New parameters:** $\Delta\lambda_i$ per forcing agent (at minimum: $\Delta\lambda_\text{aero}$, $\Delta\lambda_\text{CO2}$, $\Delta\lambda_\text{other}$).

2. **Per-timestep forcing decomposition:** Before the energy balance step, decompose $F_\text{total}(t)$ by agent to compute the fractional weights $w_i(t) = F_i(t)/F_\text{total}(t)$.

3. **Per-timestep $\lambda_\text{eff}$ update:**
   ```
   lambda_eff(t) = lambda_0 + sum_i( delta_lambda_i * w_i(t) )
   ```

4. **Per-timestep matrix reconstruction:** Rebuild $\mathbf{A}(\lambda_\text{eff}(t))$ using the updated $\lambda_\text{eff}$.

5. **Discretisation:** If using matrix exponential (`expm(A * dt)`), this must be recomputed each step (or cached with a small lookup if $\lambda_\text{eff}$ changes slowly). If using first-order Euler or RK4, just substitute the updated $A_{11}$ directly.

6. **TOA imbalance:** $N(t) = \lambda_\text{eff}(t)\,T_1(t) + k_2(T_1 - T_2) + \ldots$ and $R(t) = N(t) - F_\text{total}(t)$. The Gregory regression on $R$ vs $T$ will then naturally reflect the time-varying pattern effect.

### 5.5 Edge cases and numerical considerations

**Sign of $F_\text{total}$:** The fractional weights $w_i = F_i/F_\text{total}$ are ill-conditioned near $F_\text{total} = 0$ (pre-industrial baseline). A regularisation is needed:

$$w_i(t) = \frac{F_i(t)}{F_\text{total}(t) + \epsilon_F} \quad \text{or guard with } |F_\text{total}| > F_\text{min}$$

A practical choice is $F_\text{min} \approx 0.05$ W m$^{-2}$, holding $\lambda_\text{eff} = \lambda_0$ when total forcing is near zero.

**Cancellation between positive and negative agents:** Aerosol forcing is large and negative while greenhouse gas forcing is large and positive. During certain historical decades they nearly cancel, making $w_i$ large in magnitude but unreliable. A natural alternative is to weight by *magnitude* rather than signed fraction, and carry the sign separately:

$$\lambda_\text{eff}(t) = \lambda_0 + \sum_i \Delta\lambda_i \cdot \text{sgn}(F_i) \cdot \frac{|F_i(t)|}{\sum_j |F_j(t)|}$$

This avoids division by near-zero total forcing and is physically more sensible (each agent's local SST pattern scales with its own magnitude, not with the residual total).

**Multiple agents:** If the model tracks $k$ forcing agents, the implementation generalises straightforwardly. For CICERO-SCM, the minimum useful decomposition is likely:
- Aerosols (ARI + ACI combined, or separately)
- Well-mixed greenhouse gases (CO$_2$, CH$_4$, N$_2$O, halocarbons)
- Short-lived climate forcers (O$_3$, H$_2$O strat)
- Land use / albedo
- Natural (volcanic, solar)

Natural forcing (especially volcanic) warrants special treatment: individual eruptions produce very short, large-magnitude forcings that may not excite the same SST pattern as sustained anthropogenic forcing.

---

## 6. Calibration strategy for Tier 3

### 6.1 Targets

The existing SCM calibration targets (ECS, TCR, transient ocean heat uptake) are largely insensitive to $\Delta\lambda_i$ because they measure integrated responses over long periods dominated by CO$_2$. Additional targets are needed:

| Target | What it constrains |
|--------|-------------------|
| Time series of $\lambda_\text{diff}(t)$ from CMIP6 historical | $\Delta\lambda_\text{aero}$ primarily |
| Difference between historical ECS (energy budget method) and abrupt-4xCO$_2$ ECS | $E_\text{aero}$ or $\Delta\lambda_\text{aero}$ |
| Single-forcing AMIP Gregory slopes for aerosol vs CO$_2$ | Ratio $\Delta\lambda_\text{aero}/\Delta\lambda_\text{CO2}$ |
| Post-2000 warming acceleration (aerosol decline + CO$_2$ ramp) | Time variation of $\lambda_\text{eff}$ |

### 6.2 Prior construction

Literature values for the pattern effect contribution to historical ECS bias range from 0.1 to 0.8 K (Andrews et al. 2022, Sherwood et al. 2020). This implies:

$$\Delta T_\text{pattern} \approx \frac{\Delta\lambda_\text{aero} \cdot \bar{w}_\text{aero,hist}}{\lambda_0^2} \cdot F_\text{hist} \approx 0.3 \text{ K}$$

Working backwards from observed $\Delta\lambda_\text{diff} \approx 0.5$ W m$^{-2}$ K$^{-1}$ in AMIP runs with a mean aerosol fraction $\bar{w}_\text{aero} \approx -0.4$, gives $\Delta\lambda_\text{aero} \approx +1.2$ W m$^{-2}$ K$^{-1}$. This is a reasonable prior centre for a Bayesian calibration.

### 6.3 Identifiability

$\Delta\lambda_\text{aero}$ and $E_\text{aero}$ are not wholly independent: both make the historical warming appear to have a weaker feedback than the true $\lambda_0$. To decorrelate them, the calibration needs at least one target that distinguishes the *temperature response* change (Tier 1) from the *Gregory slope* change (Tier 3). The time profile of $\lambda_\text{diff}(t)$ provides this: Tier 1 shifts $\lambda_\text{diff}$ by a constant offset, while Tier 3 makes it time-varying in correlation with the aerosol forcing fraction.

---

## 7. Relationship to existing CICERO-SCM features

CICERO-SCM already includes:

- **Upwelling-diffusion thermocline:** This introduces a slow adjustment mode that partially mimics the Southern Ocean warming lag aspect of the pattern effect (analogous to deep-ocean efficacy $\varepsilon$ in FaIR). It is **not** a substitute for forcing-composition-dependent feedback, because it acts identically on all forcing agents.
- **Separate forcing time series by agent:** The model already disaggregates the forcing inputs, which means the $F_i(t)$ time series needed for $w_i(t)$ are already available. The main additions are (a) the $\Delta\lambda_i$ parameters and (b) the per-timestep recomputation of the feedback term.

The recommended implementation path is therefore:

1. Add $\Delta\lambda_\text{aero}$ as a single new calibrated parameter
2. Compute $w_\text{aero}(t) = |F_\text{aero}(t)| / \sum_j |F_j(t)|$ at each timestep
3. Apply $\lambda_\text{eff}(t) = \lambda_0 + \Delta\lambda_\text{aero} \cdot w_\text{aero}(t)$
4. Rebuild the upper-left element of the heat transfer matrix at each step

This single-parameter extension would capture the dominant source of time-varying feedback without requiring a full forcing decomposition across all agents.

---

## 8. Summary and recommended first step

| Tier | Parameters added | Code change | Energy-conserving | Reproduces $\lambda_\text{diff}(t)$ |
|------|-----------------|-------------|-------------------|--------------------------------------|
| 1 | $E_i$ per agent | Minimal (pre-multiply forcing) | Yes | Partially (constant offset) |
| 2 | $E_i(t)$ functions | Small (time-loop update) | Yes | Better |
| 3 | $\Delta\lambda_i$ per agent | Moderate (per-step matrix update) | Yes | Best |

The recommended first implementation is **Tier 1**, which can be validated by confirming that the Gregory regression of a CICERO-SCM run with $E_\text{aero} \approx 1.5$ reproduces the observed historical ECS bias relative to abrupt-4xCO$_2$ ECS. If this is successful and if the time-varying structure of $\lambda_\text{diff}(t)$ remains a target, **Tier 3** with a single $\Delta\lambda_\text{aero}$ parameter is the natural next step.

---

## References

- Andrews, T. et al. (2022). On the effect of historical SST patterns on radiative feedback. *Journal of Geophysical Research: Atmospheres*, 127.
- Kawaguchi, C. & Ceppi, P. (2025). Lower tropospheric stability dominates intermodel spread in the pattern effect. *Geophysical Research Letters*.
- Lewis, N. & Curry, J. (2018). The impact of recent forcing and ocean heat uptake data on estimates of climate sensitivity. *Journal of Climate*, 31, 6051-6071.
- Marvel, K. et al. (2016). Implications for climate sensitivity from the response to individual forcings. *Nature Climate Change*, 6, 386-389.
- Rugenstein, M. & Armour, K. (2021). Three flavors of radiative feedback and their implications for estimating equilibrium climate sensitivity. *Geophysical Research Letters*, 48.
- Sherwood, S. C. et al. (2020). An assessment of Earth's climate sensitivity using multiple lines of evidence. *Reviews of Geophysics*, 58.
