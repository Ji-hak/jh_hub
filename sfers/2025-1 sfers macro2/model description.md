# Modified BGG (1999) DSGE Model — Full Walkthrough (ver1)

**Author of this file:** macro2  
**Short description:** A Bernanke–Gertler–Gilchrist (1999)-style New Keynesian model with a financial accelerator. The monetary policy rule _does not_ include an explicit financial-stability term in this version.  
**Audience:** Senior undergraduate economics students; assumes core macro and basic dynamic optimization. We explain log-linearization and timing operators carefully.

---

## 0) What this model is (and is not)

This is a medium-scale New Keynesian (NK) model enriched with a **financial accelerator** à la **BGG (1999)**. The model combines sticky prices (Phillips curve), intertemporal consumption smoothing (Euler), investment with adjustment costs, and an **entrepreneurial balance-sheet channel** that makes the external finance premium (the **risk spread**) move with **leverage**/**net worth**. Monetary policy sets the **nominal policy rate** via a **Taylor-type rule** with smoothing; the **real rate** follows from the **Fisher relation**. In this _ver1_, the rule **does not** include a financial-stability term (e.g., leverage, credit growth).

Compared to the textbook BGG setup, your file keeps the essence but uses a compact, log-linear representation with explicit **measurement equations** for annualized observables, and a simple reduced-form mapping for the **spread** and **net-worth dynamics** that matches the original economic logic.

---

## 1) Variables & Notation

### 1.1 Time indexing & operators

- Curly braces denote time: `x{t-1}`, `x{t}`, `x{t+1}`.
- All model equations are **linearized around the steady state**. Unless stated otherwise, a variable is a **(log) deviation from steady state** (for real quantities and prices) or a **deviation in percent/pp** for rates.
- We use \( \text{pie}\_t \) for inflation (your mnemonic `pie`) and \( rn_t \) for the **ex-ante nominal policy rate**. The **ex-ante real rate** is \( r_t \).
- **Annualization** (measurement block) multiplies quarterly rates by 4.

### 1.2 Endogenous state/control variables

| Mnemonic | Meaning (linearized)                        | Units / interpretation |
| -------- | ------------------------------------------- | ---------------------- |
| `y`      | Output gap                                  | log-deviation          |
| `pie`    | Inflation \( \pi_t \) (quarterly)           | percent / pp           |
| `r`      | Ex-ante **real** interest rate              | percent / pp           |
| `rn`     | **Nominal** policy rate (quarterly)         | percent / pp           |
| `n`      | Entrepreneurial **net worth**               | log-deviation          |
| `a`      | Productivity                                | log-deviation          |
| `c`      | Consumption                                 | log-deviation          |
| `ce`     | Entrepreneur consumption                    | log-deviation          |
| `rk`     | **Return to capital premium** / risk spread | percent / pp           |
| `q`      | Price of capital (Tobin’s Q)                | log-deviation          |
| `k`      | Capital                                     | log-deviation          |
| `x`      | Price markup / real marginal cost gap       | log-deviation          |
| `h`      | Hours / labor                               | log-deviation          |
| `lev`    | Leverage \( \approx q + k - n \)            | log-deviation          |
| `g`      | Government spending                         | log-deviation          |
| `i`      | Investment                                  | log-deviation          |
| `nr`     | Natural real rate                           | percent / pp           |
| `z`      | Preference shock process                    | log-deviation          |

### 1.3 Observables (measurement variables)

| Mnemonic | Definition                                                           |
| -------- | -------------------------------------------------------------------- |
| `rAnn`   | Annualized nominal policy rate: \( rAnn_t = rAnnSS + 4 \cdot rn_t \) |
| `piAnn`  | Annualized inflation: \( piAnn_t = 4 \cdot pie_t \)                  |
| `rnAnn`  | Annualized natural rate: \( rnAnn_t = rAnnSS + 4 \cdot nr_t \)       |

`rAnnSS = 100 \cdot ((1/\beta)^4 - 1)` is the steady-state **annual** nominal rate implied by discount factor \( \beta \).

### 1.4 Exogenous shocks

| Shock           | Mnemonic | Law of motion                                        |
| --------------- | -------- | ---------------------------------------------------- |
| Productivity    | `ea`     | \( a*t = \rho_a a*{t-1} + s_a\, ea_t \)              |
| Policy shock    | `etamp`  | Enters Taylor rule for \( rn_t \)                    |
| Net-worth shock | `en`     | Enters net-worth law of motion                       |
| Govt spending   | `eg`     | \( g*t = \rho_g g*{t-1} + s_g\, eg_t \)              |
| Natural rate    | `etarn`  | \( nr_t = (1-\rho_a)a_t + (1-\rho_z)z_t + etarn_t \) |

Preference shock process: \( z*t = \rho_z z*{t-1} \).

---

## 2) Model blocks (big picture)

We organize the model into **four blocks**. Each block is internally coherent and the blocks interact to determine equilibrium dynamics.

### (A) Household & Demand block

- **Euler**: \( c*t = -\big(rn_t - \text{pie}\_t - nr_t\big) + c*{t+1} \).  
  _Intuition_: Higher ex-ante real rate \( rn_t - \text{pie}\_t \) (relative to the natural rate \( nr_t \)) depresses current consumption; expected future consumption raises current consumption.
- **Resource constraint**: \( y_t = c_y c_t + i_y i_t + g_y g_t + ce_y ce_t \). Output is the share-weighted sum of components.
- **Government spending** follows AR(1).

### (B) Supply & Price-setting block

- **Production**: \( y*t = a_t + \alpha k*{t-1} + (1-\alpha)\Omega h*t \).  
  *Intuition*: Cobb–Douglas in (effective) labor \( \Omega h_t \) and capital \( k*{t-1} \).
- **Labor FOC**: \( y_t = x_t + c_t + (1+\eta)h_t \). Links hours, consumption (wealth effect), and the markup \( x_t \).
- **Phillips curve (NKPC)**: \( \text{pie}_t = -\kappa x_t + \beta\, \text{pie}_{t+1} \).  
  _Intuition_: Current inflation rises with expected future inflation and falls with the **price gap** \( x_t \) being closed (price adjustment pressure).

### (C) Capital & Investment block

- **Q equation**: \( q*t = \varphi (i_t - k*{t-1}) \). Investment adjustment costs make \( q_t \) co-move with marginal installation costs.
- **Capital accumulation**: \( k*t = \delta i_t + (1-\delta)k*{t-1} \).

### (D) Financial accelerator & Policy block

- **Return to capital (gross)**:  
  \( rk*t = (1-\epsilon)\big(y_t - k*{t-1} - x*t\big) + \epsilon q_t - q*{t-1} \).  
  This is a reduced-form log-linearization: returns rise with marginal product of capital and capital gains \( q*t - q*{t-1} \).
- **Spread / premium mapping**:  
  \( rk*{t+1} = r_t - \nu \big(n_t - q_t - k_t\big) \).  
  Since \( lev_t = q_t + k_t - n_t \), this is \( rk*{t+1} = r_t + \nu \, lev_t \): **higher leverage ⇒ bigger premium**.
- **Net worth dynamics**:  
  \( n*t = \gamma\,RRks\,(rk_t - r*{t-1}) + r*{t-1} + n*{t-1} + en*t \).  
  Net worth grows with the \_excess* return on capital, carries over past net worth, and is buffeted by `en_t`.
- **Monetary policy (Taylor with smoothing)**:  
  \( rn*t = \rho\, rn*{t-1} + \zeta\, \text{pie}_{t-1} + s_{rn}\, etamp_t \).
- **Fisher relation**: \( rn*t = r_t + \text{pie}*{t+1} \). ⇒ \( r*t = rn_t - \mathbb{E}\_t[\text{pie}*{t+1}] \).

These blocks close the model with the shock processes above and the measurement equations for observables.

---

## 3) Log-linearization & timing — quick tutorial

- **Log-linearization**: For a level variable \( X*t \) with steady state \( \bar X \), define \( x_t \equiv \log X_t - \log \bar X \). For a gross rate \( R_t \) near \( \bar R \), write \( r_t \approx (R_t-\bar R) \) in percentage points. Your file uses the \_linearized* forms directly (no bars shown).
- **Timing**: \( k*{t-1} \) is predetermined (installed last period). \( q_t \) is the shadow value of capital **at t**. In the spread equation, \( rk*{t+1} \) on the LHS and \( r*t \) on the RHS reflect that the \_entrepreneur chooses at t* the portfolio that pays at \( t+1 \); the relevant safe rate is known at t (subject to policy rule and expectations).

---

## 4) Equation-by-equation (intuition + math detail)

We restate each equation as in your file and unpack its meaning and derivation sketch.

### 4.1 Resource constraint

**File:** \( y_t = c_y c_t + i_y i_t + g_y g_t + ce_y ce_t \).  
**Meaning:** Output gap is the share-weighted sum of (log) gaps in components. The weights \( c_y, i_y, g_y, ce_y \) are steady-state shares that sum to ~1. Linearization of \( Y = C+I+G+CE \) around the steady state yields the combination above.

### 4.2 Euler equation (consumption)

**File:** \( c*t = -\big(rn_t - \text{pie}\_t - nr_t\big) + c*{t+1} \).  
**Meaning:** Standard NK Euler in deviation form: higher **ex-ante real rate** relative to the **natural** real rate depresses \( c_t \). The \( nr_t \) term recovers the efficient benchmark (what consumption growth would be under flexible prices/optimal policy), so \( c_t \) is the gap from that natural benchmark.

### 4.3 Entrepreneur consumption

**File:** \( ce_t = n_t \).  
**Meaning:** In BGG, entrepreneurs consume out of net worth (up to survival). Linearization gives a tight link between \( ce_t \) and \( n_t \). It also ensures entrepreneurial consumption is procyclical when balance sheets strengthen.

### 4.4 Risk spread mapping

**File:** \( rk\_{t+1} = r_t - \nu\,(n_t - q_t - k_t) \equiv r_t + \nu\,lev_t \).  
**Meaning:** External finance premium increases with **leverage**. If net worth \( n_t \) falls relative to asset value \( q_t + k_t \), leverage rises and so must the premium to entice lenders (moral hazard/monitoring à la BGG). Coefficient \( \nu > 0 \) is the local elasticity.

### 4.5 Return to capital

**File:** \( rk*t = (1-\epsilon)(y_t - k*{t-1} - x*t) + \epsilon q_t - q*{t-1} \).  
**Meaning:** A reduced-form return combining marginal product (roughly \( y - k \) net of markups) and capital gains on \( q \). Parameter \( \epsilon \in (0,1) \) splits cash-flow vs. valuation components.

### 4.6 Price of capital (Tobin’s Q)

**File:** \( q*t = \varphi\,(i_t - k*{t-1}) \).  
**Meaning:** With quadratic adjustment costs, marginal installation cost moves with investment relative to existing capital. Higher \( i_t \) pushes up \( q_t \).

### 4.7 Production technology

**File:** \( y*t = a_t + \alpha k*{t-1} + (1-\alpha)\Omega h_t \).  
**Meaning:** Cobb–Douglas with effective labor \( \Omega h_t \). Linearization of \( y = a + \alpha k + (1-\alpha) h \) plus labor efficiency factor \( \Omega \).

### 4.8 Labor FOC

**File:** \( y_t = x_t + c_t + (1+\eta)h_t \).  
**Meaning:** Equates (log) marginal product of labor to real wage adjusted for markup and consumption (wealth effect). The Frisch elasticity \( 1/\eta \) governs hours’ sensitivity.

### 4.9 Phillips curve

**File:** \( \text{pie}_t = -\kappa x_t + \beta \text{pie}_{t+1} \).  
**Meaning:** Sticky prices with Calvo-type friction: current inflation rises with expected future inflation and falls as the markup gap is closed. \( \kappa \) collects structural parameters (Calvo probability, elasticity of substitution, etc.).

### 4.10 Capital accumulation

**File:** \( k*t = \delta i_t + (1-\delta)k*{t-1} \).  
**Meaning:** Linearization of \( K*t = (1-\delta)K*{t-1} + I_t \) under log deviations with the share \( \delta \) scaling the mapping from investment to capital growth (a standard linearization shortcut).

### 4.11 Net worth law of motion

**File:** \( n*t = \gamma\,RRks\,(rk_t - r*{t-1}) + r*{t-1} + n*{t-1} + en*t \).  
**Meaning:** Entrepreneurs survive with probability (or retention) \( \gamma \); net worth accumulates with the **excess** return \( (rk_t - r*{t-1}) \) on the installed capital and earns the safe rate on retained wealth; shock `en_t` perturbs balance sheets directly. This condenses the BGG cross-sectional contract outcome into a tractable linear process.

### 4.12 Monetary policy rule (no FS term)

**File:** \( rn*t = \rho rn*{t-1} + \zeta \text{pie}_{t-1} + s_{rn} etamp_t \).  
**Meaning:** Smoothing (\( \rho \)) plus lagged inflation response (\( \zeta \)). No explicit response to credit/leverage in this version.

### 4.13 Fisher relation

**File:** \( rn*t = r_t + \text{pie}*{t+1} \).  
**Meaning:** Definition of the ex-ante real rate. Combined with the rule it implies \( r*t = \rho rn*{t-1} + \zeta \text{pie}_{t-1} + s_{rn} etamp*t - \mathbb{E}\_t[\text{pie}*{t+1}] \).

### 4.14 Shock processes

**File:**  
\( a*t = \rho_a a*{t-1} + s*a ea_t \), \( g_t = \rho_g g*{t-1} + s*g eg_t \), \( z_t = \rho_z z*{t-1} \),  
\( nr_t = (1-\rho_a) a_t + (1-\rho_z) z_t + etarn_t \).  
**Meaning:** Stationary AR(1) drivers and a natural-rate index tied to productivity/preferences plus an iid shock.

### 4.15 Measurement equations

**File:** \( rAnn_t = rAnnSS + 4 rn_t \), \( piAnn_t = 4 pie_t \), \( rnAnn_t = rAnnSS + 4 nr_t \).  
**Meaning:** Map quarterly model objects to annualized observables. \( rAnnSS = 100((1/\beta)^4-1) \).

---

## 5) Calibration & economic roles of parameters

| Parameter                 |                  Value | Role / Interpretation                                                                  |
| ------------------------- | ---------------------: | -------------------------------------------------------------------------------------- |
| \( c_y, i_y, g_y, ce_y \) | 0.51, 0.18, 0.20, 0.12 | Steady-state shares in resource constraint                                             |
| \( \nu \)                 |                    0.2 | Elasticity of spread to leverage (financial accelerator strength)                      |
| \( \epsilon \)            |                   0.96 | Weight on capital gains in \( rk_t \) vs cash flows                                    |
| \( \rho \)                |                    0.9 | Policy-rate smoothing                                                                  |
| \( \rho_a, \rho_g \)      |             0.95, 0.95 | Persistence of productivity and government spending                                    |
| \( \zeta \)               |                    1.1 | Inflation response in rule (Hawkish if \( >1 \))                                       |
| \( \kappa \)              |                 0.0858 | NKPC slope (higher ⇒ stronger inflation-output tradeoff)                               |
| \( \Omega \)              |                   0.02 | Labor-efficiency weight in production                                                  |
| \( \eta \)                |                   0.33 | Inverse Frisch elasticity (labor supply curvature)                                     |
| \( \beta \)               |                   0.99 | Discount factor ⇒ \( rAnnSS \approx 4\cdot 100((1/0.99)^1-1) \) annualized steady rate |
| \( \delta \)              |                  0.025 | Depreciation (quarterly)                                                               |
| \( \gamma \)              |                 0.9728 | Entrepreneur survival/retention                                                        |
| \( RRks \)                |                   2.02 | Scaling for the effect of excess return on net worth                                   |
| \( s*a, s_g, s*{rn} \)    |           1, 1, 0.0625 | Shock standard deviations (scales in linear law)                                       |
| \( \varphi \)             |                   0.25 | Investment adjustment-cost scaling                                                     |
| \( \alpha \)              |                   0.35 | Capital share in production                                                            |
| \( \rho_z \)              |                    0.5 | Preference shock persistence                                                           |

**Comments:**

- \( \nu \) and \( \gamma RRks \) jointly shape the **amplification** and **persistence** of financial cycles. Larger \( \nu \) ⇒ spreads react more to leverage ⇒ stronger accelerator.
- \( \epsilon \approx 1 \) makes **valuation (Q)** movements dominate \( rk_t \); if you want returns driven more by cash flows, reduce \( \epsilon \).
- \( \zeta>1 \) helps determinacy (Taylor principle) even with lagged inflation.
- \( \delta=0.025 \) corresponds to ~10% annual depreciation.
- All AR(1) coefficients are <1 ⇒ **stationarity**.

---

## 6) How this differs from canonical BGG (1999)

1. **Spread equation (reduced-form):**  
   Canonical BGG derives the premium from a costly-state-verification contract. Locally, the premium rises with leverage. Your mapping  
   \( rk\_{t+1} = r_t + \nu\,lev_t \)  
   is a **parsimonious linear** counterpart capturing the same comparative statics without the full contract algebra.

2. **Net worth law of motion:**  
   BGG has \( N*t \) evolving with entrepreneurs’ survival rate and the realized returns on capital vs. the risk-free rate. Your  
   \( n_t = \gamma RRks (rk_t - r*{t-1}) + r*{t-1} + n*{t-1} + en_t \)  
   collapses these channels into a single linear equation with a shock `en_t`. It is **consistent in spirit** and works well in a linear solution framework.

3. **Monetary policy rule (no FS term):**  
   Some modern variants (e.g., with macroprudential flavor) include \( lev*t \) or credit growth in the rule. Your version **explicitly omits** this, making it a clean benchmark to later test \_financial-stability-augmented* rules.

4. **Return to capital spec:**  
   Your \( rk_t \) uses a convex combination of cash-flow and valuation components via \( \epsilon \). Standard BGG focuses on the stochastic return on capital from production/price dynamics; your form is an admissible linearization that is convenient when using a Q-type investment block.

**Bottom line:** The economic mechanisms align with BGG. The linear mappings you use are **valid local approximations** and are common in quantitative implementations when one wants tractability and easy estimation/simulation.

---

## 7) Internal consistency & determinacy (quick assessment)

- **Stationarity:** All AR(1) shocks are stationary; capital and net worth dynamics are stable for the stated parameter ranges.
- **No redundant equations:** Even though both **Taylor** and **Fisher** include \( rn_t \), Fisher simply **defines** \( r_t \) given \( rn_t \) and inflation expectations—this does **not** overdetermine the system.
- **Taylor principle:** With \( \zeta>1 \) (albeit on **lagged** inflation) and strong smoothing \( \rho=0.9 \), determinacy is typically satisfied in NK models. If indeterminacy appears in solution, consider adding a (small) response to **current** or **expected** inflation.
- **Financial block:** Signs are economically coherent: higher leverage \( lev*t \uparrow \Rightarrow rk*{t+1} \uparrow \Rightarrow excess return rises \Rightarrow \( n\_{t+1} \) improves (after survival), which feeds back by lowering leverage next round (**accelerator** with mean reversion). Magnitudes (\( \nu, \gamma RRks \)) control amplification/persistence.

---

## 8) Full model in compact linear state-space form (conceptual)

The model solver (e.g., MAPS) will stack endogenous variables into a state vector \( x*t \) and write:
\[
x*{t+1} = B\,x*t + \Phi\,\varepsilon*{t+1}, \qquad
y*t = F\,x_t,
\]
with \( y_t \) collecting observables (e.g., \( rAnn_t, piAnn_t, rnAnn_t \)), and \( \varepsilon*{t+1} \) the shocks (`ea`, `eg`, `etamp`, `en`, `etarn`). Your `unpack_model` returns these matrices as `B`, `PHI`, `F` for simulation and policy experiments.

---

## 9) Extending toward financial-stability policy (for later)

A natural next step is to augment the policy rule:
\[
rn*t = \rho rn*{t-1} + \zeta \text{pie}_{t-1} + \phi_{lev}\,lev*t + s*{rn} etamp_t,
\]
or use a macroprudential instrument (e.g., capital requirement or LTV cap) reacting to leverage or credit. This can be tested against your current benchmark to quantify stabilization trade-offs (inflation vs. financial volatility).

---

## 10) Quick reference — equations exactly as in the file

- **Resource constraint**: \( y_t = c_y c_t + i_y i_t + g_y g_t + ce_y ce_t \)
- **Euler**: \( c*t = -\big(rn_t - \text{pie}\_t - nr_t\big) + c*{t+1} \)
- **Entrepreneur cons.**: \( ce_t = n_t \)
- **Spread**: \( rk\_{t+1} = r_t - \nu (n_t - q_t - k_t) \)
- **Return to capital**: \( rk*t = (1-\epsilon)(y_t - k*{t-1} - x*t) + \epsilon q_t - q*{t-1} \)
- **Q**: \( q*t = \varphi (i_t - k*{t-1}) \)
- **Production**: \( y*t = a_t + \alpha k*{t-1} + (1-\alpha)\Omega h_t \)
- **Labor FOC**: \( y_t = x_t + c_t + (1+\eta)h_t \)
- **NKPC**: \( \text{pie}_t = -\kappa x_t + \beta \text{pie}_{t+1} \)
- **Capital accum.**: \( k*t = \delta i_t + (1-\delta)k*{t-1} \)
- **Net worth**: \( n*t = \gamma RRks (rk_t - r*{t-1}) + r*{t-1} + n*{t-1} + en_t \)
- **Taylor rule**: \( rn*t = \rho rn*{t-1} + \zeta \text{pie}_{t-1} + s_{rn} etamp_t \)
- **Fisher**: \( rn*t = r_t + \text{pie}*{t+1} \)
- **Shocks**: \( a*t = \rho_a a*{t-1} + s*a ea_t \), \( g_t = \rho_g g*{t-1} + s*g eg_t \), \( z_t = \rho_z z*{t-1} \), \( nr_t = (1-\rho_a)a_t + (1-\rho_z)z_t + etarn_t \)
- **Leverage**: \( lev_t = q_t + k_t - n_t \)
- **Measurement**: \( rAnn_t = rAnnSS + 4 rn_t \), \( piAnn_t = 4 pie_t \), \( rnAnn_t = rAnnSS + 4 nr_t \); \( rAnnSS = 100((1/\beta)^4-1) \).

---

## 11) Practical tips for simulation

- Keep the shock sizes moderate when simulating with a **zero lower bound** (ZLB) or an **occasionally binding** constraint—large negative natural-rate shocks can cause QP solvers to fail (in ELB routines).
- If convergence issues arise, try: slightly lower \( \nu \) (weaker accelerator), reduce \( \epsilon \) (less return volatility from Q), or increase policy inertia \( \rho \) marginally (smoother rate path).

---

알겠습니다. NK-DSGE 모형(BGG 기반)의 각 방정식에 대한 직관적 및 수학적 해석, 전체 모형의 작동 메커니즘, 그리고 자연이자율 충격 및 제로하한 정책 대응 시 반응까지 포함한 자세한 해설을 작성하여 Markdown 파일로 정리하겠습니다.

작성된 문서는 다음을 포함할 예정입니다:

- 전체 모형 개요
- 각 방정식의 경제적·수학적 의미
- 자연이자율 충격 및 제로금리 고정(선제적 지침 하) 시 동태 반응 분석
- 정리된 내용을 담은 Markdown 문서 파일

완성되면 바로 전달드리겠습니다.

# BGG 기반 신케인즈 DSGE 모형: 방정식별 구조와 작동 원리

## 개요 (Overview)

이 모형은 **Bernanke-Gertler-Gilchrist (BGG)**의 금융가속기(framework)를 신케인즈(새로운 케인즈 학파) DSGE 모델에 통합한 변형된 버전입니다. 즉, 명목경직성(가격이 한동안 고정되는 특성)으로 인해 단기적으로 산출갭과 인플레이션이 발생하며, 한편으로는 금융마찰(신용시장 마찰)을 통해 경기변동이 **금융가속기 효과**로 증폭되는 구조를 갖고 있습니다. 모형에는 가계(소비자), 기업(생산자), 기업가(사업주)라는 세 유형의 경제주체가 존재하고, 통화정책을 수행하는 중앙은행도 포함됩니다. 가계는 합리적 기대하에 소비와 저축을 결정하며 **Euler 방정식**을 만족시킵니다. 기업은 노동과 자본을 투입해 재화를 생산하고, 가격경직성 때문에 단기적으로 초과이윤(마크업)이 변동하며 **필립스 곡선** 관계를 따릅니다. 기업가(사업주)는 자신의 **순자산(net worth)**을 바탕으로 자본을 투자하고 부족한 자금을 차입하여 생산자본을 구매하며, 이때 차입여건이 순자산 규모에 영향을 받아 **위험프리미엄(외부자금 조달비용 가산금리)**이 형성됩니다. 이러한 금융부문의 작동이 BGG 모델의 핵심으로, 순자산과 레버리지(leverage)의 변동이 자본투자와 산출에 피드백되어 경기변동을 **증폭**시키는 메커니즘을 만들어냅니다. 아래에서는 모형의 각 핵심 방정식을 하나씩 살펴보며 그 **수학적 의미와 경제적 직관**을 설명하고, 나아가 이러한 방정식들이 결합되어 모형이 어떻게 작동하는지 종합적으로 정리하겠습니다. 마지막으로, **자연이자율에 큰 폭의 음(-)의 충격**이 발생한 경우 이 모형에서 어떤 동태적 반응이 나타나는지, 그리고 **명목금리를 장기간 0%로 고정하는 선제적 지침(Forward Guidance)**이 있을 때 모형의 반응이 어떻게 달라지는지를 구체적으로 논의하겠습니다.

## 주요 방정식별 설명 (Key Equations and Intuition)

### 1. 자원제약식 (Resource Constraint)

**방정식:** \(y*t = c_y \, c_t + i_y \, i_t + g_y \, g_t + c*{ey} \, ce*t\).  
이 방정식은 **국민소득의 사용처**를 나타내는 경제의 총자원 제약식입니다. \(y_t\)는 산출(output) 혹은 **산출갭**을 의미하고, \(c_t, i_t, g_t, ce_t\)는 각각 **소비**, **투자**, **정부지출**, **기업가(사업주) 소비**를 나타냅니다. \(c_y, i_y, g_y, c*{ey}\)는 이들 항목이 총산출에서 차지하는 **스테디스테이트 비중**입니다. 제시된 값에 따르면 소비가 약 51%, 투자 18%, 정부지출 20%, 기업가 소비 11% 정도의 비중을 갖으며 합이 1이 됩니다. 자원제약식의 의미는 **“산출 = 총수요의 합”**이라는 것으로, 한 경제 내에서 생산된 재화와 서비스의 양 \(y_t\)는 궁극적으로 민간소비, 투자, 정부지출, 기업가의 소비로 **배분**되어 사용됨을 나타냅니다. 이 식은 모형이 **1차 선형화(로그 선형화)**되어 있기 때문에 각 변수는 편차(예: 로그편차)로 해석되지만, 계수 \(c_y, i_y, ...\)는 해당 항목의 정상상태 분배를 보여주어, **어떤 수요항목이 산출 변동에 크게 기여하는지**를 직관적으로 알려줍니다. 예를 들어, 동일한 1%p 충격이라도 소비(\(c_t\)) 충격은 비중이 큰 만큼 \(y_t\)에 더 큰 영향(약 0.51)으로 반영되고, 투자 충격은 상대적으로 영향이 작게 반영됩니다.

### 2. 소비자 Euler 방정식 (Consumption Euler Equation)

**방정식:** \(c*t = -\big( r^n_t - \pi*{t+1} - nr*t \big) + c*{t+1}\).  
이 식은 **소비를 결정하는 Euler 방정식**의 선형화 형태로, **가계의 최적화 조건**을 나타냅니다. 한마디로 **현재 소비와 미래 소비 간 선택**에서, **실질이자율**이 어떻게 작용하는지를 보여줍니다. 좌변 \(c*t\)는 현재 시점 소비(의 편차)이고, 우변은 \(c*{t+1}\) (다음 기한의 기대 소비) **및 이자율 항**으로 구성되어 있습니다. \(r^n*t\)는 **명목 정책금리**(단기이자율)이고 \(\pi*{t+1}\)는 **다음 기의 기대 인플레이션**이며 \(nr_t\)는 **자연실질이자율**(natural real interest rate)입니다. 자연이자율 \(nr_t\)는 물가가 완전히 신축적이고 산출갭이 0인 **잠재균형**에서의 실질금리로, 경기상황과 쇼크에 따라 변동하는 **중립금리**입니다. Euler 방정식의 구조는 다음과 같습니다:

- \(r^n*t - \pi*{t+1}\)은 예상 **실질이자율**(현재 명목금리에서 기대 인플레이션을 뺀 값)입니다.
- \(r^n*t - \pi*{t+1} - nr_t\)는 **실질이자율 갭**으로, 현재 경제에서 적용되는 실질이자율과 경제의 중립적인 자연실질이자율 간의 차이를 나타냅니다.

이 방정식은 \(c*t - c*{t+1} = -(r^n*t - \pi*{t+1} - nr_t)\)와 동치이며, 이를 해석하면 **실질이자율이 자연이자율보다 높으면 현재 소비가 줄어들고**, 반대로 실질이자율이 자연이자율보다 낮으면 현재 소비가 늘어난다는 의미입니다. 이는 **소비-저축 선택**의 기본 원리를 반영한 것으로, 중앙은행이 금리를 인하하여 (명목금리를 낮추고 인플레이션 기대를 높이면) 실질이자율이 떨어지면 가계는 저축보다 소비를 늘이고, 금리를 올려 실질이자율이 높아지면 소비를 뒤로 미루어 저축을 늘린다는 설명과 일치합니다. 결과적으로 Euler 식은 **통화정책**(금리 조정)과 **가계소비** 간의 직접적인 연결 고리이며, 나아가 산출갭 \(y_t\) 변동의 핵심적 전파 경로입니다. 특히 \(nr_t\) (자연이자율) 항의 존재는 경제에 **수요충격**이나 **선호충격**이 발생하여 자연금리가 변할 때, 실제 금리가 이에 얼마나 따라가지 못하면 소비와 산출에 갭이 생기는지 보여주는 중요한 요소입니다.

### 3. 기업가 소비 방정식 (Entrepreneurial Consumption Equation)

**방정식:** \(ce_t = n_t\).  
이 식은 **기업가(사업주)의 소비** \(ce_t\)가 **기업가의 순자산** \(n_t\)와 같음을 나타냅니다. 여기서 \(n_t\) (net worth)은 기업가들이 보유한 자기자본의 실질가치로, **이전 시점까지 축적된 부(富)**를 의미합니다. 기업가 소비를 이렇게 설정한 배경에는 BGG 모형의 **기업가 행동 가정**이 자리잡고 있습니다. BGG 금융가속기 모형에서는 기업가가 무한정 자본을 축적하지 못하도록, 매 기별로 일부 기업가는 **퇴장(파산 또는 은퇴)**하고 해당 기업가의 순자산을 즉각 소비해버린다고 가정합니다. 위 식 \(ce_t = n_t\)는 **기업가 부문의 예산 제약**을 단순화한 것으로, “해당 기간에 퇴장하는 기업가들이 자기 순자산 전부를 소비한다”는 의미로 해석할 수 있습니다. 예를 들어 기업가들의 일정 비율(여기서 \(1-\gamma\), 약 2.72%)가 매기 퇴장한다고 할 때, 그들이 남긴 순자산이 모두 소비로 (가계나 새로운 기업가에게) 이전됩니다. 따라서 기업가 부문 전체로 보면, 각 시점 기업가 소비 \(ce_t\)의 증감은 **순자산 \(n_t\)의 증감과 동일**하게 움직입니다. 경제적 직관으로는, 경기 호황 등으로 기업가들의 순자산이 늘어나면 (예: 투자 수익이 커지면) 그 중 일부가 즉각 **소비**로 이어지고, 반대로 순자산이 줄어드는 상황 (예: 투자 실패나 금융손실)은 기업가들의 소비 축소로 직결됩니다. 이 관계는 금융부문의 **자기 자본 소진(decapitalization)**과 **배당/소비** 결정이 강하게 연결되어 있음을 보여주며, 나중에 나오는 \(n_t\)의 움직임이 거시경제에 피드백되는 경로 중 하나입니다. 또한 이 단순화는 기업가가 **장기 축재**하기보다 어느 정도 이윤을 즉시 소비/인출함으로써, 금융가속기 모형에서 **순자산 축적의 과도한 증가를 방지**하는 역할도 합니다 (기업가들이 너무 많은 내부자금을 갖게 되면 외부차입의 필요성이 줄어들어 금융가속기 효과가 약해지므로, 적정 수준에서 자산이 유출되게 만드는 장치입니다).

### 4. 위험스프레드 방정식 (Risk Spread Equation)

**방정식:** \(rk*{t+1} = r^n_t - \pi*{t+1} - \nu \big( n*t - q_t - k_t \big)\).  
이 식은 **기업가가 요구하는 자본 수익률** \(rk*{t+1}\)과 **안전한 실질이자율** 간의 관계를 보여주는 **위험스프레드(가산금리) 방정식**입니다. \(rk*{t+1}\)은 **다음 기회 기업가가 보유한 자본으로부터 얻는 기대수익률**을 나타내고, \(r^n_t - \pi*{t+1}\)은 위에서처럼 **실질정책이자율**입니다. 우변의 두 번째 항에는 기업가의 **순자산 대비 부채 비율**이 등장합니다. \(n_t - q_t - k_t\)는 기업가 **순자산 \(n_t\)**에서 **자본의 총가치 \(q_t + k_t\)**를 뺀 값인데, 이는 음(-)의 값을 취해 **레버리지(부채비율)**와 연결됩니다. 실제로 모형에서 **레버리지**를 \(lev_t = q_t + k_t - n_t\)로 정의하고 있으므로, \(n_t - q_t - k_t = -\,lev_t\)입니다. 따라서 위 식은 다음과 같이 다시 쓸 수 있습니다:

\[ rk*{t+1} = (r^n_t - \pi*{t+1}) + \nu \, lev_t. \]

이는 **“기업가가 요구하는 자본의 기대수익률 = 무위험 실질금리 + 레버리지에 따른 위험프리미엄”** 형태로 해석됩니다. 여기서 \(\nu\)는 **금융가속기 강도**를 나타내는 계수로, **레버리지 변화가 위험스프레드에 미치는 민감도**를 의미합니다. 주어진 \(\nu = 0.2\)은 **순자산 대비 부채비율 1%p 증가가 향후 자본수익률 요구치 약 0.2%p 상승을 가져온다**는 식으로 해석할 수 있습니다 (물론 선형화된 근사치 개념). 경제적 의미는 명확합니다: **기업가의 순자산이 부족하고 부채의존도가 높을수록(레버리지↑)** 대출자(금융기관)는 대출상환 불이행 위험이 커졌다고 보기 때문에, 기업가는 투자 프로젝트에 대해 **더 높은 기대수익률**을 내지 않으면 자금을 조달하기 어려워집니다. 즉, **외부자금 조달비용의 가산금리(External Finance Premium)**가 상승하는 것이죠. 반대로 기업가가 자기자본(내부자금)을 충분히 보유하면(레버리지↓) **위험스프레드가 축소**되어 비교적 낮은 금리로도 투자자금을 확보할 수 있게 됩니다. 이 관계를 통해 금융가속기 메커니즘이 작동하는데, 예를 들어 **경기침체기**에는 기업가 순자산 \(n*t\)가 타격을 입어 감소하고 레버리지가 올라가므로 위험프리미엄 \(\nu \, lev_t\)도 커집니다. 이는 투자에 필요한 자본조달 비용 \(rk*{t+1}\)을 높여 **투자위축**을 야기하고, 이는 다시 경기침체를 심화시켜 순자산을 더 떨어뜨리는 **악순환**이 만들어집니다. 반대로 **호황기**에는 순자산이 늘고 레버리지가 낮아져 위험스프레드가 축소되고, 낮은 차입비용이 투자와 생산을 더욱 촉진하여 순자산이 증가하는 **선순환**이 생깁니다. 이렇게 **순자산 - 레버리지 - 위험프리미엄 - 투자/산출**의 연결고리가 피드백을 형성하는 것이 금융가속기의 핵심이며, 위 방정식이 바로 그 정량적 관계를 간단히 요약한 것입니다.

### 5. 자본수익률 산정식 (Return to Capital Equation)

**방정식:** \(rk*t = (1-\epsilon_i)\big( y_t - k*{t-1} - x*t \big) + \epsilon_i \big( q_t - q*{t-1} \big).\)  
이 식은 **기업가가 보유한 자본의 현재기 수익률** \(rk*t\)를 구성하는 요소를 보여줍니다. 기업가는 이전 시점에 \(q*{t-1}\)의 가격으로 \(k\_{t-1}\)만큼의 자본을 매입하여 보유하고 있다가, 현재 시점에 그 자본을 활용해 생산에 투입하고 (임대하거나 직접 사용하고) 난 후 다시 시장에 팔게 됩니다. 그러면 두 가지 형태의 수익이 발생할 수 있습니다:

- 하나는 자본을 생산에 투입하여 얻게 되는 **임대소득 또는 생산수익**이고,
- 다른 하나는 사용 후 남은 자본을 매각할 때 발생하는 **자본가격의 변화**(자본 이득 혹은 손실)입니다.

위 방정식에서 \((1-\epsilon*i)\big( y_t - k*{t-1} - x*t \big)\) 항은 전자에 해당하고, \(\epsilon_i \big( q_t - q*{t-1} \big)\) 항은 후자에 해당합니다. 먼저 \(y*t - k*{t-1} - x*t\) 부분을 살펴보면, \(y_t\) (산출량의 편차)에서 \(k*{t-1}\) (전기 보유자본량의 편차)과 \(x*t\) (가격**마크업**의 편차)를 뺀 형태입니다. 직관적으로 이는 **1단위 자본당 얻은 실질산출(产出)** 혹은 \*\*자본의 한계생산력*(marginal product of capital)_에 대응하는 표현\*\*입니다. 마크업 \(x_t\)는 기업이 가격을 한계비용보다 얼마나 높게 책정했는지를 나타내는데, 만약 \(x_t\)가 양(+)이면 기업이 가격을 높게 받아 이윤을 많이 챙겼다는 뜻이고, 그만큼 자본소유주로서 기업가가 가져가는 임대소득(혹은 자본투입 대비 보상)은 줄어들 수 있습니다. 반대로 \(x_t\)가 음(-)이면 경쟁이 심해 가격이 한계비용 수준에 가까워졌음을 의미하고, 이때는 자본투입에 대한 대가(임대료)가 높아지게 됩니다. 따라서 \(y_t - k_{t-1} - x*t\)는 **자본 한 단위당 실질수익**이라 해석할 수 있고, \((1-\epsilon_i)\)는 그 수익이 \(rk_t\)에 기여하는 비중을 나타냅니다. \(\epsilon_i\)는 모델에서 “위험 인자” 혹은 **수익 구성비율 파라미터**로 볼 수 있는데, 값을 0.96으로 높게 잡았으므로 **자본수익의 대부분**이 두 번째 부분, 즉 자본가격 변동에서 나오도록 설정되어 있습니다. 두 번째 항 \(\epsilon_i (q_t - q*{t-1})\)는 **보유한 자본 자산의 가격 변동으로 인한 수익**입니다. 예를 들어 경기가 좋아져 현재 자본가격 \(q*t\)가 이전 대비 상승했다면 (\(q_t - q*{t-1} > 0\)), 기업가는 자본을 비싸게 팔아 그만큼 이득을 얻으므로 \(rk*t\)가 증가합니다. 반대로 자본가격이 떨어지면 \(rk_t\)는 낮아집니다. \(\epsilon_i\)가 0.96이라는 것은 **자본가격 변동이 자본수익률 변동의 핵심 요인**임을 뜻합니다. 이는 자본이 상당히 **내구재**이고 1기간만 사용하고 사라지는 투입물이 아니기 때문입니다. 실제 경제에서 기업의 자본투자로 인한 수익은 “현재 영업이익 + 자산 가치 변동”으로 나타나는데, 보통 **자산가격이 안정적일 때는 영업이익(임대료)이 주된 수입**이지만, 자산가격이 크게 움직이는 금융불안 또는 거품기에는 **자산가격 변동이 투자수익 변동의 주요인**이 될 수 있습니다. 모형은 이러한 상황을 단순화해 반영한 것입니다. 따라서 정리하면, **\(rk_t\)**는 “자본으로부터 얻는 **현재 수익**(생산참여로 얻은 소득)과 **미실현 수익**(자본자산 평가차익)의 합”을 나타내며, \(\epsilon_i\)가 그 두 가지 구성의 비율을 조절합니다. 이 방정식을 통해 기업가의 순자산에 영향을 주는 \(rk_t\)의 움직임을 설명할 수 있는데, 예를 들어 **경기침체 시**에는 산출 \(y_t\)가 감소하고 마크업 \(x_t\)는 하락(가격이 경직적이어서 생산자들이 일부 가격을 인하하거나, 수요부진으로 초과공급 상태)할 수 있습니다. \(y_t - k*{t-1} - x*t\)가 줄면 자본의 한계생산 소득이 줄어 \(rk_t\) 하락 요인이 됩니다. 게다가 침체기에는 자본재 가격 \(q_t\) 역시 떨어질 가능성이 커서 \((q_t - q*{t-1})\)도 음수가 되어 \(rk_t\)를 크게 떨어뜨립니다. 반대로 **호황 시**에는 \(y_t\)↑, \(x_t\)↑(일시적 과점력 강화) 등으로 1단위 자본생산소득이 늘고, \(q_t\) 상승까지 겹치면 \(rk_t\)가 크게 높아집니다. 이렇게 \(rk_t\)의 변동은 곧 기업가의 순자산 변동으로 이어져 (다음 식 참고) 금융가속기 경로를 통해 실물경제에 파급됩니다.

### 6. 자본 가격 방정식 (Price of Capital Equation)

**방정식:** \(q*t = \varphi \, ( i_t - k*{t-1} ).\)  
이 식은 자본의 실질가격 \(q*t\) (토빈의 Q로 해석 가능)가 어떻게 결정되는지를 보여주며, **투자와 자본 축적의 관계**를 나타내는 **투자자본 조정비용** 방정식입니다. \(\varphi\)는 **자본조정비용 계수**이고, \(i_t\)는 현재 시점의 투자, \(k*{t-1}\)은 이전 시점의 자본스톡입니다. 상식적으로 생각하면, **투자가 크게 증가**하면 (즉, 많은 신규 자본을 한꺼번에 설치하려 하면) 설치/조정 비용이 증가하여 **기존 자본의 가치**가 상대적으로 상승하게 됩니다. 반면, 투자수요가 낮을 때는 자본을 늘리는 데 큰 비용이 들지 않으므로 자본가격이 낮아집니다. 이 식을 선형화 형태로 보면, \(i*t - k*{t-1}\)은 **자본스톡 대비 투자비율의 편차** 정도로 해석할 수 있습니다. 따라서 \(i_t\)가 평소보다 많아 (자본증가 속도가 빨라) 양(+)이면 \(q_t\)도 양(+)으로 상승하고, \(i_t\)가 적으면 \(q_t\)는 음(-)으로 떨어집니다. 계수 \(\varphi=0.25\)는 **투자 변동이 자본가격에 주는 영향의 크기**를 결정하는데, 이 값이 클수록 소폭의 투자 증감에도 자본가격이 크게 움직이고, 값이 작으면 투자 증감에 대한 자본가격 반응이 작습니다. 요약하면, 이 방정식은 **토빈의 Q 이론**을 반영한 것으로, 기업이 투자 결정을 할 때 **한계투자의 비용**이 기존 자본의 가치로 나타나며, **투자가 급증하면 자본의 한계가격(한계비용)이 올라간다**는 경제적 의미를 갖습니다. 또한 이는 모형에서 **투자 지출 \(i_t\)과 자본가치 \(q_t\)**를 연결하여 투자변동이 자본가격 및 나아가 **기업가의 순자산**과 **레버리지**에도 영향을 미치는 경로를 제시합니다. 예를 들어, 기술진보나 정책 등의 영향으로 **투자 수요가 급증**하면 \(i_t\)↑로 인해 \(q_t\)↑ (자본가격 상승)이 발생합니다. 자본가격이 상승하면 기업가의 **자산 가치**가 높아져 순자산 \(n_t\)이 늘고 레버리지가 낮아지는 효과가 있으며, 이는 위험스프레드 축소를 통해 추가적인 투자촉진을 불러올 수 있습니다. 반대로 **투자 위축기**에는 \(q_t\)가 떨어져 기존 자본의 가치가 하락하고, 이는 기업가 부채비율을 높여 금융여건을 더욱 악화시킬 수 있습니다. 이렇듯 \(q_t\)는 실물투자와 금융건전성 사이를 연결하는 **매개 변수** 역할을 합니다.

### 7. 생산함수 (Production Function)

**방정식:** \(y*t = a_t + \alpha \, k*{t-1} + \Omega \, h*t.\)  
이 식은 경제의 **생산기술**을 나타내는 선형화된 **생산함수**입니다. \(a_t\)는 **총요소생산성(TFP) 쇼크**로 볼 수 있는 **생산성 수준**의 편차이며, \(k*{t-1}\)은 이전 시점 자본 투입량 (편차), \(h*t\)는 **노동 투입량**(근로시간 등)의 편차입니다. \(\alpha\)와 \(\Omega\)는 각각 **자본의 산출탄력성**과 **노동의 산출탄력성**을 나타내는 계수입니다. 주어진 값으로는 \(\alpha = 0.35\), \(\Omega = 0.64\)인데, 이는 자본과 노동의 소득분배 비중이 약 35% 대 65% 정도임을 시사합니다 (0.35 + 0.64 = 0.99로 1에 근접). 이 값들은 기업이 **규모수익 불변의 코브더글라스 생산함수** \(Y = A K^{0.35} H^{0.65}\)를 갖는다는 가정을 반영한 것으로 볼 수 있습니다 (소수점 차이는 마크업이나 기타 요소 반영에 따른 미세한 조정일 수 있습니다). 생산성 증가 \(a_t\)는 같은 자본·노동으로 더 많은 산출을 가능케 하므로 \(y_t\)를 직접 밀어올리고, 자본 투입이 많을수록 (\(k*{t-1}\)↑) 또는 노동을 많이 투입할수록 (\(h_t\)↑) 산출이 증가합니다. 이 식은 거시모형의 **공급 측면 구조**를 나타내며, 자본과 노동의 상대적 중요도를 파라미터로 명시해줍니다. 또한 **잠재산출(자연산출)**의 개념도 이 식을 통해 정의되는데, 가격이 완전히 신축적인 경제에서는 마크업이 일정(통상 0)하고 노동공급도 균형 수준에서 결정되므로, \(y_t\)는 이 생산함수로 주어지는 **최적산출량**이 됩니다. 요컨대, 이 생산함수는 **실제 산출갭 \(y_t\)**를 결정짓는 **공급능력**을 규정하며, 생산성 쇼크 \(a_t\)에 의해 **공급측 요인으로 인한 경기변동**이 발생할 수 있음을 보여줍니다. 예를 들어 \(a_t\) (기술진보)가 1% 상승하면 동일한 자본·노동으로 산출이 1% 늘어나므로 \(y_t\)가 1% 증가하고, 이는 물가하방 압력(공급충격으로 인한 인플레이션 둔화)을 주는 등 거시경제에 영향을 미칩니다.

### 8. 노동에 대한 최적조건 (FOC for Labor / 노동공급-수요 균형식)

**방정식:** \(y*t = x_t + c_t + (1+\eta) \, h_t.\)  
이 식은 **노동에 대한 1차 조건**을 나타낸 것으로, 가계의 **노동공급**과 기업의 **노동수요**가 균형을 이룰 때 성립하는 관계입니다. 이를 **노동시장 균형 조건** 혹은 **효용극대화의 한계조건**으로 이해할 수 있습니다. 직관적인 유도 과정은 다음과 같습니다. 가계는 소비 \(c_t\)와 여가의 상충관계에서 **한계대체율** 조건을 만족시키고, 기업은 노동의 한계생산물이 실질임금에 상응하게 고용을 결정합니다. 가계의 효용함수가 (일례로 \(U = \log C - \frac{1}{1+\eta}H^{1+\eta}\)) 형태라면, 노동에 대한 **한계효용 손실**과 소비의 한계효용 간 비율이 실질임금과 같아지는 조건이 성립합니다. 그 결과 **가계 측면**에서는 대략적인 관계 \(W/P = C_t \cdot H_t^{\eta}\) (실질임금 = 소비 * 노동^η)와 같은 관계가 유도되고, **기업 측면**에서는 한계생산물 조건 \(W/P = \text{MPL} \cdot MC\) (실질임금 = 한계생산물 \_ 한계비용)이나, 완전경쟁이라면 \(W/P = \text{MPL}\)이 됩니다. 여기서 \(MC\)는 **실질한계비용**, \(x_t\)는 **마크업(가격/한계비용)**이므로 \(1/x_t = MC\)로 볼 수 있습니다. 여러 식을 정리하면 (세부 유도는 생략) **노동시장 균형식**은 **산출 = 마크업 + 소비 + (1+η)\*노동투입** 형태로 나타낼 수 있는데, 그것이 바로 \(y_t = x_t + c_t + (1+\eta)h_t\)입니다. 이 식의 경제적 의미를 해석해 보면, **산출을 늘리기 위해서는 세 가지 경로가 있다는 것**을 보여줍니다:

- **마크업 조정 (\(x_t\)):** 마크업 \(x*t\)는 가격과 한계비용의 차이인데, 완전경쟁에 비해 가격이 더 높게 책정되어 있으면 (마크업↑) 생산이 억제됩니다. 따라서 \(x_t\)가 낮아지면 (즉, 기업들이 가격을 경쟁적으로 인하하거나 비용이 상승하여 마크업이 줄어들면) 산출 \(y_t\)가 늘어나는 쪽으로 유인됩니다. 식에서 \(x_t\)가 증가하면 \(y_t\)가 증가해야 등식이 맞는데, 사실 \(x_t\)는 *마크업 격차\_(desired markup과의 차이)로 이해하는 게 적절합니다. 정상상태 마크업이 일정 수준에 있고 \(x_t\)가 그 편차라면, \(x_t\)가 양(+)이라는 것은 실제 마크업이 높아 한계비용 대비 가격이 높음을 뜻해 산출이 수요 부족으로 줄어든 상황이고, \(x_t\)가 음(-)이면 마크업이 낮아 기업들이 적은 이윤을 감수하고 많이 생산하는 상황이라고 볼 수 있습니다. 따라서 \(x_t\)가 낮아지는 것은 (마크업 감소 = 한계비용 상승 또는 가격하락) 기업들이 더 많이 생산하게 만드는 효과를 냅니다. 요컨대 \(x_t\)는 **수요 측면에서 산출갭에 영향을 주는 가격변수**입니다.

- **소비 수요 (\(c_t\)):** 소비 \(c_t\)가 증가하면 (수요증가) 산출 \(y_t\)를 증가시키는 방향으로 작용합니다. 노동시장 관점에서 보면, 소비가 증가하면 가계의 노동공급 결정에 변화를 일으킵니다. 소비가 높다는 것은 **한계효용체감**으로 인해 한계효용이 낮아진 상태이므로, 가계는 여가보다 추가 소비를 얻는 것에 덜 민감해집니다. 따라서 같은 임금을 주더라도 노동공급을 늘릴 유인이 생겨, 노동투입 \(h_t\)를 증가시킬 수 있습니다. 이는 생산을 늘려 \(y_t\) 상승과 일치합니다. 식에서 \(c_t\)가 1 상승하면 \(y_t\)도 동일하게 1 상승하려면 노동투입 \(h_t\) 쪽에서 조정이 이뤄져야 합니다. 결국 소비증가는 생산증가와 일관되게 움직입니다.

- **노동투입 (\(h_t\)):** 마지막 항 \( (1+\eta) h_t \)는 노동 투입이 산출에 기여하는 부분입니다. \(\eta\)는 **프리쉬 노동공급 탄력성의 역수**(값 0.33이면 노동공급 탄력성 ~3)로, \(1+\eta\)는 **노동 투입 변화가 산출에 미치는 효과**를 나타냅니다. 노동 1 단위가 증가하면 일단 산출이 증가하는 직접효과(계수 1)가 있고, 추가로 \(\eta\)만큼 **소비-여가간 효용교환 효과**가 반영됩니다. 즉, 노동을 더 많이 할수록 가계는 여가손실로 인한 효용감소를 상쇄하기 위해 더 많은 임금(소비)을 필요로 하므로, 균형을 맞추려면 산출이 그만큼 더 커져야 한다는 뜻입니다. 결과적으로 \(h_t\)가 1% 증가하면 \(y_t\)는 1%보다 큰 \((1+\eta)\)% 증가로 연결됩니다. 이 항은 **노동공급의 비탄력성 정도**를 반영한 것으로, \(\eta\)가 클수록 (노동공급이 비탄력적일수록) 같은 노동 증가가 산출을 크게 늘려야 균형이 잡히는, 즉 실질임금이 크게 올라야 노동을 공급한다는 것을 의미합니다.

이 식은 결국 **산출, 가격마크업, 소비, 노동 간의 상호관계**를 나타낸 것으로, 신케인즈 모형에서 **총공급 곡선의 미시적 기반**이 됩니다. 특히 마크업 \(x_t\)와 산출 \(y_t\)의 반대 방향 관계는 **필립스 곡선**의 기초가 되며, 소비와 노동의 관계는 **IS 곡선**(Euler 식)과 함께 거시적 수요-공급 균형을 묘사합니다.

### 9. 필립스 곡선 (New Keynesian Phillips Curve)

**방정식:** \(\pi*t = -\kappa \, x_t + \beta \, \pi*{t+1}.\)  
이 식은 신케인즈 모형의 **물가 설정 관계**, 즉 **뉴케인즈 필립스 곡선**의 선형 형태입니다. \(\pi_t\)는 **인플레이션 (물가상승률)**이고, \(x_t\)는 앞서 나온 **마크업 갭**입니다. \(\beta\)는 가계의 할인인자(여기선 0.99로, 기대 인플레이션에 대한 가중치 역할)이고, \(\kappa\)는 **필립스 곡선의 기울기**(물가경직성 정도에 따른 계수)입니다. 이 식을 해석하면 다음과 같습니다:

- **기대항 \(\beta \pi\_{t+1}\):** 기업이 미래의 예상 인플레이션을 고려해 현재 가격을 조정하는 **앞으로형 기대** 요소입니다. \(\beta \approx 0.99\)로 1에 매우 가깝기 때문에, 미래 인플레이션 전망이 높으면 현재 인플레이션도 높아지는 경향이 큽니다. 이는 합리적 기대 하에서 **인플레이션의 관성**을 나타내는 부분입니다.

- **실물경기 압력항 \(-\kappa x_t\):** 현재 마크업 갭 \(x_t\)가 인플레이션에 미치는 영향입니다. \(\kappa\)값은 0.04로 비교적 작게 설정되어 있는데, 이는 **가격경직성이 커서** 산출이나 마크업 변화가 물가에 완만한 영향을 준다는 뜻입니다. 식의 부호가 음(-)인 점에 주목해야 합니다. \(x_t\)가 **양(+)**이라는 것은 기업들이 현재 **마크업을 높게** 가져가고 있다는 의미인데, 마크업이 높다는 것은 **한계비용에 비해 가격을 많이 올려받는 상황**이라서 비용측면의 물가 상승 압력이 낮습니다. 따라서 \(x_t\)가 높으면 (즉, 기업이 이윤을 많이 붙이고 있다는 것은 비용상승이 적거나 수요가 낮다는 신호) **인플레이션 \(\pi_t\)**는 오히려 낮아지는 방향이므로 \(-\kappa x_t\)가 음(-)으로 작용합니다. 반대로 \(x_t\)가 **음(-)**이면 기업들이 이윤을 줄이고 있다는 뜻이고, 이 경우 생산 원가 대비 가격을 낮게 책정하므로 마진이 줄어든 만큼 **한계비용 상승** 또는 **수요 과열**을 의미하여 **인플레이션을 높이는 효과**를 냅니다. 요컨대 \(x_t\)는 **경기의 과열/침체 정도**를 나타내는 **비용압력 지표**로 볼 수 있고, 경기가 좋아 마진이 줄어들고 비용이 상승하면 (\(x_t\)↓) 물가가 오르고, 경기가 나빠 마진이 늘어나면 (\(x_t\)↑) 물가가 떨어지는 식입니다. 이는 전통적인 **산출갭-인플레이션 간 반비례 관계**와 맥락을 같이 합니다. 실제로 \(x_t\)는 기업의 마진이니까, 그것을 일정 목표 대비한 격차로 보면 **(음의) 산출갭**과 비례한다고 할 수 있습니다. 산출이 잠재치를 웃돌 때 (\(y_t\) 양호) 마크업이 하락(\(x_t\)↓)해서 인플레이션이 올라가고, 산출이 부진할 때 (\(y_t\) 침체) 마크업이 높아져(\(x_t\)↑) 인플레이션이 내려가는 현상은 **필립스 곡선**의 핵심이기도 합니다.

이 필립스 곡선 식은 **현재 인플레이션**이 **미래 인플레이션 기대**와 **현재 경기상황(마크업/산출갭)**에 의해 결정됨을 보여줍니다. \(\kappa\) 값이 작다는 것은 가격 경직성이 상당히 크다는 뜻인데, 예를 들어 \(\kappa = 0.04\)이라면 산출갭(혹은 \(x_t\)) 1%의 변화가 즉각 인플레이션을 0.04%p 밖에 못 움직일 정도로 **물가 반응이 둔감**하다는 의미입니다. 이는 기업들이 가격을 자주 조정하지 못하거나, 메뉴비용 등으로 인해 경기변동이 당장의 가격결정에 큰 영향을 못 미치는 상황을 반영합니다. 따라서 통화정책 등으로 수요를 조절해도 인플레이션에 **시차를 두고 서서히** 영향이 나타나는 구조라고 할 수 있습니다. 이 식은 향후 자연이자율 충격 시나리오에서 **디플레이션 압력**이나 **인플레이션 기대**가 어떤 역할을 하는지 분석할 때 핵심적으로 작용합니다.

### 10. 자본 축적식 (Capital Accumulation Equation)

**방정식:** \(k*t = \delta \, i_t + (1-\delta) \, k*{t-1}.\)  
이 식은 **자본의 축적(dynamic)**을 나타내는 표준적인 법칙으로, 투자와 기존 자본으로부터 **다음 기의 자본스톡**이 형성되는 과정을 보여줍니다. \(\delta\)는 **자본의 감가상각률**로, 모형에서는 0.025 (분기별 약 2.5%, 연율로 약 10%)로 설정되어 있습니다. 방정식을 해석하면: 현재 시점의 신규 투자 \(i*t\) 중 \(\delta\) 비율만큼이 **자본증가로 연결**되고, 이전 시점 보유자본 \(k*{t-1}\) 중 \((1-\delta)\) 비율만큼이 **감가상각 후 잔존**하여 다음 시점 자본 \(k*t\)로 이월된다는 의미입니다. 쉽게 말해, **“다음 기의 자본 = 새로 늘린 자본 + 남은 기존 자본”**입니다. 이 식은 사실 비선형 원식 \(K*{t} = I*{t} + (1-\delta) K*{t-1}\)의 선형 근사 형태입니다 (여기 변수들은 로그편차나 비율 편차 개념으로 보시면 됩니다). 자본축적식은 경제에 **투자의 동태적 효과**를 설명하는 핵심입니다. 예를 들어 현재 어떤 충격으로 투자 \(i*t\)가 증가하면, 당장은 산출\(y_t\)에 \(i_y\) 만큼 기여하여 상승시키겠지만, 이 식에 의해 \*\*미래의 자본스톡 \(k*{t+1}\)도 증가**하게 됩니다. 자본스톡 증가분은 이후 **생산능력 확충**으로 이어져 미래 산출을 높이고, 또한 추가 자본이 기업가 재무제표에 자산으로 잡혀 순자산 증가, 레버리지 개선 등의 **금융 측면 파급효과**도 갖습니다. 감가상각률 \(\delta\)가 0.025라는 것은 **자본이 비교적 천천히 축적**됨을 의미합니다. 한 기(분기)에 자본의 2.5%만 소진되므로, 투자 증가의 효과가 누적되어 **장기간에 걸쳐 자본축적**에 영향을 줍니다. 따라서 일시적 투자쇼크라도 그 여파가 서서히 사라지며, **hump-shaped** 자본과 산출 반응을 만들 수 있습니다. 이 식은 모형의 **내부 모멘텀**을 제공하는 요소이기도 하며, 특히 금융가속기 모형에서는 자본축적과 기업가 재정상태가 얽혀 **충격의 지속성\*\*을 높여주는 방향으로 작용합니다.

### 11. 순자산 진화식 (Evolution of Net Worth)

**방정식:** \(n*t = \gamma \, RR^{ks} \big( rk_t - (r^n*{t-1} - \pi*t) \big) + \big( r^n*{t-1} - \pi*t \big) + n*{t-1} + e^n_t.\)  
이 식은 **기업가 부문의 순자산 \(n_t\) (내부자금)**가 시간에 따라 어떻게 변하는지를 나타냅니다. 쉽게 말해 **“이전 시점의 순자산이 현재 어떻게 늘어나거나 줄어드나”**를 보여주는 식입니다. 우변을 몇 부분으로 나누어 살펴보겠습니다:

- \(n\_{t-1}\): **이전 기의 순자산**이 현 시점에 이월된 부분입니다. (선형화되어 있으므로 직접 더해지는 형태)

- \(\big( r^n*{t-1} - \pi_t \big)\): **기초 순자산에 대한 안전자산 수익 부분**입니다. 기업가가 지난 기에 보유한 순자산은 (만약 은행 예금 등 안전 자산에 넣어두었을 경우) 명목이자 \(r^n*{t-1}\)만큼 불어나지만, 물가상승률 \(\pi*t\)만큼 실질가치가 깎입니다. 따라서 \(r^n*{t-1} - \pi*t\)는 **실질 안전수익률**이며, 이는 이전 순자산 \(n*{t-1}\)에 곱해진 형태로 실질가치 증가분을 주어야 하지만, 선형화된 현재 식에서는 일정한 값처럼 더해지고 있습니다. 사실 정확히 말하면 \(n*{t-1}\)과 \((r^n*{t-1}-\pi_t)\)의 곱이 들어가야 하지만, 균형부근에서 정규화된 관계로 볼 수 있습니다. **결과적으로 이 항은 기업가 순자산이 안전한 곳에 굴러 얻은 증가분** 정도로 이해할 수 있습니다.

- \(\gamma \, RR^{ks} \big( rk*t - (r^n*{t-1} - \pi*t) \big)\): **기업가 투자수익의 초과부분**이 순자산에 미치는 영향입니다. \(\gamma\)는 **기업가의 생존 확률** (또는 “다음 기까지 생존하여 재투자하는 기업가의 몫”)을 나타냅니다. 주어진 \(\gamma = 0.9728\)은 약 97.28%의 기업가가 생존하고 2.72%가 퇴장한다는 의미로 해석할 수 있습니다. \(RR^{ks}\)는 **자본수익률의 정상상태 수준**으로, 1.02 (약 1/분기 2%의 수익, 연간 약 8% 수익률)로 설정되어 있습니다. 이 둘의 곱 \(\gamma \, RR^{ks}\)는 정상상태에서 **생존한 기업가들이 유지하는 자본수익 비율** 정도로 볼 수 있습니다. 이어지는 괄호 부분 \(( rk_t - (r^n*{t-1} - \pi*t) )\)는 **기업가가 실제로 얻은 자본수익률 \(rk_t\)**에서 **안전한 실질금리 비용**을 뺀 **초과수익률**입니다. 즉, 기업가가 빚을 내 자본에 투자했을 때, 그 투자로 얻은 실질 수익률 \(rk_t\)과 차입 비용(지난 기 실질이자율) 간의 차이입니다. 만약 \(rk_t\)가 과거 대출금리보다 높았다면 그만큼 **이익**을 본 것이고, 낮았다면 **손실**을 본 것입니다. \(\gamma RR^{ks}\)를 곱한 것은 그 중 **생존한 기업가들이 실제 자기 순자산으로 회수한 몫**을 의미합니다. 예를 들어, 지난 기에 기업가가 자기자본 \(n*{t-1}\)과 부채를 합쳐 투자했는데, 실제 수익률 \(rk_t\)가 높아서 대출이자 이상으로 벌었다면 그 **차익**이 순자산에 더해집니다. 반대로 \(rk_t\)가 낮아 손실을 보면 순자산이 깎입니다. 이 메커니즘이 **금융가속기의 핵심**입니다. 산출이 좋고 자본수익률이 높으면 기업가 부가 늘어 향후 더 투자할 여력이 생기고, 불황으로 자본수익률이 낮으면 순자산이 줄어 투자능력이 약화됩니다. 특히 \(rk_t\)가 높고 인플레이션이 약간 올라간 상태에서는 부채 실질부담이 줄어 기업가가 실질이익을 더 가져가는 반면, \(rk_t\) 저조와 디플레이션 상황에서는 부채 부담이 실질적으로 늘어나 순자산에 이중 타격을 줄 수 있습니다. 본 식에서 \(\nu\)나 \(\epsilon_i\) 등의 매개로 \(rk_t\)가 결국 \(n_t\)에 영향을 미치게 연결되어 있고, 그 결과는 다시 **레버리지 증가(또는 감소)** 및 **위험스프레드 변화**로 이어져 실물경기에 피드백되는 순환이 형성됩니다.

- \(e^n_t\): **순자산 쇼크**로, 외생적으로 기업가 부에 영향을 주는 요인입니다. 이는 0에 수렴하는 충격의 편차로 생각할 수 있으며, 예를 들어 갑작스런 **금융자산 가치 하락**이나 **기업가부문 부채탕감/지원** 같은 정책 등이 모형에서 랜덤하게 **순자산 수준을 변화**시키는 역할을 합니다. \(e^n_t\)가 음(-)이면 아무런 실물변화 없이도 기업가 순자산이 증발하는 **금융위기 충격**처럼 작용하여 경기침체를 야기할 수 있고, 양(+)이면 반대로 금융여건이 개선되는 효과를 냅니다.

전체적으로 이 순자산 진화식은 **기업가 대차대조표의 축소/확대를 동태적으로 묘사**하며, **금융부분의 자체 동학(dynamics)**을 결정합니다. 정상상태에서는 \(rk = r^n - \pi\) 정도로 맞춰져 있어 순자산이 경향적으로 유지되지만, 충격이 발생하면 \(rk\)와 실질금리 차이가 벌어져 \(n*t\)에 누적적인 변화를 줍니다. 예를 들어 **금융위기적 상황**을 생각해보면, 자본가격 급락으로 \(rk_t\)가 크게 떨어지고 동시에 물가는 하락(디플레이션)하여 실질부채코스트 \(r^n*{t-1}-\pi*t\)가 높아진다면, 괄호 안 \(rk_t - (r^n*{t-1}-\pi_t)\)는 큰 폭의 음(-)의 값이 되어 기업가 순자산을 크게 까먹을 것입니다. 반면, **호황**으로 \(rk_t\)가 높고 인플레이션도 안정적이면 기업가는 초과수익을 얻어 순자산이 빠르게 회복되겠죠. 이 모든 효과가 \(\gamma\) (생존율)로 적절히 조절되어 반영됩니다. \(\gamma \approx 0.97\)이 의미하는 것은 기업가 순자산 증감의 대부분이 **생존한 기존 기업가의 몫으로 남고**, 3% 정도는 퇴장 기업가의 소비로 빠져나간다는 것입니다. 그래서 순자산 증가는 한 번에 다 쌓이지 않고 일부 소실되며 (이 부분은 앞의 \(ce_t = n_t\)와 연계), 순자산 감소도 완전히 반영되진 않고 약간 완충되는 면이 있습니다. 하지만 \(\gamma\)가 1에 매우 가까워 **순자산 충격의 자기 지속성**이 상당함을 보여주는데, 이는 금융 충격이 **장기적인 경기 여파**를 가질 수 있음을 시사합니다.

### 12. 통화정책반응 함수 (Monetary Policy Rule)

**방정식:** \(r^n*t = \rho \, r^n*{t-1} + (1-\rho)\big( \zeta \pi_t + \zeta_y y_t + \zeta_f \, lev_t \big) + s^{rn} e^{tamp}\_t.\)  
이 식은 중앙은행이 어떻게 정책금리 \(r^n_t\)를 결정하는지를 나타내는 **이자율 규칙** (Taylor rule, 테일러 준칙)의 확장 형태입니다. \(\rho\)는 **이자율 결정의 관성(이월) 정도**이며, 나머지 괄호 안은 **정책반응 함수**입니다. 수치를 대입하면 \(\rho = 0.9\)로 매우 높아, **이전 분기의 금리를 90% 비중으로 이어받고**, 나머지 10%만 새로운 목표치로 조정하는 **강한 이자율 관성(smoothing)**을 보여줍니다. 이는 중앙은행이 금리를 급격히 변동시키지 않고 점진적으로 움직이는 관찰과 부합합니다. 괄호 안을 보면, 세 가지 요소에 반응하도록 되어 있습니다:

- \(\pi_t\) (인플레이션 편차): 계수 \(\zeta = 1.5\)로, 중앙은행이 **목표 인플레이션 대비 현재 인플레이션의 괴리**에 강하게 반응함을 의미합니다. 1.5라는 값은 일반적으로 알려진 **테일러 원칙**(인플레이션 1% 상승 시 명목금리를 1%보다 더 올려 실질금리를 상승시키는 원칙)을 만족하여, 중앙은행이 물가 안정을 최우선으로 한다는 신호입니다.

- \(y_t\) (산출갭): 계수 \(\zeta_y = 0.5\)로, 중앙은행이 **경기상황(산출갭)**에도 완만하게 반응함을 의미합니다. 산출갭 1%p 증가(경제 과열) 시 약 0.5%p 금리를 높이는 식으로, 경기 부양/냉각 목적의 조절입니다.

- \(lev_t\) (레버리지): 계수 \(\zeta_f = 0.5\)로 설정되어 있는데, 이는 **중앙은행이 금융상황(레버리지 수준)**에도 반응하는 것을 보여줍니다. 일반적인 테일러룰에는 없는 항목으로, **거시건전성**을 염두에 둔 확장입니다. \(lev_t = q_t + k_t - n_t\)는 **기업가 부문의 부채 의존도** 혹은 **금융취약성 지표**라 할 수 있습니다. 이 항의 계수가 양수 0.5인 것은, 예컨대 **레버리지가 높이 상승**(기업가 부채 누증)하면 중앙은행이 금리를 추가로 0.5\*Δlev 만큼 더 올려 거품 축적을 억제하려 한다는 의미입니다. 반대로 금융 디레버리징으로 레버리지가 크게 떨어진 상황에선, 이를 보완하기 위해 금리를 낮추는 방향으로 작용합니다. 모형 설명에 “doesn't include finance stability in monetary rule”이라고 되어 있으나, 실제 파라미터는 \(\zeta_f = 0.5\)로 되어 있어 **일정 정도 금융안정에 기울인 정책**을 모형화한 것으로 보입니다. 다만 그 크기가 0.5로 인플레이션 반응과 동등하게 설정된 점은, 금융 변수를 완만히 고려하기 시작한 중앙은행의 행동을 반영한 것이라 해석할 수 있습니다. (만약 \(\zeta_f = 0\)였다면 전통적 테일러룰이 됩니다.)

마지막으로 \(s^{rn} e^{tamp}\_t\)는 **통화정책 충격**입니다. \(e^{tamp}\_t\)는 정책 당국의 **비정칙적 금리조정**이나 예기치 못한 통화정책 변동을 나타내고, \(s^{rn}=0.5\)는 그 표준편차 크기 정도입니다. 이 항이 양(+)이면 예를 들어 기준금리를 갑자기 예상보다 더 올린 깜짝 긴축, 음(-)이면 예상보다 더 내린 깜짝 완화로 볼 수 있습니다.

요약하면 이 통화정책 방정식은 **중앙은행이 인플레이션을 우선시하되 산출갭과 금융상황도 고려하여 점진적으로 금리를 조정**하는 정책반응을 형식화한 것입니다. 예를 들어, **인플레이션이 목표보다 높고 경제도 과열되고 있으며 거품까지 낀 상황**이라면 \(\pi_t\), \(y_t\), \(lev_t\) 모두 양(+)일 것이고, 이 식에 의해 상당히 큰 폭으로 금리를 인상하는 결과가 나옵니다. 반대로 **심각한 경기침체와 디플레이션 위험, 그리고 금융위기(레버리지 축소)** 상황에서는 셋 다 음(-)이어서 금리를 큰 폭 인하(혹은 0에 수렴)하는 방향을 제시합니다. \(\rho = 0.9\) 때문에 이러한 조정은 단계적으로 이뤄지지만, 장기적으로는 지정된 반응계수 비율만큼 금리가 움직여 **물가와 생산, 금융을 안정화**시키려 합니다. 이 규칙을 통해 모형에서 **명목금리 \(r^n_t\)**가 내생적으로 결정되며, 이는 다시 Euler 방정식을 통해 **실질이자율 경로**를 좌우하고, 또한 기업가 차입비용에도 영향을 줍니다. 실제 경제에 대입하면 이 식은 **테일러 준칙**이 잘 작동하던 정상시기의 정책 행동을 묘사하지만, 뒤에 논의할 **자연이자율 급락 & ZLB 시나리오**에서는 중앙은행이 이 식대로 행동하기 어려워지는 특수한 상황이 나타납니다.

### 13. 생산성 충격 (Productivity Shock Process)

**방정식:** \(a*t = \rho_a \, a*{t-1} + s_a e^a_t.\)  
이 식은 총요소생산성(TFP) 수준 \(a_t\)가 어떻게 움직이는지를 나타내는 **외생적 충격 프로세스**입니다. \(\rho_a\)는 **생산성 충격의 자기상관(지속성)**으로, 값이 0.95로 높게 설정되어 있어 **생산성 수준이 상당히 지속적으로 유지**됨을 보여줍니다. 이는 기술충격이 발생하면 천천히 소멸하고 경제에 **장기적인 영향**을 줄 수 있음을 의미합니다. \(s_a e^a_t\)는 새로운 **생산성 충격**으로, \(e^a_t\)는 평균 0의 랜덤 충격이고 \(s_a\)는 그 표준편차(영향력 크기)입니다. \(s_a = 1\)로 주어져 있지만, 이 값은 단위 조정(scale)일 뿐 실제 모형에서는 충격 분산에 맞게 조정됩니다. 예를 들어 \(e^a_t = 0.01\) (1%p 임의 상승)이라면 \(a_t\)가 1% 상승합니다. 이 충격은 생산함수를 통해 즉각 산출 \(y_t\)을 올리고, 산출 증가 -> 물가하락 압력(필립스 곡선) -> 중앙은행 금리인하 등 여러 경로로 영향을 줍니다. 또한 생산성 상승은 기업의 생산비용을 떨어뜨려 **마크업 \(x_t\)를 일정 수준 조정**하게 할 것이고, 실질임금과 여타 변수에도 파급됩니다. **요약:** \(a_t\)는 **공급측 요인**인 기술/생산성의 충격을 나타내며, 고도로 지속적(AR(1) 계수 0.95)이라서 **경제에 항구적 수준효과**(기술진보 축적)까지 줄 수 있는 요소로 취급됩니다.

### 14. 정부지출 충격 (Government Spending Shock Process)

**방정식:** \(g*t = \rho_g \, g*{t-1} + s_g e^g_t.\)  
이 식은 **정부지출 \(g_t\)**가 외생적으로 움직이는 경로를 보여줍니다. \(\rho_g = 0.95\)로 높은 지속성을 가지며, \(s_g = 1\) 크기의 랜덤충격 \(e^g_t\)에 반응합니다. 정부지출 증감은 자원제약식에서 \(g_y \approx 0.20\)의 비중으로 산출에 직접 기여하는 **총수요 충격**입니다. 예를 들어 \(e^g_t\)가 양(+)이면 정부지출이 평소보다 늘어나 \(g_t\)↑, 이는 곧바로 \(y_t\)를 증가시킵니다(정부지출乘数효과는 1차 선형 모형에서는 1보다 작거나 같지만 기본적으로 긍정). 산출갭 증가 -> 인플레이션 상승 -> 금리 인상 등 수요 충격에 준하는 파급효과가 나타납니다. \(\rho_g = 0.95\)이므로 정부지출 충격이 발생하면 여러 기간에 걸쳐 지출 수준이 서서히 이전 평균으로 돌아옵니다. 실제 정책으로 보면, 일시적인 재정지출 확대나 축소가 경제에 미치는 효과를 포착하기 위한 구성입니다. BGG 모형에서는 정부지출이 금융시장과 직접적 연관은 없지만, 산출과 인플레이션을 바꾸어 통화정책 대응, 실질이자율, 순자산 경로 등에 2차적 영향을 줄 수 있습니다.

### 15. 레버리지 정의 (Leverage Definition)

**방정식:** \(lev_t = q_t + k_t - n_t.\)  
이 정의는 앞서 여러 번 언급한 **레버리지(부채 레버리지)**의 정확한 의미를 수식으로 보여줍니다. \(lev_t\)는 기업가 부문의 **총자산 대비 자기자본의 부족분**으로, \(q_t + k_t\)는 **기업가가 보유한 총자산의 (로그) 가치**이고 \(n_t\)는 **자기자본(순자산)**입니다. 따라서 \(q_t + k_t - n_t\)는 자산에서 자기자본을 뺀 나머지, 즉 **타인자본(부채)**의 대략적인 (로그) 크기를 반영합니다. 사실 레버리지란 **총자산/자기자본 비율**을 의미하지만, 로그 선형화된 상황에서는 비율의 로그가 차이로 근사될 수 있어서 단순차 \(q_t + k_t - n_t\)로 표현하고 있습니다. \(lev_t\)가 상승하면 (자산 가치 대비 부채가 상대적으로 늘면) **재무레버리지 악화**를 뜻하고, \(lev_t\) 하락은 **레버리지 축소(디레버리지)**를 뜻합니다. 이 변수는 모형에서 **위험스프레드** 결정식에 이미 -로 나타나 있고, 통화정책반응 함수에도 직접 등장하므로 매우 중요합니다. 레버리지는 **금융위험 누적정도**를 나타내는 지표로, 예컨대 자산가격 붕괴로 \(q_t\)가 급락하거나 투자손실로 \(n_t\)가 급감하면 \(lev_t\)가 급등하고 이는 곧 신용경색과 투자위축으로 연결됩니다. 반대로 자산 가격 상승 또는 이익 증가로 자기자본이 두터워지면 \(lev_t\)가 낮아져 금융위험이 줄어듭니다. 이렇듯 \(lev_t\)는 **금융상황의 함축 요약**이라 할 수 있으며, 모형에 명시적으로 포함됨으로써 **금융안정 정책**이나 **자산시장 충격**을 분석할 수 있게 해줍니다.

### 16. 선호(소비절약 성향) 충격 (Preference Shock Process)

**방정식:** \(z*t = \rho_z \, z*{t-1}.\)  
\(z_t\)는 **가계의 시간선호나 소비성향을 나타내는 충격**으로 해석됩니다. \(\rho_z = 0.5\)이므로 꽤 빠르게 평균으로 회귀하는 충격입니다. 특이하게도 이 식에는 새로운 충격항 \(e^z_t\)가 명시되어 있지 않은데, 이는 \(z_t\) 자체보다는 곧 나오는 자연이자율 식에서 \(z_t\)의 움직임을 활용하기 위해 최소한의 AR 구조만 준 것으로 보입니다. 실제로 \(z_t\)는 Euler 방정식 등에는 직접 등장하지 않지만, **자연이자율의 변동 요인**으로 쓰입니다. 경제적 의미로 \(z_t\)가 양(+)이면 가계의 **소비성향이 높아지는** (혹은 할인율이 낮아져 미래대비 현재 소비를 선호하는) 충격이고, \(z_t\)가 음(-)이면 **소비를 억제하고 저축을 선호**하는 충격이라 볼 수 있습니다. 예를 들어 갑작스런 위험회피 성향 증대로 사람들이 소비를 줄이고자 한다면 \(z_t\)가 큰 음(-)으로 움직일 수 있습니다. 이 자체는 내생적으로 Euler 식에 영향을 줄 텐데, 본 모형에서는 직접 Euler에 등장시키는 대신, **이런 선호변화가 나타난 상황에서 자연실질이자율이 어떻게 변하는가**로 간접 나타내는 접근을 취하고 있습니다.

### 17. 자연이자율 프로세스 (Natural Rate Process)

**방정식:** \(nr*t = \rho*{nr} \, nr*{t-1} + (1-\rho*{nr}) \Big( (1-\rho*a)a_t + (1-\rho_z) z_t \Big) + e^{tar}\_t.\)  
이 식은 **자연실질이자율 \(nr_t\)**가 시간에 따라 어떻게 움직이는지를 나타냅니다. 자연이자율은 앞서 설명했듯이 **경기 균형 상태에서의 실질금리**로, 경제 내부 요인에 의해 결정되는 일종의 잠재금리입니다. 위 식은 자연이자율이 **TFP 충격**과 **선호충격**의 영향을 받으며, 자체적으로 \(\rho*{nr} = 0.9\)라는 높은 지속성을 가진 AR 과정을 따른다고 가정합니다. 자세히 보면:

- \(\rho*{nr} nr*{t-1}\): 전기의 자연이자율에 0.9의 비중을 두어 **자연이자율 자체의 관성**을 반영했습니다. 자연이자율이 한번 내려가면 꽤 오래 낮은 상태로 지속될 수 있음을 의미하는 설정입니다.

- \((1-\rho\_{nr})\big( (1-\rho_a)a_t + (1-\rho_z) z_t \big)\): 이는 자연이자율의 **장기균형 수준**이 생산성 \(a_t\)와 선호요인 \(z_t\)의 함수임을 보여줍니다. \(\rho_a = 0.95\), \(\rho_z = 0.5\)이므로 \(1-\rho_a = 0.05\), \(1-\rho_z = 0.5\입니다. 결과적으로 이 항은 \(0.1 \* (0.05 a_t + 0.5 z_t)\) = \(0.005 a_t + 0.05 z_t\)와 같습니다. 이는 **생산성 충격**과 **선호충격**이 자연이자율에 미치는 방향을 결정합니다. 부호는 둘 다 양(+)이므로:

  - \(a_t\) (생산성 상승)이 **자연이자율을 올리는 요인**입니다. 생산성이 높아져 경제가 강해지면, 잠재적 산출이 늘어나고 소비·투자 수요도 함께 늘어나 자연실질금리를 밀어올리는 것으로 이해할 수 있습니다. \(a_t\) 1%p↑ 당 \(nr_t\) 약 0.005 (0.5bp)↑로 반영됩니다.
  - \(z_t\) (현재소비 성향↑)이 **자연이자율을 올리는 요인**입니다. 선호충격 \(z_t\)는, 양(+)일 때 사람들이 현재 소비를 더 선호하므로 투자나 저축이 줄어 잠재균형에서 금리가 상승하는 효과를 냅니다. 반대로 \(z_t\)가 음(-)이면 (미래 대비 현재소비 기피, 저축선호↑) 자연금리를 끌어내리는 요인이 됩니다. \(z_t\) 1%p 변화 당 \(nr_t\) 약 0.05 (5bp) 변동이 생깁니다. \(z_t\)의 영향이 \(a_t\)보다 10배 정도 큰데, 이는 **선호(수요) 요인이 자연금리 변동에 주요**하다는 모델링입니다.

- \(e^{tar}_t\): **자연이자율 자체의 추가적인 외생충격**입니다. \(tar\)는 “target rate” 느낌이지만, 여기서는 모델링 편의를 위해 넣은 순수 자연이자율 교란항으로 보면 됩니다. \(\rho_{nr}\)이 높기에 \(e^{tar}\_t\) 충격이 발생하면 자연이자율이 상당 기간 지속적으로 바뀌는 효과를 냅니다.

이 자연이자율 프로세스는 특히 **유동성 함정**이나 **ZLB 상황**을 분석할 때 핵심적입니다. 왜냐하면 자연이자율이 큰 폭으로 하락하면, 중앙은행이 명목금리를 아무리 0 근처로 내려도 실질금리가 여전히 자연수준보다 높아져 **심각한 수요부족과 디플레이션 압력**을 초래하기 때문입니다. 여기서는 \(z*t\) 같은 **선호 충격** (예: **소비절약성향 급변** 또는 **위험회피심리 급변**)이 그런 역할을 할 수 있습니다. \(\rho*{nr}=0.9\)이므로 한 번 자연이자율이 내려가면 천천히 회복되며, 이 기간 동안 통화정책은 어려움에 처할 수 있습니다. 밑에서 구체적으로 다루겠지만, **자연이자율의 음(-)의 충격**은 모형 내에서 가계가 갑자기 소비를 줄이고 싶어하거나 투자수요가 급감하는 환경을 뜻하며, 이를 상쇄하지 못하면 경제가 침체합니다. 이 식은 바로 그런 **Wicksellian 금리** 개념을 공식화한 것으로, **실물충격들이 균형이자율을 어떻게 움직이는지**를 제시합니다.

### 18. 측정 방정식들 (Measurement Equations)

마지막으로, 모델 변수를 실제 **관측가능한 자료와 연결**시켜주는 측정식들이 있습니다:

- **연율 환산 정책금리:** \(rAnn_t = rAnn^{SS} + 4 \, r^n_t.\)  
  여기서 \(rAnn_t\)는 **연간(연율) 기준의 명목 단기금리**입니다. \(rAnn^{SS}\)는 그 **정상상태 값**인데, \(100\*((1/\beta)^4 - 1)\)로 주어져 있습니다. \(\beta = 0.99\)이면 \((1/0.99)^4 - 1 \approx 0.0407\), 즉 약 4.07%입니다. 이 값에 100을 곱했으므로 \(rAnn^{SS} \approx 4.07\%\) (연간)입니다. 다시 식을 보면, \(r^n_t\)는 분기 기준 금리의 편차(또는 실질금리 편차)인데 여기에 4를 곱하고 정상상태 연간금리를 더해줌으로써, 실제 자료의 연율 금리에 맞춥니다. 예를 들어 \(r^n_t = 0\)인 모형 상태는 연 4.07%의 정책금리에 해당합니다. \(r^n_t = -0.01\) (1%p 낮은 수준)이면 연율 4.07% - 4% = 0.07%, 즉 연 0.07%로 매우 낮은 금리가 됩니다. 이렇듯 연결해두면 모형의 결과를 중앙은행이 발표하는 기준금리나 시중단기금리 데이터와 비교할 수 있습니다.

- **연율 환산 인플레이션:** \(\piAnn_t = 4 \, \pi_t.\)  
  \(\piAnn_t\)는 연간 인플레이션 (예: CPI 연율 상승률)이고, \(\pi_t\)는 분기 인플레이션입니다. 4배를 곱해주면 단순히 분기율을 연율로 환산한 것입니다. 정상상태에서 \(\pi_t = 0\) (물가안정 목표)라면 \(\piAnn_t = 0\)로 0% 인플레이션이 기준이며, 만약 \(\pi_t = 0.005\) (분기 0.5% 상승)이라면 \(\piAnn_t = 0.02\)로 연 2% 인플레이션이 됩니다.

- **연율 환산 자연이자율:** \(rnAnn_t = rAnn^{SS} + 4 \, nr_t.\)  
  이것은 자연실질이자율 \(nr_t\)를 연율 명목금리 기준으로 환산한 것입니다. \(nr_t\)는 실질금리격 차이인데, 거기에 연환산하고 정상상태 명목금리 (4.07%)를 더해줍니다. 예컨대 \(nr_t\)가 -0.01이라면 자연실질금리가 -1%p 정도 낮다는 뜻이고, \(rnAnn_t\) = 4.07% - 4% = 0.07%로, 현재 경제가 균형을 이루려면 명목금리가 0.07% 정도 되어야 함을 의미합니다. 중앙은행 입장에서 자연금리 추정치와 실제 금리를 비교하여 정책 스탠스를 파악하는데, 이 모형에서는 그 부분을 묘사합니다.

이러한 측정식들은 모형의 이론적인 변수를 실제 데이터 (금리, 인플레이션 등)와 직접 연결하기 위한 것으로, 모형의 **추정 혹은 검증 단계**에 활용됩니다. 단순히 단위 변환 역할이므로 경제적 메커니즘에 추가되는 내용은 없습니다.

---

以上各 방정식을 통해 모형의 구조를 모두 살펴보았습니다. 다음 절에서는 이러한 방정식들이 **어떻게 맞물려 작동**하는지를 종합적으로 설명하고, 특히 **자연이자율 충격**과 **제로금리 하한(ZLB) 및 Forward Guidance** 상황에서 모형이 어떤 동태적 반응을 보이는지 상세히 다루겠습니다.

## 모형의 작동 메커니즘 (Overall Mechanism of the Model)

이 모형은 **신케인즈 경기변동 메커니즘**과 **금융가속기 메커니즘**이 결합되어 있습니다. 간략히 요약하면, **총수요와 총공급의 상호작용**으로 산출갭과 인플레이션이 결정되고, 중앙은행이 금리를 조절하여 이를 안정화하려고 하지만, 동시에 **금융부문의 상태(기업가 순자산과 레버리지)**가 투자와 산출에 영향을 미쳐 충격의 효과를 **증폭**시키는 체계입니다. 평상시(작은 충격 하)에는 중앙은행의 테일러룰이 비교적 잘 작동해서 인플레이션과 산출갭을 적절히 조정하지만, **큰 금융충격이나 자연이자율 변화**에는 정책 여력이 제약될 수 있으며, 이때 **금융가속기**가 경기 변동성을 크게 키울 수 있습니다. 메커니즘을 단계별로 정리해 보겠습니다:

- **(1) 가계측: Euler 방정식과 자연이자율 역할)** – 가계는 금리와 기대인플레이션에 반응해 소비를 조절합니다. 만약 경제에 **수요를 감소시키는 충격**이 오면 (예: 선호충격 \(z*t\) 음(-) 또는 정부지출 감소), **자연실질이자율 \(nr_t\)**가 떨어집니다. 그러면 균형을 맞추기 위해서는 실제 실질금리 \(r^n_t - \pi*{t+1}\)도 내려와야 하지만, 중앙은행이 금리를 충분히 낮추지 못하거나 인플레이션이 충분히 올라주지 않으면 (특히 ZLB 근처에서는 금리인하 한계) Euler 방정식에 의해 **현재소비 \(c_t\)**가 크게 위축됩니다. 이는 곧 산출갭 \(y_t\)의 감소로 이어지죠. 반대로 수요가 과열되는 충격(예: 정부지출↑, 선호충격 \(z_t\) 양(+))에는 자연금리가 상승하고, 금리가 그만큼 못 따라올라가면 \(c_t\)가 과열되어 \(y_t\)↑ 압력을 줍니다. **자연이자율**은 이렇게 **수요충격의 크기와 지속성**을 대변하는데, 중앙은행은 이를 쫓아 정책금리를 조정함으로써 \(c_t\)와 \(y_t\) 변동을 완화하려 합니다.

- **(2) 기업 및 가격측: 필립스 곡선과 마크업 \(x_t\) 역할)** – 산출갭 \(y_t\)의 변화는 곧 **노동수요와 마진**에 영향을 줍니다. 경제가 침체되어 \(y_t\)가 떨어지면 기업들은 생산을 줄이고 노동 \(h_t\)를 감축하며, 수요 부진으로 가격을 공격적으로 올리지 못해 **마크업 \(x_t\)**가 높아지는 경향이 있습니다 (한계비용이 내려가거나 수요 부족으로 낮은 가격인상을 통해 이윤을 지키려는 행동). 이때 필립스 곡선에 의해 **인플레이션 \(\pi_t\)**가 하락하게 됩니다 (디플레이션 압력). 반대로 경제가 과열되어 \(y_t\)↑, \(h_t\)↑이면 한계비용이 올라가서 마크업 \(x_t\)는 낮아지고, 기업들은 가격을 올리기 쉬운 환경이 되어 \(\pi_t\)가 상승합니다. 따라서 \(y_t\) 증감 -> \(x_t\) 반대로 -> \(\pi_t\) 증감의 연결이 생깁니다. 중앙은행은 이를 감안하여 금리 \(r^n_t\)를 조정하는데, 예를 들어 인플레이션이 목표보다 높으면 금리를 인상하여 수요를 누르고(실질금리↑ -> \(c_t\)↓), 낮으면 금리인하로 부양하려 합니다. 이 **통화정책 피드백**은 테일러룰에 의해 형성되어, \(\pi_t\)와 \(y_t\)의 지나친 변동을 억제하는 방향으로 작동합니다.

- **(3) 금융부문: 순자산 \(n_t\)와 위험스프레드 \(rk - (r^n-\pi)\) 역할** – 여기서 BGG 금융가속기가 개입합니다. 예를 들어 **경기침체 초기**를 가정해 보겠습니다. \(y*t\)↓에 따라 기업 수익이 줄고 자본가격 \(q_t\)도 하락할 것입니다. 그러면 기업가의 **자본수익률 \(rk_t\)**가 낮아지고, 심지어 디플레이션이 발생하면 실질부채이자 부담은 오히려 늘어나 \(rk_t < (r^n*{t-1}-\pi*t)\)가 될 수 있습니다. 순자산 진화식에 따르면 이런 상황에서 기업가 **순자산 \(n_t\)**는 크게 감소합니다. 한편 경기침체는 대개 금리인하로 대응하지만, **만일 금리 인하가 한계에 부딪히거나(0% 하한)** 혹은 금리를 내려도 디플레이션으로 실질금리가 충분히 낮아지지 않는다면, 실물경기는 더 악화되고 \(n_t\) 감소 폭은 더 커집니다. \(n_t\)가 줄면 레버리지 \(lev_t\)는 상승하여 **금융위험**이 커지고, 위험스프레드 식에 의해 \*\*\(rk*{t+1}\) (미래 투자에 대한 요구수익률)**가 상승합니다. 이는 투자자들이 **투자에 더욱 신중**해짐을 뜻합니다. 구체적으로 말하면, 은행 등 대출자가 볼 때 기업가의 재무건전성이 나빠졌으니 대출금리에 가산금리를 더 붙이거나 대출을 꺼리게 됩니다. 결과적으로 기업가들은 투자비용 상승에 직면하고 실제 **투자 \(i_t\)**를 크게 줄일 수밖에 없습니다. 투자 감소는 다시 **자본수요 축소 -> 자본가격 \(q_t\) 추가 하락 -> 기업가 순자산 악화**로 피드백됩니다. 이 일련의 과정이 **금융가속기: “불황 시 악순환, 호황 시 선순환”**이며, 경기침체를 더욱 심화시키고 경기호황을 더욱 부추기는 힘으로 작용합니다. 모형에서는 작은 \(\nu\) (0.2)와 \(\kappa\) (0.04) 등으로 이 효과가 천천히 진행되지만, 충격이 지속되면 상당한 **진폭 증폭과 지속성 증가**를 가져옵니다. 금융가속기가 없다면 예컨대 총수요 충격의 영향으로 \(y_t\)가 1분기만 -1% 줄고 끝날 것을, 금융가속기 존재 시에는 -1% → -1.5% → -1% → -0.5% 이런 식으로 **더 크게, 더 오래** 지속시키는 양상입니다. 이는 모형 내에서 기업가 순자산 \(n_t\)라는 **상태변수\*\*가 생김으로써 발생하는 중요한 동태적 특징입니다.

- **(4) 정책 및 기대역할: Forward Guidance 효과)** – 중앙은행이 충격에 대응해 금리를 조정하면 Euler 방정식과 필립스 곡선을 통해 향후 \(c*{t+1}, \pi*{t+1}\) 등에 대한 **기대**를 변경시킵니다. 특히 **Forward Guidance(선제적 안내)**란 중앙은행이 미래의 금리경로에 대한 신호를 주어 **민간의 기대를 조작**함으로써 현재 경제에 영향을 주는 정책입니다. 모형의 Euler 방정식을 다시 보면 \(c*t = c*{t+1} - (r^n*t - \pi*{t+1} - nr*t)\)로 정리할 수 있는데, 여기서 \(c*{t+1}\)나 \(\pi*{t+1}\)에 대한 기대가 올라가거나, 미래 금리 경로 기대 \(r^n_t\)가 낮아지면 현재 소비 \(c_t\)를 높이는 힘이 생깁니다. 예를 들어 충격으로 인해 자연실질금리가 크게 떨어져 **ZLB 상황**이 온다면 (즉 필요한 금리는 음수인데 0으로 묶여있는 상황), 그냥 두면 현재 소비와 산출이 크게 위축됩니다. 그러나 중앙은행이 **“경기가 회복돼도 한동안 금리를 0에 묶어두겠다”**고 약속(Forward Guidance)하면, 민간은 미래 어느 시점에 원래라면 금리가 올라야 할 때도 낮은 금리가 유지되리라고 예상하게 됩니다. 이는 **미래의 실질이자율 경로를 낮추는 효과**를 가져오고, Euler 방정식에 의해 **현재 소비를 증가**시키는 요인으로 작용합니다. 또한 중앙은행이 향후 완화적 정책을 시사하면 \*\*미래의 인플레이션 기대 \(\pi*{t+1}\)**도 높아질 수 있습니다 (왜냐하면 중앙은행이 과열을 용인할 거라는 시그널로 이해). 그러면 현재 실질금리 \(r^n*t - \pi*{t+1}\)가 더 떨어지는 효과가 있어 역시 소비와 투자에 긍정적 영향을 줍니다. 요컨대 Forward Guidance는 **기대 경로 조정을 통해 현재의 총수요를 자극**하는 수단인 것입니다. 다만 이러한 효과는 모형에서 가정한 완전정보, 합리적 기대 하에서 상당히 크게 나타날 수 있는데, 문헌에서는 이를 **Forward Guidance Puzzle**이라고도 합니다 (먼 미래의 금리 약속도 현재 소비를 비현실적으로 크게 늘릴 수 있다는 문제). 현실적으로는 민간의 기대가 완전하지 않거나 할인되기 때문에 어느 정도 제약은 있지만, **질적인 효과\*\*는 분명합니다. 후술할 구체 시나리오에서 확인하겠습니다.

- **(5) 거시건전성 정책의 역할)** – 이 모형에는 테일러룰에 \(\zeta_f lev_t\) 항이 포함되어 있어 **중앙은행이 레버리지 동태를 주시하고 개입**하는 구조입니다. 예를 들어 자산가격 상승 등으로 민간 부채가 과도하게 늘어 레버리지가 상승하는 국면에서는 중앙은행이 금리를 추가로 인상하여 거품을 억제하려고 합니다. 이때 금리인상은 단기에 소비와 투자를 둔화시켜 산출갭을 줄이고 물가를 낮출 수 있지만, 금융위험 축적을 막음으로써 향후 **더 큰 버블 붕괴**를 예방하는 의도가 있습니다. 반대로 위기 시에는 레버리지 급락(\(lev_t\)↓)에 대응해 금리를 평소보다 더 낮춰주는 것이 금융완화에 도움이 됩니다. 그러나 금리 하나로 물가안정과 금융안정을 동시에 달성하기는 어려운 트레이드오프도 있습니다. \(\zeta_f = 0.5\)라는 임의의 가중치를 통해 모형에서는 **일정 정도만 금융안정에 기울이는** 것으로 타협하고 있습니다. 예컨대 금융위기가 발생해 인플레이션은 낮고 산출갭도 마이너스지만 레버리지는 너무 높다면, \(\pi_t, y_t\)는 금리인하를 요구하지만 \(lev_t\)는 금리인상을 요구하는 딜레마가 생깁니다. \(\zeta_f\)값에 따라 어느 쪽을 더 고려할지가 결정되고, 이 모형에서는 전자(경제활성화)에 더 방점을 두면서도 약간의 **속도 조절**을 하는 수준으로 볼 수 있습니다.

이와 같이 모형의 다양한 요소들은 **서로 연결**되어 작동합니다. **공급측 충격**(생산성 \(a_t\))이 온다면 필립스곡선을 통해 물가에 직접 영향, 산출에 직접 영향, Euler를 통한 간접 영향 등 다경로로 파급되고, **수요측 충격**(정부지출 \(g_t\)나 선호 \(z_t\))은 자연이자율과 Euler를 통해 1차 영향, 이어 필립스곡선, 금융가속기 등을 거쳐 2차 3차 효과를 냅니다. **금융측 충격**(순자산 \(e^n_t\) 혹은 자산가격 변동)도 마찬가지로 투자경로와 소비경로에 모두 효과를 미쳐 경기변동을 만들어냅니다. 이제 이러한 메커니즘이 두 가지 중요한 시나리오 – **(A) 자연이자율에 큰 음의 충격 발생**, **(B) 그 상황에서 중앙은행이 금리를 장기간 0으로 유지하는 선제지침 시행** – 에서 구체적으로 어떻게 전개되는지 살펴보겠습니다.

### 자연이자율의 급락 충격: 경기침체와 모형의 반응

**시나리오 (A):** 어떤 요인으로 **자연실질이자율 \(nr_t\)**이 갑작스럽게 크게 하락하는 충격을 가정합시다. 이는 모형상으로는 \(e^{tar}\_t\)가 큰 음(-) 또는 선호충격 \(z_t\)가 큰 음(-) 값으로 나타날 수 있습니다. 예를 들어 **심각한 소비절약 성향 증가**(가계가 갑자기 지출을 줄이고 저축을 늘리려는 성향, 공포심리)로 \(z_t = -0.5\) 정도의 충격이 왔다고 합시다. \(\rho_z = 0.5\)이므로 이 충격은 상당 기간 지속되며, 자연이자율 식에서 \(nr_t\)를 하락시킵니다 (대략 \(\Delta nr_t \approx 0.05 \* (-0.5) = -0.025\) 즉 2.5%p 정도 하락 일시적으로 가능). **자연이자율의 급락**은 말하자면 “경제의 중립금리가 마이너스로 떨어지는” 상황입니다. 이는 **Wicksell의 견해**로 보면 강한 **과잉저축/수요부족** 상태로, 중앙은행이 명목금리를 이만큼 내리지 않으면 경제가 불균형을 겪게 된다는 뜻입니다.

모형에서 이 충격이 발생하면 연쇄적인 반응이 일어납니다:

- **초기 충격과 금리 제약:** 자연이자율이 -2%까지 떨어졌다고 가정하죠 (연율 기준으론 -8%p 정도 필요). 그런데 중앙은행은 현실적으로 명목금리를 영(0)에 가까이밖에 내릴 수 없습니다. 모형의 테일러룰에 따르면 인플레이션이 떨어지고 산출갭이 크게 마이너스일 것이므로 금리를 낮추려 할 것입니다. 예컨대 \(\pi*t\) 목표대비 -2%, \(y_t\) -5% 상황이면 테일러룰 원하는 금리는 \(\approx 1.5*(-2) + 0.5*(-5) = -3.5\)% (연율 -14%p에 해당!)일 수도 있습니다. 하지만 **제로 하한(ZLB)**이 존재하면 실제 \(r^n_t\)는 0까지만 내려갈 수 있습니다 (모형에는 명시적 ZLB 제약은 없지만, 시나리오 상정). 그러면 **실질금리** \(r^n_t - \pi*{t+1}\)는 여전히 높게 남습니다. 인플레이션 기대 \(\pi*{t+1}\)도 낮을 테니 실질금리는 상당히 상승한 상태가 됩니다. Euler 방정식으로 보면 이때 \(r^n_t - \pi*{t+1} - nr_t\)는 크게 **양(+)**이어서 **현재소비 \(c_t\)**가 급격히 위축되는 결과를 초래합니다. 가계는 “금리가 충분히 낮지 않아 앞으로의 소비 대비 현재소비가 비싸졌다”고 느껴 소비를 미루고 저축을 늘립니다. 이로 인해 \(c_t\)가 큰 폭 감소하고, 자원제약에 따라 산출 \(y_t\) 역시 급격히 위축됩니다. 한편 금리가 0까지 떨어졌음에도 불구하고 디플레이션이 예상되므로 명목금리가 실질금리로 전환되는 효과(피셔 효과)까지 고려하면 실질금리는 오히려 상승한 셈입니다. 이것을 흔히 **디플레이션 압력 하에서의 실질금리 상승 효과**라고 하며, **유동성 함정**의 특징입니다. 중앙은행이 아무리 완화정책을 펴도 물가하락 기대 때문에 실질 차입비용은 내려가지 않는 함정에 빠진 것이죠.

- **실물경기 급락과 2차 파급:** 산출갭이 크게 음(-)으로 벌어지면 (예: -5% 이상 침체), **실업률 상승, 투자 급감** 등이 발생합니다. 투자 \(i*t\)도 Euler류의 투자방정식(사실 q식, Euler식 결합)으로 보면 실질금리↑, 수요↓로 함께 크게 감소합니다. 자본가격 \(q_t\)가 떨어지고, 미래 수익 전망이 나빠 자본의 기대수익률 \(rk*{t+1}\)도 낮아집니다. 동시에 현재의 \(\pi*t\)는 디플레이션(-)을 보이고, 마크업 \(x_t\)는 수요부진으로 상승해 있습니다. 이런 상황에서 기업가들은 보유자산 가치가 하락하고 사업수익도 악화되어, **순자산 \(n_t\)**에 큰 손실을 입습니다. 또한 디플레이션 때문에 과거에 빌린 부채의 실질부담이 늘어나 \(n_t\)를 추가로 깎아먹습니다 (BGG 모형의 **부채 디플레이션 효과**). \(n_t\) 감소 -> 레버리지 \(lev_t\) 폭증 -> 위험스프레드 확대(\(rk*{t+1}\) 요구치↑) -> **신규 투자 위축**이라는 **금융가속기 악순환**이 강하게 돌기 시작합니다. 예를 들어 순자산이 10% 증발하고 레버리지가 10%p 높아졌다면 \(\nu lev \approx 0.2\*10 = 2\) %p의 추가 이자상승 효과, 이것이 투자위축으로 \(y_t\)를 더 떨어뜨리고, 다음 기에도 순자산이 더 줄어드는 연쇄입니다. **은행 등 금융중개 역할**도 경색되어, 기업이 아무리 낮은 정책금리를 적용받는다 해도 실제 대출금리는 위험프리미엄 때문에 크게 높아진다 볼 수 있습니다. 이런 환경에서는 통화정책 전파경로도 막혀버립니다. 즉, 금리를 0까지 내렸지만 민간에 체감되는 차입금리는 높아 “완화정책이 먹히지 않는” 현상이 발생합니다.

- **심화되는 디플레이션과 부채부담:** 필립스 곡선에 의해 \(x_t\) 높은 상황에서 \(\pi_t\)는 계속 하락압력을 받습니다. 실제 1930년대 대공황이나 1990년대 일본의 사례를 보면, 이러한 충격 하에서 **디플레이션-경기침체-디플레이션 심화**의 악순환이 나타납니다. 물가가 떨어지면 실질금리 하락 대신 오히려 상승, 실질부채 상승, 수요 더 위축... 모형에서도 \(\pi_t\)가 음의 지속으로 가면 Euler식의 부담을 가중시킵니다. 중앙은행은 금리를 더 내릴 수 없으니 유일한 완충은 **기대인플레이션을 높이는 것**인데, 기대인플레이션이 높아지려면 민간이 미래 정책이 상당히 완화적으로 유지되어 경제를 과열시킬 거라 믿어야 합니다. 하지만 모두가 어렵다고 생각하면 오히려 **기대 인플레이션이 더 낮아지는** 악순환도 가능합니다 (자기실현적 디플레이션).

- **정리:** 자연이자율 급락 충격에 통화정책 대응이 제약되면, 모형은 **심각한 경기침체와 디플레이션, 금융위기의 복합국면**을 보입니다. 산출은 추락하고 실업 증가, 인플레이션 급락, 정책금리 0, 그러나 실질금리 높음, 기업가 순자산 대량 파괴, 은행 부실과 신용경색 (모형엔 은행직접 표현 없지만 위험스프레드↑로 나타남). **재정정책** 등 다른 대응이 없는 한 모형 내에서 이 상황은 \(\rho\_{nr}=0.9\) 때문에 상당 기간 지속됩니다. 산출갭은 서서히 개선되겠지만(자연이자율 서서히 정상화), 그동안 잃어버린 생산과 후생 손실은 막대합니다. 이 서술은 2008 금융위기나 2020 팬데믹 초기의 상황과 흡사한 면이 있습니다. BGG 모형은 이러한 **실물-금융 교착 상태**를 제대로 잡아내는 데 유용하며, 여기선 더 나아가 **Forward Guidance 정책**의 효과까지 살펴볼 수 있습니다.

### ZLB 하에서 Forward Guidance: 장기간 0금리 고정의 효과

**시나리오 (B):** 위와 비슷한 자연이자율 급락 상황에서, 중앙은행이 **명목금리를 향후 상당 기간 0%에 묶어두겠다고 공개적으로 약속**했다고 가정합시다. 예컨대 충격 발생 시점에 “앞으로 8분기(2년) 동안 정책금리를 0%로 유지”하겠다고 선제적 지침을 발표하는 것입니다. 이는 통화정책의 **비전형적 수단** 중 하나로, 미래 정책경로에 대한 강력한 신호를 줌으로써 현재 기대를 관리하는 전략입니다. 모형의 agents(가계, 기업)는 합리적 기대를 가지므로, 중앙은행 발표를 믿고 금리가 8분기째 0일 것을 예상합니다. 그렇다면 이 상황이 (A)와 어떻게 달라지는지 단계별로 설명하겠습니다:

- **기대경로 변화:** Forward Guidance가 나오자마자, 민간의 **미래 인플레이션 기대 \(\pi*{t+1}, \pi*{t+2}, ...\)**와 **미래 금리 기대 \(r^n*{t+1}, r^n*{t+2}, ...\)**가 변합니다. 중앙은행이 금리를 올리지 않겠다고 하니, 자연이자율이 추후 회복되더라도 한동안 실질금리가 **자연수준보다 낮게** 유지될 것이라는 믿음이 형성됩니다. 이는 **미래의 수요를 자극**할 전망으로 이어집니다. 특히, 먼 장래보다는 비교적 **가까운 미래 (1~2년 내)**의 정책이 확실히 완화적이라 정해졌기 때문에 Euler 방정식에 큰 영향이 있습니다. 가계는 “앞으로 2년간 돈 빌릴 비용이 매우 낮겠군, 그리고 중앙은행이 일부러 경기를 뜨겁게 할 테니 인플레이션도 오를 거야”라고 기대할 수 있습니다. 그러면 **향후 소비** \(c*{t+1}, c*{t+2},...\)에 대한 기대 경로가 상향조정되고, 심지어 기대인플레이션도 약간씩 목표치에 근접하도록 올라갈 겁니다. Euler 방정식 \(c*t = c*{t+1} - (r^n*t - \pi*{t+1} - nr*t)\)를 보면, 이때 \(r^n_t - \pi*{t+1} - nr*t\)는 여전히 양(+)일지라도 \(c*{t+1}\) 증가 효과가 있습니다. Forward Guidance 없이는 \(c\_{t+1}\)도 침체라 \(c_t\)를 더 줄였겠지만, 이제 “미래엔 좋아질테니 지금 너무 움츠리지 말자”는 기대가 생겨 **현재소비 감소 폭이 줄어듭니다**. **투자**의 경우에도 비슷합니다. 가령 기업은 금리가 0이 오래 가면 투자프로젝트의 할인율이 낮아지니 채산성이 좋아진다고 판단하여, 완전히 투자를 멈추지 않고 일부 유지할 수 있습니다. 또한 기대인플레이션 상승은 현재 실질이자율을 떨어뜨려, 위험스프레드가 주는 압박도 일부 희석될 수 있습니다. (디플레이션 예상 -> 실질부채부담 증가 효과가 경감됨)

- **실제 전개:** Forward Guidance 덕분에 경제추락이 완화되었다고 해도, 초기 충격이 워낙 크다면 여전히 산출갭은 마이너스, 인플레이션은 낮겠지만, **그 깊이가 얕아지게** 됩니다. 예컨대 FG 없이는 \(y_t^{min} = -5\%\)였다면, FG로 -3% 정도에서 바닥치고 반등할 수 있습니다. 인플레이션도 -2%까지 떨어질 것을 -0.8%에서 막는다든지의 차이입니다. 중앙은행은 금리를 어차피 0으로 유지 중이라 단기 대응은 같지만, **민간의 기대가 더 낙관적**이라 경제주체들의 행동이 덜 위축된 것입니다. 이로써 금융부문의 악순환도 덜 심화됩니다. \(n_t\) 감소가 덜하고, \(lev_t\) 상승도 조금 억제되어 \(\nu lev\) 부담이 낮아집니다. 위험스프레드 상승이 완만하면 투자 하락폭도 줄고, \(q_t\)가 너무 내려가지 않아 \(n_t\)를 더 지켜주는 긍정피드백이 일어납니다.

- **정책 의도:** Forward Guidance의 핵심은 **“미래에 약간 인플레이션 과열을 용인하여 현재의 디플레이션을 만회한다”**는 약속입니다. 중앙은행 입장에서는 나중에 자연이자율이 회복되면 사실 금리를 올려야 하지만, FG로 인해 일정 기간 올리지 못하게 스스로 구속한 것이죠. 그 기간에는 경제가 잠시 과열될 수 있습니다. 그러나 이는 지금 심각한 침체를 완화하기 위한 의도적 조치로, **“일종의 후속 보상”**이라 할 수 있습니다. Eggertsson & Woodford가 설명한 바와 같이, ZLB 상황에서는 **미래에 물가가 약간 더 오르도록(인플레 기대 상승)** 함으로써 현재 실질금리를 떨어뜨려 수요를 살리는 정책이 효과적입니다. FG는 바로 그런 메커니즘을 노린 것입니다. 모형에서 보면 FG로 미래 금리를 누르니 미래 \(y, \pi\)가 평소보다 높아지고 (약간 과열), 그 기대가 현재를 떠받쳐 현재 \(y, \pi\)의 최저치를 덜 낮게 만들며, 전체 경기궤적의 **U자 곡선이 덜 깊게** 됩니다.

- **한계와 위험:** 물론 Forward Guidance에는 **신뢰성 문제**가 따라붙습니다. 중앙은행이 나중에 약속을 어기고 일찍 금리를 올려버리면, 민간 기대가 배신당해 경제에 충격을 줄 수 있습니다. 그래서 FG가 효과를 보려면 **강력한 커뮤니케이션과 제도적 장치** (예: “일정 기간 또는 인플레이션 일정 수준 도달까지 금리인상 안한다”와 같은 조건부 약속)가 필요합니다. 모형 자체는 이러한 기대이행 문제를 다루지 않지만, 현실에서는 매우 중요합니다. 또한 FG로 경기를 너무 띄우면 나중에 **인플레이션이 목표치를 초과**하거나 자산버블이 생길 위험도 있습니다. 본 모형에서는 \(\zeta_f\) 같은 항을 통해 금융안정을 조금 고려하긴 하지만, FG 기간에는 오히려 레버리지 높아지는 것을 용인하게 될 수도 있습니다 (저금리가 오래 지속되니 부채가 늘 가능성). 따라서 FG는 **양날의 검**으로, 적정 기간과 강도를 잘 설정해야 합니다.

- **결과 비교:** 결국 Forward Guidance를 실시한 시나리오(B)에서는, **같은 자연이자율 충격에 대해 경기침체의 깊이와 지속기간이 완화**됩니다. 인플레이션도 깊은 디플레이션에 빠지지 않고 비교적 빨리 목표 수준으로 회귀할 수 있습니다. 산출갭도 ZLB 탈출 후(2년 후 금리인상 재개 시점)에는 플러스 영역으로 오르며, FG 없었다면 잃어버렸을 생산을 어느 정도 만회하는 **메이크업 효과**가 나타납니다. 금융가속기의 증폭효과도 줄어 전체 변동성이 낮아집니다. 이처럼 **명확한 정책 약속이 기대를 안정시키고 자금흐름을 유지시켜 위기 상황에서 유용한 대응책이 될 수 있음**을 모형이 보여줍니다. 다만 모형은 민간이 이 신호를 100% 신뢰하고 합리적으로 행동한다는 이상 조건이므로, 실제 정책에서는 FG의 효과가 불확실하고 때로는 예상만큼 크지 않을 수 있다는 점을 언급해둘 필요가 있습니다.

**요약:** 자연이자율에 큰 음의 충격이 발생하면 모형은 **심각한 침체와 금융위기** 국면을 재현하며, Forward Guidance로 장기간 0금리를 약속하면 **민간의 미래에 대한 낙관 기대를 높여 현재 침체를 완화**시킬 수 있음을 보여줍니다. 이는 현실에서 중앙은행들이 ZLB 상황에서 **커뮤니케이션 강화**와 **기대관리**에 힘쓰는 이유를 뒷받침합니다. BGG 금융가속기 모형에서 Forward Guidance의 도입은 흥미로운 상호작용을 시사하는데, 저금리가 오래 유지되면 금융시장에 안도감을 주어 위험스프레드가 덜 오르지만, 한편으론 부채 축적을 부추겨 **미래 취약성**을 키울 수도 있습니다. 이런 트레이드오프까지 종합적으로 고려하는 것은 확장된 모형(예: 거시건전성 정책 규제 병행 등)의 영역일 것입니다.

## 결론 (맺음말)

以上을 통해 **BGG 기반 신케인즈 DSGE 모형**의 핵심 방정식과 경제적 의미를 모두 살펴보았습니다. 각 방정식은 가계의 소비 결정, 기업의 가격 결정, 기업가의 투자 및 차입 결정, 그리고 정책당국의 금리 결정 등 **경제의 다양한 의사결정과 제약**을 수량화하고 있음을 알 수 있습니다. 특히 이 모형은 **금융부문의 역할**을 강조하여, 순자산과 레버리지 변동이 어떻게 **실물경제에 파급**되고 **충격을 증폭**시키는지 명확히 보여줍니다. 마지막으로 다룬 시나리오들은 이러한 메커니즘이 현실적인 상황에서 어떤 모습을 보이는지 가늠케 해주었습니다. **자연이자율 충격**은 경기침체의 근원을 이해하게 하고, **명목금리 제약과 Forward Guidance**는 비전통적 정책의 중요성을 일깨워줍니다. 이 모형은 학부 수준에서 이해하기엔 방정식이 많고 복잡해 보일 수 있지만, 각각의 의미를 직관적으로 풀어보면 우리가 거시경제에서 접하는 원리들과 크게 다르지 않음을 알 수 있습니다. 결국 **총수요-총공급, 통화정책, 그리고 금융시장** – 이 세 축의 상호작용이 경제를 움직이며, 모델은 이를 체계적으로 연결짓고 있을 뿐입니다. 이러한 DSGE 모형 분석을 통해 학생들은 **이론과 현실의 연결고리**를 보다 입체적으로 이해할 수 있고, 정책 충격이나 외부 충격에 대한 경제의 반응을 시뮬레이션해봄으로써 거시경제정책 수립에 필요한 통찰을 얻게 될 것입니다.

참고 문헌: Bernanke, Gertler, Gilchrist (1999), Gertler & Kiyotaki (2010), Christiano et al. (2014) 등 금융마찰을 포함한 DSGE 연구, 및 Eggertsson & Woodford (2003) 등의 ZLB 상황 정책 연구.

**End of document.**
