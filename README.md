This repository contains MATLAB implementations of **constant-step** and **adaptive** Generalized-Alpha (Chung–Hulbert) time-integration schemes for a 2-DOF structural system. It includes error estimation (Zienkiewicz–Xie), adaptive time-stepping, and visualization of displacement, velocity, acceleration, and error evolution.

---

## 🔧 Features

### ✔ 1. Generalized-Alpha Time Integrator  
Implements the Chung–Hulbert method with parameters:
- ρ∞ (spectral radius at infinity)  
- α₁, α₂ (Rayleigh damping)  
- αₘ, α_f, β, γ computed internally  

### ✔ 2. Constant Time-Step Solver  
`Const_Chung_Hulbert_Solver`  
- Marches solution forward with fixed Δt  
- Precomputes inverse effective stiffness for efficiency  
- Computes Zienkiewicz–Xie error at each step  

### ✔ 3. Adaptive Chung–Hulbert Solver  
`Adaptive_Chung_Hulbert_Solver`  
- Uses ZX relative error to adapt Δt  
- Time step update:  
  $\[
  \Delta t_{\text{new}} = \Delta t_{\text{old}} \sqrt{\frac{\eta}{\eta_{\text{rel}}}}
  \]$
- Automatically expands storage arrays during runtime  
- Achieves accuracy with fewer total steps  

### ✔ 4. Zienkiewicz–Xie Error Estimator  
`ZX_Error`  
Computes:
- Error vector  
- Error norm  
- Relative error norm  
- Accumulated error  

### ✔ 5. Visualization Tools  
`CH_Plots` generates:
- Displacement, velocity, acceleration histories  
- Time-step evolution  
- ZX relative error evolution (with bounds in adaptive mode)  
- Cumulative error norm  

---


---

## ▶ Usage

### Run the full simulation:
```matlab
run main.m
Included examples:
Case 1: Constant Δt, ρ∞ = 1

Case 2: Constant Δt, ρ∞ = 0.1

Case 3: Adaptive time-stepping, ρ∞ = 1

Each case outputs:

Final displacements/rotations

Cumulative error

Number of time steps

History plots and error evolution figures

📊 Outputs
The script automatically prints:
vbnet
Copy code
Value of theta and u after 5 seconds...
Constant-step: cumulative error = ...
Adaptive-step: cumulative error = ...
