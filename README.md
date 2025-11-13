# DMA-Multicarrier-WPT

This repository accompanies the open-access IEEE Transactions on Wireless Communications article "Waveform Optimization and Beam Focusing for Near-Field Wireless Power Transfer With Dynamic Metasurface Antennas and Non-Linear Energy Harvesters".

## Repository layout

| File | Description |
| --- | --- |
| `DMA_opt_run.m` | Main experiment script. Sweeps user counts, DMA apertures, and sub-carrier sets while alternately optimizing the transmit waveform and DMA tuning vector to satisfy rectenna DC targets with minimum RF/baseband power consumption. Saves intermediate `.mat` files per configuration. 【F:DMA_opt_run.m†L1-L169】 |
| `DMA_deploy.m` | Places the DMA aperture in space, computes the per-element attenuation matrix, and returns element coordinates based on the desired operating frequency and panel dimensions. 【F:DMA_deploy.m†L1-L35】 |
| `Channel_comp.m` | Builds the frequency-selective channel tensor between each DMA element and user by combining free-space path loss and a boresight gain model. 【F:Channel_comp.m†L1-L24】 |
| `DMA_scenario.m` | Given current waveform and DMA settings, synthesizes the multi-sine transmit samples, evaluates rectified DC voltages, and accounts for high-power amplifier efficiency and consumption. 【F:DMA_scenario.m†L1-L53】 |
| `fixed_W_SCA.m`, `fixed_Q_SCA.m` | Successive convex approximation (SCA) sub-problems that, respectively, refine the DMA Lorentzian weights and the digital beamforming waveform while enforcing hardware and rectifier constraints. 【F:fixed_W_SCA.m†L1-L31】【F:fixed_Q_SCA.m†L1-L45】 |
| `phase_init_func.m` | Heuristic initializer that allocates RF chains to users and seeds the DMA phase profile and waveform amplitudes before SCA iterations. 【F:phase_init_func.m†L1-L116】 |
| `H_DMA.m`, `Q_DMA.m`, `Rect_K.m` | Helper routines for the DMA propagation loss, tunable element response, and rectifier Taylor-series coefficients. 【F:H_DMA.m†L1-L29】【F:Q_DMA.m†L1-L4】【F:Rect_K.m†L1-L4】 |
| `taylor_coef_q.m`, `taylor_coef_w.m` | Generate first-order Taylor expansions of the rectifier output with respect to the DMA weights (`q`) and beamforming waveform (`w`), used by the SCA solvers. 【F:taylor_coef_q.m†L1-L37】【F:taylor_coef_w.m†L1-L36】 |

## Requirements

- MATLAB R2023a or later (tested with the Communications Toolbox for `physconst`).
- [CVX](http://cvxr.com/cvx/) with the MOSEK solver configured, required by the SCA sub-problems in `fixed_W_SCA.m` and `fixed_Q_SCA.m`. 【F:fixed_W_SCA.m†L12-L26】【F:fixed_Q_SCA.m†L13-L45】
- The user location datasets referenced in `DMA_opt_run.m` (`user_location.mat`) placed in the repository root; the scripts expect `user_loc1`, `user_loc2`, and `user_loc4` tensors inside the MAT-file. 【F:DMA_opt_run.m†L33-L47】

## Reproducing the paper's experiments

1. Install MATLAB dependencies and ensure CVX/MOSEK are on the MATLAB path.
2. Place `user_location.mat` (provided with the publication's supplementary material) in the project root.
3. Open `DMA_opt_run.m` and adjust the search grids if desired:
   - `L_set` – DMA aperture dimensions (meters). 【F:DMA_opt_run.m†L12-L13】
   - `M_set` – Number of simultaneous WPT users. 【F:DMA_opt_run.m†L11-L12】【F:DMA_opt_run.m†L56-L69】
   - `N_set` – Sub-carrier counts for the multi-sine waveform. 【F:DMA_opt_run.m†L11-L13】【F:DMA_opt_run.m†L69-L75】
   - `DC_Threshold_main` – Rectifier power target per user. 【F:DMA_opt_run.m†L24-L25】
4. Run the script from MATLAB: `>> DMA_opt_run`.
5. For each configuration and user index, the script reports progress and stores a `.mat` snapshot (for example, `DMA_OverN_User_1_Length_1_User_num_2_sub_carrier_4.mat`) with the optimized variables and power accounting. 【F:DMA_opt_run.m†L155-L166】

The alternating procedure iterates between:

- Optimizing the DMA Lorentzian weight vector (`Q_dma`) with `fixed_W_SCA` while holding the digital waveform fixed.
- Optimizing the digital waveform (`w_beam`) with `phase_init_func` and `fixed_Q_SCA` while keeping the DMA profile constant.
- Re-evaluating the time-domain scenario with `DMA_scenario` to update the rectifier outputs and power consumption.

These steps mirror the methodology in the IEEE TWC article, where the DMA is leveraged to sculpt both spatial and frequency responses that maximize end-user DC power under transmitter efficiency constraints.

## Tips for extending the code

- **Alternative channel models:** Replace `Channel_comp.m` with a more detailed propagation model (e.g., clustered multipath) while preserving the output tensor dimensions `(users × subcarriers × RFC × elements)`. 【F:Channel_comp.m†L1-L24】
- **Hardware constraints:** The Lorentzian response ball constraint in `fixed_W_SCA.m` encodes the feasible DMA tuning manifold. Modify the quadratic inequality if a different metasurface tuning range is required. 【F:fixed_W_SCA.m†L22-L24】
- **Power accounting:** `DMA_scenario.m` centralizes rectifier modeling and high-power amplifier efficiency. Adjust rectifier parameters (`i_s`, `eta_0`, `V_o`, `R_L`) or the amplifier efficiency curve (`eff_max`, `As`) to match alternative hardware. 【F:DMA_scenario.m†L4-L53】

## Reference

If you use this code, please cite:

> O. Al-Khateeb, A. M. Elbir, N. Shlezinger, H. V. Poor, and Y. C. Eldar, "Dynamic Metasurface Antenna Enabled Multi-Carrier Wireless Power Transfer," *IEEE Transactions on Wireless Communications*, Early Access, 2024.

---

For questions or clarifications, please refer to the corresponding author details in the journal article.
