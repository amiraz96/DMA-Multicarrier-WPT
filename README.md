# DMA-Multicarrier-WPT

This repository accompanies the open-access IEEE Transactions on Wireless Communications article "Waveform Optimization and Beam Focusing for Near-Field Wireless Power Transfer With Dynamic Metasurface Antennas and Non-Linear Energy Harvesters".

## Repository layout

| File | Description |
| --- | --- |
| `DMA_opt_run.m` | Main experiment script. Sweeps user counts, DMA apertures, and sub-carrier sets while alternately optimizing the transmit waveform and DMA tuning vector to satisfy rectenna DC targets with minimum RF/baseband power consumption. Saves intermediate `.mat` files per configuration. 【F:DMA_opt_run.m†L1-L169】 |
| `DMA_deploy.m` | Places the DMA aperture in space, computes the per-element attenuation matrix, and returns element coordinates based on the desired operating frequency and panel dimensions. 
| `Channel_comp.m` | Builds the frequency-selective channel tensor between each DMA element and user by combining free-space path loss and a boresight gain model. 
| `DMA_scenario.m` | Given current waveform and DMA settings, synthesizes the multi-sine transmit samples, evaluates rectified DC voltages, and accounts for high-power amplifier efficiency and consumption. 
| `fixed_W_SCA.m`, `fixed_Q_SCA.m` | Successive convex approximation (SCA) sub-problems that, respectively, refine the DMA Lorentzian weights and the digital beamforming waveform while enforcing hardware and rectifier constraints. 
| `phase_init_func.m` | Heuristic initializer that allocates RF chains to users and seeds the DMA phase profile and waveform amplitudes before SCA iterations. 
| `H_DMA.m`, `Q_DMA.m`, `Rect_K.m` | Helper routines for the DMA propagation loss, tunable element response, and rectifier Taylor-series coefficients.
| `taylor_coef_q.m`, `taylor_coef_w.m` | Generate first-order Taylor expansions of the rectifier output with respect to the DMA weights (`q`) and beamforming waveform (`w`), used by the SCA solvers.

## Requirements

- MATLAB R2021a or later.
- [CVX](http://cvxr.com/cvx/) with the MOSEK solver configured, required by the SCA sub-problems in `fixed_W_SCA.m` and `fixed_Q_SCA.m`. 
- The user location datasets referenced in `DMA_opt_run.m` (`user_location.mat`); the scripts expect `user_loc` tensors inside the MAT-file. 

## Reproducing the paper's experiments

1. Install MATLAB dependencies and ensure CVX/MOSEK are on the MATLAB path.
2. Place `user_location.mat` (provided with the publication's supplementary material) in the project root.
3. Open `DMA_opt_run.m` and adjust the search grids if desired:
   - `L_set` – DMA aperture dimensions (meters). 
   - `M_set` – Number of simultaneous WPT users.
   - `N_set` – Sub-carrier counts for the multi-sine waveform.
   - `DC_Threshold_main` – Rectifier power target per user. 
4. Run the script from MATLAB: `>> DMA_opt_run`.
5. For each configuration and user index, the script reports progress and stores a `.mat` snapshot (for example, `DMA_OverN_User_1_Length_1_User_num_2_sub_carrier_4.mat`) with the optimized variables and power accounting. 【F:DMA_opt_run.m†L155-L166】

The alternating procedure iterates between:

- Optimizing the DMA Lorentzian weight vector (`Q_dma`) with `fixed_W_SCA` while holding the digital waveform fixed.
- Optimizing the digital waveform (`w_beam`) with `phase_init_func` and `fixed_Q_SCA` while keeping the DMA profile constant.
- Re-evaluating the time-domain scenario with `DMA_scenario` to update the rectifier outputs and power consumption.

These steps mirror the methodology in the IEEE TWC article, where the DMA is leveraged to sculpt both spatial and frequency responses that maximize end-user DC power under transmitter efficiency constraints.

## Tips for extending the code

- **Alternative channel models:** Replace `Channel_comp.m` with a more detailed propagation model (e.g., clustered multipath) while preserving the output tensor dimensions `(users × subcarriers × RFC × elements)`. 
- **Hardware constraints:** The Lorentzian response ball constraint in `fixed_W_SCA.m` encodes the feasible DMA tuning manifold. Modify the quadratic inequality if a different metasurface tuning range is required. 
- **Power accounting:** `DMA_scenario.m` centralizes rectifier modeling and high-power amplifier efficiency. Adjust rectifier parameters (`i_s`, `eta_0`, `V_o`, `R_L`) or the amplifier efficiency curve (`eff_max`, `As`) to match alternative hardware. 

## Citing

If you use this code, please cite the associated article. A BibTeX entry is provided below.

```bibtex
@ARTICLE{DMA_azarbahram,
  author={Azarbahram, Amirhossein and López, Onel L. A. and Latva-Aho, Matti},
  journal={IEEE Transactions on Wireless Communications}, 
  title={Waveform Optimization and Beam Focusing for Near-Field Wireless Power Transfer With Dynamic Metasurface Antennas and Non-Linear Energy Harvesters}, 
  year={2025},
  volume={24},
  number={2},
  pages={1031-1045},
  keywords={Array signal processing;Radio frequency;Wireless communication;Optimization;Power demand;Complexity theory;Rectifiers;Radio transmitters;Costs;Transmitting antennas;Radio frequency wireless power transfer;waveform design;beamforming;dynamic metasurface antennas;near-field channels},
  doi={10.1109/TWC.2024.3503908}}
```

## License

The repository inherits the usage rights granted by the original authors. Refer to the paper or contact the authors for explicit licensing terms.
