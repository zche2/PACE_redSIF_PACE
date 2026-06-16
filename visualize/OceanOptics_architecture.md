# Ocean Optics Module — Architecture

```
┌─────────────────────────────────────────────────────────┐
│                  GEOPHYSICAL INPUTS                     │
│                                                         │
│   λ (wavelength arrays)                                 │
│   Chl-a · CDOM · SPM  (concentration profiles)         │
│   dz  (layer discretisation)                            │
└───────────────────────┬─────────────────────────────────┘
                        │
                        ▼
┌─────────────────────────────────────────────────────────┐
│          LAYER A — Type Hierarchy & Dispatch            │
│                                                         │
│   AbstractWaterConstituent                              │
│   ├── AbstractAbsorptionModel                           │
│   │       └── PureWater · CDOM · Phyto types           │
│   └── AbstractScatteringModel                           │
│           └── Phase functions · VSF types              │
└───────────────────────┬─────────────────────────────────┘
                        │
                        ▼
┌─────────────────────────────────────────────────────────┐
│          LAYER B — Bio-Optical Numerical Engine         │
│                                                         │
│   Absorption     a_tot = a_w + a_cdom + a_ph            │
│   Scattering     b_tot  and  b_bp (backscattering)      │
│   Angular VSF    Fournier–Forand · Petzold              │
└───────────────────────┬─────────────────────────────────┘
                        │
                        ▼
┌─────────────────────────────────────────────────────────┐
│          LAYER C — Radiative Transfer Bridge            │
│                                                         │
│   Optical thickness    τ  = (a_bulk + b_bulk) × dz      │
│   Single-scat. albedo  ω₀ = b_bulk / (a_bulk + b_bulk)  │
│   Dimensionless scaling & matrix packing                │
└───────────────────────┬─────────────────────────────────┘
                        │
                        ▼
┌─────────────────────────────────────────────────────────┐
│                   OceanState struct                     │
│          (single hand-off to the RTM layer)             │
└──────────┬────────────────────────────┬─────────────────┘
           │                            │
           ▼                            ▼
┌──────────────────────┐   ┌────────────────────────────┐
│  vSmartMOM.jl        │   │  Fresnel Boundary          │
│  Core RTM Solver     │──▶│  Condition Matcher         │
│                      │   │  (air–sea interface)       │
└──────────────────────┘   └────────────────────────────┘
```
