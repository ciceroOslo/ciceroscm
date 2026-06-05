# AMMONIA Updates From WP2-4 Model Means

This note extracts the quantities needed from the June 3, 2026 AMMONIA summary using the `Anthropogenic` rows.

The requested derived quantities are:

- O3 burden change per unit emission change
- Methane lifetime change per unit emission change
- Combined aerosol RF change, defined as `NO3 RF + SO4 RF`, per unit emission change

For `NH3`, the normalization is per `Tg NH3`.

For `NOx`, the AMMONIA table is first interpreted as `Tg NOx`, then converted to `Tg N` because the gas parameter file uses `NOx` in `Mt_N`.

## Summary Table

| Species | Base perturbation | Normalization used here | O3 burden change | CH4 lifetime change | Combined aerosol RF |
| --- | --- | --- | --- | --- | --- |
| NH3 | `10.00 Tg NH3` | per `Tg NH3` | `-0.00282 DU Tg(NH3)-1` | `0.01544 % Tg(NH3)-1` | `-1.01722 mW m-2 Tg(NH3)-1` |
| NOx | `10.00 Tg NOx` | per `Tg N` | `0.13143 DU Tg(N)-1` | `-0.46652 % Tg(N)-1` | `-0.81486 mW m-2 Tg(N)-1` |

## NH3 Calculations

AMMONIA anthropogenic NH3 values:

- Emission perturbation: `10.00 Tg NH3`
- O3 burden change: `-0.0282 DU`
- CH4 relative lifetime change: `0.1544 %`
- NO3 RF: `-10.0698 mW m-2`
- SO4 RF: `-0.1024 mW m-2`

O3 burden change per Tg NH3:

```text
-0.0282 / 10.00 = -0.00282 DU Tg(NH3)-1
```

Methane lifetime change per Tg NH3:

```text
0.1544 / 10.00 = 0.01544 % Tg(NH3)-1
```

Combined aerosol RF per Tg NH3:

```text
(-10.0698 - 0.1024) / 10.00 = -1.01722 mW m-2 Tg(NH3)-1
```

In W m-2 per Tg NH3:

```text
-1.01722 mW m-2 Tg(NH3)-1 = -0.00101722 W m-2 Tg(NH3)-1
```

Because the current gas parameter file stores NH3 emissions in `Mt`, and `1 Tg = 1 Mt`, the numeric values are unchanged when interpreted per `Mt NH3`.

## NOx Calculations

AMMONIA anthropogenic NOx values:

- Emission perturbation: `10.00 Tg NOx`
- O3 burden change: `0.40 DU`
- CH4 relative lifetime change: `-1.4204 %`
- NO3 RF: `-2.39 mW m-2`
- SO4 RF: `-0.09 mW m-2`

The gas parameter file uses `NOx` in `Mt_N`, so the emission perturbation must be converted from `Tg NOx` to `Tg N`.

Assuming NOx is reported as NO2-equivalent mass:

```text
1 Tg NOx = 14 / 46 Tg N
10.00 Tg NOx = 10.00 x 14 / 46 = 3.04348 Tg N
```

O3 burden change per Tg N:

```text
0.40 / 3.04348 = 0.13143 DU Tg(N)-1
```

Methane lifetime change per Tg N:

```text
-1.4204 / 3.04348 = -0.46652 % Tg(N)-1
```

Combined aerosol RF per Tg N:

```text
(-2.39 - 0.09) / 3.04348 = -0.81486 mW m-2 Tg(N)-1
```

In W m-2 per Tg N:

```text
-0.81486 mW m-2 Tg(N)-1 = -0.00081486 W m-2 Tg(N)-1
```

## N2O Reference Table

The N2O block can stay as reference material:

| Source | Direct ERF vs. SARF (%) | Indirect O3+aerosol (mW m-2 ppb(N2O)-1) | Indirect CH4 (mW m-2 ppb(N2O)-1) |
| --- | --- | --- | --- |
| IPCC AR6 | `7 +/- 13` | `0.55 +/- 0.04` | `-0.94` |
| AMMONIA | `-7.7` | `0.40` | `-0.39` |
| CESM2 | `-10` | `0.99` | `-0.21` |
| GISS | `-9` | `-0.54` | `-0.82` |
| GFDL | `-4` | `0.75` | `-0.15` |