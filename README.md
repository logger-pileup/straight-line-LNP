# Proof sizes for proving MLWE sample

## Overview

These scripts can be used used to predict the proof sizes of a NIZK proof for proving an MLWE sample in some lattice-based ZKP schemes. Concretely, we want to prove knowledge of a vector s such that the ||(s,Ds-u)|| is at most B.

The scripts accept configurable parameters (ring dimension, modulus, etc.) and have methods for setting the other parameters (based on the security of lattice assumptions), verifying the required conditions on the parameters hold, and calculating the proof size. 

We use the lattice estimator 
> Martin R. Albrecht, Rachel Player and Sam Scott. On the concrete hardness of Learning with Errors. Journal of Mathematical Cryptology. Volume 9, Issue 3, Pages 169–203, ISSN (Online) 1862-2984, ISSN (Print) 1862-2976 DOI: 10.1515/jmc-2015-0016, October 2015.

to predict the security of LWE and SIS. Since this can take some time to run, we provide an option to set the MLWE dimension manually and continue with the rest of the script.

## Run the scripts

```
git clone https://github.com/malb/lattice-estimator.git 
cd script
sage sl_lnp.sage # etc.
```

Note: when running a script with `sage`, it executes code in the `if __name__ === "__main__"` of any sripts imported with `load`, which may lead to undesirable behaviour.

## Schemes
- `lnp.sage` - original LNP.
- `katsumata.sage` - Katsumata LinHC scheme.
- `lnp_katsumata.sage` - LNP with Katsumata transform.
- `sl_lnp.sage` - straight-line LNP.

## Results

| Ring dimension | Modulus | LNP / KB | LNP + Katsumata / KB | Straight-line LNP / KB |
| - | - | - | - | - |
| 128 | 32 | 14.96 | -     | - |
| 128 | 64 | 19.58 | 120.7 | - |
| 128 | 86 | 23.78 | 172.6 | 87.17 |
| 64  | 64 | 17.95 | 114.1 | - |
| 64  | 90 | 22.07 | 173.4 | 86.67 |
