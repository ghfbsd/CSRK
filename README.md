# CSRK
Implementation of the CSRK (Carnahan-Starling Redlich-Kwong) fluid equation of state (EOS) in [R](https://www.r-project.org/about.html).

This is a function (CSRK) that will return a function that will give volumes
and log(fugacity coefficient) for any pure fluid species at a given pressure
(bars) and temperature (Kelvin).

The code defines the function and runs a simple test showing how to use it
to get volumes and fugacity coefficients for water.  From R, you can try it
directly by typing,
```

source('https://raw.githubusercontent.com/ghfbsd/CSRK/main/CSRK.R')

```
(provided your R release has URL accessing built into it; most do).

Citation:

Helffrich, G. and Connolly, J. A. D. (2024).  A fluid equation of state for use in planetary interiors.  _American Mineralogist_ (submitted).
