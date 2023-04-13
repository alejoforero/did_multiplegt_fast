# did_multiplegt_fast
Parallel version of did_multiplegt. Parallelizes the bootstrap calculation, everything else is the same and credit goes the creators of original command.


# Requirements
`ssc install matsave`

Install `ftools` and from github, not SSC

```
cap ado uninstall ftools

net install ftools, from("https://raw.githubusercontent.com/sergiocorreia/ftools/master/src/")
```
Check his website for more details on installation.
http://scorreia.com/software/reghdfe/install.html

Install `parallel` from github; don't install from SSC:
```
cap ado uninstall parallel

net install parallel, from(https://raw.github.com/gvegayon/parallel/stable/) replace

mata mata mlib index
```

# Installation

`net install did_multiplegt_fast, from("https://raw.githubusercontent.com/alejoforero89/did_multiplegt_fast/master/"`

# Validation
You can compare the results of the fast versus the original command with the dofile: 'did_multiplegt_fast_validation.do'. Once the seed is set, the results are identical.



# References: 

Clément de Chaisemartin & Xavier D'Haultfoeuille & Yannick Guyonvarch, 2019. "DID_MULTIPLEGT: Stata module to estimate sharp Difference-in-Difference designs with multiple groups and periods," Statistical Software Components S458643, Boston College Department of Economics, revised 08 Feb 2023
https://ideas.repec.org/c/boc/bocode/s458643.html
