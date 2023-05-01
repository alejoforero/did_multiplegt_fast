# did_multiplegt_fast
Parallel version of did_multiplegt. Parallelizes the bootstrap calculation, everything else is the same and credit goes to the creators of the original command.

Using Stata MP4 and a 16 core CPU I achieved estimation time reduction of up to 90%.

Set maxprocessors(#) to a number close to the cores of your computer such that it minimizes mod(breps,#). For example if you set breps to 50 and have a 16 core processor, the most efficient parameter is either 13 or 17. Your mileage may vary.
Also note that this will open # Stata sessions each with a copy of your data. Make sure you have enough RAM memory for that, and if not, minimize the estimation dataset with only the needed variables.

This command adaptation accelerates the boostrap estimations, and will work better with a higher number of breps. With a low number of breps, (e.g. <5) it will be ->slower<- than the original command, due to the overhead costs of setting up the parallel Stata instances.


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

`net install did_multiplegt_fast, from("https://raw.github.com/alejoforero/did_multiplegt_fast/master/")`

# Validation
You can compare the results of the fast versus the original command with the dofile: `did_multiplegt_fast_validation.do`. Once the seed is set, the results are identical.



# References: 

Clément de Chaisemartin & Xavier D'Haultfoeuille & Yannick Guyonvarch, 2019. "DID_MULTIPLEGT: Stata module to estimate sharp Difference-in-Difference designs with multiple groups and periods," Statistical Software Components S458643, Boston College Department of Economics, revised 08 Feb 2023
https://ideas.repec.org/c/boc/bocode/s458643.html
