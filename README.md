
# Fast Multivariate Multi-point Evaluation

This is an implementation of the paper "Fast multivariate multi-point evaluation revisited" by Joris van der Hoeven and Grégoire Lecerf.
The purpose of this project is to provide an open-source implementation of the Kedlaya-Umans algorithm for multivariate multi-point evaluation over finite fields.

## Milestones

### Bivariate case

- [x] Whole field evaluation (base case when $t=0$) in $\tilde O(p d_Y + d_Y d_X + p^2 + n)$ bivariate case
- [x] Naive evaluation in $\tilde O(d_Y d_X n)$ bivariate case
- [x] Main recursion / multimodular bivariate case
- [x] Timing/bounds heuristics
- [ ] CRT Tree instead of iterative
- [x] **KU in the bivariate case**

### Multivariate case

- [x] Dense arithmetic in $\mathbb F_p[X_1, \dots, X_k]$
- [x] Whole field evaluation
- [x] Naive evaluation
- [x] Main recursion / multimodular
- [ ] CRT Tree instead of iterative
- [x] **KU in the multivariate case**

### Reduction to modular composition

The modular composition as described in the paper is implemented in gives the following timing on my laptop :

```text
Degree      Brent-Kung NTL    Naive Precomp     Naive             KU                Equal (1 vs 2-4) 
-----------------------------------------------------------------------------------------------------
100         0.000394          0.001805          0.001795          0.077471          1 | 1 | 1        
200         0.001074          0.010998          0.013413          0.364637          1 | 1 | 1        
400         0.003127          0.054707          0.097664          1.919801          1 | 1 | 1        
800         0.009534          0.211145          0.649146          9.475037          1 | 1 | 1        
1600        0.028280          0.925718          2.754807          47.084018         1 | 1 | 1        
3200        0.084810          3.881357          10.914877         255.298157        1 | 1 | 1        
6400        0.252576          16.453398         46.066293         1212.349421       1 | 1 | 1   
```

## Installation

First clone the repo with the submodule :

```bash
git clone --recurse-submodules https://github.com/maelhos/FMME.git
```

**Note that if you downloaded this repo as a zip or didn't use `--recurse-submodules` you must manually run `git clone https://github.com/libntl/ntl.git` in the main directory as well as `git clone https://github.com/vneiger/pml.git`.**

Then, patch and build the NTL/PML :

```bash
chmod +x patchNTL.sh
./patchNTL.sh
cd ntl/src
./configure
make clean && make
cd ../..
chmod +x patch_pml_static.sh
./patch_pml_static.sh
cd pml/ntl-extras
make clean && make
cd ../..
```

*Patching the NTL is required for all instances above a few rounds. If you get a "Polynomial too big for FFT" error, check that the following lines are correctly set in `ntl/include/FFT.h` :*

```c++
#if (36 <= NTL_FFTMaxRootBnd)
#define NTL_FFTMaxRoot (36)
```

Note that we use the NTL with static linking.
After building the NTL you should have it as a static library in `ntl/src/ntl.a`

Then you can just run

```bash
make
```

And you should have the current test file named `ku`.
