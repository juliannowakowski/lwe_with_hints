# LWE with Hints
An efficient Python library for attacking LWE with hints, accompanying the paper

**Too Many Hints – When LLL Breaks LWE.**


## Requirements

* [NumPy](https://numpy.org/)
* [fpylll](https://github.com/fplll/fpylll)
* [SciPy](https://scipy.org/)

We note that all required libraries ship with [Sage](https://www.sagemath.org/). (However, fpylll can be very slow in Sage, when not upgrading to the latest fpylll version, see [https://github.com/fplll/fpylll/pull/239](https://github.com/fplll/fpylll/pull/239).)

If you want to run the library in Sage, then you have run Sage in its Python mode.
For instance, to run `tutorial.py` in Sage, simply execute the following command:
```console
sage --python tutorial.py
```
The code is written for Python 3.x / Sage 9.x.

## Tutorial

### Getting started

We import the library and load an LWE instance from `tutorial/lwe_instance.json`.
```py
>>> from lwe_with_hints import *
>>> A,b,q = loadLWEInstanceFromFile("tutorial/lwe_instance.json")
```

As the following code snippet shows, our LWE instance has parameters n = 70, m = 80, q = 521.
```py
>>> A.shape
(70,80)
>>> q
521
```


Let us first attack the LWE instance *without* hints. To this end, we create an `LWELattice` object and run the (progressive) BKZ lattice reduction algorithm on it.
```py
>>> lattice = LWELattice(A,b,q)
>>> lattice.reduce()
```
Executing the above code should take roughly 1 - 2 minutes on a laptop.

When execution is finished, the LWE secret is stored in `lattice.s`.
```py
>>> lattice.s
array([ 0,  0,  2,  1,  2, -1, -1,  0, -1,  1,  1,  1,  0,  1,  1, -3, -2,
       -2, -3, -1,  1,  1, -1,  2, -1,  0, -2,  0,  1,  0,  0,  0,  1,  0,
        1,  0,  0,  0,  0,  0,  1,  0,  0,  0, -1,  0, -1, -1, -1,  1,  1,
        1,  0,  0,  0,  0,  2,  1,  0, -1,  1, -2,  0,  0, -2,  0,  1,  0,
       -1,  0])
```

The shortest vector, that was recovered by lattice reduction, is stored in `lattice.shortestVector`.
```py
>>> lattice.shortestVector
array([ 2, -1,  1,  1,  0, -2,  2,  2,  0, -1,  1, -1,  0,  0,  0,  0,  1,
       -2,  0, -2,  2, -3,  0, -1,  0,  0, -1, -1, -1, -1, -1,  1, -1, -2,
        2,  1,  1,  1,  0,  0, -1,  0, -1,  2,  2,  1,  0, -1,  0,  0, -1,
       -1, -1, -1, -1,  1,  1,  1,  1, -1,  1,  1,  1,  1,  0,  3,  0,  1,
        1,  0,  0, -2, -1,  0,  0,  1,  2, -2, -1,  1,  0,  0, -2, -1, -2,
        1,  1,  0,  1, -1, -1, -1,  0, -1, -1,  3,  2,  2,  3,  1, -1, -1,
        1, -2,  1,  0,  2,  0, -1,  0,  0,  0, -1,  0, -1,  0,  0,  0,  0,
        0, -1,  0,  0,  0,  1,  0,  1,  1,  1, -1, -1, -1,  0,  0,  0,  0,
       -2, -1,  0,  1, -1,  2,  0,  0,  2,  0, -1,  0,  1,  0,  1])
```

The BKZ blocksize, at which the LWE secret was successfully recovered, is stored in `lattice.successBlocksize`.
```py
>>> lattice.successBlocksize
5
```

### Integrating hints
Now let us attack the LWE instance again, but this time *with* hints.

Let us first reset the `LWELattice` object.
```py
>>> lattice = LWELattice(A,b,q)
```

Suppose we obtain the following hints via some side channel.

1. The first four coordinates of the LWE secret  are 0, 0, 2 and 1.
2. The inner product between the LWE secret and the following two vectors equals -1670 and 2381, respectively:
```py
(-459, -441, 107, -207, 30, 358, -221, -483, 457, 96, 118, -241, 400, -478, 374, -46, -376, 415, 213, 476, -195, 25, -486, 444, 228, 313, -252, -182, -314, 105, -248, 163, 489, -388, 222, 110, -493, -491, 378, 213, 493, 48, 497, 138, 441, 140, 351, 135, -123, 414, -7, -344, -320, 54, 400, 230, -80, -85, -76, -475, 342, 276, 340, 1, 477, 158, -378, 146, 274, -355)
(-315, 212, 432, 236, 423, -389, 67, -313, 365, 416, -180, -121, -472, 56, -468, -234, -305, 173, 444, -348, 261, -120, 249, -466, 491, -349, -140, 49, 99, 383, -321, 127, 447, -199, 443, -89, -384, -102, -106, 253, 495, 500, 8, -354, 115, 488, -88, 471, -291, 323, 349, 68, -253, -35, 34, -175, 172, -162, 123, -266, -6, 321, 481, -116, -60, 227, 163, -24, 56, 114)
```
3. The sum over all coordinates of the LWE secret is an even number.

We integrate the hints as follows. (Note: `vec` is simply an alias for `numpy.array`)
```py
>>> n = 70
>>> v_1 = vec( [1] + [0]*(n-1) )
>>> v_2 = vec( [0] + [1] + [0]*(n-2) )
>>> v_3 = vec( [0]*2 + [1] + [0]*(n-3) )
>>> v_4 = vec( [0]*3 + [1] + [0]*(n-4) )
>>> v_5 = vec( [-459, -441, 107, -207, 30, 358, -221, -483, 457, 96, 118, -241, 400, -478, 374, -46, -376, 415, 213, 476, -195, 25, -486, 444, 228, 313, -252, -182, -314, 105, -248, 163, 489, -388, 222, 110, -493, -491, 378, 213, 493, 48, 497, 138, 441, 140, 351, 135, -123, 414, -7, -344, -320, 54, 400, 230, -80, -85, -76, -475, 342, 276, 340, 1, 477, 158, -378, 146, 274, -355] )
>>> v_6 = vec( [-315, 212, 432, 236, 423, -389, 67, -313, 365, 416, -180, -121, -472, 56, -468, -234, -305, 173, 444, -348, 261, -120, 249, -466, 491, -349, -140, 49, 99, 383, -321, 127, 447, -199, 443, -89, -384, -102, -106, 253, 495, 500, 8, -354, 115, 488, -88, 471, -291, 323, 349, 68, -253, -35, 34, -175, 172, -162, 123, -266, -6, 321, 481, -116, -60, 227, 163, -24, 56, 114] )
>>> v_7 = vec( [1]*n )
>>> 
>>> lattice.integratePerfectHint( v_1, 0 ) #Coordinates
>>> lattice.integratePerfectHint( v_2, 0 )
>>> lattice.integratePerfectHint( v_3, 2 )
>>> lattice.integratePerfectHint( v_4, 1 )
>>> lattice.integratePerfectHint( v_5, -1670 ) #Inner products
>>> lattice.integratePerfectHint( v_6, 2381 )
>>> lattice.integrateModularHint( v_7, 0, 2 ) #Sum over all coordinates is even

```
Integrating the hints reduces the security of our LWE instance. When now running BKZ on the lattice, the secret will be recovered significantly faster. The required runtime on a laptop is roughly half a minute, and the required BKZ blocksize is 3.
```py
>>> lattice.reduce()
>>> lattice.successBlocksize
3
>>> lattice.s
array([ 0,  0,  2,  1,  2, -1, -1,  0, -1,  1,  1,  1,  0,  1,  1, -3, -2,
       -2, -3, -1,  1,  1, -1,  2, -1,  0, -2,  0,  1,  0,  0,  0,  1,  0,
        1,  0,  0,  0,  0,  0,  1,  0,  0,  0, -1,  0, -1, -1, -1,  1,  1,
        1,  0,  0,  0,  0,  2,  1,  0, -1,  1, -2,  0,  0, -2,  0,  1,  0,
       -1,  0])
```

### Generating LWE instances
Our library implements key generation algrotihms for various LWE-/NTRU-based schemes. To generate an LWE instance `(A,b,q)` with secret `s` and error `e`, simply run
```py
>>> A,b,q,s,e = generateLWEInstance(scheme)
```
where `scheme` can be any of the following strings:
* `"Kyber512"`, `"Kyber768"`, `"Kyber1024"`,
* `"Dilithium2"`, `"Dilithium3"`, `"Dilithium5"`,
* `"Falcon2`", `"Falcon4`", `"Falcon8`", `"Falcon16`", `"Falcon32`", `"Falcon64`", `"Falcon128`", `"Falcon256`", `"Falcon512`", `"Falcon1024`",
* `"NTRU-HPS-509`", `"NTRU-HPS-677`", `"NTRU-HPS-821`", `"NTRU-HRSS`".

Additionally, one can generate Kyber-like toy instances, where both secret and error follow a binomial distribtuion with parameter `eta`, as follows:
```py
A,b,q,s,e = generateToyInstance(n,m,q,eta)
```

**Disclaimer:** Our key generation algorithms are not suitable for production enviroments!

## Reproducing experiments from the paper

To recreate our experiments for perfect hints, please run the following commands:
```console
python3 experiments.py Kyber512 -hints="192:256:1" -trials=32
python3 experiments.py Kyber512 -hints="220:256:1" -trials=32 -centered
python3 experiments.py Kyber768 -hints="359:390:1" -trials=32
python3 experiments.py Falcon512 -hints="220:256:1" -trials=32
python3 experiments.py NTRU-HRSS -hints="310:350:1" -trials=32
python3 experiments.py Dilithium2 -hints="440:503:4" -trials=16
```
To recreate our experiments for modular hints, please run the following commands:
```console
python3 experiments.py Kyber512 -hints="440:455:1" -trials=16 -modular
python3 experiments.py Falcon512 -hints="440:455:1" -trials=16 -modular
python3 experiments.py NTRU-HRSS -hints="610:625:1" -trials=16 -modular
python3 experiments.py Kyber768 -hints="690:705:1" -trials=16 -modular
python3 experiments.py Dilithium2 -hints="870:885:1" -trials=16 -modular

```
If you want to run the experiments in verbose mode, simply add the flag `-verbose` to the above commands.

## Acknowledgments

For generating Falcon keys we use Thomas Prest's great [falcon.py](https://github.com/tprest/falcon.py) library.