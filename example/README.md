# Example

This directory contains an example using `fhelib` as well as some additional explanations on design choices.


## Parameter Configuration
Currently, the only supported scheme is the leveled BGV scheme where multiple configurations take place during build time.
Most importantly, we specify the ring layer in `TIIFHE_RING` (default: `hexl`) and the message layer in `TIIFHE_MSG` (default: `mpz`).
The library can easily be extended with custom layers for both layer types as long as a new layer provides the required functions.

Additionally, the following scheme parameters are configured via `CMakeLists.txt` or via command line arguments to `cmake`:

- the polynomial degree `TIIFHE_N`;
- the plaintext modulus `TIIFHE_T`;
- the ciphertext modulus `TIIFHE_Q`;
- the key switching modulus `TIIFHE_P`; and
- the key switching decomposition number `TIIFHE_OMEGA`.

For the currently provided layers, only a power-of-two degree is supported and we require `p == 1 % (2 * TIIFHE_N)` for all primes `p`.
This includes the plaintext modulus `TIIFHE_T` if packing is to be used and all primes have to be co-prime.
The sage script `params.py` provides a simple prime generation and security estimate.
For more involved parameter settings, a powerful parameter generator such as [fhegen](https://github.com/Crypto-TII/fhegen) can be used.


## Project Integration
The example assumes that a library copy exists in the subdirectory `fhelib/` (here, it simply is a symbolic link to the parent directory).
Execute the following to build and run the example:
```
cmake -S . -B build
cmake --build build
./build/example
```


## High-Level API
Due to the build time configuration, all functions have an implicit global context with the scheme parameters.
Within this context, there exists a high-level API (`tiifhe_bgv_XXX`) and a low-level API (`tiifhe_bgv_XXX_idx`).
For our example, we use the high-level API providing further information in the comments along the way.

## Low-Level API
In our tests (see also `tests/bgv.c`), we mostly use the more efficient but more difficult to use low-level API.
Due to the implicit context, we can access individual primes of the ciphertext and key switching modulus using index-based functions where, for index `i`, `0 <= i < TIIFHE_QLEN` operates on the ciphertext modulus and `TIIFHE_QLEN <= i < TIIFHE_QPLEN` operates on the key switching modulus.
For key switching, we decompose into `TIIFHE_OMEGA` groups where the prime at index `i` belongs to group `i % TIIFHE_OMEGA`.

To reduce common bugs, the low-level API manages the NTT state of each ciphertext (`ntt`) and tracks the current message scaling for each individual prime with the `drop` value.
However, we expose all data structures and polynomial functions transparently in the header files; experts using the low-level API are thus encouraged to modify data to their heart's content.

## Error Logging
With the build option `TIIFHE_LOG_ERROR`, each decryption appends the ciphertext error to a file `error.csv` in the current directory.

## Seed Expansion
With the build option `TIIFHE_EXPAND_SEED`, the library expands the random values instead of storing seeds for key material.

## Seed Fixing
With the build option `TIIFHE_FIX_SEED`, the library generates seeds deterministically, simplifying debugging or generating reference values.
