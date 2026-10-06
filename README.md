# Green's boundary condition solver for the nonautonomous Fisher-KPP equation

This is a C++ solver for propagating fronts in the nonautonomous Fisher-KPP equation, using the **Green's function boundary condition (GBC)** method to simulate an infinite domain on a finite one. The GBC method applies a dynamic boundary condition at the leading edge of the front to account for the part of the solution that lies outside the finite domain.

This code accompanies the paper:

> T. Tsubota, S. Mahajan, A. van Kan, and E. Knobloch, *Accurate simulation of pulled and pushed fronts in the nonautonomous  Fisher-Kolmogorov-Petrovsky-Piskunov equation* (2026).

The solver integrates, in a frame moving with speed $c(t)$,

$$u_t = a(t)\,u(1-u) + d(t)\,u_{xx} + c(t)\,u_x \qquad \text{(fisher-basic)}$$

$$u_t = a(t)\,u + d(t)\,u_{xx} + c(t)\,u_x + u^2 - (1 + a(t))\,u^3 \qquad \text{(fisher23)}$$

where the frame speed is tied to the diffusion coefficient by $c(t) = \gamma\, d(t)$ for a constant $\gamma$.

**Note:** In this codebase, including the various descriptions below, the spatial coordinate $x$ is the comoving coordinate $z = x_{\text{lab}} - \int_0^t c(t')\,dt'$ in the text.

## Installation

Requirements:

- A C++14 compiler
- [CMake](https://cmake.org) ≥ 3.19
- [HDF5](https://www.hdfgroup.org/solutions/hdf5/), including the C++ library
- [TCLAP](https://tclap.sourceforge.net) (header-only command-line parser)
- [SuperLU](https://github.com/xiaoyeli/superlu)
- [Armadillo](https://arma.sourceforge.net) **built with SuperLU support**. The code was developed and tested with Armadillo 15.0.1.

Armadillo must be built with `-DARMA_USE_SUPERLU=ON`, which prebuilt packages may not do. The steps below build it from source and install it to `~/opt/armadillo`, which is where the `Makefile` looks by default.

### macOS

1. Install the Xcode Command Line Tools: `xcode-select --install`.
2. Install [Homebrew](https://brew.sh).
3. Run `./setup.sh`. This installs CMake, HDF5, TCLAP, and SuperLU with Homebrew, then clones Armadillo 15.0.1, builds it with SuperLU, and installs it to `~/opt/armadillo`.

Other operating systems are likely to work but are not tested. Examine `setup.sh` to see what it does, and adapt it to your system as necessary.

## Building

From the repository root:

```sh
make
```

You can ignore this compiler warning, which comes from Armadillo:

```
warning: definition of implicit copy constructor for 'SpMat_MapMat_val<double>' is deprecated ...
```

## Quick start

```sh
./build/fronts -p fisher-basic -P "a const 1 d const 1 gamma manual 2" -i step \
    -T 100 -d 0.1 -x -30 -X 30 -D 0.01 -o 10 -b 10
```

This simulates the Fisher equation $u_t = u(1-u) + u_{xx} + 2u_x$ from a step initial condition, in a frame moving at the asymptotic speed $c = 2$, on $x \in [-30, 30]$ with $dx = 0.01$, up to $T = 100$ with $dt = 0.1$, saving every 10th time step, using a GBC buffer of width 10. It takes a few seconds.

Output is written to `sims/` in the current working directory, which is created if needed.

## Command-line reference

| Flag | Long form      | Default                                  | Meaning |
| ---- | -------------- | ---------------------------------------- | ------- |
| `-p` | `--pde`        | `fisher-basic`                           | PDE: `fisher-basic` or `fisher23` |
| `-P` | `--pdeParams`  | `"a const 1 d const 1 gamma manual 2"`   | Coefficients $a(t)$, $d(t)$ and frame parameter $\gamma$ (see below) |
| `-i` | `--ic`         | *(required)*                             | Initial condition (see below) |
| `-t` | `--t0`         | `0`                                      | Initial time. Only `0` is supported. |
| `-T` | `--t1`         | `10`                                     | Final time |
| `-d` | `--dt`         | `0.1`                                    | Time step |
| `-x` | `--x0`         | `-20`                                    | Left boundary of the domain |
| `-X` | `--x1`         | `20`                                     | Right boundary of the domain |
| `-D` | `--dx`         | `0.01`                                   | Grid spacing |
| `-o` | `--tout`       | `1`                                      | Save output every `tout` time steps |
| `-b` | `--buf`        | `0`                                      | Right boundary condition: buffer width $\delta > 0$ uses GBC; `0` uses zero Dirichlet; `-1` uses zero Neumann |
| `-m` | `--impulse`    | *(none)*                                 | Optional HDF5 cache for the Green's function integrals (see below) |
| `-n` | `--icintegral` | *(none)*                                 | Optional HDF5 cache for the initial-condition integral (see below) |
| `-q` | `--quiet`      |                                          | Less progress output |

### Coefficients (`-P`)

The string lists `a`, `d`, and `gamma`, each followed by a type and its parameters, e.g. `"a linear 1 0.01 d const 1 gamma manual 2.5"`. Omitted coefficients default to `a const 1`, `d const 1`, `gamma manual 2`. The frame speed is always $c(t) = \gamma\, d(t)$.

Types for `a` and `d`:

| Type        | Form                         | Parameters                                                        |
| ----------- | ---------------------------- | ----------------------------------------------------------------- |
| `const`     | $c$                          | `const c`; `const` → $1$                                          |
| `linear`    | $a + bt$                     | `linear a b`; `linear b` → $1 + bt$; `linear` → $1 + t$           |
| `quadratic` | $a + bt + ct^2$              | `quadratic a b c`; `quadratic a c` → $a + ct^2$; `quadratic c` → $1 + ct^2$; `quadratic` → $1 + t^2$ |
| `exp`       | $a e^{bt}$                   | `exp a b`; `exp b` → $e^{bt}$; `exp` → $e^t$. Requires $b \neq 0$. |
| `algebraic` | $a + bt^p$                   | `algebraic a b p`; `algebraic b p` → $1 + bt^p$; `algebraic p` → $1 + t^p$; `algebraic` → $1 + t^{1/2}$ |
| `tanh`      | $a\tanh(t - b) + c$          | `tanh a b c`; `tanh a b` → $c = 0$; `tanh a` → $b = c = 0$; `tanh` → $\tanh t$ |

### Initial conditions (`-i`)

| Name                     | Form                                                     | Default parameter |
| ------------------------ | -------------------------------------------------------- | ----------------- |
| `step`                   | $1$ for $x < 0$, $0$ otherwise                           |                   |
| `exp λ` (or `exponential`) | $e^{-\lambda \lvert x \rvert}$                         | $\lambda = 1$     |
| `tanh a`                 | $\tfrac12 \tanh(-ax) + \tfrac12$                         | $a = 1$           |
| `gaussian σ`             | $\frac{1}{\sqrt{2\pi}\,\sigma} e^{-x^2/(2\sigma^2)}$     | $\sigma = 1$      |
| `line`                   | $1$ everywhere                                           |                   |

### Green's function caches (`-m`, `-n`)

By default, the Green's function integrals needed by the GBC are computed inline at every step. Alternatively, pass a file path to `-m` (and `-n` for the initial-condition integral):

- if the file exists, the integrals are loaded from it;
- otherwise, they are precomputed for the whole run, saved to that path, and then loaded.

This lets you reuse the integrals across runs with the same coefficients, `dt`, `T`, and buffer width. The cache stores three $N \times N$ matrices with $N = T/dt$, so it can become very large.

## Examples

### Naïve boundary condition

```sh
./build/fronts -p fisher-basic -P "a const 1 d const 1 gamma manual 2" -i step \
    -T 100 -d 0.1 -x -30 -X 30 -D 0.01 -o 10 -b 0
```

The same run as the [quick start](#quick-start), but with a zero Dirichlet condition at the right boundary instead of the GBC (`-b 0`; use `-b -1` for zero Neumann).

### Shallow initial condition

```sh
./build/fronts -p fisher-basic -P "a const 1 d const 1 gamma manual 2.5" -i "exp 0.5" \
    -T 100 -d 0.1 -x -30 -X 30 -D 0.01 -o 10 -b 10
```

Starts from $u(x, 0) = e^{-|x|/2}$. A shallow initial condition with steepness $\lambda < 1$ propagates at speed $\lambda + 1/\lambda$, which is $2.5$ here, faster than the pulled speed $2$. The frame moves at $c = 2.5$, so the front stays in place. The GBC accounts for the part of the initial condition beyond the right boundary.

### Time-dependent growth rate

```sh
./build/fronts -p fisher-basic -P "a linear 1 0.01 d const 1 gamma manual 2.5" -i step \
    -T 100 -d 0.02 -x -40 -X 30 -D 0.01 -o 50 -b 10
```

Simulates $u_t = a(t)\,u(1-u) + u_{xx} + 2.5\,u_x$ with $a(t) = 1 + 0.01t$. The front speeds up from $2$, so in this frame it first drifts left and then turns around.

### Pushed front

```sh
./build/fronts -p fisher23 -P "a const 0.5 d const 1 gamma manual 1.4434" -i step \
    -T 200 -d 0.05 -x -30 -X 30 -D 0.01 -o 20 -b 10
```

Simulates the quadratic-cubic equation $u_t = \tfrac12 u + u_{xx} + u^2 - \tfrac32 u^3 + \gamma u_x$. For $a = 1/2$, the front is pushed, with speed $v^\dagger = 3\sqrt{3/4} - \sqrt{4/3} \approx 1.4434$, faster than the pulled speed $\sqrt{2} \approx 1.414$.

## Output

Each run writes `sims/<timestamp>_sim.h5`, where the timestamp is in UTC, and appends a row with the run parameters to `sims/sims.csv`.

The HDF5 file contains:

| Dataset | Shape      | Contents |
| ------- | ---------- | -------- |
| `t`     | `(Nt,)`    | Output times, from `t0` to `t1 - tout*dt` |
| `x`     | `(Nx,)`    | Grid points in the comoving frame |
| `u`     | `(Nt, Nx)` | Solution; `u[i]` is the profile at time `t[i]` |

and the attributes `pde`, `pde_params`, `ic`, `t0`, `t1`, `dt`, `x0`, `x1`, `dx`, `tout`, and `buf`.

Reading it in Python:

```python
import h5py

with h5py.File("sims/20261005T160259Z_sim.h5", "r") as f:
    t, x, u = f["t"][:], f["x"][:], f["u"][:]
    params = dict(f.attrs)
```

To convert to lab-frame positions, add $\int_0^t c(t')\,dt' = \gamma \int_0^t d(t')\,dt'$.

## Notes and limitations

- **Initial conditions outside the domain.** The GBC accounts for the part of the initial condition lying in the linear region only for `exp` initial conditions. For the others, that part is treated as zero, which is a good approximation when the initial condition is steep (e.g. `step`, `gaussian`).
- **Frame speed.** The closed-form Green's function requires $c(t) = \gamma\, d(t)$, so the frame cannot follow an arbitrary $c(t)$. See the paper's appendix on choosing $\gamma$.
- **Initial time.** Only `t0 = 0` is supported.
- **Time step.** For time-dependent coefficients, use a small `dt` to resolve the time dependence; large `dt` gives poor agreement with theory.
- **Memory.** The full solution is held in memory and written at the end of the run, and nothing is saved if the run is interrupted. Use `-o` to reduce the output size for long runs.
- **Output names.** Output files are named by the second at which the run finished, so two runs finishing in the same second in the same directory overwrite each other.

## Extending the code

**New PDE.** Subclass `FisherBasicSolver` and override `update_rhs` with the new nonlinearity (see `Fisher23Solver`), then add it to the solver selection in `main.cxx`. The GBC assumes the linearization about $u = 0$ is $u_t = a(t)u + d(t)u_{xx} + c(t)u_x$, so new nonlinearities must preserve that.

**New coefficient type.** Add a `Coefficient` subclass to `options.h` and `options.cxx` implementing both `func(t)` and its antiderivative `integral(t)`, which the Green's function integrals use. Then add it to `make_coefficient` and `scale_params` in `main.cxx`. `scale_params` describes how to multiply the coefficient by a constant: for $a + bt$ both parameters scale, but for $ae^{bt}$ only $a$ does. It is used to build $c(t) = \gamma\, d(t)$ from `d`.

**New initial condition.** Add an `InitialCondition` subclass to `options.h` and `options.cxx`, and add it to the initial-condition selection in `main.cxx`. Its contribution from outside the domain will be treated as zero.

## License

MIT. See [`LICENSE`](LICENSE).
