# fluid_pde
Fun filled flowy fluid dynamics partial differential equation solver with 3rd
order accuracy (PPM SSPRK4) and minimal dependencies.

Cell center vs cell average is not distinguished so this is only 3rd order
accurate even with PPM and SSPRK4.

Supersonic wind:

![supersonic_wind](supersonic.png)


Explosive lens:

![explosive_lens](explosive_lens.png)


Random:

![post_blast](post_blast.png)

## Methods
Multithreaded with shared memory.

Time integrators:
- Euler
- RK2
- SSPRK3 (3rd order strong-stability preserving Runge-Kutta)
- SSPRK4 (4th order strong-stability preserving Runge-Kutta)

Riemann solvers:
- HLLC
- HLLE

Reconstruction methods (with van Leer limiter):
- Piecewise linear (PLM)
- Piecewise parabolic (PPM)

## Usage
Build dependencies: g++ or similar for C++17, make, pthreads

```
cd src
make
./fluid
```

View in a browser while running: `cd viewer && python -m http.server` and
open browser to `http://localhost:8000/` (via
[websocket_ctube](https://github.com/bryance-oyang/websocket_ctube))

In `src/init_cond` make a initial condition file from `template_init_cond.h` and
include in `init_cond.cc`.

Settings in `src/config.h`
