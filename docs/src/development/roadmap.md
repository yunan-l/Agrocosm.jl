# Roadmap

## Phase 1: completed foundation

The rainfed single-crop C3/C4 pathway, daily process order, CPU/GPU kernels,
precision support, balance diagnostics, high-level simulation API, lifecycle
state architecture, and v2 checkpoints are implemented.

## Phase 2: modular and differentiable runtime

Near-term work includes:

1. a one-day transition over `ProcessModules` and `ModelState`;
2. explicit active/inactive parameter boundaries and Enzyme CPU tests;
3. soil/ecosystem spin-up and spin-up-to-transient continuity;
4. stable interfaces for alternative photosynthesis and stomatal models;
5. complete soil/climate output metadata;
6. multi-crop stands, rotations, and broader management.

## Later phases

Later work targets gradient-based calibration, data assimilation, hybrid
process–machine-learning models, broader validation, and large-domain
applications.

Detailed implementation notes remain in `docs/roadmap.md` in the source tree.
