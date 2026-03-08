# AGENTS.md — ICON-Land (JSBACH4+)

Fortran 90/2003 land surface model, part of the ICON Earth system model.
Repository hosted on GitLab at DKRZ. License: BSD-3-Clause.
All source code is in `src/` (216 `.f90` files across 27 process directories).

## Build and Test

There is **no standalone build system** in this repository. JSBACH is compiled
as a component of the ICON model using ICON's build infrastructure.

- **Documentation build:** `./build-doc.sh [-q] [-f] [-v] [-c]` (uses FORD + uv)
- **DSL preprocessor:** `python3 scripts/dsl4jsb/dsl4jsb.py -d src -t output_dir -k`

There are **no unit tests**. Testing is done via ICON integration experiments:
`land_jsbach_R2B4_test` (with/without QUINCY, CPU/GPU), `atm_nwp_jsbach_test`,
`aes_amip_test`, `exp.land_quincy_canopy_R2B4_test`,
`exp.land_quincy_standalone_test_1d`/`_2m`/`_2y`.
Compilers: **NAG** (debug mode) and **nvhpc** (NVIDIA, CPU/GPU).

## CI Checks (GitLab CI)

Two scripts run on merge request pipelines:

- `./scripts/ci_scripts/check_file_format.sh` — enforces:
  valid UTF-8 encoding, ASCII-only characters, no trailing whitespace.
  Checked patterns: `*.f90`, `*.F90`, `*.inc`, `*.incf`, `*.h`, `*.c`, `*.cu`.
  Checked directories: `src/`, `scripts/`.
- `./scripts/ci_scripts/check_license.sh` — all `.f90`/`.F90`/`.inc`/`.incf` files
  in `src/` must contain the required license header (see below).

**Forbidden Doxygen commands:** `@par Revision history` and `@author` are rejected by CI.

## Required License Header

Every Fortran source file in `src/` must begin with this header:

```fortran
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
```

## File and Directory Conventions

File naming pattern: `mo_<domain>_<role>.f90` (`mo_` prefix is standard MPI-M convention).

| Suffix | Purpose |
|---|---|
| `_interface` | Task registration, `Integrate` and `Aggregate` methods |
| `_config_class` | Configuration type and namelist reading |
| `_memory_class` | Memory (state variable) declarations |
| `_process` | Pure computational / scientific routines |
| `_init` | Initialization routines |
| `_constants` | Physical and model constants |
| `_parameters` | Tunable parameters |
| `_util` | Utility / helper functions |

Each directory under `src/` corresponds to one physical process (e.g. `radiation/`,
`hydrology/`, `flake/`). The framework lives in `src/base/` and `src/shared/`.

## Module Structure

Every module follows this pattern without exception:

```fortran
IMPLICIT NONE
PRIVATE
PUBLIC :: exported_symbol1, exported_symbol2
CHARACTER(len=*), PARAMETER :: modname = 'mo_module_name'
```

Every procedure defines a local routine name for error reporting:
```fortran
CHARACTER(len=*), PARAMETER :: routine = modname//':procedure_name'
```

## Naming Conventions

- **Types:** `t_` prefix, snake_case — `t_jsb_model`, `t_rad_config`, `t_jsb_var_real2d`
- **Task types:** `tsk_` prefix — `tsk_surface_radiation`, `tsk_albedo`
- **OOP methods:** PascalCase — `Get_model`, `Register_task`, `Init_rad_memory`
- **Compute routines:** snake_case with `calc_` prefix — `calc_radiation_surface_net`
- **Variables:** snake_case — `sw_srf_net`, `no_of_models`
- **Temporaries:** `z` prefix — `zswvis_net`, `zsw_net`
- **Logicals:** `l_` prefix — `l_changes_fractions`, `l_ice_on_lakes`
- **Enum constants:** UPPER_CASE with trailing underscore — `RAD_`, `SEB_`, `FLAKE_`
- **Physical constants:** PascalCase — `LatentHeatVaporization`, `molarMassCO2`
- **Implementation constants:** UPPER_CASE — `SHORT_NAME_LEN`, `WB_ERROR`

## Import Style

Always use explicit imports: `USE module_name, ONLY: symbol1, symbol2`.

Ordering convention:
1. External / host model modules (`mo_kind`, `mo_exception`)
2. Framework base modules (`mo_jsb_control`, `mo_jsb_class`)
3. Process-specific modules
4. DSL macro imports (`dsl4jsb_Use_processes`, `dsl4jsb_Use_memory`)

## Error Handling

- **Fatal errors:** `CALL finish(routine, 'error message')` — halts execution
- **Diagnostics:** `CALL message(routine, 'info message')` — non-fatal logging
- **Debug guard:** `IF (debug_on()) CALL message(...)` — conditional output
- Every procedure must define a `routine` parameter for traceability (see Module Structure)

## Function and Subroutine Style

- All functions use an explicit `RESULT(return_value)` clause
- All arguments have explicit `INTENT(in)`, `INTENT(out)`, or `INTENT(inout)`
- Use `ELEMENTAL` or `PURE` where applicable
- Mark GPU-eligible scalar routines with `!$ACC ROUTINE SEQ`

## Comment Style

- `!>` for module-level and procedure-level documentation (Doxygen-compatible)
- `!!` for continuation of doc blocks; `!<` for trailing member documentation
- `!` for inline comments; include units in brackets: `! [W m-2]`, `! [kg m-2 s-1]`
- Mark DSL variable access with `! in` or `! out` comments:
  `dsl4jsb_Get_var2D_onChunk(A2L_, swvis_srf_down)   ! in`

## DSL4JSB Macros

The codebase uses a custom DSL preprocessor (`scripts/dsl4jsb/dsl4jsb.py`) that
translates `dsl4jsb_*` macro calls into standard Fortran. Key macros:

- `dsl4jsb_Use_processes(PROC1_, PROC2_)` — USE statements for process modules
- `dsl4jsb_Use_config(PROC_)` / `dsl4jsb_Use_memory(PROC_)` — USE config/memory class
- `dsl4jsb_Def_config(PROC_)` / `dsl4jsb_Def_memory(PROC_)` — declare pointer to config/memory
- `dsl4jsb_Get_config(PROC_)` / `dsl4jsb_Get_memory(PROC_)` — retrieve from tile
- `dsl4jsb_Get_var2D_onChunk(PROC_, varname)` / `dsl4jsb_Get_var3D_onChunk` — get variable slice
- `dsl4jsb_Aggregate_onChunk(PROC_, var, aggregator)` — aggregate variable
- `dsl4jsb_Lctlib_param(param)` — access land cover type library parameter

## Process Module Pattern

Each process follows a standard file structure (templates in `templates/`):

1. `mo_<name>_interface.f90` — task registration + Integrate/Aggregate
2. `mo_<name>_config_class.f90` — configuration type + namelist
3. `mo_<name>_memory_class.f90` — state variable declarations
4. `mo_<name>_process.f90` — scientific computation (optional)
5. `mo_<name>_init.f90` — initialization (optional)

To register a new process, add entries in these 4 files:
`src/base/mo_jsb_process_class.f90` (process ID enum),
`src/base/mo_jsb_model_usecases.f90` (usecase configuration),
`src/base/mo_jsb_model_init.f90` (task queue and init calls),
`src/base/mo_jsb_process_factory.f90` (interface, config, memory).

## Conditional Compilation and GPU

Feature guards wrap module content with preprocessor directives:
`#ifndef __NO_JSBACH__` (all JSBACH code), `#ifndef __NO_QUINCY__` (QUINCY features),
`#ifdef __QUINCY_STANDALONE__` (standalone QUINCY), `#ifndef __NO_JSBACH_HD__`
(hydrological discharge), `#ifdef _OPENACC` (GPU acceleration).

OpenACC compute loops use `!$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)`.

## Branch and Commit Conventions

- Feature branches: `feature/[issue-number]-brief-feature-info`
- Commit messages: short imperative style (e.g. "Add flake process", "Fix albedo calculation")
- Optional prefix for subsystem: `IQ:` for QUINCY changes
