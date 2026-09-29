# AGENTS.md

Guidance for AI coding agents working in the ICON-Land (JSBACH) repository. The [ICON-Land documentation](https://jsbach.gitlab-pages.dkrz.de/jsbach) (FORD-generated from inline comments plus `doc/guide/`) and the [GitLab wiki](https://gitlab.dkrz.de/jsbach/jsbach/-/wikis/home) hold the full technical/scientific docs. The wiki is the richest source; where it and the source tree disagree, prefer the source tree.

This directory is normally checked out as `externals/jsbach` inside an ICON working copy. See the parent `AGENTS.md` for ICON-wide guidance.

## What this is

ICON-Land — historically and still in repo/namelist naming "JSBACH" — is the land component of ICON: a Fortran framework for bio-geophysical and bio-geochemical land-surface processes (energy/water balance, soil, snow, hydrology, radiation, phenology, photosynthesis, carbon/nitrogen/phosphorus cycling, land cover change, hydrological discharge).

It is a **framework plus interchangeable processes**, not a monolithic land model. Currently two core sets of processes are used, these sets share several common processes but also use custom processes:

- **JSBACH4** — the classic MPI-M land scheme (custom processes: `ASSIMI_`, `PHENO_`, `CARBON_`, `DISTURB_`, `NLCC_`, ...).
- **QUINCY** (a.k.a. IQJ, "ICON-Land with QUINCY configuration and JSBACH physics") — the MPI-BGC terrestrial biosphere model (custom processes: `VEG_`, `SB_`, `Q_ASSIMI_`, `Q_PHENO_`, `Q_SYL_`, `Q_AGR_`).

ICON-Land **cannot be compiled or run on its own**. It always needs a *host model* that provides I/O, domain decomposition, parallelization and time control — currently ICON, the only supported host, via the adapter file `src/shared/mo_util_jsbach.f90` on the ICON side. See "Adapters" below.

Development happens in the standalone repo `gitlab.dkrz.de:jsbach/jsbach`, but building and testing happens inside a host-model clone (normally `gitlab.dkrz.de:icon/icon-mpim`).

## Core concepts

You cannot usefully edit this code without these. They are all wiki-documented; the authoritative code lives in `src/base/`.

### Model → tiles → process memory

A **model** (`t_jsb_model`, `src/base/mo_jsb_model_class.f90`) owns a **hierarchical tile tree**, a list of **processes**, and a **task queue**. It extends a hierarchical state machine (`t_Hsm`, `src/base/mo_hsm_class.f90`).

**Tiles** (`t_jsb_tile_abstract`, `src/base/mo_jsb_tile_class.f90`, `src/base/mo_jsb_tile.f90`) subdivide each grid cell by surface type, e.g. `box → {land, lake}`, `land → {veg, glacier}`, `veg → {pft01 … pftNN}`. Each tile carries a fraction `tile%fract`; children sum to their parent, and the `box` fraction is the land(+lake) fraction of the ICON grid cell. Each tile holds **its own memory per process**: `tile%mem(PROC_)%p`.

**Process memory** (`t_jsb_memory`, `src/base/mo_jsb_memory_class.f90`) is a list of variables of type `t_jsb_var_real1d` (1D - scalar/monitoring), `t_jsb_var_real2d` (2D - horizontal layer) and `t_jsb_var_real3d` (3D - horizontal and vertical layer), registered with `mem%Add_var(...)`. The process memories of the two QUINCY processes `VEG_` and `SB_` further define biogeochemical materials (BGCMs), which are collections of 2D or 3D variables (see below). Note: Output variable names are composed as `PROC_variable_tile` so grepping for an output field name will fail; grep for the bare `variable` part inside the owning process directory.

### Processes

Each scientific process lives in its own `src/<name>/` directory and is usually identified by an enumerator ending in `_`, whereby some processes have subdivisions represented by additional enumerators. The enumerators are defined in `src/base/mo_jsb_process_class.f90`.

- Commonly used processes: `srf_energy_bal/` (`SEB_`), `turbulence/` (`TURB_`), `soil_snow_energy/` (`SSE_`), `hydrology/` (`HYDRO_`), `hd/` (`HD_`, hydrological discharge + ocean coupling), `radiation/` (`RAD_`), `anthropogenic_lcc/` (`ALCC_`), `pre_and_post_lcc/` (`PPLCC_`), `atmland/` (`A2L_`, `L2A_`).
- JSBACH custom processes: `assimilation/` (`ASSIMI_`), `phenology/` (`PHENO_`), `carbon/` (`CARBON_`), `climate_buffer/` (`CLIMBUF_`), `disturbance/` (`DISTURB_`, `FUEL_`, `FLCC_`, `WLCC_`), `natural_lcc/` (`NLCC_`), `forest_age/` (`FAGE_`).
- QUINCY custom processes: `q_vegetation/` (`VEG_`), `q_soil_biogeochemistry/` (`SB_`), `q_assimilation/` (`Q_ASSIMI_`), `q_phenology/` (`Q_PHENO_`), `q_sylviculture/` (`Q_SYL_`), `q_agriculture/` (`Q_AGR_`).

Conventional module set per process `PROC` (see `templates/` and `templates/README.md`. Note: the templates might not be up-to-date and copying a comparable existing process is the recommended route):

| module | required | role |
| --- | --- | --- |
| `mo_PROC_interface.f90` | yes | defines the process **tasks**: their `Integrate` (`update_*`) and `Aggregate` (`aggregate_*`) procedures |
| `mo_PROC_config_class.f90` | yes | namelist parsing (`jsb_PROC_nml` / `lnd_*_nml`), the `active` switch, and creation of process-specific vertical grids; for QUINCY processes with BGCMs: ENUM with compartments |
| `mo_PROC_memory_class.f90` | if it has state | variable declarations + `Init_PROC_memory` with the `Add_var` calls |
| `mo_PROC_init.f90` | often | reads ic/bc files and sets initial/boundary conditions; called from `jsbach_init` |
| `mo_PROC_process.f90` | optional | the actual formulas, infrastructure-free |
| `mo_PROC_constants.f90` | optional | process parameters or constants, `PARAMETER` or `SAVE` + `init_PROC_constants()`;  |
| `mo_PROC_parameters.f90` | optional | process parameters — in addition to or instead of `mo_PROC_constants.f90`, `SAVE` + `init_PROC_parameters()` |
| `mo_PROC_util.f90` | optional | helper routines |
| `mo_PROC_init_lcc.f90` | lcc only | declares active/passive conserved quantities and tiles subject to changing land cover |

Some specifics:
- As opposed to JSBACH processes, several custom QUINCY processes keep their science in task specific modules (`mo_q_veg_growth.f90`, `mo_q_veg_turnover.f90`, ...) rather than in a single `mo_PROC_process.f90`. In addition, these task specific modules are not infrastructure-free.
- `A2L_` and `L2A_` share `src/atmland/mo_atmland_interface.f90`, have no config class, and their update routines are called directly from `mo_jsb_interface.f90`, not via the task queue.
- `PPLCC` is a process required by the lcc framework with a near minimal set of modules. The two tasks of the process call the same routines for the integrate and aggregate procedure.
- The file naming scheme of `VEG_` and `SB_` is — historically and still — inconsistent with some files having the current process enumerator as tag (`_veg_`, `_sb_`) while others carry an extended tag (`_q_veg_`, `_q_sb_`).

**Registering a new process** touches several infrastructure modules — miss one and it either crashes or silently does nothing:

1. add the ID + name mapping in `src/base/mo_jsb_process_class.f90` (`Get_process_name`, `Get_process_id`),
2. associate it with tiles in `src/base/mo_jsb_model_usecases.f90`,
3. queue its tasks in `jsbach_init` (`src/base/mo_jsb_model_init.f90`) — order matters, see below,
4. add interface + config class (+ memory class) to `src/base/mo_jsb_process_factory.f90` and bump `max_no_of_processes`; set `has_memory = .FALSE.` in `Create_process` if it has none,
5. call its `PROC_init` (and `PROC_init_lcc`) from `jsbach_init`,
6. call `init_PROC_constants()` and / or `init_PROC_parameters()` from `jsbach_setup_models` (`src/base/mo_jsb_base.f90`) if applicable.

### Process actions

Where a process runs on the tile tree (see associations in `src/base/mo_jsb_model_usecases.f90`) is an *action* per (process, tile), enumerated in `src/base/mo_jsb_process_class.f90`: `SKIP_`, `INHERIT_`, `AGGREGATE_`, `ON_LEAFS_`, `ON_TILE_`, `ON_SUBTREE_`. `ON_TILE_` on a tile gives descendants `INHERIT_` (their memory points at the parent's); `ON_LEAFS_` gives leaves `ON_TILE_` and intermediate tiles `AGGREGATE_`.

Query helpers (see FAQ): `model%Is_process_enabled` (configured at all), `tile%Is_process_active` (integrated *or* aggregated here), `tile%Is_process_calculated` (integrated here), `tile%Has_process_memory`. In addition, the config of a process `PROC` is often directly used to test if a process is active — using `dsl4jsb_Config(PROC_)%active` or (mostly in the infrastructure) `model%processes(PROC_)%p%config%active`.

### Task queue

The task queue (`src/base/mo_jsb_task_class.f90`, built in `jsbach_init` in `src/base/mo_jsb_model_init.f90`) fixes the order of `update_*` calls within a timestep. Tiles are walked **leaves → top** so that `aggregate_*` works. `update_`/`aggregate_` are ultimately dispatched from `task%Do_it` inside `Process_task_on_tile` (`src/base/mo_jsb_tile.f90`) — deep inside the state machine, which is why call sites are hard to grep for.

Ordering rules:

1. diagnostics needed by lcc processes but not changing state,
2. lcc tasks, sandwiched between `pplcc_pre_lcc` and `pplcc_post_lcc`,
3. everything else by scientific dependency. 

Currently two different task queues are defined depending on the model scheme (`model%config%model_scheme` — `MODEL_JSBACH` or `MODEL_QUINCY`). For the model scheme `MODEL_QUINCY` additional modularity on task queueing is added by the selected mode (`model%config%qmodel_id`) with `QCANOPY` only adding tasks required for bio-geophysical land-surface processes, `QPLANT` additionally to `QCANOPY` adding tasks required for plant bio-geochemistry, `QSOIL` additionally to `QCANOPY` adding tasks required for soil bio-geochemistry (whereby the implementation of `QSOIL` is not completed yet), and `QLAND` enabling all tasks.

Tasks are queued regardless of the process's `active` namelist flag, so **every `update_`/`aggregate_` routine should start with `IF (.NOT. tile%Is_process_active(XYZ_)) RETURN`**.

### Aggregation

`Aggregate_weighted_by_fract_2d`/`_3d` (`src/base/mo_jsb_tile_class.f90`) average child values weighted by tile fraction. The `l_aggregate_all` flag on `Add_var` decides whether children *lacking* the variable still contribute their fraction to the denominator (`.TRUE.`, correct for fluxes that must be per whole grid box, e.g. what the atmosphere receives) or are ignored entirely (`.FALSE.`, default, correct for intensive quantities like temperature).

A variable only reaches `_box` output if it is actually listed in the process's `aggregate_` routine — otherwise `_box` output will stay the initial value `initval_r` as specified in the `Add_var` call in the process memory module — which is often but not always zero. Note: technically it is possible to only add variables for certain tiles — aggregation routines can of course only be called for variables added along the tile tree.

### Land cover change (lcc) infrastructure

See [lcc-infrastructure](https://gitlab.dkrz.de/jsbach/jsbach/-/wikis/Documentation/ICON-Land-framework/lcc-infrastructure) in the GitLab wiki. Type definitions and main routines are located in `src/base/mo_jsb_lcc.f90`, `src/base/mo_jsb_lcc_class.f90`, `src/base/mo_jsb_cqt_class.f90`. Cover-fraction changes among tiles, as induced by lcc processes, often not only need to conserve total fractions but also (tracked) matter, aka conserved quantities (CQs). Variables that describe such to CQs get a `cons_quan_type_id` in their `Add_var` call, grouping them into **conserved quantity types (CQTs)** (same unit, same rank, same treatment). Currently `CARBON_` (JSBACH),  `VEG_` and `SB_` (QUINCY) declare CQs.

Three conservation mechanisms: **implicit scaling** (rescale in place — default behaviour), **passive relocation** (relocated to same variable, shrinking tile → extending tile), **active relocation** (relocated to different variable, shrinking tile → extending tile, e.g. biomass → litter/product pools). Provided proper definition, implicit scaling and passive relocation are handled automatically by the lcc infrastructure. Active relocation requires process specific tailored routines. For QUINCY, some BGCM variables are explicitly handled inside the relocating `ALCC` task.

The `PPLCC_` process brackets all lcc tasks: `pplcc_pre_lcc` converts CQs to per-grid-cell units, `pplcc_post_lcc` converts back — which is also what makes implicit scaling happen automatically.

lcc processes doing passive/active relocation must provide `mo_PROC_init_lcc.f90` and follow the protocol in their update routine:

```fortran
dsl4jsb_Get_lcc_relocations(PROC_NAME, lcc_relocations)
CALL init_lcc_reloc(lcc_relocations, options, tile, initial_area)
CALL start_lcc_reloc(lcc_relocations, options, lost_area, initial_area)
! ... process-specific handling of active CQs ...
CALL end_lcc_reloc(lcc_relocations, options, gained_area)
```

Current lcc processes: `NLCC_`, `WLCC_`, `FLCC_`, `FAGE_` (JSBACH) and `ALCC_` (JSBACH and QUINCY).

### Usecases

`src/base/mo_jsb_model_usecases.f90` — selected via `usecase` in `jsb_model_nml`. Each defines the tile hierarchy and the process/action association:

- `jsbach_lite` — minimal physics-only setup
- `jsbach_pfts` — full JSBACH with plant functional type (PFT) tiles
- `jsbach_forest_age_classes` — PFTs plus forest age classes (custom process: `FAGE_`)
- `quincy_13_pfts` — QUINCY (aka IQJ) with 13 PFTs (uses `data/lctlib_quincy_nlct14.def`)

`jsbach_lite` and `jsbach_pfts` have separate `_tmx` variants selected by `use_tmx` (explicit surface coupling in the AES/Sapphire atmosphere).

Tile construction passes name, `lct` (land cover type, `src/base/mo_jsb_lct_class.f90`), parent, `fract_varname` (which field in the fractions file, `bc_land_fract.nc`, supplies the initial fraction; `'equal'` means "split evenly, read nothing"), and the process/action arrays. Processes can also be attached conditionally afterwards via `tile%Set_process_action(...)`.

PFT parameters come from the land-cover-type library in `data/` (`lctlib_nlct21.def`, `lctlib_nlct21_cmip7.def`, `lctlib_quincy_nlct14.def`).

### Decomposition and grids

Horizontal: the global cell vector is split across MPI domains, each domain into blocks of `nproma` cells. 2D variables are `var(nproma, nblks)`, 3D are `var(nproma, nlev, nblks)`. Inside a task you work on one chunk, addressed by `options%iblk`, `options%ics`, `options%ice` — i.e. `var(options%ics:options%ice, options%iblk)`.

Vertical: `t_jsb_vgrid` (`src/base/mo_jsb_grid_class.f90`), created with `new_vgrid(...)` in a process config class, published with `Register_vgrid`, retrieved anywhere with `Get_vgrid('name')`. Convention: **`ubounds > lbounds` always**, and layer 1 is always the one at the ground with `lbounds(1) = 0.0_wp`, independent of whether the axis points down (soil) or up (canopy).

### Adapters (host-model abstraction)

ICON-Land should not `USE` host-model modules directly. Instead the host provides a fixed set of generic interface modules — `mo_jsb_parallel_iface`, `mo_jsb_domain_iface`, `mo_jsb_namelist_iface`, `mo_read_namelists_iface`, `mo_jsb_time_iface`, `mo_jsb_io_iface`, `mo_jsb_io_netcdf_iface`, `mo_jsb_grid_iface`, `mo_jsb_utils_iface`, `mo_jsb_varlist_iface`, `mo_jsb_vertical_axes_iface`, `mo_physical_constants_iface` — all implemented in a single host-side file, `src/shared/mo_util_jsbach.f90` in ICON. ICON-Land *requests* a feature through the generic interface; the host decides how to provide it.

**Rule: only files in `src/base/` may `USE` the adapter modules.** Everything outside `base` must go through a `base` module. If you need new host functionality you must define the host-agnostic interface here *and* implement it in the host adapter file (paired MRs — see "Repository layout and merge workflow").

### BGCMs (QUINCY)

Biogeochemical materials (BGCMs) are collections of 2D or 3D variables. A BGCM can either have one layer, which is a set of variables representing chemical elements of a pool or a flux (currently carbon, nitrogen, phosphorus and different isotopes), or a two layered set of biological compartments, which themselves each contain such a set of chemical elements.
`mo_lnd_bgcm_class.f90`, `mo_lnd_bgcm_store*.f90`, `mo_lnd_bgcm_idx.f90` hold type definitions for and routines handling QUINCY's BGCMs as multi-element structures and routines to synchronise the actual variables with their matrix representation. (All at once procedures `Store_bgcm_vars_in_matrices(...)` and `Write_matrices_to_bgcm_vars(...)`, as well as such handling selected matrices, e.g. `Write_matrix_to_1l_2d_bgcm_vars(...)`, `Update_matrix_with_2d_bgcm_element_vars(...)` or `Get_matrix_from_store_1l_3d(...)`.) Scientific routines usually work on the matrix representation accessed via the `dsl4jsb_*_mt1L2D`/`mt2L3D`-style DSL macros instead of using plain variable pointers as done for conventional variables. Note: As mentioned above output variable names are composed as `PROC_variable_tile`, with `_variable_` being composed as `_bgcm-prefix(_compartment)_element_` in case of BGCMs.

## The DSL (`dsl4jsb`)

Because everything is an object behind pointers, process code is written against a small macro language, expanded by `scripts/dsl4jsb/dsl4jsb.py` (regex substitution) **before compilation**.

```fortran
dsl4jsb_Use_processes HYDRO_
dsl4jsb_Use_memory(HYDRO_)
dsl4jsb_Def_memory(HYDRO_)
dsl4jsb_Real2D_onChunk :: evapotrans
dsl4jsb_Get_memory(HYDRO_)
dsl4jsb_Get_var2D_onChunk(HYDRO_, evapotrans)
```

expands to `USE` statements, a `TYPE(t_hydro_memory), POINTER :: hydro__mem`, a `SELECT TYPE` on `tile%mem(hydro_)%p`, and `evapotrans => hydro__mem%evapotrans%ptr(ics:ice,iblk)`.

The macros assume that the names `tile`, `options`, `iblk`, `ics`, `ice` exist and are assigned, which can lead to segmentation faults at runtime if this is not the case. In a task routine they therefore need to be declared and set from `options` first. `_tile` variants (`dsl4jsb_Get_memory_tile`, `dsl4jsb_Get_var2D_onChunk_tile`, ...) let you name a different tile, which mainly lcc code needs. `_onDomain` variants point at the whole array and are used in `mo_PROC_init` modules, which run before decomposition on the full domain.

Full macro inventory: grep `dsl4jsb_` in `scripts/dsl4jsb/dsl4jsb.py`. Families are `Use_*`/`Def_*`/`Get_*` for processes, memory, config, lctlib parameters, lcc relocations, plain variables (`var2D/var3D` × `onChunk/onDomain`), and the BGCM matrices (`mt1L2D`, `mt1L3D`, `mt2L2D`, `mt2L3D`).

**Consequences you will hit:**

- Preprocessed files are `<name>.pp-jsb.f90` under `<build-dir>/externals/jsbach/src/...`. **Never edit those** — edits are lost on the next build.
- Compiler and runtime **line numbers are those of the original `.f90`**: `dsl4jsb.py` writes a `#line 1 "<original file>"` directive after its 5-line header, so the offset is compensated.
- If a rebuild seems not to pick up your change, check that the corresponding `.pp-jsb.f90` was regenerated. `make clean` in the ICON root removes the preprocessed files and fixes stale-DSL problems.

## Source layout (`src/`)

Infrastructure:

- `base/` — the framework: model/tile/process/task/memory/variable classes, HSM, usecases, task-queue construction, lcc infrastructure, grids, time, I/O, parallel, lctlib, CQT classes, BGCM store, `mo_jsb_version.f90`
- `interfaces/` — `mo_jsb_interface.f90`, the entry points the host calls: `jsbach_interface`, `jsbach_start_timestep`, `jsbach_finish_timestep`, `jsbach_get_var`, `invalidate_cuda_graphs`
- `drivers/` — `mo_jsbach_model.f90` (standalone driver, invoked from ICON's `src/drivers/icon.f90`) and `mo_jsb4_forcing.f90` (offline forcing)
- `shared/` — cross-process helpers: physical/math constants, `mo_phy_schemes.f90`, `mo_orbit_solar.f90`, isotope and time-averaging utilities, OpenACC helpers

Process directories: see the inventory under "Processes" above.

Other top-level directories: `data/` (lctlib files), `doc/` (FORD guide pages, media, CSS), `scripts/` (DSL, CI checks, preprocessing, restart helpers), `templates/` (new-process skeletons), `LICENSES/`.

## Configuration (namelists)

- `jsb_control_nml` and one `jsb_model_nml` **per model instance** live in the *master* namelist file: technical settings (model id/name/namelist filename, debug and timer levels, `is_standalone`, `restart_jsbach`), parsed in `mo_jsb_control.f90`.
- A second, differently-scoped `jsb_model_nml` in the *model* namelist file carries the scientific configuration (`model_scheme` (`jsbach`/`quincy`), `usecase`, `npft`, `use_tmx`, `quincy_model` (which translates to `qmodel_id`: `QLAND`, `QPLANT`, `QCANOPY`, `QSOIL`), element-inclusion switches, ...), as does `jsb_grid_nml`; grep `mo_jsb_config_class.f90` for the full list.
- Per-process namelists: **`jsb_<proc>_nml`** historically for JSBACH processes (`jsb_hydro_nml`, `jsb_carbon_nml`, `jsb_sse_nml`, ...) and **`lnd_<proc>_nml`** for the QUINCY ones (`lnd_veg_nml`, `lnd_sb_nml`, `lnd_q_assimi_nml`, `lnd_q_pheno_nml`, `lnd_q_syl_nml`, `lnd_q_agr_nml`). Each carries at least `active`, plus `bc_filename`/`ic_filename`.
- In coupled runs ICON-Land output namelists belong to `NAMELIST_atm`, not `NAMELIST_lnd` as in standalone runs.

`__NO_JSBACH__` and `__NO_JSBACH_HD__` guard the whole component and the HD process respectively; other cpp keys in use are `_OPENACC`, `_OPENMP`, `HAVE_CDI_PIO`, `HAVE_F2003`, `HAVE_FC_ATTRIBUTE_CONTIGUOUS`. Keep preprocessor use to technical concerns only, and mind capitalization.

## Build and run (via the host model; here for levante)

From an ICON(-MPIM) clone — ICON-Land is built as part of ICON, there is no build system in ICON-Land. The `config/dkrz/levante.*` wrappers pin their whole toolchain by absolute path, so **no module needs to be loaded to build**:

```bash
mkdir build && cd build
../config/dkrz/levante.intel     # or another config/dkrz/levante.* wrapper
make -j32
```

`--enable-jsbach` is on by default; `--enable-jsbach-hd` defaults to following it. **`--enable-quincy` is deprecated** — QUINCY code is always built when JSBACH is enabled. ICON's `icon.mk.in` runs `scripts/dsl4jsb/dsl4jsb.py` over every `externals/jsbach/src/**/*.f90` before compiling. If running ICON-Land standalone, i.e. without ICON atmosphere and with forcing from a file, more code can be disabled (`--disable-les --disable-upatmo --disable-ocean --disable-coupling --disable-yaxt --disable-atmo --disable-comin --disable-art --disable-waves`).

Running:

- **mkexp** (preferred for real experiments): work in ICON's `run/` directory, copy an example config, then run mkexp *by path* setting `ACCOUNT=<account-id>` where `account-id` is the account to be used for accounting compute time, if required, e.g. on levante. In many cases, `account-id` will be your "effective group ID" (the first listed in the terminal with command `groups`). Then submit the generated `run_start`:

  ```bash
  cd run
  cp examples/JSB-standalone-R2B4.config myexp.config
  ../utils/mkexp/mkexp myexp.config ACCOUNT=<account-id>
  cd ../experiments/myexp/scripts && sbatch myexp.run_start
  ```

  Nothing needs to be on `PATH` and no module is needed. For an out-of-source build, `export MKEXP_PATH=.:<build-dir>/run` so mkexp finds the build's `run/SETUP.config`.

  ICON-Land examples are `run/examples/JSB-standalone-R2B4.config` and `run/examples/IQJ-standalone-R2B4.config`, with short test variants in `run/checksuite.icon-land/`. Their `EXP_TYPE` resolves under `run/mkexp/types/` (`land/icon-land-standalone-R02B04.config`). Frequently used `EXP_OPTIONS` are defined under `run/mkexp/options/`, in particular the options `land/quincy/use_IQJ` and `land/quincy/use_IQJ_standalone` used with `run/examples/IQJ-standalone-R2B4.config`. mkexp also still searches legacy locations for experiment types (`run/standard_experiments/`) and options (`run/standard_options/`).

- **`exp.*` runscript templates** (short tests, and the only form Buildbot accepts): `run/exp.land_jsbach_R2B4_test`, `exp.land_lite_R2B4_test`, `exp.land_quincy_canopy_R2B4_test`, `exp.land_iqj_sitelevel_test`, `exp.atm_nwp_jsbach_test`, ... The actual run script is generated for the current machine with `../make_runscripts <name>` where <name> is the name of the template without "exp.". This will generate the runscript `run/exp.<name>.run. If using SLURM, update the <account-id> in the header of the runscript, and then submit.

## Documentation

The [ICON-Land documentation](https://jsbach.gitlab-pages.dkrz.de/jsbach) is generated by FORD from the `!>`/`!<` inline comments plus the guide pages in `doc/guide/`, and published by the scheduled `pages` job — which runs only on the default branch, from a schedule, never on a merge request.

Build it locally with `./build-doc.sh`: `-f` skips the search index and is much faster, `-c` removes the output again. The result is `doc/html/index.html` (`output_dir` in `project.md`).

**The script needs outbound network access** — it bootstraps `uv` from `astral.sh` into `3rdparty/` and then installs a FORD fork from GitHub, so it cannot run offline or in a sandbox. `graphviz` is needed in addition if `project.md` sets `graph: true` (currently `false`).

Nothing verifies the rendered result automatically: the `Test build doc` CI job is `manual`. So when you can run the build, do, and check that the annotations you added come out as intended — FORD markup mistakes stay invisible until someone renders the pages.

## Conventions

- **License header**: every source file starts with a one-line module summary, then the ICON-Land block using `!>` comment markers — `ICON-Land`, `Copyright (C) 2013-2026, MPI-M, MPI-BGC`, `Contact: icon-model.org`, `Authors: AUTHORS.md`, `See LICENSES/ for license information`, `SPDX-License-Identifier: BSD-3-Clause`. **This differs from ICON's own header** (different years, institutions and comment marker) — copy from a neighbouring file in this repository, not from ICON. See `templates/module_header_with_BSD-3-license.template`; enforced by `scripts/ci_scripts/check_license.sh`.
- **Model name in prose**: always write "ICON-Land" — never "ICON-L" or "ICON-LAND" — in documentation, wiki pages, commit messages and merge-request text. (Code identifiers, paths and namelist names keep their historical "jsbach" spelling.)
- **File hygiene**: ASCII only, no non-UTF-8 bytes, no trailing blanks (`scripts/ci_scripts/check_file_format.sh`). Fortran line limit 132 chars, but keep to ≤110 for readable FORD pages. Two-space indentation, no tabs.
- **Style**: spaces around operators; real literals as `_wp` (`2._wp`); cast integers explicitly (`REAL(iexp, wp)`) rather than mixing types; `==`/`<`/`<=` rather than `.EQ.`/`.LT.`/`.LE.`.
- **Naming**: variables `small_caps`; 3D variables carry a suffix naming their third dimension (`_cl` canopy layer, `_sl` soil layer, `_snl` snow layer); enumerators `ALL_CAPS` (`ASSIMI_`, `MODEL_QUINCY`); type-bound procedures `Capitalised_small_caps`, usually verb-first (`Init_vars`, `Configure_model_processes`).
- **Argument order** for scientific routines: dimensions → lctlib parameters → config variables → `in` → `inout` → `out`, with bgc_material structures after the plain 2D/3D variables in each group.
- **Inline documentation is part of the deliverable.** FORD uses `!>` before a code line and `!<` after it. Every memory variable and every constant needs a long name and `[unit]` at declaration; every constant needs a literature reference (or an explicit "tuned against X") at its definition. See "Documentation" for checking how it renders.
- **Only `update_*` (task Integrate) routines should touch model infrastructure**; where possible keep formulas in infrastructure-free routines so they stay testable and portable.

### GPU (OpenACC)

GPU support is OpenACC. `Add_var` already emits the `!$ACC ENTER DATA COPYIN` / `!$ACC UPDATE DEVICE` for registered variables — anything not created via `Add_var` must be placed on the device explicitly. No array/`WHERE` notation inside accelerated regions: write loops out with `!$ACC PARALLEL DEFAULT(PRESENT)` + `!$ACC LOOP GANG VECTOR COLLAPSE(n)`, cell index innermost. Never issue `!$ACC UPDATE HOST/DEVICE` inside the time loop — it is very expensive. Copy an existing accelerated routine rather than inventing a pattern.

## CI and testing

This repository's own GitLab CI (`.gitlab-ci.yml`) is deliberately thin: a `Check Style` job running the two `scripts/ci_scripts/` checks on merge-request pipelines, plus the `Test build doc` and `pages` jobs covered under "Documentation". **Run both check scripts locally before pushing** — they are the most common avoidable pipeline failure:

```bash
./scripts/ci_scripts/check_file_format.sh
./scripts/ci_scripts/check_license.sh
```

Nothing here builds or runs on its own, so **all functional and regression testing happens in the paired host-model merge request** (normally in `icon-mpim`). There it is covered by *two independent systems* — do not treat either as a substitute for the other:

- **ICON's GitLab CI**, which runs automatically on the merge-request pipeline (stage/job inventory: ICON's own `AGENTS.md`). Its `unit` stage compiles ICON including ICON-Land with every supported toolchain, so a build break in this repository surfaces there without you doing anything. The `system` stage runs experiments and compares checksums, but all of its builds are `manual` — a green pipeline therefore does **not** mean your change was run. Start the system jobs your change can affect.
- **Buildbot**, which is pull-based and independent of GitLab CI; there is no Buildbot for this repository. Trigger it by commenting `runBB(<builder>)` on the host-model MR — e.g. `runBB(levante;atm)` or `runBB(merge2rc)`. New tests must be `exp.*` runscripts (not mkexp configs) placed under ICON's `run/` and registered with an `$addexp` line in `scripts/buildbot_scripts/` (usually `create_list_icon_dev`). If results change, reference data must be regenerated via `create_list_update` in agreement with the gatekeeper; GPU/CPU tolerance testing uses `create_list_tolerance-update`.

Two `basic` jobs bite early and are unrelated to the code: `cla-check` (contributor agreement — see "Commit messages") and `release-notes`, which fails unless `RELEASE_NOTES.md` contains an entry quoting the merge request's own URL, or the MR description contains the literal marker `[no_release_notes]`. Pick one deliberately: user-visible ICON-Land changes belong in the release notes; a pure `externals/jsbach` pointer bump for internal work does not. The job must be retried by hand after you edit the description.

Minimum manual testing before proposing a merge: levante with **both NAG (debugging) and Intel**, using `run/exp.land_quincy_canopy_R2B4_test` and `run/exp.land_jsbach_R2B4_test`, or the `run/checksuite.icon-land/test_jsb_standalone.config`, `run/checksuite.icon-land/test_iqj_canopy_standalone.config` and `amip-aes.config` mkexp setups from ICON's `run/` directory. If a new feature is off by default, demonstrate bit-identical results; if it is on by default, document and explain the changes.

## Repository layout and merge workflow

Three repositories are involved:

- **jsbach** (this one) — ICON-Land development. `land4icon-mpim` is the default and main development branch; `master` is what ICON release candidates point at. Branches: `feature/...`, `bugfix/...` (older commits use `hotfix/`), `project/...` (formerly `application/...`), lowercase, hyphen-separated.
- **icon-mpim** — MPI-M/MPI-BGC ICON development; carries this repo as `externals/jsbach`. Land branches there are named `feature-land-...` (hyphens, no slash).
- **icon** — the shared/official repository with the release candidates and releases.

Work cycle: issue → branch from `land4icon-mpim` → develop → merge `land4icon-mpim` in frequently → test → merge request into `land4icon-mpim` with a verbose description (it becomes the squashed commit message) → review by a gatekeeper.

Backmerging a feature needs **two coordinated merge requests** — one here into `land4icon-mpim`, one in `icon-mpim` into `master` — merged in that order by a gatekeeper. While the pair is under test, the `icon-mpim` MR has `externals/jsbach` pointing at your feature branch; once the jsbach MR is merged, **update the pointer to the resulting commit on `land4icon-mpim`** before the `icon-mpim` MR goes in. Merges here are squashed and the source branch is deleted afterwards, so your feature-branch commit never becomes part of `land4icon-mpim` — leaving it in place would point `icon-mpim` `master` at a commit no branch contains. When updating an ICON release candidate, `externals/jsbach` must point at a commit on jsbach **`master`**, never `land4icon-mpim`; `master` is updated from `land4icon-mpim` first, then the version string `icon_land_version` in `src/base/mo_jsb_version.f90` is bumped (`YYYY.MM` / `YYYY.MM-N`) and tagged `icon-land-YYYY.MM`.

If you change infrastructure, remember the host adapter: a matching change to `src/shared/mo_util_jsbach.f90` on the ICON side must go into the paired `icon-mpim` merge request, and the dependency order between the two MRs must be stated so the dependent one is not merged first.

**Always open merge requests as Draft**, in this repository and in `icon-mpim` alike. Never merge one yourself — a gatekeeper does that, and undrafting is their signal that the MR is ready.

**Cross-link the pair in both directions.** The `icon-mpim` MR names the jsbach one from the start (it is created second, so the number is known). As soon as both exist, go back and append a reference to the `icon-mpim` MR — full cross-project form, `icon/icon-mpim!<number>` — to the end of the jsbach MR description, so a reviewer landing on either one finds the other.

## Known pitfalls

- **Reading another process's memory** is technically allowed and often wrong. Reading a variable from a process that runs further down the hierarchy may hand you the initialisation value if the variable is not in that process's `aggregate_` routine. Writing another process's variables can be silently discarded (overwritten by a later aggregation) or leak to siblings and parent through `INHERIT_` pointers. Do it only in exceptional, documented cases.
- `l_aggregate_all` does **not** mean "aggregate over all child tiles" in the naive sense — see Aggregation above.
- `isteptype=TSTEP_CONSTANT` variables appear to require an output file of their own; mixing them with time-dependent variables in one ICON output file produces corrupted fields.
- Stale preprocessed files: if a rebuild does not reflect your change, see the DSL consequences above.
- Checking out a jsbach branch that does not match ICON's `src/shared/mo_util_jsbach.f90` produces confusing compile errors — verify the adapter implementations match the ICON-Land version before debugging further.
- Wiki drift: the wiki was largely re-synced with the source tree in August 2026; known remaining drift: the usecases page misses `jsbach_forest_age_classes` and the `_tmx` variants. MPI-ESM/ECHAM as a second host is described *on purpose*, marked as historical — older ICON-Land versions can be run with it; do not "fix" those sections away.

## Commit messages

Do not add a `Co-Authored-By: Claude` (or other AI) trailer to commits. This repository feeds into ICON, whose contribution process gates merge requests on contributor-agreement status, and no CLA covering AI co-authorship is in place.

When a commit message references a merge request or an issue, write the **full URL** — `https://gitlab.dkrz.de/jsbach/jsbach/-/merge_requests/450` — not GitLab's `project!number` shorthand. Commit messages are read with `git log` and in mail, outside GitLab, where the shorthand resolves to nothing. In merge request descriptions the shorthand is the better choice and should be used instead.

In the `icon-mpim` repository, **never prepend `[mpim]` to commit messages on a feature branch.** That prefix marks commits that have landed on `icon-mpim`'s `master`, and it is added automatically when the merge request is merged. Writing it yourself on a feature branch produces a doubled or misleading prefix.
