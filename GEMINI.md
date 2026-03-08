# ICON-Land Architectural Overview

This document summarizes the core architecture and process management of the ICON-Land model as of March 2026.

## Core Design Principles
The ICON-Land model uses a modular, task-based architecture governed by a **Hierarchical State Machine (HSM)**. This design decouples physical processes from the spatial and temporal orchestration of the model.

## Key Architectural Components

### 1. Process Management (`t_jsb_process`)
- **Definition:** Every physical process (e.g., `carbon`, `hydrology`, `phenology`, `flake`) is an instance of `t_jsb_process`.
- **Registry:** Processes are registered in the global `t_jsb_model`.
- **Components:** Each process maintains its own configuration (`t_jsb_config`) and memory (`t_jsb_memory`).

### 2. Task-Based Execution (`t_jsb_process_task`)
- **Granularity:** Processes do not execute as monolithic blocks. Instead, they register specific **tasks** (e.g., `update_flake_temp`).
- **Interfaces:** Tasks implement abstract methods:
    - `Integrate`: For local physical updates (top-down or point-wise).
    - `Aggregate`: For bottom-up information flow (e.g., summing sub-tile fluxes).
- **Task Queue:** The model (`t_jsb_model`) maintains a `task_queue` that defines the exact execution sequence of all enabled tasks.

### 3. Orchestration via HSM
- **Execution Flow:** The `Run_tasks` method in `mo_jsb_model_class.f90` iterates through the task queue.
- **Message Passing:** For each task, the model sends an 'INTEGRATE' or 'AGGREGATE' message to the HSM.
- **HSM (`mo_hsm_class.f90`):** Manages the traversal of the **Tile Tree**.

### 4. Spatial Structure: Tiles and Grids
- **Tile Tree:** Represents the hierarchical sub-structure of the land surface (e.g., Grid -> Land -> Soil/Vegetation/Lake).
- **Tile Handler:** Each tile's `Handler` (implemented as `Process_task_on_tile` in `mo_jsb_tile.f90`) receives the task message and calls `task%Do_it` for that specific tile.
- **Memory Access:** Interaction between processes occurs primarily through shared access to the tile's memory list.

## Core File Map (`src/base/`)
- `mo_jsb_process_class.f90`: Definition of the process structure and IDs.
- `mo_jsb_process_factory_core.f90`: Central factory linking processes to their task registration logic.
- `mo_jsb_model_class.f90`: Orchestrates the task queue and HSM triggers.
- `mo_jsb_task_class.f90`: Abstract interface for physical tasks.
- `mo_jsb_tile.f90`: Logic for executing tasks on specific spatial tiles.
- `mo_hsm_class.f90`: Underlying Hierarchical State Machine logic.

## Process Case Study: FLAKE (`src/flake/`)
The FLAKE model implements freshwater lake physics.

### 1. Task Registration (`mo_flake_interface.f90`)
- **Task:** `flake_wbphysics` (Water Body Physics).
- **Methods:**
    - `Integrate`: `update_flake_wbphysics` - Executes the physical time-stepping.
    - `Aggregate`: `aggregate_flake_wbphysics` - Aggregates lake variables to parent tiles.

### 2. Physical Logic (`mo_flake_process.f90` & `mo_flake_sfcflx_process.f90`)
- `calc_radflux`: Solar radiation penetration.
- `calc_momsenlat`: Surface fluxes (momentum, sensible, latent heat).
- `calc_lake_ice`: Ice and snow cover dynamics.
- `calc_lake_temp_stratification`: Vertical temperature profile using a self-similarity approach (mixed-layer and thermocline).

### 3. State Variables (`mo_flake_memory_class.f90`)
Key state variables managed in `t_flake_memory`:
- `h_ice`, `h_snow`: Ice and snow thicknesses.
- `t_wml`, `h_ml`: Mixed-layer temperature and depth.
- `t_bot`: Bottom temperature.
- `c_t`: Shape factor for the thermocline temperature profile.
- `t_lice`: Lake ice temperature.

## Integration of FLAKE with RAD and SEB
As of March 2026, the FLAKE model is being better integrated into the core ICON-Land processes (RAD and SEB).

### 1. Porting of Albedo Computations to RAD Process
- **Goal:** Centralize surface albedo calculations in the `RAD` process to maintain architectural consistency.
- **Changes in `src/radiation/`:**
    - **Constants (`mo_rad_constants.f90`):** Added FLAKE-specific parameters for water and ice albedo (`AlbedoLakeWaterFlake`, `AlbedoLakeWhiteIce`, etc.).
    - **Logic (`mo_rad_process.f90`):** Implemented `calc_alb_lwtr_flake` and `calc_alb_lice_flake` using Mironov and Ritter (2004) empirical formulas.
    - **Interface (`mo_rad_interface.f90`):** The `update_albedo` task now detects if `FLAKE_` is active on a lake tile. It uses `t_lwtr` (lake water temperature) from `SEB` memory as input for these calculations.
- **Changes in `src/flake/`:**
    - **Interface (`mo_flake_interface.f90`):** Removed local albedo calculation logic. The `albedo_lwtr` and `albedo_lice` variables are now retrieved as inputs from the `RAD` process memory.
