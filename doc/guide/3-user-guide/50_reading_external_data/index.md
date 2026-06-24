---
title: Reading External Boundary-Condition Data
---

## Overview

Many JSBach processes need to prescribe time-varying boundary-condition (BC) fields that
are not computed by the model itself — observed land-use change, external LAI forcing,
prescribed albedo, and similar datasets. These fields are typically stored in per-year
NetCDF files following a common naming convention.

The module `src/shared/mo_jsb_extvar_reader.f90` provides a general-purpose API for
reading, buffering, and interpolating such fields, so that each new use does not have to
re-implement the underlying file I/O and slot-management logic.

---

## Contents

| Page | Description |
|---|---|
| [Using `mo_jsb_extvar_reader`](extvar_reader.md) | API reference, buffer management, time interpolation, and step-by-step implementation recipe |

---

## When to use this module

Use `mo_jsb_extvar_reader` whenever a process needs to:

- read per-year NetCDF files following the `{prefix}{YYYY}.nc` naming pattern,
- maintain a monthly interpolation buffer across year boundaries (monthly mode), or
- read a single 2D target field from a yearly file and apply it as a ramp over the year
  (annual mode).

The `external_lai` phenology scheme (`scheme = 'external_lai'` in `jsb_pheno_nml`) is the
reference implementation; its source in `src/phenology/mo_pheno_init.f90` illustrates
the recommended usage pattern.
