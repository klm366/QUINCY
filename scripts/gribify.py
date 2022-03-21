#!/usr/bin/env python3

import re

from pathlib import Path

from numpy import short

memory_class_files = {
    '_tile_' : Path('src/base/mo_jsb_tile.f90'),
    'alcc'   : Path('src/anthropogenic_lcc/mo_alcc_memory_class.f90'),
    'assimi' : Path('src/assimilation/mo_assimi_memory_class.f90'),
    'a2l'    : Path('src/atmland/mo_atm2land_memory_class.f90'),
    'l2a'    : Path('src/atmland/mo_land2atm_memory_class.f90'),
    'carbon' : Path('src/carbon/mo_carbon_memory_class.f90'),
    'disturb': Path('src/disturbance/mo_disturbance_memory_class.f90'),
    'fuel'   : Path('src/disturbance/mo_fuel_memory_class.f90'),
    'hd'     : Path('src/hd/mo_hd_memory_class.f90'),
    'hydro'  : Path('src/hydrology/mo_hydro_memory_class.f90'),
    'pheno'  : Path('src/phenology/mo_pheno_memory_class.f90'),
    'pplcc'  : Path('src/pre_and_post_lcc/mo_pplcc_memory_class.f90'),
    'rad'    : Path('src/radiation/mo_rad_memory_class.f90'),
    'sse'    : Path('src/soil_snow_energy/mo_sse_memory_class.f90'),
    'seb'    : Path('src/srf_energy_bal/mo_seb_memory_class.f90'),
    'tcq'    : Path('src/test_conserved_quantities/mo_tcq_memory_class.f90'),
    'tlcc'   : Path('src/test_lcc/mo_tlcc_memory_class.f90'),
    'turb'   : Path('src/turbulence/mo_turb_memory_class.f90'),
    'veg'    : Path('src/vegetation/mo_veg_memory_class.f90'),
}

tiles = [ 'box', 'land', 'lake', 'glac', 'veg', 'pft01', 'pft02', 'pft03', 'pft04', 'pft05', 'pft06', 'pft07', 'pft08', 'pft09', 'pft10', 'pft11' ]

jsbach_discipline = 2 # land-surface products
jsbach_category_start = 198

force = False
write_files = False
write_shortnames = True
show_mapping = False

basedir = Path('..')

whitespace_pattern = re.compile(r'^\s+(?:!.*)?$|^\s+|\s+(?:!.*)?$')
add_var_pattern = re.compile(r"(?i)mem%add_var\(\s*'(?P<var>[^']+)'")
restart_pattern = re.compile(r'(?i)lrestart\s*=\s*\.false\.')
owner_tile_pattern = re.compile(r'(?i)IF\s*\(TRIM\(mem%owner_tile_name\)')
grib2_pattern = re.compile(r'''(?ix)
    t_grib2\(
        \s*
        (?P<discipline>(?:\+|-)?[0-9]{1,3})
        \s*,\s*
        (?P<category>(?:\+|-)?[0-9]{1,3})
        \s*,\s*
        (?P<code>(?:\+|-)?[0-9]{1,3})
        \s*,\s*
        [^)]+
    \)
    ''')
jsb_add_var_pattern = re.compile(r"(?i)jsb_add_var\s*\(\s*[^,]+,\s*'(?P<var>[^']+)'")

process_vars = {}
for process, memfile in memory_class_files.items():
    with open(basedir / memfile) as f:
        vars = {}
        last_var_end = -1
        contents = list(iter(f))
        for i, line in enumerate(contents):
            if ((match := add_var_pattern.search(line)) is not None
               or (process == '_tile_' and (match := jsb_add_var_pattern.search(line)) is not None)):
                if '!' in line[:match.start()]:
                    continue

                var = match.group('var')

                add_var_call = []
                j = i
                while whitespace_pattern.sub('', contents[j]).endswith('&'):
                    add_var_call.append(contents[j])
                    j += 1
                add_var_call.append(contents[j])

                add_var_call_str = ' '.join(whitespace_pattern.sub('', s).strip('&') for s in add_var_call)

                # This assumes that lrestart is TRUE by default.
                is_restart = process == '_tile_' or restart_pattern.search(add_var_call_str) is None

                if not is_restart:
                    continue

                for k, avline in enumerate(add_var_call):
                    grib2 = grib2_pattern.search(avline)
                    if grib2 is not None:
                        break

                grib2_templ = None

                k += i
                if process != '_tile_' and grib2 is None:
                    k = max(i-10, last_var_end+1)
                    for k in range(max(i-10, last_var_end+1), i):
                        if owner_tile_pattern.search(contents[k]) is not None:
                            grib2_templ = grib2_pattern.search(contents[k+1])

                    k = i
                    while grib2 is None and k > last_var_end + 1:
                        k -= 1
                        grib2 = grib2_pattern.search(contents[k])

                if grib2 is None:
                    print(f'Warning: {var}: no GRIB2 description found. Skipping...')
                    continue

                if abs(i - k) > 5:
                    print(f'Warning: {var}: GRIB2 description at line {k} more than 5 lines away from add_var call (line {i}).')

                tpl_triple = None
                if grib2_templ is not None:
                    tpl_triple = tuple(map(int, grib2_templ.groups()))
                    if any(255 == t for t in tpl_triple):
                        tpl_triple = None
                        print(f'Warning: {var}: Found template triple but one of the entries is 255.')

                grib2_triple = tuple(map(int, grib2.groups()))
                vi = [(i, match.start()), j, (k, grib2.start(), grib2.end()), grib2_triple, tpl_triple]
                vars.setdefault(var, []).append(vi)

                #if process == '_tile_':
                #    print(f'{var}: add_var @ {vi[0]} to line {vi[1]}. GRIB2 info @ {vi[2]}\n    GRIB2 info: {vi[3]}\n    Template: {vi[4]}')
                #    print(contents[vi[2][0]][vi[2][1]:vi[2][2]])

                last_var_end = j
    process_vars[process] = vars

print('Found variables:')
total = 0
for process, vars in process_vars.items():
    total += len(vars)
    print('{}: {}'.format(process, len(vars)))
print('Total: {}'.format(total))


grib_map = {}
jsbach_category = jsbach_category_start
jsbach_parameter = 0
for process, vars in process_vars.items():
    for varname, vis in vars.items():
        vi = vis[0]
        if not force and not any(255 == t for t in vi[3]):
            if vi[3] in grib_map:
                print('Warning: {}/{}: Code triple {!s} already taken by {}! Codes will not be unique.'.format(process, varname, vi[3], '/'.join(grib_map[vi[3]])))
            else:
                grib_map[vi[3]] = (process, varname)
            continue

        if any(vi[4] is not None for vi in vis):
            for vi in vis:
                if vi[4] is not None:
                    break

            vis[0][3] = vi[4]
        else:
            vis[0][3] = (jsbach_discipline, jsbach_category, jsbach_parameter)
            jsbach_parameter += 1
            if jsbach_parameter == 255:
                jsbach_parameter = 0
                jsbach_category += 1

        if vis[0][3] in grib_map:
            print('Warning: {}/{}: Code triple {!s} already taken by {}! Codes will not be unique.'.format(process, varname, vis[0][3], '/'.join(grib_map[vis[0][3]])))
        else:
            grib_map[vis[0][3]] = (process, varname)

if show_mapping:
    for trip, pv in grib_map.items():
        print('{!s}: {}'.format(trip, '/'.join(pv)))

for process, memfile in memory_class_files.items():
    if not process_vars[process]:
        continue

    with open(basedir / memfile) as f:
        contents = list(iter(f))

    for vis in process_vars[process].values():
        dis, cat, par = vis[0][3]

        for vi in vis:
            k, start, stop = vi[2]

            line = contents[k]
            line = ''.join((line[:start], f't_grib2({dis:3d}, {cat:3d}, {par:3d}, grib_bits)', line[stop:]))
            contents[k] = line

    if write_files:
        print(f'Patching {memfile}...')
        with open(basedir / memfile, 'w') as f:
            f.writelines(contents)

shortnamedef = ['# auto-generated by JSBACHs gribify.py']

for trip, (process, varname) in grib_map.items():
    for ti, tile in enumerate(tiles):
        sd = [
            "",
            f"'{process}_{varname}_{tile}' = {{" if process != '_tile_' else f"'{varname}_{tile}' = {{",
            f"    discipline = {trip[0]} ;",
            f"    parameterCategory = {trip[1]} ;",
            f"    parameterNumber = {trip[2]} ;",
        ]

        # box tile is special
        if ti != 0:
            sd += [
                f"    tileAttribute = 0 ;",
                f"    tileIndex = {ti} ;"
            ]
        sd.append(f"}}")

        shortnamedef += sd

if write_shortnames:
    print('Writing shortName.def')
    with open('shortName.def', 'w') as f:
        f.write('\n'.join(shortnamedef))