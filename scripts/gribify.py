#!/usr/bin/env python3

import argparse
import collections
import re

from pathlib import Path

memory_class_files = {
    '_tile_' : Path('src/base/mo_jsb_tile.f90'),
    'alcc'   : Path('src/anthropogenic_lcc/mo_alcc_memory_class.f90'),
    'assimi' : Path('src/assimilation/mo_assimi_memory_class.f90'),
    'a2l'    : Path('src/atmland/mo_atm2land_memory_class.f90'),
    'l2a'    : Path('src/atmland/mo_land2atm_memory_class.f90'),
    'carbon' : Path('src/carbon/mo_carbon_memory_class.f90'),
    'disturb': Path('src/disturbance/mo_disturbance_memory_class.f90'),
    'fuel'   : Path('src/disturbance/mo_fuel_memory_class.f90'),
    'fage'   : Path('src/forest_age/mo_fage_memory_class.f90'),
    'hd'     : Path('src/hd/mo_hd_memory_class.f90'),
    'hydro'  : Path('src/hydrology/mo_hydro_memory_class.f90'),
    'nlcc'   : Path('src/natural_lcc/mo_nlcc_memory_class.f90'),
    'pheno'  : Path('src/phenology/mo_pheno_memory_class.f90'),
    'pplcc'  : Path('src/pre_and_post_lcc/mo_pplcc_memory_class.f90'),
    'rad'    : Path('src/radiation/mo_rad_memory_class.f90'),
    'sse'    : Path('src/soil_snow_energy/mo_sse_memory_class.f90'),
    'seb'    : Path('src/srf_energy_bal/mo_seb_memory_class.f90'),
    'turb'   : Path('src/turbulence/mo_turb_memory_class.f90'),
}

tiles = [ 'box', 'land', 'lake', 'glac', 'veg', 'pft01', 'pft02', 'pft03', 'pft04', 'pft05', 'pft06', 'pft07', 'pft08', 'pft09', 'pft10', 'pft11' ]

jsbach_discipline = 2 # land-surface products
jsbach_category_start = 198

parser = argparse.ArgumentParser(description='GRIB code generation and analysis tool.')
parser._get_formatter = lambda: parser.formatter_class(prog=parser.prog, max_help_position=26)
parser.add_argument('-f', '--force', action='store_true', help='Reassign codes for variables with already valid GRIB codes.')
parser.add_argument('-s', '--show', action='store_true', help='Print a table of code mappings.')
parser.add_argument('-W', '--write-files', action='store_true', help='Write generated codes to the memory files.')
parser.add_argument('-S', '--write-shortnames', action='store_true', help='Write shortName.def')
parser.add_argument('-D', '--dry-run', action='store_true', help='Do not actually write files (for debugging).')

args = parser.parse_args()

force = args.force
write_files = args.write_files
write_shortnames = args.write_shortnames
show_mapping = args.show
dry_run = args.dry_run

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
cf_pattern = re.compile(r'''(?ix)
    t_cf\(
        \s*
        '(?P<standard_name>[^']*)'
        (?:
            \s*,\s*
            '(?P<units>[^']*)'
            \s*,\s*
            '(?P<long_name>[^']*)'
        )?
        [^)]*
    \)
    ''')
jsb_add_var_pattern = re.compile(r"(?i)jsb_add_var\s*\(\s*[^,]+,\s*'(?P<var>[^']+)'")

VariableInfo = collections.namedtuple('VariableInfo', [
    'start', 'endline', 'grib_location', 'grib_triple', 'template_triple', 'cf_info', 'vgrid'])

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

                # Grab the full add_var call.
                add_var_call = []
                j = i
                while whitespace_pattern.sub('', contents[j]).endswith('&'):
                    add_var_call.append(contents[j])
                    j += 1
                add_var_call.append(contents[j])

                add_var_call_str = ' '.join(whitespace_pattern.sub('', s).strip('&') for s in add_var_call)

                vgrid = add_var_call_str.split(',')[4 if process == '_tile_' else 3].strip()

                # This assumes that lrestart is TRUE by default.
                is_restart = process == '_tile_' or restart_pattern.search(add_var_call_str) is None

                if not is_restart:
                    continue

                cf_info = None
                cf_match = cf_pattern.search(add_var_call_str)
                if cf_match is not None:
                    cf_info = cf_match.groupdict()

                # Find the GRIB2 metadata inside the add_var call, remembering its line offset in `k`.
                for k, avline in enumerate(add_var_call):
                    grib2 = grib2_pattern.search(avline)
                    if grib2 is not None:
                        break

                k += i

                grib2_templ = None

                # The GRIB2 metadata might have been defined a few lines before the call, usually in an IF construct.
                # Look at up to 10 lines before the call but make sure that we don't accidentally grab the previous
                # variable's metadata. Also match `owner_tile_pattern` to try and find a non-dummy GRIB triple.
                if process != '_tile_' and grib2 is None:
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
                vi = VariableInfo((i, match.start()), j, (k, grib2.start(), grib2.end()), grib2_triple, tpl_triple, cf_info, vgrid)
                vars.setdefault(var, []).append(vi)

                last_var_end = j
    process_vars[process] = vars

print('Found variables:')
total = 0
for process, vars in process_vars.items():
    total += len(vars)
    print('{}: {}'.format(process, len(vars)))
print('Total: {}'.format(total))


grib_map = {}

for process, vars in process_vars.items():
    for varname, vis in vars.items():
        if force and vis[0].grib_triple >= (2, jsbach_category_start, 0):
            # Var has a generated triple and --force is in effect.
            # Reset triple so it will be filled below.
            vis[0] = vis[0]._replace(grib_triple=(255, 255, 255))

        # If no GRIB triplet but a template triple is present, take it.
        if any(255 == t for t in vis[0].grib_triple) and any(vi.template_triple is not None for vi in vis):
            for vi in vis:
                if vi.template_triple is not None:
                    break

            vis[0] = vis[0]._replace(grib_triple=vi.template_triple)

        # Check if variable already has a non-dummy GRIB triple (either from template or from field definition).
        if not any(255 == t for t in vis[0].grib_triple):
            vi = vis[0]

            if vi.grib_triple in grib_map and grib_map[vi.grib_triple] != (process, varname):
                print('Warning: {}/{}: Code triple {!s} already taken by {}! Codes will not be unique.'.format(
                    process, varname, vi.grib_triple, '/'.join(grib_map[vi.grib_triple])))

            grib_map[vi.grib_triple] = (process, varname)


# At this point, grib_map has been prepopulated with all known non-dummy triplets.
# Populate everything else with unused triplets from our namespace.
jsbach_category = jsbach_category_start
jsbach_parameter = 0
for process, vars in process_vars.items():
    for varname, vis in vars.items():
        # Check if variable has a dummy GRIB triple.
        if any(255 == t for t in vis[0].grib_triple):
            while (jsbach_discipline, jsbach_category, jsbach_parameter) in grib_map:
                jsbach_parameter += 1
                if jsbach_parameter == 255:
                    jsbach_parameter = 0
                    jsbach_category += 1

            vis[0] = vis[0]._replace(grib_triple=(jsbach_discipline, jsbach_category, jsbach_parameter))
            grib_map[vis[0].grib_triple] = (process, varname)

if show_mapping:
    print()
    print('|   GRIB    |  proc   |             short              | CF standard | CF long | vgrid |')
    print('|-----------|---------|--------------------------------|-------------|---------|-------|')
    for trip, (process, varname) in sorted(grib_map.items()):
        vi = process_vars[process][varname][0]
        print('| {:<9s} | {:7s} | {:30s} | {} | {} | {} |'.format(
            '{:1d}.{:3d}.{:3d}'.format(*trip),
            process,
            varname,
            vi.cf_info['standard_name'] if vi.cf_info['standard_name'].casefold() != varname.casefold() else '<-',
            vi.cf_info['long_name'],
            vi.vgrid))

if write_files:
    for process, memfile in memory_class_files.items():
        if not process_vars[process]:
            continue

        with open(basedir / memfile) as f:
            contents = list(iter(f))

        for vis in process_vars[process].values():
            dis, cat, par = vis[0].grib_triple

            for vi in vis:
                k, start, stop = vi.grib_location

                line = contents[k]
                line = ''.join((line[:start], f't_grib2({dis:3d}, {cat:3d}, {par:3d}, grib_bits)', line[stop:]))
                contents[k] = line

        print(f'Patching {memfile}...' + (' (DRY RUN)' if dry_run else ''))
        if not dry_run:
            with open(basedir / memfile, 'w') as f:
                f.writelines(contents)


if write_shortnames:
    shortnamedef = ['# auto-generated by JSBACHs gribify.py']

    for trip, (process, varname) in sorted(grib_map.items()):
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

    shortnamedef.append('')
    print('Writing shortName.def...' + (' (DRY RUN)' if dry_run else ''))
    if not dry_run:
        with open('shortName.def', 'w') as f:
            f.write('\n'.join(shortnamedef))
