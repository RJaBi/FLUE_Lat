project('FLUE', 'c',
  version : '0.1.0',
  license: 'MIT',
  meson_version: '>=0.64.0',
  default_options : ['warning_level=2'],
)

add_languages('fortran', native: false)

meson_path = meson.current_build_dir()

py_mod = import('python')
py = py_mod.find_installation(pure: false)
py_dep = py.dependency()


incdir_numpy = run_command(py,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()


incdir_f2py = run_command(py,
    ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()


# fpm
fflags = '-g -fPIC -funroll-loops -fstack-arrays -Ofast -frepack-arrays -faggressive-function-elimination -fopenmp'
cflags = '-fPIC'
fpm_build_dir = meson_path / 'fpm_build'

run_command(
    'fpm',
    'install',
    '--profile',
    'release',
    '--flag',
    fflags,
    '--c-flag',
    cflags,
    '--prefix',
    fpm_build_dir,
    check : true
)

# f2py
source_path = meson.global_source_root() / '..' / 'c_wrapper' / 'FLUE_c.f90'
FLUE_c_f90 = files(source_path)

FLUE_source = custom_target('FLUE_compiled.c',
  input: FLUE_c_f90,
  output: ['FLUE_compiledmodule.c', 'FLUE_compiled-f2pywrappers2.f90'],
  command: [
    py,
    '-m',
    'numpy.f2py',
    '@INPUT@',
    '-m',
    'FLUE_compiled',
    '--lower'
  ]
)

# Install Python sources
f_sources = ['FLUE/compiled/__init__.py']
py.install_sources(f_sources, subdir:'FLUE/compiled')

main_sources = ['FLUE/__init__.py']
py.install_sources(main_sources, subdir:'FLUE')

gluonAna_sources = [
  'FLUE/gluonAna/__init__.py',
  'FLUE/gluonAna/gpManip.py',
  'FLUE/gluonAna/momentum.py'
]
py.install_sources(gluonAna_sources, subdir:'FLUE/gluonAna')

gluonCalc_sources = [
  'FLUE/gluonCalc/__init__.py',
  'FLUE/gluonCalc/su3matrixops.py',
  'FLUE/gluonCalc/gluonprop.py'
]
py.install_sources(gluonCalc_sources, subdir:'FLUE/gluonCalc')


# Extension fortran module
inc_np = include_directories(incdir_numpy, incdir_f2py, 'fpm_build' / 'include')

py.extension_module('FLUE_compiled',
  [FLUE_c_f90, FLUE_source],
  incdir_f2py / 'fortranobject.c',
  include_directories: [inc_np],
  dependencies : py_dep,
  link_args: ['-L' + 'fpm_build' / 'lib', '-lFLUE'],
  subdir: 'FLUE' / 'compiled',
  install : true,
)
