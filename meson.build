project('plasMD', 'cpp')

blasdep = dependency('openblas')
lapackdep = dependency('lapack')
fftwdep = dependency('fftw3')
ompdep = dependency('openmp')

global_deps = [blasdep, lapackdep, fftwdep, ompdep]

subdir('src')
subdir('tests')

