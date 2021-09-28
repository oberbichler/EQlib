import sys, re, os

try:
    from skbuild import setup
except ImportError:
    print(
        'Please update pip, you need pip 10 or greater,\n'
        ' or you need to install the PEP 518 requirements in pyproject.toml yourself',
        file=sys.stderr,
    )
    raise

VERSION_REGEX = re.compile(
    r'^\s*#\s*define\s+EQLIB_VERSION_([A-Z]+)\s+(.*)$', re.MULTILINE)

this_directory = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join('include/eqlib/common.h')) as f:
    matches = dict(VERSION_REGEX.findall(f.read()))
    eqlib_version = '{MAJOR}.{MINOR}.{PATCH}'.format(**matches)

with open(os.path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

long_description = long_description#[long_description.find('## Introduction'):]

setup(
    name='eqlib',
    version=eqlib_version,
    author='Thomas Oberbichler',
    author_email='thomas.oberbichler@gmail.com',
    description='',
    url='https://github.com/oberbichler/EQlib',
    license='MIT',
    long_description=long_description,
    long_description_content_type='text/markdown',
    cmake_args=[
        '-DCMAKE_INSTALL_LIBDIR=eqlib',
        '-DCMAKE_INSTALL_BINDIR=eqlib',
        '-DCMAKE_INSTALL_INCLUDEDIR=eqlib/include',
        '-DCMAKE_INSTALL_DATAROOTDIR=eqlib/share',
    ],
    packages=['eqlib'],
)