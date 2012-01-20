from distutils.core import setup

setup(name='fasta',
      version='1.0',
      author='Mike DeLaurentis',
      author_email='delaurentis@gmail.com',
      url='https://github.com/mdelaurentis/fasta-parser',
      scripts=['parse_fasta'],
      py_modules=['fasta'],
      packages=['fasta'],
      package_data={'fasta' : ['test.fna']})
