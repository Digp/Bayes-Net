Detailed instructions for installing and running NAST
can be found in the documents:
  NastC2A-installing-guide.pdf
  NastTutorial.v4.pdf

Both of which can be found on the NAST project site:
  simtk.org/home/nast
Or in the docs directory of this package.

Quick installing instructions.

1. Untar the package and move into the directory:
  tar xzf nast-1.0.tar.gz
  cd nast-1.0

2. Build
  python setup.py build

3. Install
3A. global:
  sudo python setup.py install
3B. or local:
  python setup.py install --prefix=$HOME


Quick Troubleshooting (see documentation for more details)

Mac OSX Snow Leopard
You MUST use python version 2.5
One way to do this it to explicitely type "python2.5" 
  instead of "python" for every command.

Local installations
You MUST supply the correct PYTHONPATH to the local installation.
You can do this by adding the following to your .bash_profile file:
Mac OSX:
  PYTHONPATH=$PYTHONPATH:$HOME/lib/python2.x/site-packages
  export PYTHONPATH
Unix:
  PYTHONPATH=$PYTHONPATH:$HOME/lib/python2.x/site-packages
  export PYTHONPATH
Replace x with the subversion of python that you are using (5 or 6)
Make sure to source the new file before continuing:
  source .bash_profile
