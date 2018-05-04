% cd to right dir before running compile!
% if gsl is missing:
% sudo port install gsl

function compile
  mex -O -I/opt/local/include/ -L/opt/local/lib/ -lgsl -o Lattmin Lattmin.cpp energy.cpp Graph.cpp %-argcheck
    % To check for alloc: ./matlab -debug -check_malloc
