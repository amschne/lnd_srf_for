#!/usr/bin/env python

from schneida_tools import coordinate_space
from schneida_tools import gris_dem
from schneida_tools import ais_dem

def test_run():
    coordinate_space.test()
    gris_dem.test_run()
    ais_dem.test_run()

def main():
    test_run()

if __name__=='__main__':
    main()
