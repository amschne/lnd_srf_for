#!/usr/bin/env python

from schneida_tools import coordinate_space
from schneida_tools import gris_dem

def run():
    coordinate_space.test()
    gris_dem.test()
    ais_dem.test()

def main():
    run()

if __name__=='__main__:'
    main()