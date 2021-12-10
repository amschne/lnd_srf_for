#!/usr/bin/env python

from schneida_tools import cruncep
from schneida_tools import wfde5
from schneida_tools import gswp3

def run():
    #cruncep.clean_data()
    gswp3.clean_data()
    #wfde5.clean_data()
    
def main():
    run()

if __name__=='__main__':
    main()
