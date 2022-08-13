#!/usr/bin/env python

""" Create quarter-page figure that shows coarse resolution grid imprinting
"""

import analysis_era5

def run():
    axes = analysis_era5.run_grl()
    

def main():
    run()

if __name__=='__main__':
    main()
