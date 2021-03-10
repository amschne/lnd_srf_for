#!/usr/bin/env python

""" Use subprocess.run() to call "ncks --dfl_lvl=0 --mk_rec_dmn time
    [INPUT_FILE] [OUTPUT_FILE]" to change the time dimension to the
    record dimension.
"""

import subprocess

def run(input_file, output_file, time_var='time'):
    ncks_command = ['ncks', '--dfl_lvl=0', '--mk_rec_dmn', time_var,
                    input_file, output_file]
    print('%s %s %s %s %s %s' % (ncks_command[0], ncks_command[1],
                                 ncks_command[2], ncks_command[3],
                                 ncks_command[4], ncks_command[5]))
    subprocess.run(ncks_command)

def main():
    run('in.nc', 'out.nc')

if __name__=='__main__':
    main()