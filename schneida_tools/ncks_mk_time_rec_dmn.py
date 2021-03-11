#!/usr/bin/env python

""" Use subprocess.run() to call "ncks --dfl_lvl=0 --mk_rec_dmn time
    [INPUT_FILE] [OUTPUT_FILE]" to change the time dimension to the
    record dimension.

    The netCDF Operator (NCO) site is http://nco.sourceforge.net.
"""

from os import path
import subprocess

def call_ncks(input_path, input_files, output_path, time_var='time'):
    for i, input_file in enumerate(input_files):
        run(time_var, path.abspath(path.join(input_path, input_file)),
            path.abspath(path.join(output_path, input_file)))

def run(time_var, input_file, output_file):
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
