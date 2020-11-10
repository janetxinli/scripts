#!/usr/bin/env python3
"""
Read time in 'hh:mm:ss' format and convert to a given unit.
For help use: python3 read_time.py -h
"""

import sys
import argparse
import datetime

def convert_unit(time_string, desired_unit):
    """
    Convert a time string representing an amount of time and ending in a unit
    to another unit.
    """
    from_unit = time_string[-1]
    time = float(time_string[:-1])
    if from_unit == "s":
        if desired_unit == "s":
            return time
        elif desired_unit == "m":
            return time / 60
        else:
            return time / 3600
    elif from_unit == "m":
        if desired_unit == "s":
            return time * 60
        elif desired_unit == "m":
            return time
        else:
            return time / 60
    else:
        if desired_unit == "s":
            return time * 3600
        elif desired_unit == "m":
            return time * 60
        else:
            return time

def convert_clock(time_string, desired_unit):
    """
    Given a time string in the format 'hh:mm:ss', converts the time to a desired unit
    and prints the result.
    """
    time_split = time_string.split(":")
    num_pos = len(time_split)
    i = 0
    hours = 0
    mins = 0
    secs = 0
    if (num_pos - i) > 2:  # Hours
        hours += float(time_split[i].strip("h\n"))
        i += 1
    if (num_pos - i) > 1:  # Minutes
        mins += float(time_split[i].strip("m\n"))
        i += 1
    if (num_pos - i) > 0:  # Seconds
        secs += float(time_split[i].strip("s\n"))
        i += 1
    
    if desired_unit == "s":
        time = secs + (60 * mins) + (3600 * hours)
    elif desired_unit == "m":
        time = (secs / 60) + mins + (60 * hours)
    else:
        time = (secs / 3600) + (mins / 60) + hours
    
    print(time)

def elapsed_time(start_time, end_time, desired_unit):
    """
    Calculate the amount of time (in a specified unit) elapsed between two
    given time points.
    """
    FMT = "%Y-%m-%d %H:%M:%S"
    start_time = start_time.split(".")[0]
    end_time = end_time.split(".")[0]
    start_time_obj = datetime.datetime.strptime(start_time, FMT)
    end_time_obj = datetime.datetime.strptime(end_time, FMT)
    time_delta = end_time_obj - start_time_obj
    elapsed_seconds = time_delta.total_seconds()
    if desired_unit == "s":
        print(elapsed_seconds)
    elif desired_unit == "m":
        print(elapsed_seconds / 60)
    else:
        print(elapsed_seconds / 3600)

def main():
    parser = argparse.ArgumentParser(description="Read and convert time strings")
    parser.add_argument("-m", "--mode", \
        metavar="mode for time calculation", dest="mode", type=str, choices=["convert_clock", "elapsed", "convert_unit"], \
        help="Mode of calculation (choose between 'convert_clock', 'convert_unit' and 'elapsed'")
    parser.add_argument("-i", "--input", \
        metavar="input time (hh:mm:ss) or - for stdin", dest="input", type=str, \
        help="Input time to be converted from")
    parser.add_argument("-u", "--unit", \
        metavar="unit", dest="to_unit", default="m", type=str, choices=["s", "m", "h"], \
        help="Desired unit for time to be converted to ['s', 'm', 'h']")
    parser.add_argument("-s", "--start", \
        metavar="start", dest="start_time", type=str, \
        help="Start time for elapsed time calculation")
    parser.add_argument("-e", "--end", \
        metavar="end", dest="end_time", type=str, \
        help="End time for elapsed time calculation")
    
    arguments = parser.parse_args()
    input_time = arguments.input
    output_units = arguments.to_unit
    mode = arguments.mode
    start = arguments.start_time
    end = arguments.end_time

    if input_time == "-":
        input_time = sys.stdin.readline()

    if mode == "convert_clock":
        if input_time != None:
            convert_clock(input_time, output_units)
        else:
            raise InputError("Must include an input time in ##h:##m:##s format and \
                optionally a desired unit for time to be converted to.")
    
    elif mode == "elapsed":
        if start != None or end != None:
            elapsed_time(start, end, output_units)
        else:
            raise InputError("Must include a start and end time in datetime format \
                for elapsed time to be calculated, as well as desired units.")
    
    elif mode == "convert_unit":
        if input_time[-1] in ("HhMmSs"):
            print(str(convert_unit(input_time, output_units)))
        else:
            raise InputError("Must include unit at end of time string (e.g. 350m)")


if __name__ == "__main__":
    main()
