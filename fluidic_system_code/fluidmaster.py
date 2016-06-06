"""
Contains the functions for different fluid transfer schemes.

Name: fluidmaster
Author: Clayton Little
"""

## We need some modules
import time
import datetime
import drive

## Set up the serial connection
drive.setup(3,9600) # (COM port number,baud rate)


def pulse(wait,stims,stim_len,stim_btwn,vol1,spd1,vol2,spd2,bmedia):
    """
    Creates a simple pulsing stimulation scheme, where the concentration
    changes instantaneously (a square wave).
    The scheme starts out with blank media, and ends with blank media.
    m1 is ligand, m2 is blank media.
    Usage:
    wait is the seconds before the first pulse starts
    stims is the number of pulses
    stim_len is the pulse length in seconds
    stim_btwn is the time between pulses in seconds
    vol1 is the volume of m1 microliters
    spd1 is the speed of m1 microliters/second
    vol2 is the volume of m2 microliters
    spd2 is the speed of m2 in microliters/second
    bmedia is the initial blank media in microliters
    """
    # CALCULATE USEFUL VALUES
    trans_time1 = (round(vol1/spd1) + 1)  # how long the fluid transfer takes
    trans_time2 = (round(vol2/spd2) + 1)  # how long the fluid transfer takes
    # ESTABLISHING BASELINE
    print(timestamp(),"Waiting to establish baseline.",sep=" ")
    time.sleep(wait)
    # REMOVE INITIAL MEDIA
    print(timestamp(),"Removing initial media.",sep=" ")
    drive.m2(0,bmedia,100);time.sleep(bmedia/100)
    # START THE PULSES
    for i in range(stims):
        # ligand pulse
        drive.m1(1,vol1,spd1);time.sleep(trans_time1)
        print(timestamp(),"Pulse",i+1,"START",sep=" ")
        time.sleep(stim_len)
        print(timestamp(),"Pulse",i+1,"END",sep=" ")
        drive.m1(0,vol1,spd1);time.sleep(trans_time1)
        # blank media
        drive.m2(1,vol2,spd2);time.sleep(trans_time2)
        if i != (stims-1):
            print(timestamp(),"Blank",i+1,"START",sep=" ")
            time.sleep(stim_btwn)
            print(timestamp(),"Blank",i+1,"END",sep=" ")
            drive.m2(0,vol2,spd2);time.sleep(trans_time2)
        # if this is the last pulse, don't remove the blank media
        else:
            print(timestamp(),"DONE",sep=" ")


def timestamp():
    ts = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
    return ts
