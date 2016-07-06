"""
Contains functions for different fluid transfer schemes.
Functions:
-pulse2                two-motor ligand pulse scheme
-pulse3                three-motor ligand pulse scheme2
-updatelog             writes to console and external file
-timestamp             returns nicely formatted data+time

Name: fluidMaster
Author: Clayton Little
"""

## We need some modules
import time
import datetime
import drive

## Set up the serial connection
drive.setup(4,9600) # (COM port number,baud rate)
## And the print files
logfilename = "TransferRecord.txt"



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
    updatelog("PULSE EXPERIMENT")
    # CALCULATE USEFUL VALUES
    trans_time1 = (round(vol1/spd1) + 1)  # how long the fluid transfer takes
    trans_time2 = (round(vol2/spd2) + 1)  # how long the fluid transfer takes
    # ESTABLISHING BASELINE
    updatelog("Waiting to establish baseline.")
    time.sleep(wait)
    # REMOVE INITIAL MEDIA
    updatelog("Removing initial media.")
    drive.m2(0,bmedia,100);time.sleep(bmedia/100 + 1)
    # START THE PULSES
    for i in range(stims):
        # ligand pulse
        drive.m1(1,vol1,spd1);time.sleep(trans_time1)
        updatelog("Pulse " + str(i+1) + " START")
        time.sleep(stim_len)
        updatelog("Pulse " + str(i+1) + " END")
        drive.m1(0,vol1,spd1);time.sleep(trans_time1)
        # blank media
        drive.m2(1,vol2,spd2);time.sleep(trans_time2)
        if i != (stims-1):
            updatelog("Blank " + str(i+1) + " START")
            time.sleep(stim_btwn)
            updatelog("Blank " + str(i+1) + " END")
            drive.m2(0,vol2,spd2);time.sleep(trans_time2)
        # if this is the last pulse, don't remove the blank media
        else:
            updatelog("EXPERIMENT COMPLETE\n")


def updatelog(msg):
    """
    Prints logs to the screen, as well as an external text file.
    msg is the message associated with a particular timestamp.
    """
    f = open(logfilename,"a")
    f.write(timestamp() + " " + msg + "\n")
    f.close()
    print(timestamp() + " " + msg)
    
    
def timestamp():
    ts = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
    return ts