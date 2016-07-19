"""
Contains functions for different fluid transfer schemes.
Functions:
-pulse2                two-motor ligand pulse scheme
-pulse3                three-motor ligand pulse scheme2
-updatelog             writes to console and external files
-timestamp             returns nicely formatted data+time

Name: fluidMaster
Author: Clayton Little
"""

## We need these modules
import time
import datetime
import drive

## Set up the serial connection
drive.setup(4,9600) # (COM port number,baud rate)



def pulse2(wait,bmedia,numpulses,pulse,blank,
           vol1,spd1,vol2,spd2):
    """
    USES TWO SYRINGE PUMPS
    Creates a simple pulsing stimulation scheme, where the concentration
    changes instantaneously (a square wave).
    The scheme starts out with blank media, and ends with blank media.
    m1 is ligand, m2 is blank media.
    Usage:
    wait is the seconds before the first pulse starts
    bmedia is the initial blank media in the dish in microliters
    numpulses is the number of pulses
    pulse is the pulse length in seconds
    blank is the time between pulses in seconds
    vol1 is the volume of m1 microliters
    spd1 is the speed of m1 microliters/second
    vol2 is the volume of m2 microliters
    spd2 is the speed of m2 in microliters/second
    """
    # SET UP LOGS
    global transferLog
    global timeLog
    transferLog = "transferLog.txt"
    timeLog = datetime.datetime.fromtimestamp(time.time()).strftime('%y%m%d_%H%M%S') + "_timeLog.txt"
    updatelog("PULSE EXPERIMENT")
    # CALCULATE USEFUL VALUES
    trans_time1 = (round(vol1/spd1) + 1)  # how long the fluid transfer takes
    trans_time2 = (round(vol2/spd2) + 1)  # how long the fluid transfer takes
    # ESTABLISH BASELINE
    time.sleep(wait)
    # REMOVE INITIAL MEDIA
    updatelog("Removing initial media.")
    drive.m2(0,bmedia,100);time.sleep(bmedia/100 + 1)
    # START THE PULSES
    for i in range(numpulses):
        # ligand pulse
        drive.m1(1,vol1,spd1);time.sleep(trans_time1)
        updatelog("Pulse " + str(i+1) + " START")
        time.sleep(pulse)
        updatelog("Pulse " + str(i+1) + " END")
        drive.m1(0,vol1,spd1);time.sleep(trans_time1)
        # blank media
        drive.m2(1,vol2,spd2);time.sleep(trans_time2)
        if i != (numpulses-1):
            updatelog("Blank " + str(i+1) + " START")
            time.sleep(blank)
            updatelog("Blank " + str(i+1) + " END")
            drive.m2(0,vol2,spd2);time.sleep(trans_time2)
        # if this is the last pulse, don't remove the blank media
        else:
            updatelog("EXPERIMENT COMPLETE\n")
            
            
def pulse3(wait,bmedia,numpulses,pulse,blank,
           bVOL,wVOL,lVOL,
           bSPD,wSPD,lSPD):
    """
    USES THREE SYRINGE PUMPS
    Creates a simple pulsing stimulation scheme, where the concentration
    changes instantaneously (a square wave).
    The scheme starts out with blank media, and ends with blank media.
    m1 is blank media, m2 is waste media, m3 is concentrated ligand
    Usage:
    ...
    wait is the seconds before the first pulse starts
    bmedia is the initial blank media in the dish in microliters
    numpulses is the number of pulses
    pulse is the pulse length in seconds
    blank is the time between pulses in seconds
    ...
    bVOL is the volume of blank media in microliters
    wVOL is the volume of waste media in microliters
    lVOL is the volume of concentrated ligand in microliters
    ...
    bSPD is the speed of blank media in microliters/second
    wSPD is the speed of waste media in microliters/second
    lSPD is the speed of concentrated ligand in microliters/second
    """
    # SET UP LOGS
    global transferLog
    global timeLog
    transferLog = "transferLog.txt"
    timeLog = datetime.datetime.fromtimestamp(time.time()).strftime('%y%m%d_%H%M%S') + "_timeLog.txt"
    updatelog("PULSE EXPERIMENT")
    # CALCULATE USEFUL VALUES
    bTT = (round(bVOL/bSPD) + 0.5)  # Transfer Time of blank media
    wTT = (round(bVOL/bSPD) + 0.5)  # Transfer Time of waste media
    lTT = (round(bVOL/bSPD) + 0.5)  # Transfer Time of concentrated ligand
    # ESTABLISH BASELINE
    time.sleep(wait)
    # REMOVE INITIAL MEDIA
    updatelog("Removing initial media.")
    drive.m2(0,bmedia,wSPD);time.sleep(bmedia/bSPD + 0.5)
    print("SWAG")
    # START THE PULSES
    for i in range(numpulses):
        # ligand pulse
        drive.m3(1,lVOL,lSPD);time.sleep(lTT)
        drive.m1(1,bVOL,bSPD);time.sleep(bTT)
        updatelog("Pulse " + str(i+1) + " START")
        time.sleep(pulse)
        updatelog("Pulse " + str(i+1) + " END")
        drive.m2(0,wVOL,wSPD);time.sleep(wTT)
        # blank media
        drive.m1(1,bVOL,bSPD);time.sleep(bTT)
        if i != (numpulses-1):
            updatelog("Blank " + str(i+1) + " START")
            time.sleep(blank)
            updatelog("Blank " + str(i+1) + " END")
            drive.m2(0,wVOL,wSPD);time.sleep(wTT)
        # if this is the last pulse, don't remove the blank media
        else:
            updatelog("EXPERIMENT COMPLETE\n")


def updatelog(msg):
    """
    Prints logs to the screen, as well as an external text file.
    msg is the message associated with a particular timestamp.
    """
    ## first the readable file
    f = open(transferLog,"a")
    f.write(timestamp() + " " + msg + "\n")
    f.close()
    print(timestamp() + " " + msg)
    ## then the file for MatLab to read
    f = open(timeLog,"a")
    f.write(str(round(time.time())) + ',')
    f.close()
    
    
def timestamp():
    ## Returns full date and time
    ts = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
    return ts