"""
    Name: fluidMasterClean
    collects cleaned up version of code used in experiments
    """

## We need these modules
import time
import datetime
import drive
import tkinter as tk

## Set up the serial connection
drive.setup("/dev/cu.usbmodem1411",9600) # (port number,baud rate), e.g. COM11 on Windows or /dev/cu.usbmodem1411 on Mac

def wash(volume, inspd, outspd):
    drive.m2(0,volume,int(outspd))
    print("wash-sleeping")
    time.sleep(15)
    print("wash-done-sleeping")
    drive.m1(1,volume,int(inspd))

def countdown(waittime):
    """
    countdown(waittime) 

    counts down and displays time

    waittime : time in seconds
    """
    waittime = int(waittime)
    for t in range(waittime):
        print (str(waittime-t)+" seconds remaining    ", end="\r")
        time.sleep(1)

# Ibidi 4-well setup
def pulsesetup(numpulses, pulselength, waitlength, numwashes = 5, washinterval = 5, speedup = 1):
    """
    pulsesetup(npulses, pulselength, waitlength)    
    pulsesetup(npulses, pulselength, waitlength, numwashes, testing)

    numpulses     : number of pulses
    pulselength : length of pulses in minutes
    waitlength  : length of wait after pulse in minutes
    numwashes   : number of washes after pulse, default 5
    washinterval: interval between washes in minutes, default 5
    speedup     : speedup factor for testing (e.g 60 to make pulselength 1 last 1 minute)

    example: 6 pulses, 2 hours with 5 hour wait: pulsesetup(6, 120, 300)

    """

    # define motor directions: fluid in or out of dish?
    OUT = 0; 
    IN = 1;

    # volume parameters in ul
    OUT_SPD = 50
    IN_SPD = 100
    LIG_VOL = 60
    MEDIA_VOL = 600

    MINUTE = 60/speedup; # minutes = 60 seconds unless testing is provided
    WASHINTERVAL = 5; # time between washes in in minutes
    DRAINWAIT = 15; # wait after pumping out, in seconds 
                    # (to let pressures equilibrate so fluid is done flowing out)

    # SHOW ESTIMATED TIME AND VOLUME, GIVE 30 SECONDS TO ABORT

    # total experiment volume in microliter
    Vpercycle = (numwashes + 2)*MEDIA_VOL + LIG_VOL;
    Vtot = Vpercycle*numpulses;

    # total experiment time in hours
    Ttot = (pulselength + waitlength)*numpulses/60; 

    print("Total media volume required : " + str(Vtot) + " ul")
    print("Total experiment time : " + str(Ttot) + " hours")
    print("Countdown before start")
    countdown(30)
    print("Starting experiment")

    # START EXPERIMENT

    for i in range(numpulses):

        print ("\nPulse "+str(i+1))

        #initial drain
        print ("Draining")
        drive.m2(OUT, MEDIA_VOL + 50, OUT_SPD)
        time.sleep(DRAINWAIT)
        
        # ligand pulse
        print ("Ligand dilution in")
        drive.m3(IN, LIG_VOL,   OUT_SPD) 
        drive.m1(IN, MEDIA_VOL, IN_SPD)
        countdown(pulselength*MINUTE)

        print ("Pulse done, starting washes")
        #initial quick wash after pulse
        wash(MEDIA_VOL, IN_SPD, OUT_SPD)
        time.sleep(DRAINWAIT)
        #slower washes
        for k in range(numwashes):

            print ("Wash " + str(k+1) + " of " + str(numwashes) + "              ")
            wash(MEDIA_VOL, IN_SPD, OUT_SPD)
            countdown(WASHINTERVAL*MINUTE)

        # wait until next pulse
        print ("Washes done, waiting for next pulse")
        remainingtime = (waitlength - WASHINTERVAL*numwashes)*MINUTE
        countdown(remainingtime)
        
        print ("Pulse cycle " + str(i+1) + " of " + str(numpulses) + " complete")

    print ("Experiment done")

