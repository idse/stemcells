"""
    Contains functions for different fluid transfer schemes.
    Functions:
    -pulse2                two-motor ligand pulse scheme
    -pulse3wash            three-motor ligand pulse scheme w/ media replacement washes
    -pulse3flow            three-motor ligand pulse scheme w/ continuous flow washes
    -updatelog             writes to console and external log files
    -timestamp             returns nicely formatted data+time
    
    Name: fluidMaster
    Author: Clayton Little
    """

## We need these modules
import time
import datetime
import drive
import tkinter as tk

# Set up the serial connection
drive.setup("/dev/cu.usbmodem1411",9600) # (port number,baud rate)


def wash(volume, inspd, outspd):
    drive.m2(0,volume,int(outspd))
    print("wash-sleeping")
    time.sleep(30)
    print("wash-done-sleeping")
    drive.m1(1,volume,int(inspd))

def swish(times, volume, spd):
    """
        uses the waste pump to successively remove and add liquid three times in a row.
    """
    drive.m2(0,volume, spd)
    drive.m2(1,volume, spd)
    drive.m2(0,volume, spd)
    drive.m2(1,volume, spd)
    drive.m2(0,volume, spd)
    drive.m2(1,volume, spd)

def ramp_step_setup(duration, initial_volume, steps, final_concentration, testing=60):
    """
    We have not been able to get this code to work because of the problem with mixing
    Swishing did not fix the issue because air and bubbles would get trapped in the swish
    line and the water used to swish would sometimes get mixed in with the cells
    """

    MEDIA_OUT_SPD = 500
    LIGAND_OUT_SPD = 20
    MEDIA_VOL = 80
    SWISH_SPD = 100

    #---------------------
    total_volume = initial_volume
    past_lig = 0
    for i in range(steps):
        print ("Step "+str(i+1))
        partial_concentration = (final_concentration/steps)*(i+1)
        print(partial_concentration)
        lig_vol = (partial_concentration*(total_volume+MEDIA_VOL)-past_lig)/(1-partial_concentration)
        past_lig+=lig_vol
        total_volume += (lig_vol+MEDIA_VOL)
        print ("lig_vol " + str(lig_vol))
        print ("total_vol " + str(total_volume))
        print ("================")
        #pump ligand
        drive.m3(1,int(lig_vol), LIGAND_OUT_SPD)
        time.sleep(5)
        # #pump media
        drive.m1(1,MEDIA_VOL, MEDIA_OUT_SPD)
        time.sleep(2)
        #swish
        drive.m2(0,int(total_volume),SWISH_SPD)
        time.sleep(5)
        #800 is the estimated tube volume so 600 to be safe
        drive.m2(1,int(int(total_volume)+600/steps), SWISH_SPD)
        #sleep
        print("Sleeping for " + str(int(duration/steps)) + " seconds")
        time.sleep(int(duration/steps))
    print ("p= "+str(past_lig))


def ramp_wash_setup(duration, initial_volume ,steps, final_concentration):  

    MEDIA_IN_SPD = 200
    LIGAND_OUT_SPD = 20
    MEDIA_OUT_SPD = 50
    OUT_OFFSET = 40
    #Assumed that initial volume is constant total volume
    #-------------------
    for i in range(steps):
        print ("Step"+str(i))
        #Take off old media
        drive.m2(0, int(initial_volume)+OUT_OFFSET, MEDIA_OUT_SPD)
        #delay for continued suction
        time.sleep(int(int(initial_volume)*1/float(MEDIA_OUT_SPD)))
        print ("wait for old media")
        time.sleep(10)
        #Calculate new concentration
        partial_concentration = (i+1)*float(final_concentration)/(steps)
        in_lig = int(partial_concentration*initial_volume)
        print (in_lig)
        in_med = initial_volume-in_lig
        #Pump out ligand
        drive.m3(1, in_lig, LIGAND_OUT_SPD)
        time.sleep(5)
        #Pump out media
        drive.m1(1, in_med, MEDIA_IN_SPD)
        print ("------End of step------")
        print ("Wating "+str(int(duration/float(steps)))+"s for the next step")
        time.sleep(int(duration/float(steps)))




# def ramp_setup(iterations, duration, volume, testing=60):
    """
    This function did not work because there was no mixing involved
    """
#     OUT_SPD = 50
#     IN_SPD = 500
#     MEDIA_VOL = 500
#     #The amount of liquid used to flush each mini-pulse (40 minimum)
#     RAMP_WASH_VOL = 40
    
#     #------------------------
#     for i in range(iterations):
#         #initial drain
#         drive.m2(0,MEDIA_VOL+100, OUT_SPD)
#         time.sleep(30)
#         drive.m1(1,MEDIA_VOL-volume*2, OUT_SPD)
#         time.sleep(30)
#         for k in range(duration):
#             inc_vol = volume/float(duration)
#             drive.m3(1, int(inc_vol), int(inc_vol/2.0))
#             drive.m1(1, 40, int(inc_vol/2.0))
#             time.sleep(1)
#         #washes
#         for k in range(6):
#             print ("\nSwish/Wash" + str(k))
#             # swish(3,500,250)
#             wash(MEDIA_VOL,IN_SPD, OUT_SPD)
#             for t in range(testing*10):
#                 print (str(testing*10-t)+" Sec Remaining", end = "\r")
#                 time.sleep(1)
#         for k in range (testing*testing*5):
#             print (str(testing*testing*5-t)+" Sec Remaining", end = "\r")
#             time.sleep(1)
#         print ("PULSE " + str(pulses) + " COMPLETE")

#4 well setup
def alt_pulsesetup(pulses, testing = 60):
    
    
    
    OUT_SPD = 50
    IN_SPD = 500
    LIG_VOL = 50
    MEDIA_VOL = 500
    
    
    #Do not edit below this line
    #----------------------------
    for i in range(pulses):
        #initial drain
        drive.m2(0,MEDIA_VOL+100, OUT_SPD)
        time.sleep(30)
        
        print ("\nPulse "+str(i))
        drive.m3(1,LIG_VOL,OUT_SPD)
        time.sleep(5)
        drive.m1(1,MEDIA_VOL,IN_SPD)
        for t in range(60*testing):
            print (str(60*testing-t)+" Sec Remaining", end="\r")
            time.sleep(1)
        #initial wash for double wash
        # swish(3,250,250)
        wash(MEDIA_VOL,IN_SPD, OUT_SPD)
        time.sleep(30)
        #washes
        for k in range(6):
            print ("\nSwish/Wash" + str(k))
            # swish(3,500,250)
            wash(MEDIA_VOL,IN_SPD, OUT_SPD)
            for t in range(testing*10):
                print (str(testing*10-t)+" Sec Remaining", end = "\r")
                time.sleep(1)
        for k in range (testing*testing*5):
            print (str(testing*testing*5-t)+" Sec Remaining", end = "\r")
            time.sleep(1)
        print ("PULSE " + str(pulses) + " COMPLETE")
    print ("EXPERIMENT ENDED")



#2000 mcl setup
def pulsesetup(pulses):
    for i in range(pulses):
        drive.m2(0,2000,500)
        #pulse
        time.sleep(5)
        print ("\nPulse "+str(i))
        drive.m3(1,80,80)
        drive.m1(1,2000, 50)
        for t in range(60*60):
            print (str(60*60-t)+" Sec Remaining", end="\r")
            time.sleep(1)
        #initial wash for double wash
        swish(3,1000,1000)
        wash(2000,50)
        #washes
        for k in range(6):
            print ("\nSwish/Wash" + str(k))
            swish(3,1000,1000)
            wash(2000,50)
            for t in range(60*60):
                print (str(60*60-t)+" Sec Remaining", end = "\r")
                time.sleep(1)
        
    print ("EXPERIMENT ENDED")


def pulse2(wait,bmedia,numpulses,pulse,blank,
           vol1,spd1,vol2,spd2):
    """
        USES TWO SYRINGE PUMPS
        Creates a simple pulsing stimulation scheme, where the concentration
        changes instantaneously (a square wave).
        The scheme starts out with blank media, and ends with blank media.
        m1 is ligand, m2 is blank media.
        Usage:
        ...
        wait is the seconds before the first pulse starts
        bmedia is the initial blank media in the dish in microliters
        numpulses is the number of pulses
        pulse is the pulse length in seconds
        blank is the time between pulses in seconds
        ...
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


def pulse3wash(wait,pulse,blank,wash,numpulses,numwashes,bmedia,
               bVOL,wVOL,lVOL,
               bSPD,wSPD,lSPD):
    """
        USES THREE SYRINGE PUMPS WITH MEDIA REPLACEMENT WASH
        Creates a simple pulsing stimulation scheme, where the concentration
        changes instantaneously (a square wave).
        The scheme starts out with blank media, and ends with blank media.
        m1 is blank media, m2 is waste media, m3 is concentrated ligand
        Usage:
        ...
        wait is the seconds before the first pulse starts
        pulse is the pulse length in seconds
        blank is the time between pulses in seconds
        wash is the time between washes in seconds
        numpulses is the number of pulses
        numwashes is the number of washes
        bmedia is the initial blank media in the dish in microliters
        ...
        bVOL is the volume of blank media in microliters
        wVOL is the volume of waste media in microliters
        lVOL is the volume of concentrated ligand in microliters
        ...
        bSPD is the speed of blank media in microliters/second
        wSPD is the speed of waste media in microliters/second
        lSPD is the speed of concentrated ligand in microliters/second
        """
    params = [wait,pulse,blank,wash,numpulses,numwashes,bmedia,bVOL,wVOL,lVOL,bSPD,wSPD,lSPD]
    ### ADJUST PARAMETERS
    # we change some of the times so that the total intervals are accurate.
    # for example, if the "wash" argument is ten minutes, that does not account
    # for the time of the actual washing. We substract that time here.
    bTT = (bVOL/bSPD + 0.5)  # Transfer Time of blank media
    wTT = (wVOL/wSPD + 0.5)  # Transfer Time of waste media
    lTT = (lVOL/lSPD + 0.5)  # Transfer Time of concentrated ligand
    pulse = pulse - lTT - bTT
    blank = blank - wash*(numwashes+1) - wTT
    wash = wash - wTT - bTT
    
    # CHECK PARAMETERS
    if (pulse < 0) | (blank < 0) | (wash < 0):
        print("One or more time intervals are too short.")
        print('Pulse:',pulse,' Blank:',blank,' Wash:',wash)
        return
    
    # PROJECT FLUID USAGE
    print('Projected fluid usage in mcl---',
          'Blank:',bVOL*(numpulses + (numpulses-1)*(numwashes+1) + 1),
          ' Waste:',wVOL*(numpulses + (numpulses-1)*(numwashes+1) + 1),
          ' Ligand:',lVOL*numpulses)

    # SET UP LOGS
    global transferLog
    global timeLog
    transferLog = "transferLog.txt"
    timeLog = datetime.datetime.fromtimestamp(time.time()).strftime('%y%m%d_%H%M%S') + "_timeLog.txt"
    updatelog("PULSE EXPERIMENT with params: " + str(params).strip('[]'))
    app.destroy() #close GUI
    
    # ESTABLISH BASELINE
    time.sleep(wait - (bmedia/bSPD + 0.5))
    
    # REMOVE INITIAL MEDIA
    updatelog("Removing initial media.")
    drive.m2(0,bmedia,wSPD);time.sleep(bmedia/wSPD + 0.5)
    
    # START THE PULSES
    for p in range(numpulses):
        # ligand pulse
        drive.m3(1,lVOL,lSPD);time.sleep(lTT)
        drive.m1(1,bVOL,bSPD);time.sleep(bTT)
        updatelog("Pulse " + str(p+1) + " START")
        time.sleep(pulse)
        updatelog("Pulse " + str(p+1) + " END")
        # not last pulse
        if p != (numpulses-1):
            # end of ligand, now washing
            for j in range(numwashes+1): # initial removing of ligand is not a "wash"
                drive.m2(0,wVOL,wSPD);time.sleep(wTT)
                drive.m1(1,bVOL,bSPD);time.sleep(bTT)
                time.sleep(wash);
                if j == 0:
                    updatelog("Blank " + str(p+1) + " START")
            time.sleep(blank)
            updatelog("Blank " + str(p+1) + " END")
            drive.m2(0,wVOL,wSPD);time.sleep(wTT)
        # last pulse, don't remove the blank media
        else:
            drive.m2(0,wVOL,wSPD);time.sleep(wTT)
            drive.m1(1,bVOL,bSPD);time.sleep(bTT)
            updatelog("EXPERIMENT COMPLETE\n")


def pulse3flow(wait,pulse,blank,wash,numpulses,numwashes,bmedia,
               washVOL,washSPD,
               bVOL,wVOL,lVOL,
               bSPD,wSPD,lSPD):
    """
        USES THREE SYRINGE PUMPS WITH CONTINUOUS FLOW WASH
        Creates a simple pulsing stimulation scheme, where the concentration
        changes instantaneously (a square wave).
        The scheme starts out with blank media, and ends with blank media.
        m1 is blank media, m2 is waste media, m3 is concentrated ligand
        Usage:
        ...
        wait is the seconds before the first pulse starts
        pulse is the pulse length in seconds
        blank is the time between pulses in seconds
        wash is the time between washes in seconds
        numpulses is the number of pulses
        numwashes is the number of washes
        bmedia is the initial blank media in the dish in microliters
        ...
        washVOL is the volume transfer of a single wash in microliters
        washSPD is the transfer speed of a wash in microliters/s
        ...
        bVOL is the volume of blank media in microliters
        wVOL is the volume of waste media in microliters
        lVOL is the volume of concentrated ligand in microliters
        ...
        bSPD is the speed of blank media in microliters/second
        wSPD is the speed of waste media in microliters/second
        lSPD is the speed of concentrated ligand in microliters/second
        """
    params = [wait,pulse,blank,wash,numpulses,numwashes,bmedia,washVOL,washSPD,bVOL,wVOL,lVOL,bSPD,wSPD,lSPD]
    ### ADJUST PARAMETERS
    # we change some of the times so that the total intervals are accurate.
    # for example, if the "wash" argument is ten minutes, that does not account
    # for the time of the actual washing. We substract that time here.
    bTT = (bVOL/bSPD + 0.5)  # Transfer Time of blank media
    wTT = (wVOL/wSPD + 0.5)  # Transfer Time of waste media
    lTT = (lVOL/lSPD + 0.5)  # Transfer Time of concentrated ligand
    washTT = (washVOL/washSPD + 0.5) #Transfer Time of single wash
    pulse = pulse - lTT - bTT
    blank = blank - wash*(numwashes) - wTT
    wash = wash - washTT
    
    # CHECK PARAMETERS
    if (pulse < 0) | (blank < 0) | (wash < 0):
        print("One or more time intervals are too short.")
        print('Pulse:',pulse,' Blank:',blank,' Wash:',wash)
        return
    
    # PROJECT FLUID USAGE
    print('Projected fluid usage in mcl---',
          'Blank:',bVOL*2*numpulses + washVOL*(numpulses-1)*numwashes,
          ' Waste:',wVOL*(2*numpulses-1) + bmedia + washVOL*(numpulses-1)*numwashes,
          ' Ligand:',lVOL*numpulses)

    # SET UP LOGS
    global transferLog
    global timeLog
    transferLog = "transferLog.txt"
    timeLog = datetime.datetime.fromtimestamp(time.time()).strftime('%y%m%d_%H%M%S') + "_timeLog.txt"
    updatelog("PULSE EXPERIMENT with params: " + str(params).strip('[]'))
    app.destroy() #close GUI
    
    # ESTABLISH BASELINE
    time.sleep(wait - (bmedia/bSPD + 0.5))
    
    # REMOVE INITIAL MEDIA
    updatelog("Removing initial media.")
    drive.m2(0,bmedia,wSPD);time.sleep(bmedia/wSPD + 0.5)
    
    # START THE PULSES
    for p in range(numpulses):
        # ligand pulse
        drive.m3(1,lVOL,lSPD);time.sleep(lTT)
        drive.m1(1,bVOL,bSPD);time.sleep(bTT)
        updatelog("Pulse " + str(p+1) + " START")
        time.sleep(pulse)
        updatelog("Pulse " + str(p+1) + " END")
        if p != (numpulses-1): # not last pulse
            # end of ligand
            drive.m2(0,wVOL,wSPD);time.sleep(wTT)
            drive.m1(1,bVOL,bSPD);time.sleep(bTT)
            for j in range(numwashes):
                time.sleep(wash);
                drive.m4(1,washVOL,washSPD);time.sleep(washTT)
                if j == 0:
                    updatelog("Blank " + str(p+1) + " START")
            time.sleep(blank)
            updatelog("Blank " + str(p+1) + " END")
            drive.m2(0,wVOL,wSPD);time.sleep(wTT)
        else: # last pulse, don't remove the blank media
            drive.m2(0,wVOL,wSPD);time.sleep(wTT)
            drive.m1(1,bVOL,bSPD);time.sleep(bTT)
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










## Set up GUI

LARGE_FONT = ("Verdana", 15)
SMALL_FONT = ("Verdana", 11)

class SyringeGUI(tk.Tk):
    
    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)
        self.resizable(True,False)
        
        container = tk.Frame(self)
        container.pack(side="top", fill="both", expand = True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)
        
        self.frames = {}
        for F in (StartPage, pulse2Page, pulse3washPage, pulse3flowPage):
            frame = F(container, self)
            self.frames[F] = frame
            frame.grid(row=0, column=0, sticky="nsew")
        
        self.show_frame(StartPage)
    
    def show_frame(self, cont):
        frame = self.frames[cont]
        frame.tkraise()
    
    def createentries(self, child, fields):
        entries = []
        for field in fields:
            fi = fields.index(field)
            child.lab = tk.Label(child,text=field,anchor="e",font=SMALL_FONT)
            child.ent = tk.Entry(child,font=SMALL_FONT)
            child.lab.grid(row=(3+fi),column=1,sticky='EW')
            child.ent.grid(row=(3+fi),column=2,sticky='EW')
            entries.append(child.ent)
        return entries
    
    def getparams(self, entries):
        params = []
        for entry in entries:
            params.append(int(entry.get()))
        return params

class StartPage(tk.Frame):
    
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="SYRINGE SYSTEM INTERFACE",font=LARGE_FONT)
        label.pack(pady=10,padx=10)
        label = tk.Label(self, text="select your pulse configuration:",font=SMALL_FONT)
        label.pack()
        label = tk.Label(self, text=" ")
        label.pack()
        
        button1 = tk.Button(self, text="TWO SYRINGE PUMPS", font=SMALL_FONT, command=lambda: controller.show_frame(pulse2Page))
        button1.pack()           
        button2 = tk.Button(self, text="THREE SYRINGE PUMPS W/ EXCHANGE WASH", font=SMALL_FONT, command=lambda: controller.show_frame(pulse3washPage))
        button2.pack()            
        button3 = tk.Button(self, text="THREE SYRINGE PUMPS W/ FLOW WASH", font=SMALL_FONT, command=lambda: controller.show_frame(pulse3flowPage))
        button3.pack()


class pulse2Page(tk.Frame):
    
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="SYRINGE SYSTEM INTERFACE",font=LARGE_FONT)
        label.grid(row=0,column=0,columnspan=4,pady=10,padx=10)
        label = tk.Label(self, text="Two-Syringe Control",font=SMALL_FONT)
        label.grid(row=1,column=0,columnspan=4)
        #self.grid_columnconfigure(1, weight=1)
        self.grid_columnconfigure(2, weight=1)
        
        button1 = tk.Button(self, text="Back", font=SMALL_FONT,command=lambda: controller.show_frame(StartPage))
        button1.grid(row=2,column=0,pady=0,padx=10,)

        ## The fields:
        fields = ('Initial wait (s)','Initial blank media (mcl)','Number of pulses','Pulse length (s)','Blank length (s)','Syringe 1 transfer volume (mcl)','Syringe 1 transfer speed(mcl/s)','Syringe 2 transfer volume (mcl)','Syringe 2 transfer speed(mcl/s)')
        self.entries = controller.createentries(self, fields)
        button2 = tk.Button(self, text="Run", font=SMALL_FONT, command=lambda : self.runfunc(controller))
        button2.grid(row=(3+len(fields)),column=3,pady=10,padx=10,)
    
    def runfunc(self, controller):
        params = controller.getparams(self.entries)
        pulse2(*params)

class pulse3washPage(tk.Frame):
    
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="SYRINGE SYSTEM INTERFACE",font=LARGE_FONT)
        label.grid(row=0,column=0,columnspan=4,pady=10,padx=10)
        label = tk.Label(self, text="Wash by removing THEN replacing the media",font=SMALL_FONT)
        label.grid(row=1,column=0,columnspan=4)
        #self.grid_columnconfigure(1, weight=1)
        self.grid_columnconfigure(2, weight=1)
        
        button1 = tk.Button(self, text="Back", font=SMALL_FONT,command=lambda: controller.show_frame(StartPage))
        button1.grid(row=2,column=0,pady=0,padx=10,)
                            
        ## The fields:
        fields = ('Initial wait (s)','Pulse length (s)','Blank length (s)','Wash interval (s)','Number of pulses','Number of washes','Initial blank media (mcl)','Blank transfer volume (mcl)','Waste transfer volume (mcl)','Ligand transfer volume (mcl)','Blank transfer speed (mcl/s)','Waste transfer speed (mcl/s)','Ligand transfer speed (mcl/s)')
        self.entries = controller.createentries(self, fields)
        button2 = tk.Button(self, text="Run", font=SMALL_FONT,command=lambda : self.runfunc(controller))
        button2.grid(row=(3+len(fields)),column=3,pady=10,padx=10,)
    
    def runfunc(self, controller):
        params = controller.getparams(self.entries)
        pulse3wash(*params)

class pulse3flowPage(tk.Frame):
    
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="SYRINGE SYSTEM INTERFACE",font=LARGE_FONT)
        label.grid(row=0,column=0,columnspan=4,pady=10,padx=10)
        label = tk.Label(self, text="Wash by continuous flow of specified volume",font=SMALL_FONT)
        label.grid(row=1,column=0,columnspan=4)
        #self.grid_columnconfigure(1, weight=1)
        self.grid_columnconfigure(2, weight=1)
        button1 = tk.Button(self, text="Back", font=SMALL_FONT,command=lambda: controller.show_frame(StartPage))
        button1.grid(row=2,column=0,pady=0,padx=10,)
        ## The fields:
        fields = ('Initial wait (s)','Pulse length (s)','Blank length (s)','Wash interval (s)','Number of pulses','Number of washes','Initial blank media (mcl)','Wash transfer volume (mcl)','Wash transfer speed (mcl/s)','Blank transfer volume (mcl)','Waste transfer volume (mcl)','Ligand transfer volume (mcl)','Blank transfer speed (mcl/s)','Waste transfer speed (mcl/s)','Ligand transfer speed (mcl/s)')
        self.entries = controller.createentries(self, fields)                 
        button2 = tk.Button(self, text="Run", font=SMALL_FONT, command=lambda : self.runfunc(controller))
        button2.grid(row=(3+len(fields)),column=3,pady=10,padx=10,)
    
    def runfunc(self, controller):
        params = controller.getparams(self.entries)
        pulse3flow(*params)

# global app
# app = SyringeGUI()
# app.title('GUI')
# app.mainloop()
