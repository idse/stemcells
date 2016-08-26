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

## Set up the serial connection
drive.setup(10,9600) # (COM port number,baud rate)





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
        
        button1 = tk.Button(self, text="TWO SYRINGE PUMPS", font=SMALL_FONT,
                            command=lambda: controller.show_frame(pulse2Page))
        button1.pack()

        button2 = tk.Button(self, text="THREE SYRINGE PUMPS W/ EXCHANGE WASH", font=SMALL_FONT,
                            command=lambda: controller.show_frame(pulse3washPage))
        button2.pack()
        
        button3 = tk.Button(self, text="THREE SYRINGE PUMPS W/ FLOW WASH", font=SMALL_FONT,
                            command=lambda: controller.show_frame(pulse3flowPage))
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

        button1 = tk.Button(self, text="Back", font=SMALL_FONT,
                            command=lambda: controller.show_frame(StartPage))
        button1.grid(row=2,column=0,pady=0,padx=10,)
        
        ## The fields:
        fields = (
        'Initial wait (s)','Initial blank media (mcl)','Number of pulses',
        'Pulse length (s)','Blank length (s)',
        'Syringe 1 transfer volume (mcl)','Syringe 1 transfer speed(mcl/s)',
        'Syringe 2 transfer volume (mcl)','Syringe 2 transfer speed(mcl/s)')
        self.entries = controller.createentries(self, fields)
            
        button2 = tk.Button(self, text="Run", font=SMALL_FONT,
                            command=lambda : self.runfunc(controller))
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

        button1 = tk.Button(self, text="Back", font=SMALL_FONT,
                            command=lambda: controller.show_frame(StartPage))
        button1.grid(row=2,column=0,pady=0,padx=10,)
        
        ## The fields:
        fields = (
        'Initial wait (s)','Pulse length (s)','Blank length (s)','Wash interval (s)',
        'Number of pulses','Number of washes','Initial blank media (mcl)',
        'Blank transfer volume (mcl)','Waste transfer volume (mcl)','Ligand transfer volume (mcl)',
        'Blank transfer speed (mcl/s)','Waste transfer speed (mcl/s)','Ligand transfer speed (mcl/s)')
        self.entries = controller.createentries(self, fields)
            
        button2 = tk.Button(self, text="Run", font=SMALL_FONT,
                            command=lambda : self.runfunc(controller))
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

        button1 = tk.Button(self, text="Back", font=SMALL_FONT,
                            command=lambda: controller.show_frame(StartPage))
        button1.grid(row=2,column=0,pady=0,padx=10,)
        
        ## The fields:
        fields = (
        'Initial wait (s)','Pulse length (s)','Blank length (s)','Wash interval (s)',
        'Number of pulses','Number of washes','Initial blank media (mcl)',
        'Wash transfer volume (mcl)','Wash transfer speed (mcl/s)',
        'Blank transfer volume (mcl)','Waste transfer volume (mcl)','Ligand transfer volume (mcl)',
        'Blank transfer speed (mcl/s)','Waste transfer speed (mcl/s)','Ligand transfer speed (mcl/s)')
        self.entries = controller.createentries(self, fields)
            
        button2 = tk.Button(self, text="Run", font=SMALL_FONT,
                            command=lambda : self.runfunc(controller))
        button2.grid(row=(3+len(fields)),column=3,pady=10,padx=10,)
    
    def runfunc(self, controller):
        params = controller.getparams(self.entries)
        pulse3flow(*params)
        
global app
app = SyringeGUI()
app.title('GUI')
app.mainloop()