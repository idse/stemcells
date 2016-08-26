"""
Contains functions which perform serial communication with arduino.
Functions:
-setup                set up the serial connection
-disconnect           close the serial connection
-m1                   drive first motor (fresh media)
-m2                   drive second motor (waste media)
-m3                   drive third motor (ligand )

Name: drive
Author: Clayton Little
"""

ser = 0  # set up serial object

def setup(port,br):
    """
    Sets up connection to the Arduino
    port is int of USB port. Look at arduino console to find it.
    br is the baud rate. This should be 9600.
    """
    import serial        # the PySerial package
    global ser           # so other functions can access
    ser = serial.Serial('COM'+str(port),br) # open serial connection
    connected = False    # set up Boolean indicating connection status
    ## Because the serial connection is not immediate, we want to
    ## wait until it is available before sending commands:
    while not connected:
        serin = ser.read()
        connected = True
    print("READY")
    print("Serial information:")
    print(ser)
    return ser  
    
def disconnect():
    ser.close()
    
    
def m1(dr,vol,spd=50):
    """
    Sends a single command to the first motor.
    
    dr is an int specifying direction. 1 is inject. 0 is remove.
    vol is the volume transfer in microLiters.
    spd is the speed of fluid transfer in microLiters per seconds.
    """
    ser.write(bytes(','.join(['1',str(dr),str(vol),str(spd)])+'*',"ascii"))
    

def m2(dr,vol,spd=50):
    """
    Sends a single command to the first motor.
    
    dr is an int specifying direction. 1 is inject. 0 is remove.
    vol is the volume transfer in microLiters.
    spd is the speed of fluid transfer in microLiters per seconds.
    """
    ser.write(bytes(','.join(['2',str(dr),str(vol),str(spd)])+'*',"ascii"))
    
    
def m3(dr,vol,spd=50):
    """
    Sends a single command to the first motor.
    
    dr is an int specifying direction. 1 is inject. 0 is remove.
    vol is the volume transfer in microLiters.
    spd is the speed of fluid transfer in microLiters per seconds.
    """
    ser.write(bytes(','.join(['3',str(dr),str(vol),str(spd)])+'*',"ascii"))

def m4(dr,vol,spd=50):
    """
    Sends a single command to the first motor.
    
    dr is an int specifying direction. 1 is inject. 0 is remove.
    vol is the volume transfer in microLiters.
    spd is the speed of fluid transfer in microLiters per seconds.
    """
    ser.write(bytes(','.join(['4',str(dr),str(vol),str(spd)])+'*',"ascii"))
