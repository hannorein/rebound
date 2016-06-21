# -*- coding: utf-8 -*-

"""
Pull data from HORIZONS and format it for use as a REBOUND particle. 

"""
from __future__ import print_function
import telnetlib
import datetime
import re
import sys
import math

__all__ = ["getParticle"]

# Default date for orbital elements is the current time when first particle added, if no date is passed.
# Cached at the beginning to ensure that all particles are synchronized.
# If a date is passed, the same date is used for all subsequent particle adds (that don't themselves pass a date).
INITDATE = None

def getParticle(particle=None, m=None, x=None, y=None, z=None, vx=None, vy=None, vz=None, primary=None, a=None, anom=None, e=None, omega=None, inc=None, Omega=None, MEAN=None, date=None, plane="ecliptic"):   
    if plane not in ["ecliptic","frame"]:
        raise AttributeError("Reference plane needs to be either 'ecliptic' or 'frame'. See Horizons for a definition of these coordinate systems.")
    if date is not None:
        if type(date) is datetime.datetime:
            pass
        elif type(date) is str:
            date = datetime.datetime.strptime(date,"%Y-%m-%d %H:%M")
        else:
            raise AttributeError("Unknown date format.")
    # set the cached initialization time if it's not set
    global INITDATE
    if INITDATE is None:
        INITDATE = date if date is not None else datetime.datetime.utcnow()

    if date is None: # if no date passed, used cached value
        date = INITDATE
    print("Searching NASA Horizons for '%s'... "%(particle),end="")
    sys.stdout.flush()

    t = telnetlib.Telnet()
    t.open('horizons.jpl.nasa.gov', 6775, 20)
    expect = ( ( b'Horizons>', particle+'\n' ),
               ( b'Continue.*:', 'y\n' ),
               ( b'Select.*E.phemeris.*:', 'E\n'),
               ( b'Observe.*:', 'v\n' ),
               ( b'Coordinate center.*:', '@0\n' ),
               ( b'Reference plane.*:', plane+'\n' ),
               ( b'Starting.* :', date.strftime("%Y-%m-%d %H:%M:%S")+'\n' ),
               ( b'Ending.* :', (date + datetime.timedelta(minutes=1)).strftime("%Y-%m-%d %H:%M")+'\n' ),
               ( b'Output interval.*:', '2\n' ),
               ( b'Accept default output \[.*:', 'n\n' ),
               ( b'Output reference frame \[.*:', 'J2000\n' ),
               ( b'Corrections \[.*:', 'NONE\n' ),
               ( b'Output units \[.*:', '1\n' ),
               ( b'Spreadsheet CSV format.*\[.*:', 'NO\n' ),
               ( b'Label cartesian output.*\[.*:', 'NO\n' ),
               ( b'Select output table type.*\[.*:', '3\n' ),
               ( b'Scroll . Page: .*%', ' '),
               ( b'Select\.\.\. .A.gain.* :', 'X\n' ),
               ( b'Select \.\.\. .F.tp.*:', 'selectID' )
    )
    p = Particle()
    startdata = 0
    message = ""
    idn = None
    while True:
        try:
            answer = t.expect(list(i[0] for i in expect), 8)
        except EOFError:
            break
        a = answer[2].decode()
        if "$$SOE" in a:
            lines = a.split("\n")
            for line in lines:
                if line.strip() == "$$EOE":
                    break
                if startdata == 2:
                    pos = [float(i) for i in line.split()]
                    p.x = pos[0]
                    p.y = pos[1]
                    p.z = pos[2]
                if startdata == 3:
                    vel = [float(i) for i in line.split()]
                    p.vx = vel[0]
                    p.vy = vel[1]
                    p.vz = vel[2]
                if startdata > 0:
                    startdata += 1
                if "Target body name:" in line:
                    print("Found: %s." % line[18:50].strip())
                    try:
                        idn = re.search(r"\(([0-9]+)\)", line).group(1)
                    except:
                        pass
                if line.strip() == "$$SOE":
                    startdata = 1
        message += a.replace(chr(27)+"[H","")
        if answer[0] != -1:
            if expect[answer[0]][1] == "selectID":
                try:
                    idn = a.split("ID#")[1].split("\n")[2].split()[0]
                except:
                    try:
                        idn = a.split("Record #")[1].split("\n")[2].split()[0]
                    except:
                        raise Exception("Error while trying to find object.")
                t.write((idn+"\n").encode())
            else:
                t.write(expect[answer[0]][1].encode())
        else:
            pass
            #print "NOT FOUND!!"
            #print answer
    if startdata == 0:
        print(message)
        raise SyntaxError("Object not found. See above output from HORIZONS. Please try different identifier or look up JPL Body Number.")
    if m is not None:
        p.m = m
    elif idn is not None:
        try:
            p.m = float(re.search(r"BODY%d\_GM .* \( *([\.E\+\-0-9]+ *)\)"%int(idn), HORIZONS_MASS_DATA).group(1))
            p.m /= Gkmkgs # divide by G (horizons masses give GM)
        except:
            print("Warning: Mass cannot be retrieved from NASA HORIZONS. Set to 0.")
            p.m = 0
    else:
        print("Warning: Mass cannot be retrieved from NASA HORIZONS. Set to 0.")
        p.m = 0
    return p


# There is currently no way to get mass data from HORIZONS.
# The following data was provided by Jon Giorgini (10 May 2015)
# Units: km^3/s^2

Gkmkgs = 6.674e-20 # units of km^3/kg/s^2

HORIZONS_MASS_DATA = """
    BODY1_GM       = ( 2.2031780000000021E+04 )
    BODY2_GM       = ( 3.2485859200000006E+05 )
    BODY3_GM       = ( 4.0350323550225981E+05 )
    BODY4_GM       = ( 4.2828375214000022E+04 )
    BODY5_GM       = ( 1.2671276480000021E+08 )
    BODY6_GM       = ( 3.7940585200000003E+07 )
    BODY7_GM       = ( 5.7945486000000080E+06 )
    BODY8_GM       = ( 6.8365271005800236E+06 )
    BODY9_GM       = ( 9.7700000000000068E+02 )
    BODY10_GM      = ( 1.3271244004193938E+11 )

    BODY199_GM     = ( 2.2031780000000021E+04 )
    BODY299_GM     = ( 3.2485859200000006E+05 )
    BODY399_GM     = ( 3.9860043543609598E+05 )
    BODY499_GM     = ( 4.282837362069909E+04  )
    BODY599_GM     = ( 1.266865349218008E+08  )
    BODY699_GM     = ( 3.793120749865224E+07  )
    BODY799_GM     = ( 5.793951322279009E+06  )
    BODY899_GM     = ( 6.835099502439672E+06  )
    BODY999_GM     = ( 8.696138177608748E+02  )

    BODY301_GM     = ( 4.9028000661637961E+03 )

    BODY401_GM     = ( 7.087546066894452E-04 )
    BODY402_GM     = ( 9.615569648120313E-05 )

    BODY501_GM     = ( 5.959916033410404E+03 )
    BODY502_GM     = ( 3.202738774922892E+03 )
    BODY503_GM     = ( 9.887834453334144E+03 )
    BODY504_GM     = ( 7.179289361397270E+03 )
    BODY505_GM     = ( 1.378480571202615E-01 )

    BODY601_GM     = ( 2.503522884661795E+00 )
    BODY602_GM     = ( 7.211292085479989E+00 )
    BODY603_GM     = ( 4.121117207701302E+01 )
    BODY604_GM     = ( 7.311635322923193E+01 )
    BODY605_GM     = ( 1.539422045545342E+02 )
    BODY606_GM     = ( 8.978138845307376E+03 )
    BODY607_GM     = ( 3.718791714191668E-01 )
    BODY608_GM     = ( 1.205134781724041E+02 )
    BODY609_GM     = ( 5.531110414633374E-01 )
    BODY610_GM     = ( 1.266231296945636E-01 )
    BODY611_GM     = ( 3.513977490568457E-02 )
    BODY615_GM     = ( 3.759718886965353E-04 )
    BODY616_GM     = ( 1.066368426666134E-02 )
    BODY617_GM     = ( 9.103768311054300E-03 )

    BODY701_GM     = ( 8.346344431770477E+01 )
    BODY702_GM     = ( 8.509338094489388E+01 )
    BODY703_GM     = ( 2.269437003741248E+02 )
    BODY704_GM     = ( 2.053234302535623E+02 )
    BODY705_GM     = ( 4.319516899232100E+00 )

    BODY801_GM     = ( 1.427598140725034E+03 )

    BODY901_GM     = ( 1.058799888601881E+02 )
    BODY902_GM     = ( 3.048175648169760E-03 )
    BODY903_GM     = ( 3.211039206155255E-03 )
    BODY904_GM     = ( 1.110040850536676E-03 )

    BODY2000001_GM = ( 6.3130000000000003E+01 )
    BODY2000002_GM = ( 1.3730000000000000E+01 )
    BODY2000003_GM = ( 1.8200000000000001E+00 )
    BODY2000004_GM = ( 1.7289999999999999E+01 )
    BODY2000006_GM = ( 9.3000000000000005E-01 )
    BODY2000007_GM = ( 8.5999999999999999E-01 )
    BODY2000010_GM = ( 5.7800000000000002E+00 )
    BODY2000015_GM = ( 2.1000000000000001E+00 )
    BODY2000016_GM = ( 1.8100000000000001E+00 )
    BODY2000029_GM = ( 8.5999999999999999E-01 )
    BODY2000052_GM = ( 1.5900000000000001E+00 )
    BODY2000065_GM = ( 9.1000000000000003E-01 )
    BODY2000087_GM = ( 9.8999999999999999E-01 )
    BODY2000088_GM = ( 1.0200000000000000E+00 )
    BODY2000433_GM = (  4.463E-4 )
    BODY2000511_GM = ( 2.2599999999999998E+00 )
    BODY2000704_GM = ( 2.1899999999999999E+00 )
"""


# Import at the end to avoid circular dependence
from .particle import *
