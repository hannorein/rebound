import telnetlib
import datetime

__all__ = ["getParticle"]

INITDATE = datetime.datetime.today()

def getParticle(particle=None, m=None, x=None, y=None, z=None, vx=None, vy=None, vz=None, primary=None, a=None, anom=None, e=None, omega=None, inc=None, Omega=None, MEAN=None, date=None):   
    if date is not None:
        date = datetime.datetime.strptime(date, "%Y-%m-%d %H:%M")
    else:
        date = INITDATE
    print("Trying to get data from NASA Horizons for body '%s' at date '%s'..."%(particle,date.strftime("%Y-%m-%d %H:%M")))

    t = telnetlib.Telnet()
    t.open('horizons.jpl.nasa.gov', 6775)
    expect = ( ( r'Horizons>', particle+'\n' ),
               ( r'Continue.*:', 'y\n' ),
               ( r'Select.*E.phemeris.*:', 'E\n'),
               ( r'Observe.*:', 'v\n' ),
               ( r'Coordinate center.*:', '@Sun\n' ),
               ( r'Reference plane.*:', 'eclip\n' ),
               ( r'Starting.* :', date.strftime("%Y-%m-%d %H:%M")+'\n' ),
               ( r'Ending.* :', (date + datetime.timedelta(hours=1)).strftime("%Y-%m-%d %H:%M")+'\n' ),
               ( r'Output interval.*:', '1d\n' ),
               ( r'Accept default output.*:', 'y\n' ),
               ( r'Scroll . Page: .*%', ' '),
               ( r'Select\.\.\. .A.gain.* :', 'X\n' )
    )
    p = Particle()
    startdata = 0
    message = ""
    while True:
        try:
            answer = t.expect(list(i[0] for i in expect), 4)
        except EOFError:
            break
        if "$$SOE" in answer[2]:
            lines = answer[2].split("\n")
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
                    auperday2auperyeartwopi = 58.130101
                    p.vx = vel[0]*auperday2auperyeartwopi 
                    p.vy = vel[1]*auperday2auperyeartwopi 
                    p.vz = vel[2]*auperday2auperyeartwopi 
                if startdata > 0:
                    startdata += 1
                if line[0:18] == "Target body name: ":
                    print("Found: %s" % line[18:50].strip())
                if line.strip() == "$$SOE":
                    startdata = 1
        message += answer[2].replace(chr(27)+"[H","")
        t.write(expect[answer[0]][1])
    if startdata == 0:
        print(message)
        raise SyntaxError("Object not found. See above output from HORIZONS. Please try different identifier or look up JPL Body Number.")
    if m is not None:
        p.m = m
    else:
        print("Mass cannot be retrieved from NASA HORIZONS. Please add manually.")
    return p

# Import at the end to avoid circular dependence
from .particle import *
