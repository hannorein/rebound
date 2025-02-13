# -*- coding: utf-8 -*-
"""
Pull data from HORIZONS and format it for use as a REBOUND particle. 

"""
import datetime
import re
import warnings
import sys
from .units import convert_mass

HORIZONSBASEURL = "https://ssd.jpl.nasa.gov/api/horizons.api?"

if "pyodide" in sys.modules:
    from urllib.parse import urlencode
    from pyodide.http import open_url as urlopen
    # Use CORS proxy
    HORIZONSBASEURL = "https://rebound.hanno-rein.de/api/horizons.api?"
else:
    try:
        from urllib.parse import urlencode
        from urllib.request import urlopen
    except ImportError:
        from urllib import urlencode
        from urllib2 import urlopen

# Default date for orbital elements is the current time when first particle added, if no date is passed.
# Cached at the beginning to ensure that all particles are synchronized.
# If a date is passed, the same date is used for all subsequent particle adds (that don't themselves pass a date).

INITDATE = None

# If this variable is set to "unverified" then the SSL context for API requests does not check certificates.
SSL_CONTEXT = None

def quote(text):
    return "'{}'".format(text)


def api_request(particle, datestart, dateend, plane):
    get_params = {
        "format": "text",
        "COMMAND": quote(particle),
        "START_TIME": quote(str(datestart)),
        "STOP_TIME": quote(str(dateend)),
        "MAKE_EPHEM": quote("YES"),
        "EPHEM_TYPE": quote("VECTORS"),
        "CENTER": quote("@0"),
        "REF_PLANE": quote(plane),
        "STEP_SIZE": quote("2"),  # seconds
        "REF_SYSTEM": quote("J2000"),
        "VEC_CORR": quote("NONE"),
        "OUT_UNITS": quote("KM-S"),
        "CSV_FORMAT": quote("NO"),
        "VEC_DELTA_T": quote("NO"),
        "VEC_TABLE": quote("3"),
        "VEC_LABELS": quote("NO")

    }
    url =  HORIZONSBASEURL + urlencode(get_params)
    # don't use a context manager for python2 compatibility

    if SSL_CONTEXT == "unverified":
        import ssl
        ssl_context = ssl._create_unverified_context()
    else:
        ssl_context = None
    try:
        f = urlopen(url,context=ssl_context)
    except Exception as e:
        raise RuntimeError("An error occured while accessing NASA HORIZONS. If this is a SSL certificate issue, you can try disabling the certificate verification by setting rebound.horizons.SSL_CONTEXT = 'unverified'.") from e

    if "pyodide" in sys.modules:
        body = f.read()
    else:
        body = f.read().decode()
    f.close()
    return body


def query_horizons_for_particle(mass_unit=None, particle=None, m=None, x=None, y=None, z=None, vx=None, vy=None, vz=None, primary=None, a=None,
                anom=None, e=None, omega=None, inc=None, Omega=None, MEAN=None, date=None, plane="ecliptic", hash=0):
    if plane not in ["ecliptic", "frame"]:
        raise AttributeError(
            "Reference plane needs to be either 'ecliptic' or 'frame'. See Horizons for a definition of these coordinate systems.")
    if date is not None:
        if isinstance(date, datetime.datetime):
            pass
        elif isinstance(date, str):
            if date[0:2] != "JD":
                formats = ["%Y-%m-%d", "%Y-%m-%d %H:%M", "%Y-%m-%d %H:%M:%S"] # allowed formats
                found_match = False
                for f in formats:
                    try:
                        date = datetime.datetime.strptime(date, f)
                        found_match = True
                    except:
                        continue
                if found_match == False:
                    raise AttributeError("An error occured while calculating the date. Use one "+" or ".join(formats) + " or JDxxxxxxx.xxxxxx")
    # set the cached initialization time if it's not set
    global INITDATE
    if INITDATE is None:
        INITDATE = date if date is not None else datetime.datetime.utcnow()

    if date is None:  # if no date passed, used cached value
        date = INITDATE

    if isinstance(date, datetime.datetime):
        # date is a datetime object
        datestart = date.strftime("%Y-%m-%d %H:%M:%S")
        dateend = (date + datetime.timedelta(minutes=1)).strftime("%Y-%m-%d %H:%M:%S")
    else:
        # Assume date is in JD with format JDxxxxxx.xxxx
        datestart = date
        date_f = float(re.sub("[^0-9\\.]","",date))
        dateend = "JD%.8f"%(date_f+0.1)

    print("Searching NASA Horizons for '{}'... ".format(particle))
    idn = None
    body = api_request(particle, datestart, dateend, plane)
    made_choice = False
    if "Multiple major-bodies match string" in body:
        try:
            idn = body.split("ID#")[1].split("\n")[2].split()[0]
        except KeyError:
            try:
                idn = body.split("Record #")[1].split("\n")[2].split()[0]
            except:
                raise Exception("Error while trying to find object.")

        made_choice = True
        body = api_request(idn, datestart, dateend, plane)
    elif "Matching small-bodies" in body:
        for line in body.split("\n"):
            try:
                first_word = line.split()[0]
            except IndexError:
                continue
            if first_word.isdecimal():
                idn = first_word
                break
        if not idn:
            raise Exception("Error while trying to find object.")
        made_choice = True
        body = api_request(idn, datestart, dateend, plane)

    lines = body.split("$$SOE")[-1].split("\n")
    p = Particle()

    p.x, p.y, p.z = [float(i) for i in lines[2].split()]
    p.vx, p.vy, p.vz = [float(i) for i in lines[3].split()]

    match = re.search(r"Target body name: (.+) \(([0-9]+)\)", body)
    if match:
        bodyname = match.group(1).strip()
        idn = match.group(2)
        print("Found: {} ({})".format(bodyname, idn), "(chosen from query '{}')".format(particle) if made_choice else "")
    else:
        # fall back to more general regex
        match = re.search(r"Target body name: (.+) {", body)
        if match:
            bodyname = match.group(1).strip()
            print("Found: {}".format(bodyname), "(chosen from query '{}')".format(particle) if made_choice else "")
        else:
            print("Found body (Name could not be detected)")
    if m is not None:
        if mass_unit is not None:
            p.m = convert_mass(m, mass_unit, "kg")
        else:
            ## Assume kg
            p.m = m
    elif idn is not None:
        try:
            p.m = float(
                re.search(r"BODY{:d}\_GM .* \( *([\.DE\+\-0-9]+ *)\)".format(int(idn)), HORIZONS_MASS_DATA)
                    .group(1).replace("D+", "E+")
            )
            p.m /= Gkmkgs  # divide by G (horizons masses give GM)
        except AttributeError:
            warnings.warn("Warning: Mass cannot be retrieved from NASA HORIZONS. Set to 0.", RuntimeWarning)
            p.m = 0
    else:
        warnings.warn("Warning: Mass cannot be retrieved from NASA HORIZONS. Set to 0.", RuntimeWarning)
        p.m = 0
    p.hash = hash
    return p


# There is currently no way to get mass data from HORIZONS.
# The following data was provided by Jon Giorgini (10 May 2015)
# Last updated: Sep 15 2021.
# Source: ftp://ssd.jpl.nasa.gov/pub/xfr/gm_Horizons.pck
# Units: km^3/s^2

Gkmkgs = 6.67408e-20  # units of km^3/kg/s^2

HORIZONS_MASS_DATA = """
    BODY1_GM       = ( 2.2031868551400003D+04 )
    BODY2_GM       = ( 3.2485859200000000D+05 )
    BODY3_GM       = ( 4.0350323562548019D+05 )
    BODY4_GM       = ( 4.2828375815756102D+04 )
    BODY5_GM       = ( 1.2671276409999998D+08 )
    BODY6_GM       = ( 3.7940584841799997D+07 )
    BODY7_GM       = ( 5.7945563999999985D+06 )
    BODY8_GM       = ( 6.8365271005803989D+06 )
    BODY9_GM       = ( 9.7550000000000000D+02 )
    BODY10_GM      = ( 1.3271244004127942D+11 )

    BODY199_GM     = ( 2.2031868551400003D+04 )
    BODY299_GM     = ( 3.2485859200000000D+05 )
    BODY399_GM     = ( 3.9860043550702266D+05 )
    BODY499_GM     = ( 4.282837362069909E+04  )
    BODY599_GM     = ( 1.266865319003704E+08  )
    BODY699_GM     = ( 3.793120615901047E+07  )
    BODY799_GM     = ( 5.793951322279009E+06  )
    BODY899_GM     = ( 6.835099968446816E+06  )
    BODY999_GM     = ( 8.699633756209835E+02  )
 
    BODY301_GM     = ( 4.9028001184575496D+03 )

    BODY401_GM     = ( 7.087546066894452E-04 )
    BODY402_GM     = ( 9.615569648120313E-05 )

    BODY501_GM     = ( 5.959915466180539E+03 )
    BODY502_GM     = ( 3.202712099607295E+03 )
    BODY503_GM     = ( 9.887832752719638E+03 )
    BODY504_GM     = ( 7.179283402579837E+03 )
    BODY505_GM     = ( 1.645634534798259E-01 )
    BODY506_GM     = ( 1.515524299611265E-01 )
    BODY514_GM     = ( 3.014800000000000E-02 )
    BODY515_GM     = ( 1.390000000000000E-04 )
    BODY516_GM     = ( 2.501000000000000E-03 )
 
    BODY601_GM     = ( 2.503617062809250E+00 )
    BODY602_GM     = ( 7.210497553340731E+00 )
    BODY603_GM     = ( 4.121405263872402E+01 )
    BODY604_GM     = ( 7.311617801921636E+01 )
    BODY605_GM     = ( 1.539409077211430E+02 )
    BODY606_GM     = ( 8.978137369591670E+03 )
    BODY607_GM     = ( 3.704182596063880E-01 )
    BODY608_GM     = ( 1.205081845217891E+02 )
    BODY609_GM     = ( 5.581081743011904E-01 )
    BODY610_GM     = ( 1.265765099012197E-01 )
    BODY611_GM     = ( 3.512333288208074E-02 )
    BODY612_GM     = ( 4.551624250415933E-04 )
    BODY615_GM     = ( 3.718871247516475E-04 )
    BODY616_GM     = ( 1.075208001007610E-02 )
    BODY617_GM     = ( 9.290325122028795E-03 )

    BODY701_GM     = ( 8.346344431770477E+01 )
    BODY702_GM     = ( 8.509338094489388E+01 )
    BODY703_GM     = ( 2.269437003741248E+02 )
    BODY704_GM     = ( 2.053234302535623E+02 )
    BODY705_GM     = ( 4.319516899232100E+00 )

    BODY801_GM     = ( 1.428495462910464E+03 )
    BODY803_GM     = ( 8.530281246540886E-03 )
    BODY804_GM     = ( 2.358873197992170E-02 )
    BODY805_GM     = ( 1.167318403814998E-01 )
    BODY806_GM     = ( 1.898985039060690E-01 )
    BODY807_GM     = ( 2.548437405693583E-01 )
    BODY808_GM     = ( 2.583422379120727E+00 )
  
    BODY901_GM     = ( 1.061744232879427E+02 )
    BODY902_GM     = ( 1.800000000000000E-03 )
    BODY903_GM     = ( 2.249146225742025E-03 )
    BODY904_GM     = ( 9.000000000000001E-05 )
    BODY905_GM     = ( 2.000000000000000E-06 )

    BODY2000001_GM = ( 6.2628888644409933D+01 )
    BODY2000002_GM = ( 1.3665878145967422D+01 )
    BODY2000003_GM = ( 1.9205707002025889D+00 )
    BODY2000004_GM = ( 1.7288232879171513D+01 )
    BODY2000007_GM = ( 1.1398723232184107D+00 )
    BODY2000010_GM = ( 5.6251476453852289D+00 )
    BODY2000015_GM = ( 2.0230209871098284D+00 )
    BODY2000016_GM = ( 1.5896582441709424D+00 )
    BODY2000031_GM = ( 1.0793714577033560D+00 )
    BODY2000052_GM = ( 2.6830359242821795D+00 )
    BODY2000065_GM = ( 9.3810575639151328D-01 )
    BODY2000087_GM = ( 2.1682320736996910D+00 )
    BODY2000088_GM = ( 1.1898077088121908D+00 )
    BODY2000107_GM = ( 1.4437384031866001D+00 )
    BODY2000433_GM = ( 4.463E-4 )
    BODY2000511_GM = ( 3.8944831481705644D+00 )
    BODY2000704_GM = ( 2.8304096393299849D+00 )
"""

# Import at the end to avoid circular dependence
from .particle import *
