# -*- coding: utf-8 -*-

"""
Pull data from HORIZONS and format it for use as a REBOUND particle. 

"""
from __future__ import print_function

import datetime
import re
import urllib
import urllib.parse
import urllib.request
import warnings

__all__ = ["getParticle"]

# Default date for orbital elements is the current time when first particle added, if no date is passed.
# Cached at the beginning to ensure that all particles are synchronized.
# If a date is passed, the same date is used for all subsequent particle adds (that don't themselves pass a date).

INITDATE = None


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
        "STEP_SIZE": quote("2"),  # not sure what's the unit
        "REF_SYSTEM": quote("J2000"),
        "VEC_CORR": quote("NONE"),
        "OUT_UNITS": quote("KM-S"),
        "CSV_FORMAT": quote("NO"),
        "VEC_DELTA_T": quote("NO"),
        "VEC_TABLE": quote("3"),
        "VEC_LABELS": quote("NO")

    }
    url = "https://ssd.jpl.nasa.gov/api/horizons.api?" + urllib.parse.urlencode(get_params)
    with urllib.request.urlopen(url) as f:
        return f.read().decode()


def getParticle(particle=None, m=None, x=None, y=None, z=None, vx=None, vy=None, vz=None, primary=None, a=None,
                anom=None, e=None, omega=None, inc=None, Omega=None, MEAN=None, date=None, plane="ecliptic", hash=0):
    if plane not in ["ecliptic", "frame"]:
        raise AttributeError(
            "Reference plane needs to be either 'ecliptic' or 'frame'. See Horizons for a definition of these coordinate systems.")
    jd = False
    if date is not None:
        if isinstance(date, datetime.datetime):
            pass
        elif isinstance(date, str):
            try:
                if date[0:2] == "JD":
                    jd = True
                else:
                    date = datetime.datetime.strptime(date, "%Y-%m-%d %H:%M")
            except:
                raise AttributeError(
                    "An error occured while calculating the date. Use either YYYY-MM-DD HH:MM or JDxxxxxxx.xxxxxx")
    # set the cached initialization time if it's not set
    global INITDATE
    if INITDATE is None:
        INITDATE = date if date is not None else datetime.datetime.utcnow()

    if date is None:  # if no date passed, used cached value
        date = INITDATE

    if isinstance(date, datetime.datetime):
        # date is a datetime object
        datestart = date.strftime("%Y-%m-%d %H:%M:%S")
        dateend = (date + datetime.timedelta(minutes=1)).strftime("%Y-%m-%d %H:%M")
    else:
        # Assume date is in JD
        datestart = date
        if "." in date:  # end date slightly later
            dateend = date + "1"
        else:
            dateend = date + ".1"

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
    lines = body.split("$$SOE")[-1].split("\n")
    p = Particle()

    p.x, p.y, p.z = [float(i) for i in lines[2].split()]
    p.vx, p.vy, p.vz = [float(i) for i in lines[3].split()]

    match = re.search(r"Target body name: (.+) \(([0-9]+)\)", body)
    if match:
        bodyname = match.group(1).strip()
        idn = match.group(2)
        print("Found: {}.".format(bodyname), "(chosen from query '{}')".format(particle) if made_choice else "")

    if m is not None:
        p.m = m
    elif idn is not None:
        try:
            p.m = float(
                re.search(r"BODY{:d}\_GM .* \( *([\.E\+\-0-9]+ *)\)".format(int(idn)), HORIZONS_MASS_DATA).group(1))
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
# Last updated: May 31 2018. 
# Source: ftp://ssd.jpl.nasa.gov/pub/xfr/gm_Horizons.pck
# Units: km^3/s^2

Gkmkgs = 6.67408e-20  # units of km^3/kg/s^2

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
    BODY599_GM     = ( 1.266865349115908E+08  )
    BODY699_GM     = ( 3.793120723493890E+07  )
    BODY799_GM     = ( 5.793951322279009E+06  )
    BODY899_GM     = ( 6.835099502439672E+06  )
    BODY999_GM     = ( 8.693390780381926E+02  )
 
    BODY301_GM     = ( 4.9028000661637961E+03 )

    BODY401_GM     = ( 7.087546066894452E-04 )
    BODY402_GM     = ( 9.615569648120313E-05 )

    BODY501_GM     = ( 5.959924010272514E+03 )
    BODY502_GM     = ( 3.202739815114734E+03 )
    BODY503_GM     = ( 9.887819980080976E+03 )
    BODY504_GM     = ( 7.179304867611079E+03 )
    BODY505_GM     = ( 1.487604677404272E-01 )
    BODY506_GM     = ( 1.380080696078966E-01 )
 
    BODY601_GM     = ( 2.503458199931431E+00 )
    BODY602_GM     = ( 7.211185066509890E+00 )
    BODY603_GM     = ( 4.120856508658532E+01 )
    BODY604_GM     = ( 7.311574218947423E+01 )
    BODY605_GM     = ( 1.539419035933117E+02 )
    BODY606_GM     = ( 8.978137030983542E+03 )
    BODY607_GM     = ( 3.712085754472412E-01 )
    BODY608_GM     = ( 1.205095752388872E+02 )
    BODY609_GM     = ( 5.532371285376407E-01 )
    BODY610_GM     = ( 1.265765099012197E-01 )
    BODY611_GM     = ( 3.512333288208074E-02 )
    BODY612_GM     = ( 3.424829447502984E-04 )
    BODY615_GM     = ( 3.718871247516475E-04 )
    BODY616_GM     = ( 1.075208001007610E-02 )
    BODY617_GM     = ( 9.290325122028795E-03 )

    BODY701_GM     = ( 8.346344431770477E+01 )
    BODY702_GM     = ( 8.509338094489388E+01 )
    BODY703_GM     = ( 2.269437003741248E+02 )
    BODY704_GM     = ( 2.053234302535623E+02 )
    BODY705_GM     = ( 4.319516899232100E+00 )

    BODY801_GM     = ( 1.427598140725034E+03 )
  
    BODY901_GM     = ( 1.062509269522026E+02 )
    BODY902_GM     = ( 2.150552267969335E-03 )
    BODY903_GM     = ( 3.663917107480563E-03 )
    BODY904_GM     = ( 4.540734312735987E-04 )
    BODY905_GM     = ( 2.000000000000000E-20 )

    BODY2000001_GM = ( 6.2809393000000000E+01 )
    BODY2000002_GM = ( 1.3923011000000001E+01 )
    BODY2000003_GM = ( 1.6224149999999999E+00 )
    BODY2000004_GM = ( 1.7288008999999999E+01 )
    BODY2000010_GM = ( 5.5423920000000004E+00 )
    BODY2000015_GM = ( 2.0981550000000002E+00 )
    BODY2000016_GM = ( 1.5300480000000001E+00 )
    BODY2000031_GM = ( 2.8448720000000001E+00 )
    BODY2000048_GM = ( 1.1351590000000000E+00 )
    BODY2000052_GM = ( 1.1108039999999999E+00 )
    BODY2000065_GM = ( 1.4264810000000001E+00 )
    BODY2000087_GM = ( 9.8635300000000004E-01 )
    BODY2000088_GM = ( 1.1557990000000000E+00 )
    BODY2000433_GM = ( 4.463E-4 )
    BODY2000451_GM = ( 1.0295259999999999E+00 )
    BODY2000511_GM = ( 2.3312860000000000E+00 )
    BODY2000704_GM = ( 2.3573170000000001E+00 )
"""

# Import at the end to avoid circular dependence
from .particle import *
