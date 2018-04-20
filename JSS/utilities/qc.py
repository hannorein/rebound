import re
import rebound as rb

MASS_SOLAR		= float(re.search(r"BODY%d\_GM .* \( *([\.E\+\-0-9]+ *)\)"%int("10"), rb.horizons.HORIZONS_MASS_DATA).group(1))
MASS_SOLAR     /= rb.horizons.Gkmkgs # divide by G (horizons masses give GM. units of km^3/kg/s^2)

MASS_MERCURY  	= float(re.search(r"BODY%d\_GM .* \( *([\.E\+\-0-9]+ *)\)"%int("199"), rb.horizons.HORIZONS_MASS_DATA).group(1))
MASS_MERCURY   /= rb.horizons.Gkmkgs # divide by G (horizons masses give GM. units of km^3/kg/s^2)

MASS_VENUS	 	= float(re.search(r"BODY%d\_GM .* \( *([\.E\+\-0-9]+ *)\)"%int("2"), rb.horizons.HORIZONS_MASS_DATA).group(1))
MASS_VENUS	   /= rb.horizons.Gkmkgs # divide by G (horizons masses give GM. units of km^3/kg/s^2)

MASS_EARTH_BC	= float(re.search(r"BODY%d\_GM .* \( *([\.E\+\-0-9]+ *)\)"%int("3"), rb.horizons.HORIZONS_MASS_DATA).group(1))
MASS_EARTH_BC  /= rb.horizons.Gkmkgs # divide by G (horizons masses give GM. units of km^3/kg/s^2)

MASS_MARS_BC	= float(re.search(r"BODY%d\_GM .* \( *([\.E\+\-0-9]+ *)\)"%int("4"), rb.horizons.HORIZONS_MASS_DATA).group(1))
MASS_MARS_BC   /= rb.horizons.Gkmkgs # divide by G (horizons masses give GM. units of km^3/kg/s^2)
