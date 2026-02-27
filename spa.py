from dataclasses import dataclass
import numpy as np
import spa_tables

class SpaData:
    # INPUT VALUES
    year: int = 0         # 4-digit year,      valid range: -2000 to 6000, error code: 1
    month: int = 0        # 2-digit month,         valid range: 1 to  12,  error code: 2
    day: int = 0          # 2-digit day,           valid range: 1 to  31,  error code: 3
    hour: int = 0         # Observer local hour,   valid range: 0 to  24,  error code: 4
    minute: int = 0       # Observer local minute, valid range: 0 to  59,  error code: 5
    second: float = 0     # Observer local second, valid range: 0 to <60,  error code: 6

    delta_ut1: float = 0  # Fractional second difference between UTC and UT which is used
                          # to adjust UTC for earth's irregular rotation rate and is derived
                          # from observation only and is reported in this bulletin:
                          # http:#maia.usno.navy.mil/ser7/ser7.dat,
                          # where delta_ut1 = DUT1
                          # valid range: -1 to 1 second (exclusive), error code 17

    delta_t: float = 0   # Difference between earth rotation time and terrestrial time
                         # It is derived from observation only and is reported in this
                         # bulletin: http:#maia.usno.navy.mil/ser7/ser7.dat,
                         # where delta_t = 32.184 + (TAI-UTC) - DUT1
                         # valid range: -8000 to 8000 seconds, error code: 7

    timezone: float = 0  # Observer time zone (negative west of Greenwich)
                         # valid range: -18   to   18 hours,   error code: 8

    longitude: float = 0 # Observer longitude (negative west of Greenwich)
                         # valid range: -180  to  180 degrees, error code: 9

    latitude: float = 0  # Observer latitude (negative south of equator)
                         # valid range: -90   to   90 degrees, error code: 10

    elevation: float = 0 # Observer elevation [meters]
                         # valid range: -6500000 or higher meters,    error code: 11

    pressure: float = 0  # Annual average local pressure [millibars]
                         # valid range:    0 to 5000 millibars,       error code: 12

    temperature: float   # Annual average local temperature [degrees Celsius]
                         # valid range: -273 to 6000 degrees Celsius, error code 13

    slope: float = 0     # Surface slope (measured from the horizontal plane)
                         # valid range: -360 to 360 degrees, error code: 14

    azm_rotation: float = 0 # Surface azimuth rotation (measured from south to projection of
                            # surface normal on horizontal plane, negative east)
                            # valid range: -360 to 360 degrees, error code: 15

    atmos_refract: float = 0.5667 # Atmospheric refraction at sunrise and sunset (0.5667 deg is typical)
                             # valid range: -5   to   5 degrees, error code: 16

    function: int = 0        # Switch to choose functions for desired output (from enumeration)

    # Intermediate VALUES
    jd: float = 0       #Julian day
    jc: float = 0       #Julian century

    jde: float = 0      #Julian ephemeris day
    jce: float = 0      #Julian ephemeris century
    jme: float = 0      #Julian ephemeris millennium

    l: float = 0        #earth heliocentric longitude [degrees]
    b: float = 0        #earth heliocentric latitude [degrees]
    r: float = 0        #earth radius vector [Astronomical Units, AU]

    theta: float = 0    #geocentric longitude [degrees]
    beta: float = 0     #geocentric latitude [degrees]

    x0: float = 0       #mean elongation (moon-sun) [degrees]
    x1: float = 0       #mean anomaly (sun) [degrees]
    x2: float = 0       #mean anomaly (moon) [degrees]
    x3: float = 0       #argument latitude (moon) [degrees]
    x4: float = 0       #ascending longitude (moon) [degrees]

    delta_psi: float = 0     #nutation longitude [degrees]
    delta_epsilon: float = 0 #nutation obliquity [degrees]
    epsilon0: float = 0      #ecliptic mean obliquity [arc seconds]
    epsilon: float = 0       #ecliptic true obliquity  [degrees]

    delta_tau: float = 0    #aberration correction [degrees]
    lamda: float = 0        #apparent sun longitude [degrees]
    nu0: float = 0          #Greenwich mean sidereal time [degrees]
    nu: float = 0           #Greenwich sidereal time [degrees]

    alpha: float = 0        #geocentric sun right ascension [degrees]
    delta: float = 0        #geocentric sun declination [degrees]

    h: float = 0            #observer hour angle [degrees]
    xi: float = 0           #sun equatorial horizontal parallax [degrees]
    delta_alpha: float = 0  #sun right ascension parallax [degrees]
    delta_prime: float = 0  #topocentric sun declination [degrees]
    alpha_prime: float = 0  #topocentric sun right ascension [degrees]
    h_prime: float = 0      #topocentric local hour angle [degrees]

    e0: float = 0           #topocentric elevation angle (uncorrected) [degrees]
    delta_e: float = 0      #atmospheric refraction correction [degrees]
    e: float = 0            #topocentric elevation angle (corrected) [degrees]

    eot: float = 0          #equation of time [minutes]
    srha: float = 0         #sunrise hour angle [degrees]
    ssha: float = 0         #sunset hour angle [degrees]
    sta: float = 0          #sun transit altitude [degrees]

    # Final OUTPUT VALUES
    zenith: float = 0       #topocentric zenith angle [degrees]
    azimuth_astro: float = 0#topocentric azimuth angle (westward from south) [for astronomers]
    azimuth: float = 0      #topocentric azimuth angle (eastward from north) [for navigators and solar radiation]
    incidence: float = 0    #surface incidence angle [degrees]

    suntransit: float = 0   #local sun transit time (or solar noon) [fractional hour]
    sunrise: float = 0      #local sunrise time (+/- 30 seconds) [fractional hour]
    sunset: float = 0       #local sunset time (+/- 30 seconds) [fractional hour]

# Calculate SPA output values (in structure) based on input values passed in structure
def spa_calculate(spa: SpaData) -> SpaData:
    spa.jd = julian_day(spa.year, spa.month, spa.day, spa.hour, spa.minute, spa.second, spa.delta_ut1, spa.timezone)

    spa = calculate_geocentric_sun_right_ascension_and_declination(spa)
    
    spa.h  = observer_hour_angle(spa.nu, spa.longitude, spa.alpha)
    spa.xi = sun_equatorial_horizontal_parallax(spa.r)

    spa.delta_prime, spa.delta_alpha = right_ascension_parallax_and_topocentric_dec(
        spa.latitude, spa.elevation, spa.xi, spa.h, spa.delta)

    spa.alpha_prime = topocentric_right_ascension(spa.alpha, spa.delta_alpha)
    spa.h_prime     = topocentric_local_hour_angle(spa.h, spa.delta_alpha)

    spa.e0      = topocentric_elevation_angle(spa.latitude, spa.delta_prime, spa.h_prime)
    spa.delta_e   = atmospheric_refraction_correction(spa.pressure, spa.temperature,
                                                        spa.atmos_refract, spa.e0)
    spa.e       = topocentric_elevation_angle_corrected(spa.e0, spa.delta_e)

    spa.zenith        = topocentric_zenith_angle(spa.e)
    spa.azimuth_astro = topocentric_azimuth_angle_astro(spa.h_prime, spa.latitude,
                                                                        spa.delta_prime)
    spa.azimuth       = topocentric_azimuth_angle(spa.azimuth_astro)

    return spa

def calculate_geocentric_sun_right_ascension_and_declination(spa: SpaData) -> SpaData:
    spa.jd = julian_day(spa.year, spa.month, spa.day, spa.hour, spa.minute, spa.second, spa.delta_ut1, spa.timezone)
    spa.jde = julian_ephemeris_day(spa.jd, 67)
    spa.jc = julian_century(spa.jd)
    spa.jce = julian_ephemeris_century(spa.jde)
    spa.jme = julian_ephemeris_millenium(spa.jce)

    spa.l = earth_heliocentric_longitude(spa.jme)
    spa.b = earth_heliocentric_latitude(spa.jme)
    spa.r = earth_radius_vector(spa.jme)

    spa.theta = geocentric_longitude(spa.l)
    spa.beta = geocentric_latitude(spa.b)

    spa.x0 = mean_elongation_moon_sun(spa.jce)
    spa.x1 = mean_anomaly_sun(spa.jce)
    spa.x2 = mean_anomaly_moon(spa.jce)
    spa.x3 = argument_latitude_moon(spa.jce)
    spa.x4 = ascending_longitude_moon(spa.jce)
    spa.delta_psi, spa.delta_epsilon = nutation_longitude_and_obliquity(spa.jce, [spa.x0, spa.x1, spa.x2, spa.x3, spa.x4])

    spa.epsilon0 = ecliptic_mean_obliquity(spa.jme)
    spa.epsilon = ecliptic_true_obliquity(spa.delta_epsilon, spa.epsilon0)

    spa.delta_tau = aberration_correction(spa.r)
    spa.lamda = apparent_sun_longitude(spa.theta, spa.delta_psi, spa.delta_tau)
    spa.nu0 = greenwich_mean_sidereal_time(spa.jd, spa.jc)
    spa.nu = greenwich_sidereal_time(spa.nu0, spa.delta_psi, spa.epsilon)

    spa.alpha = geocentric_right_ascension(spa.lamda, spa.epsilon, spa.beta)
    spa.delta = geocentric_declination(spa.beta, spa.epsilon, spa.lamda)
    
    return spa

def julian_day(year: int, month: int, day: int, hour: int, minute: int, second: float, dut1: float, tz: float) -> float:
    day_decimal = day + (hour - tz + (minute + (second + dut1) / 60.0) / 60.0) / 24.0

    if (month < 3):
        month = month + 12
        year = year - 1

    julian_day = int(365.25 * (year + 4716)) + int(30.6001 * (month + 1)) + day_decimal - 1524.5

    if (julian_day > 2299160.0):
        a = int(year/100)
        julian_day += (2 - a + int(a/4))

    return julian_day

def julian_ephemeris_day(julian_day: float, delta_t: float) -> float:
    return julian_day + delta_t/86400

def julian_century(jd: float) -> float:
    return (jd - 2451545) / 36525

def julian_ephemeris_century(jde: float) -> float:
    return (jde - 2451545) / 36525

def julian_ephemeris_millenium(jce: float) -> float:
    return jce / 10

def earth_heliocentric_longitude(jme: float) -> float:
    longitude = 0
    for i, terms in enumerate(spa_tables.L_TERMS):
        temp_longitude = 0

        for a, b, c in terms:
            temp_longitude = temp_longitude + a * np.cos(b + c * jme)

        longitude = longitude + temp_longitude * jme**i

    longitude = longitude / 10**8
    longitude = angle_to_0_360_range(np.rad2deg(longitude))
    return longitude

def earth_heliocentric_latitude(jme: float) -> float:
    latitude = 0
    for i, terms in enumerate(spa_tables.B_TERMS):
        temp_latitude = 0

        for a, b, c in terms:
            temp_latitude = temp_latitude + a * np.cos(b + c * jme)

        latitude = latitude + temp_latitude * jme**i

    latitude = latitude / 10**8
    latitude = np.rad2deg(latitude)
    return latitude 

def earth_radius_vector(jme: float) -> float:
    radius = 0
    for i, terms in enumerate(spa_tables.R_TERMS):
        temp_radius = 0

        for a, b, c in terms:
            temp_radius = temp_radius + a * np.cos(b + c * jme)

        radius = radius + temp_radius * jme**i

    radius = radius / 10**8
    return radius 

def geocentric_longitude(longitude: float) -> float:
    return longitude + 180

def geocentric_latitude(latitude: float) -> float:
    return -latitude

def mean_elongation_moon_sun(jce: float) -> float:
    return 297.85036 + 445267.111480 * jce - 0.0019142 * jce**2 + jce**3 / 189474

def mean_anomaly_sun(jce: float) -> float:
    return 357.52772 + 35999.050340 * jce - 0.0001603 * jce**2 + jce**3 / 300000

def mean_anomaly_moon(jce: float) -> float:
    return 134.96298 + 477198.867398 * jce + 0.0086972 * jce**2 + jce**3 / 56250

def argument_latitude_moon(jce: float) -> float:
    return 93.27191 + 483202.017538 * jce - 0.0036825 * jce**2 + jce**3 / 327270

def ascending_longitude_moon(jce: float) -> float:
    return 125.04452 - 1934.136261 * jce + 0.0020708 * jce**2 + jce**3 / 450000

def nutation_longitude_and_obliquity(jce: float, x: list[float]) -> tuple[float, float]:
    nutation_longitude = 0
    nutation_obliquity = 0
    for y_terms, pe_terms in zip(spa_tables.Y_TERMS, spa_tables.PE_TERMS):
        xy_sum = np.dot(x, y_terms)
        sin_xy_sum = np.sin(np.deg2rad(xy_sum))
        cos_xy_sum = np.cos(np.deg2rad(xy_sum))
        a, b, c, d = pe_terms
        nutation_longitude = nutation_longitude + (a + b * jce) * sin_xy_sum
        nutation_obliquity = nutation_obliquity + (c + d * jce) * cos_xy_sum

    nutation_longitude = nutation_longitude / 36000000
    nutation_obliquity = nutation_obliquity / 36000000

    return nutation_longitude, nutation_obliquity

def ecliptic_mean_obliquity(jme: float) -> float:
    u = jme / 10
    epsilon0 = (84381.448 - 4680.93 * u - 1.55 * u**2 + 1999.25 * u**3 - 51.38 * u**4
                - 249.67 * u**5 - 39.05 * u**6 + 7.12 * u**7 + 27.87 * u**8 + 5.79 * u**9 + 2.45 * u**10)
    return epsilon0

def ecliptic_true_obliquity(delta_epsilon: float, epsilon0: float) -> float:
    return epsilon0 / 3600 + delta_epsilon

def aberration_correction(r: float) -> float:
    return -20.4898 / (3600 * r)

def apparent_sun_longitude(theta: float, delta_psi: float, delta_tau: float) -> float:
    return theta + delta_psi + delta_tau

def greenwich_mean_sidereal_time(jd: float, jc: float) -> float:
    nu0 = 280.46061837 + 360.98564736629 * (jd - 2451545.0) + 0.000387933 * jc**2 - jc**3 / 38710000.0
    return angle_to_0_360_range(nu0)

def greenwich_sidereal_time(nu0: float, delta_psi: float, epsilon: float) -> float:
    return nu0 + delta_psi * np.cos(np.deg2rad(epsilon))

def geocentric_right_ascension(lamda: float, epsilon: float, beta: float) -> float: 
    lamda_rad = np.deg2rad(lamda)
    eps_rad = np.deg2rad(epsilon)
    beta_rad = np.deg2rad(beta)
    
    alpha = np.arctan2(np.sin(lamda_rad) * np.cos(eps_rad) - np.tan(beta_rad) * np.sin(eps_rad), np.cos(lamda_rad))
    return angle_to_0_360_range(np.rad2deg(alpha))

def geocentric_declination(beta: float, epsilon: float, lamda: float) -> float:
    lamda_rad = np.deg2rad(lamda)
    eps_rad = np.deg2rad(epsilon)
    beta_rad = np.deg2rad(beta)

    delta = np.arcsin(np.sin(beta_rad) * np.cos(eps_rad) + np.cos(beta_rad) * np.sin(eps_rad) * np.sin(lamda_rad))
    return np.rad2deg(delta)

def observer_hour_angle(nu: float, longitude: float, alpha: float) -> float:
    return angle_to_0_360_range(nu + longitude - alpha)

def sun_equatorial_horizontal_parallax(r: float) -> float:
    xi = 8.794 / (3600.0 * r)
    return xi 

def right_ascension_parallax_and_topocentric_dec(lat: float, elevation: float, xi: float, h:float, delta: float) -> tuple[float, float]:
    lat_rad = np.deg2rad(lat)
    xi_rad = np.deg2rad(xi)
    h_rad = np.deg2rad(h)
    delta_rad = np.deg2rad(delta)
    one_minus_f = 0.99664719 
    elev_term = elevation/6378140.0
    u_term = np.arctan(one_minus_f * np.tan(lat_rad))
    x_term = np.cos(u_term) + elev_term*np.cos(lat_rad)
    y_term = one_minus_f * np.sin(u_term) + elev_term * np.sin(lat_rad)

    delta_alpha_rad = np.arctan2( -x_term*np.sin(xi_rad)*np.sin(h_rad), np.cos(delta_rad) - x_term * np.sin(xi_rad) * np.cos(h_rad))
    delta_alpha = np.rad2deg(delta_alpha_rad)

    delta_prime = np.rad2deg(np.arctan2((np.sin(delta_rad) - y_term*np.sin(xi_rad)) * np.cos(delta_alpha_rad), 
                           np.cos(delta_rad) - x_term * np.sin(xi_rad) * np.cos(h_rad)))

    return delta_prime, delta_alpha

def topocentric_right_ascension(alpha_deg: float, delta_alpha: float) -> float:
    return alpha_deg + delta_alpha

def topocentric_local_hour_angle(h: float, delta_alpha: float) -> float:
    return h - delta_alpha

def topocentric_elevation_angle(lat: float, delta_prime: float, h_prime: float) -> float:
    lat_rad = np.deg2rad(lat)
    delta_prime_rad = np.deg2rad(delta_prime)
    h_prime_rad = np.deg2rad(h_prime)

    e0 = np.arcsin( np.sin(lat_rad) * np.sin(delta_prime_rad) + np.cos(lat_rad)*np.cos(delta_prime_rad)*np.cos(h_prime_rad))
    return np.rad2deg(e0)

def atmospheric_refraction_correction(pressure: float, temperature: float, atmos_refract: float, e0: float) -> float:
    delta_e = 0 
    if (e0 >= -(spa_tables.SUN_RADIUS + atmos_refract)):
        delta_e = (pressure/1010.0) * (283/(273+temperature)) * 1.02/(60*np.tan(np.deg2rad(e0+10.3/(e0+5.11))))
    return delta_e

def topocentric_elevation_angle_corrected(e0: float, delta_e: float) -> float:
    return e0 + delta_e

def topocentric_zenith_angle(e: float) -> float:
    return 90.0 - e

def topocentric_azimuth_angle_astro(h_prime: float, lat: float, delta_prime: float) -> float:
    lat_rad = np.deg2rad(lat)
    h_prime_rad = np.deg2rad(h_prime)
    delta_prime_rad = np.deg2rad(delta_prime)

    tau = np.arctan2( np.sin(h_prime_rad), np.cos(h_prime_rad)*np.sin(lat_rad) - np.tan(delta_prime_rad)*np.cos(lat_rad))

    return angle_to_0_360_range(np.rad2deg(tau))

def topocentric_azimuth_angle(azimuth_astro: float) -> float:
    return angle_to_0_360_range(azimuth_astro + 180)

def angle_to_0_360_range(degrees: float) -> float:
    limited_degrees = degrees % 360             # Limit from -360 to 360
    if limited_degrees < 0:         
        limited_degrees = limited_degrees + 360 # Limit from 0 to 360
    return limited_degrees

if __name__== "__main__":
    spa = SpaData()
    spa.year, spa.month, spa.day, spa.hour, spa.minute, spa.second, spa.delta_ut1, spa.timezone = [2003, 10, 17, 12, 30, 30, 0, -7]

    spa.temperature = 11.0
    spa.pressure = 820
    spa.elevation = 1830.14
    spa.delta_t = 67.0    
    spa.longitude = -105.1786
    spa.latitude = 39.742476 
    spa.slope = 30.0 
    spa.azm_rotation = -10.0
    spa = spa_calculate(spa)

    #print(spa.h, spa.h_prime, spa.alpha_prime, spa.delta_prime, spa.zenith, spa.azimuth)
    print("alpha: ", spa.alpha)
    print("delta: ", spa.delta)
    print("H: ", spa.h)
    print("H': ", spa.h_prime)
    print("alpha_prime: ", spa.alpha_prime)
    print("delta_prime: ", spa.delta_prime)
    print("zenith: ", spa.zenith)
    print("azimuth: ", spa.azimuth)