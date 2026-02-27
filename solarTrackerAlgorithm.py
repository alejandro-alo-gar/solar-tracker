import pvlib
import pandas as pd
import matplotlib.pyplot as plt
from pvlib.pvsystem import PVSystem, Array, SingleAxisTrackerMount
from pvlib.location import Location
from pvlib.modelchain import ModelChain
from pvlib.tracking import singleaxis
from dataclasses import dataclass

# Data class to collect the location, weather, PV module, inverter and temperature parameters data for the model
@dataclass
class modelContext:
    location: Location
    weather: object
    pv_module: dict
    inverter: dict
    temp_params: dict

def compute_energy(context : modelContext, surface_tilt, surface_azimuth):

    # Compute solar position (zenith & azimuth) for each timestamp
    solpos = context.location.get_solarposition(context.weather.index)
    
    # Compute extraterrestrial DNI (solar constant adjusted for day of year)
    dni_extra = pvlib.irradiance.get_extra_radiation(context.weather.index)

    # Transpose horizontal irradiance (GHI, DNI, DHI) to plane-of-array irradiance
    total_irrad = pvlib.irradiance.get_total_irradiance(surface_tilt, 
                                                        surface_azimuth,
                                                        solpos['apparent_zenith'], 
                                                        solpos['azimuth'], 
                                                        context.weather['dni'], 
                                                        context.weather['ghi'], 
                                                        context.weather['dhi'], 
                                                        dni_extra = dni_extra, 
                                                        model='haydavies')
    
    # Estimate PV cell temperature using SAPM (Sandia Array Performance Model) temperature model
    cell_temp = pvlib.temperature.sapm_cell(total_irrad['poa_global'], 
                                            context.weather['temp_air'], 
                                            context.weather['wind_speed'], 
                                            **context.temp_params)
    
    # Compute relative airmass from solar zenith angle
    airmass = pvlib.atmosphere.get_relative_airmass(solpos['apparent_zenith'])

    # Convert site altitude to atmospheric pressure
    pressure = pvlib.atmosphere.alt2pres(context.location.altitude)

    # Compute absolute airmass (corrected for local pressure)
    am_abs = pvlib.atmosphere.get_absolute_airmass(airmass, pressure)

    # Compute angle of incidence between sun rays and module surface
    aoi = pvlib.irradiance.aoi(surface_tilt, surface_azimuth, solpos['apparent_zenith'], solpos['azimuth'])

    # Adjust irradiance for spectral effects and AOI losses (SAPM model)
    effective_irrad = pvlib.pvsystem.sapm_effective_irradiance(total_irrad['poa_direct'], 
                                                               total_irrad['poa_diffuse'], 
                                                               am_abs, 
                                                               aoi, 
                                                               context.pv_module)

    # Compute DC module output (Vmp, Pmp) using Sandia PV model
    dc = pvlib.pvsystem.sapm(effective_irrad, cell_temp, context.pv_module)
    
    # Convert DC power to AC power using Sandia inverter model
    ac = pvlib.inverter.sandia(dc['v_mp'], dc['p_mp'], context.inverter)

    return ac.fillna(0) # Replace NaN values (e.g. nighttime) with zero power

def run_automatic_benchmark(context: modelContext):
    mount = SingleAxisTrackerMount(axis_tilt = 0, axis_azimuth=180, max_angle=60, cross_axis_tilt=0, racking_model='open_rack')
    array = Array( mount=mount, module_parameters=context.pv_module, temperature_model_parameters=context.temp_params)
    system = PVSystem(arrays=[array], inverter_parameters=context.inverter)
    modelChain = ModelChain(system, context.location)
    modelChain.run_model(context.weather)
    return modelChain.results.ac.fillna(0)

def run_pvlib_tracking_benchmark(context : modelContext):
    solpos = context.location.get_solarposition(context.weather.index)

    tracking_data = singleaxis(apparent_zenith=solpos['apparent_zenith'],
                               solar_azimuth=solpos['azimuth'],
                               axis_tilt=0,
                               axis_azimuth=180,
                               max_angle=60,
                               backtrack=True,
                               gcr=0.3,
                               cross_axis_tilt=0)

    surface_tilt = tracking_data['surface_tilt']
    surface_azimuth = tracking_data['surface_azimuth']

    return compute_energy(context, surface_tilt, surface_azimuth)

def run_my_tracking(context: modelContext):
    solpos = context.location.get_solarposition(context.weather.index)

    # Example: naive perfect east-west tracker
    tracker_theta = 90 - solpos['apparent_zenith']

    # Limit to mechanical range
    tracker_theta = tracker_theta.clip(-60, 60)

    surface_tilt = tracker_theta.abs()
    surface_azimuth = 180  # simplify first

    return compute_energy(context,surface_tilt, surface_azimuth)

# Select location and date
loc = Location(57.05, 9.93, name='Aalborg', altitude=10, tz='Etc/GMT-1')
date = '2018-05-15'
year = date.split('-')[0]

ctx = modelContext(location=loc, 
                   weather=pvlib.iotools.get_pvgis_tmy(loc.latitude,loc.longitude,coerce_year=int(year))[0].loc[date], 
                   pv_module=pvlib.pvsystem.retrieve_sam('SandiaMod')['Canadian_Solar_CS6X_300M__2013_'], 
                   inverter=pvlib.pvsystem.retrieve_sam('cecinverter')['ABB__PVI_3_6_OUTD_US__240V_'],
                   temp_params=pvlib.temperature.TEMPERATURE_MODEL_PARAMETERS['sapm']['open_rack_glass_glass'])
print(" Weather, PV Module, Inverter and Temperature models loaded ")

ac_auto = run_automatic_benchmark(ctx)
ac_pvlib = run_pvlib_tracking_benchmark(ctx)
ac_mine = run_my_tracking(ctx)

plt.figure()
plt.title("AC Power comparison")
plt.ylabel("AC Power (W)")
plt.plot(ac_auto, label='Auto')
plt.plot(ac_pvlib, label='pvlib singleaxis', linestyle='--')
plt.plot(ac_mine, label='my tracking', linestyle='dashdot')
plt.legend()
plt.grid()
plt.show()

'''
# tracker_angle: rotation of the tracker at each timestep
solpos = location.get_solarposition(weather.index)
tracker_angle = singleaxis(
    apparent_zenith=solpos['apparent_zenith'],
    solar_azimuth=solpos['azimuth'],
    axis_tilt=mount.axis_tilt,
    axis_azimuth=mount.axis_azimuth,
    max_angle=mount.max_angle,
    backtrack=True,           # if you want backtracking
    gcr=0.3,                  # ground coverage ratio, adjust if known
    cross_axis_tilt=mount.cross_axis_tilt
)

tracker_theta = tracker_angle['tracker_theta'].bfill().ffill()
aoi = tracker_angle['aoi'].bfill().ffill()
surface_tilt = tracker_angle['surface_tilt'].bfill().ffill()
surface_azimuth = tracker_angle['surface_azimuth'].bfill().ffill()

# Plot tracker rotation
plt.figure()
plt.plot(tracker_theta, label='tracker_theta')
plt.plot(aoi, label='aoi')
plt.plot(surface_tilt, label='surface_tilt')
plt.plot(surface_azimuth, label='surface_azimuth')
plt.plot(solpos['apparent_elevation'].loc[date], label='solar_elevation', linestyle='--')
plt.legend()
plt.show()
'''