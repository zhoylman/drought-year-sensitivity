#Hargreaves-Samani reference evapotranspiration

#The FAO56 method recommendsthat, in the absence of all 
#other meteorological data required for the Penman-Monteith equation, 
#the temperature-based estimate of ETo proposed by Hargreaves and Samani should be used.

#temperature is in degrees C

hs_eto = function(Tmax, Tmin, Julian_day, lat){
  #derived following http://www.fao.org/3/x0490e/x0490e07.htm#TopOfPage
  #https://doi.org/10.1371/journal.pone.0174045 See as reference
  #define solar constant
  Gsc = 0.0820 #MJ m-2 min-1
  #derive inverse relative distance
  dr = 1 + 0.033*cos(((2*pi)/365)*Julian_day)
  # derive solar declination
  d = 0.409*sin((((2*pi)/365)*Julian_day)-1.39)
  #lat in radians
  j = (pi/180)*lat
  #derive sunset hour angle
  omega = acos(-tan(j) * tan(d))
  #derive extra terrestrial solar radiation MJ m-2 d-1
  Ra = (((24*60)/pi)*Gsc*dr)*((omega*(sin(j)*sin(d)))+(cos(j)*cos(d)*(sin(omega))))
  #convert to depth (mm d-1)
  Ra_depth = 0.408*Ra
  #now for the ETo calculation (Hargreaves-Samani reference evapotranspiration)
  #define constants
  a = 0.0023
  ETo = a*Ra_depth*((Tmax - Tmin)^0.5)*((((Tmax+Tmin)/2))+17.8)
  return(ETo)
}