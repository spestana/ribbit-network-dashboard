
import numpy as np

def get_terminator(datetime):
	terminator = {"type":"FeatureCollection","features":[{"type":"Feature","properties":{},"geometry":{"coordinates":[[[-106.76819661251528,49.51129792824884],[-116.15268948499772,46.34267698384144],[-121.61598228901042,40.7192418896266],[-110.30211492574841,31.418378479764357],[-97.82997319420195,29.111582552975122],[-91.81128823303374,35.458684030347854],[-94.78424958305635,44.4134732619591],[-98.99366831628987,48.352280555504535],[-106.76819661251528,49.51129792824884]]],"type":"Polygon"}}]}
	return terminator
	
#def julian(datetime):
#	'''Calculate the present UTC Julian Date. Function is valid after
#	the beginning of the UNIX epoch 1970-01-01 and ignores leap
#	seconds.'''
#	return (datetime / 86400000) + 2440587.5
#
#def GMST(julianDay):
#	''' Calculate Greenwich Mean Sidereal Time according to 
#		 http://aa.usno.navy.mil/faq/docs/GAST.php '''
#	d = julianDay - 2451545.0
#	# Low precision equation is good enough for our purposes.
#	return (18.697374558 + 24.06570982441908 * d) % 24
#	
#def sunEclipticPosition(julianDay)
#		''' Compute the position of the Sun in ecliptic coordinates at
#			 julianDay.  Following
#			 http://en.wikipedia.org/wiki/Position_of_the_Sun '''
#		# Days since start of J2000.0
#		n = julianDay - 2451545.0;
#		# mean longitude of the Sun
#		L = 280.460 + 0.9856474 * n;
#		L %= 360;
#		# mean anomaly of the Sun
#		g = 357.528 + 0.9856003 * n;
#		g %= 360;
#		# ecliptic longitude of Sun
#		l = L + 1.915 * np.sin(g * this._D2R) + 0.02 * np.sin(2 * g * this._D2R);
#		# distance from Sun in AU
#		R = 1.00014 - 0.01671 * np.cos(g * this._D2R) - 0.0014 * np.cos(2 * g * this._D2R);
#		return {l: l, R: R}
#		
#		
#def eclipticObliquity: function (julianDay) {
#		// Following the short term expression in
#		// http://en.wikipedia.org/wiki/Axial_tilt#Obliquity_of_the_ecliptic_.28Earth.27s_axial_tilt.29
#		var n = julianDay - 2451545.0;
#		// Julian centuries since J2000.0
#		var T = n / 36525;
#		var epsilon = 23.43929111 -
#			T * (46.836769 / 3600
#				- T * (0.0001831 / 3600
#					+ T * (0.00200340 / 3600
#						- T * (0.576e-6 / 3600
#							- T * 4.34e-8 / 3600))));
#		return epsilon
#		
#def sunEquatorialPosition: function (sunEclLng, eclObliq) {
#		/* Compute the Sun's equatorial position from its ecliptic
#		 * position. Inputs are expected in degrees. Outputs are in
#		 * degrees as well. */
#		var alpha = Math.atan(Math.cos(eclObliq * this._D2R)
#			* Math.tan(sunEclLng * this._D2R)) * this._R2D;
#		var delta = Math.asin(Math.sin(eclObliq * this._D2R)
#			* Math.sin(sunEclLng * this._D2R)) * this._R2D;
#
#		var lQuadrant = Math.floor(sunEclLng / 90) * 90;
#		var raQuadrant = Math.floor(alpha / 90) * 90;
#		alpha = alpha + (lQuadrant - raQuadrant);
#
#		return {alpha: alpha, delta: delta}
#
#def hourAngle: function (lng, sunPos, gst) {
#		/* Compute the hour angle of the sun for a longitude on
#		 * Earth. Return the hour angle in degrees. */
#		var lst = gst + lng / 15;
#		return lst * 15 - sunPos.alpha;
#	},
#
#def latitude: function (ha, sunPos) {
#		/* For a given hour angle and sun position, compute the
#		 * latitude of the terminator in degrees. */
#		var lat = Math.atan(-Math.cos(ha * this._D2R) /
#			Math.tan(sunPos.delta * this._D2R)) * this._R2D;
#		return lat;
#	},
#
#def compute: function (time) {
#		var today = time ? new Date(time) : new Date();
#		var julianDay = julian(today);
#		var gst = GMST(julianDay);
#		var latLng = [];
#
#		var sunEclPos = this._sunEclipticPosition(julianDay);
#		var eclObliq = this._eclipticObliquity(julianDay);
#		var sunEqPos = this._sunEquatorialPosition(sunEclPos.lambda, eclObliq);
#		for (var i = 0; i <= 720 * this.options.resolution; i++) {
#			var lng = -360 + i / this.options.resolution;
#			var ha = this._hourAngle(lng, sunEqPos, gst);
#			latLng[i + 1] = [this._latitude(ha, sunEqPos), lng];
#		}
#		if (sunEqPos.delta < 0) {
#			latLng[0] = [90, -360];
#			latLng[latLng.length] = [90, 360];
#		} else {
#			latLng[0] = [-90, -360];
#			latLng[latLng.length] = [-90, 360];
#		}
#		return latLng;