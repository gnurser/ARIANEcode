#From: http://www.movable-type.co.uk/scripts/latlong.html
import math


def distance(lon1, lat1, lon2, lat2):
    # http://code.activestate.com/recipes/576779-calculating-distance-between-two-geographic-points/
    # http://en.wikipedia.org/wiki/Haversine_formula
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    q = math.sin(dlat/2)**2 + (math.cos(lat1) * math.cos(lat2) * (math.sin(dlon/2)**2))
    return 2. * math.atan2(math.sqrt(q), math.sqrt(1-q))


def small_circle(radius,clon,clat,llon,llat,meters):
    # clon = circle central longitude (deg)
    # clat = circle central latitude (deg)
    # llon = longitude of point on circle line (deg)
    # llat = latitude of point on circle line (deg)
    # meters = point spacing (meters). This is not the exact spacing, it will be used to calculate nrpoints.
    # radius = radius of sphere

    dlat = math.radians(llat-clat)
    dlon = math.radians(llon-clon)
    clat = math.radians(clat)
    clon = math.radians(clon)
    llat = math.radians(llat)
    llon = math.radians(llon)

    a = math.sin(dlat/2) * math.sin(dlat/2) + math.sin(dlon/2) * math.sin(dlon/2) * math.cos(clat) * math.cos(llat)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    d = radius * c

    # bearing
    y = math.sin(dlon) * math.cos(llat)
    x = math.cos(clat)*math.sin(llat) - math.sin(clat)*math.cos(llat)*math.cos(dlon)
    brng0 = math.atan2(y, x)

    # circle
    angle = distance(clon,clat,llon,llat)
    circleradius = math.sin(angle) * radius
    circlelength = 2*math.pi*circleradius

    # print(angle, circleradius, circlelength)

    # loop circle
    list = []
    list.append([math.degrees(llon),math.degrees(llat)])
    nrpoints = int(circlelength / meters)
    print(nrpoints)
    for i in range(1, nrpoints):
        brng = brng0 + (float(i) / float(nrpoints))*2.*math.pi
        lat = math.asin( math.sin(clat)*math.cos(d/radius) + math.cos(clat)*math.sin(d/radius)*math.cos(brng) )
        lon = clon + math.atan2(math.sin(brng)*math.sin(d/radius)*math.cos(clat),
                                math.cos(d/radius)-math.sin(clat)*math.sin(lat))
        list.append([math.degrees(lon),math.degrees(lat)])

    return list
