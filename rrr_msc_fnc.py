#******************************************************************************
#rrr_msc_fnc.py
#******************************************************************************

# Purpose:
# Miscellaneous useful functions
# Authors:
# Alex Christopher Lim, Cedric H. David, 2022-2022

import numpy as np


class DistFncs(object):
    """
    Functions that are helpful when using coordinate distances
    """
    @staticmethod
    def great_circle_dist(long1, long2, lat1, lat2, is_rad=False):
        """
        Calculates the great circle distance between two coordinates

        :param long1: A longitude in degrees or radians
        :type long1: float
        :param long2: A longitude in degrees or radians
        :type long2: float
        :param lat1: A latitude in degrees or radians
        :type lat1: float
        :param lat2: A latitude in degrees or radians
        :type lat2: float
        :param is_rad: Whether the coordinates are given in degrees or radians,
            defaults to False (degrees)
        :type is_rad: bool
        :return: The great circle distance in degrees or radians
        :rtype: float
        """
        if not is_rad:
            long1 = np.deg2rad(long1)
            long2 = np.deg2rad(long2)
            lat1 = np.deg2rad(lat1)
            lat2 = np.deg2rad(lat2)
        d = 2 * np.arcsin(np.sqrt(
            (np.sin(lat1 - lat2) / 2) ** 2
            + np.cos(lat1) * np.cos(lat2) * (np.sin((long1 - long2) / 2)) ** 2
        ))
        if not is_rad:
            return np.rad2deg(d)
        return d

    @staticmethod
    def m2km(dist):
        """
        Converts coordinate distance from meters to kilometers

        :param dist: The coordinate distance in meters
        :type dist: float
        :return: The coordinate distance in kilometers
        """
        return dist / 1000

    @staticmethod
    def km2m(dist):
        """
        Converts coordinate distance from kilometers to meters

        :param dist: The coordinate distance in kilometers
        :type dist: float
        :return: The coordinate distance in meters
        """
        return dist * 1000

    @staticmethod
    def deg2M(dist):
        """
        Converts coordinate distance from degrees to nautical miles

        :param dist: The coordinate distance in degrees
        :type dist: float
        :return: The coordinate distance in nautical miles
        :rtype: float
        """
        return dist * 60

    @staticmethod
    def rad2M(dist):
        """
        Converts coordinate distance from radians to nautical miles

        :param dist: The coordinate distance in radians
        :type dist: float
        :return: The coordinate distance in nautical miles
        :rtype: float
        """
        return DistFncs.deg2M(np.rad2deg(dist))

    @staticmethod
    def M2deg(dist):
        """
        Converts coordinate distance from nautical miles to degrees

        :param dist: The coordinate distance in nautical miles
        :type dist: float
        :return: The coordinate distance in degrees
        :rtype: float
        """
        return dist / 60

    @staticmethod
    def M2rad(dist):
        """
        Converts coordinate distance from nautical miles to radians

        :param dist: The coordinate distance in nautical miles
        :type dist: float
        :return: The coordinate distance in radians
        :rtype: float
        """
        return np.deg2rad(DistFncs.M2deg(dist))

    @staticmethod
    def M2m(dist):
        """
        Converts coordinate distance from nautical miles to meters

        :param dist: The coordinate distance in nautical miles
        :type dist: float
        :return: The coordinate distance in meters
        :rtype: float
        """
        return dist * 1852

    @staticmethod
    def m2M(dist):
        """
        Converts coordinate distance from meters to nautical miles

        :param dist: The coordinate distance in meters
        :type dist: float
        :return: The coordinate distance in nautical miles
        :rtype: float
        """
        return dist / 1852

    @staticmethod
    def M2km(dist):
        """
        Converts coordinate distance from nautical miles to kilometers

        :param dist: The coordinate distance in nautical miles
        :type dist: float
        :return: The coordinate distance in kilometers
        :rtype: float
        """
        return DistFncs.m2km(DistFncs.M2m(dist))

    @staticmethod
    def km2M(dist):
        """
        Converts coordinate distance from kilometers to nautical miles

        :param dist: The coordinate distance in kilometers
        :type dist: float
        :return: The coordinate distance in nautical miles
        :rtype: float
        """
        return DistFncs.m2M(DistFncs.km2m(dist))

    @staticmethod
    def deg2km(dist):
        """
        Converts coordinate distance from degrees to kilometers

        :param dist: The coordinate distance in degrees
        :type dist: float
        :return: The coordinate distance in kilometers
        :rtype: float
        """
        return DistFncs.M2km(DistFncs.deg2M(dist))

    @staticmethod
    def km2deg(dist):
        """
        Converts coordinate distance from kilometers to degrees

        :param dist: The coordinate distance in kilometers
        :type dist: float
        :return: The coordinate distance in degrees
        :rtype: float
        """
        return DistFncs.M2deg(DistFncs.km2M(dist))

    @staticmethod
    def rad2km(dist):
        """
        Converts coordinate distance from radians to kilometers

        :param dist: The coordinate distance in radians
        :type dist: float
        :return: The coordinate distance in kilometers
        :rtype: float
        """
        return DistFncs.M2km(DistFncs.rad2M(dist))

    @staticmethod
    def km2rad(dist):
        """
        Converts coordinate distance from kilometers to radians

        :param dist: The coordinate distance in kilometers
        :type dist: float
        :return: The coordinate distance in radians
        :rtype: float
        """
        return DistFncs.M2rad(DistFncs.km2M(dist))

    @staticmethod
    def mi2km(dist):
        """
        Converts coordinate distance from miles to kilometers

        :param dist: The coordinate distance in miles
        :type dist: float
        :return: The coordinate distance in kilometers
        :rtype: float
        """
        return dist * 25146 / 15625

    @staticmethod
    def km2mi(dist):
        """
        Converts coordinate distance from kilometers to miles

        :param dist: The coordinate distance in kilometers
        :type dist: float
        :return: The coordinate distance in miles
        :rtype: float
        """
        return dist * 15625 / 25146

    @staticmethod
    def mi2deg(dist):
        """
        Converts coordinate distance from miles to degrees

        :param dist: The coordinate distance in miles
        :type dist: float
        :return: The coordinate distance in degrees
        :rtype: float
        """
        return DistFncs.km2deg(DistFncs.mi2km(dist))

    @staticmethod
    def deg2mi(dist):
        """
        Converts coordinate distance from degrees to miles

        :param dist: The coordinate distance in degrees
        :type dist: float
        :return: The coordinate distance in miles
        :rtype: float
        """
        return DistFncs.km2mi(DistFncs.deg2km(dist))
