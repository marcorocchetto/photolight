
import logging
etaslog = logging.getLogger('etaslog')


class Observatory:
    def __init__(self, kwargs):

        if 'name' in kwargs:
            self.name = kwargs['name']
        else:
            self.name = 'Unknown'
            etaslog.warning('Observatory name not specified')

        if 'elevation' in kwargs:
            self.elevation = kwargs['elevation']
        else:
            etaslog.error('You have to specify the elevation of the observatory')
            self.exit()
        if 'longitud' in kwargs:
            self.elevation = kwargs['longitud']
        else:
            etaslog.error('You have to specify the longitude of the observatory')
            self.exit()
        if 'latitud' in kwargs:
            self.elevation = kwargs['latitud']
        else:
            etaslog.error('You have to specify the latitude of the observatory')
            self.exit()

