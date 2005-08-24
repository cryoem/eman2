from libpyEMData2 import *

class EMData(EMData_cpp):
    def __getitem__(self, key):
        try:
            size = len(key)
            if 3 == size:
                ix = key[0]
                iy = key[1]
                iz = key[2]
                return self(ix, iy, iz)
            elif 2 == size:
                ix = key[0]
                iy = key[1]
                return self(ix, iy)
            elif 1 == size:
                ix = key[0]
                return self(ix)
            else:
                raise IndexError, "Invalid index: " + str(key)
        except TypeError:
            # A number, not a tuple, hopefully
            ix = key
            return self(ix)
    def __setitem__(self, key, val):
        xoff,yoff,zoff = self.get_array_offsets()
        try:
            size = len(key)
            if 3 == size:
                ix = key[0] - xoff
                iy = key[1] - yoff
                iz = key[2] - zoff
                self.set_value_at(ix, iy, iz, val)
            elif 2 == size:
                ix = key[0] - xoff
                iy = key[1] - yoff
                self.set_value_at(ix, iy, val)
            elif 1 == size:
                ix = key[0] - xoff
                self.set_value_at(ix, val)
            else:
                raise IndexError, "Invalid index: " + str(key)
        except TypeError:
            # A number, not a tuple, hopefully
            ix = key - xoff
            self.set_value_at(ix, val)


