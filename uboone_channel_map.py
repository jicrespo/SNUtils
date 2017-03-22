# Python module defining the ChannelMap class for converting from crate, FEM, channel to plane and wire, and vice versa

class ChannelMap:
    def __init__(self, mapfile):
        self.mapfile = mapfile
        self.table_readout2larsoft = {} # The keys are the trio crate, FEM, channels
        self.table_larsoft2readout = {} # The keys are plane, wire number
        with open(mapfile) as inf:
            for line in inf.readlines():
                crate, fem, ch, plane, wire = line.split()
                self.table_readout2larsoft[(int(crate), int(fem), int(ch))] = [plane, int(wire)]
                self.table_larsoft2readout[(plane, int(wire))] = [int(crate), int(fem), int(ch)]

    def CrateFEMCh2PlaneWire(self, crate, fem, ch):
        try:
            return self.table_readout2larsoft[(int(crate), int(fem), int(ch))]
        except KeyError:
            return ["", -1]

    def PlaneWire2CrateFEMCh(self, plane, wire):
        try:
            return self.table_larsoft2readout[(plane, int(wire))]
        except KeyError:
            return [-1, -1, -1]
