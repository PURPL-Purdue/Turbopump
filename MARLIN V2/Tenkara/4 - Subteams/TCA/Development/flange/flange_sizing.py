## Gasketed Flange Sizing code
# Author: Louis DeSano
# Date: 08/29/2026

# Sources: BPVC Appendix 2
import numpy as np
import pandas as pd

# Unit conversions
IN2M = 0.0254
PSI2PA = 6894.76
N2LBF = 0.224809
boltNames = ["#0", "#2", "#4", "#5", "#6", "#8", "#10", "1/4", "5/16", "3/8", "7/16", "1/2", 
                "9/16", "5/8", "3/4", "7/8", "1", "1-1/8", "1-1/4", "1-3/8", "1-1/2"] 
bolt_data = [
    # name, diameter, UNC TPI, UNC tensile area, UNC minor area,
    #                  UNF TPI, UNF tensile area, UNF minor area
    ("#0",     0.0600,  None, None,   None, 80,   0.0018, 0.0015),
    ("#2",     0.0860,  56,   0.0037, 0.0031, 64, 0.0039, 0.0034),
    ("#4",     0.1120,  40,   0.0060, 0.0050, 48, 0.0066, 0.0057),
    ("#5",     0.1250,  40,   0.0080, 0.0067, 44, 0.0083, 0.0072),
    ("#6",     0.1380,  32,   0.0091, 0.0075, 40, 0.0102, 0.0087),
    ("#8",     0.1640,  32,   0.0140, 0.0120, 36, 0.0147, 0.0129),
    ("#10",    0.1900,  24,   0.0175, 0.0145, 32, 0.0200, 0.0175),
    ("1/4",    0.2500,  20,   0.0318, 0.0269, 28, 0.0364, 0.0326),
    ("5/16",   0.3125,  18,   0.0524, 0.0454, 24, 0.0580, 0.0524),
    ("3/8",    0.3750,  16,   0.0775, 0.0678, 24, 0.0878, 0.0809),
    ("7/16",   0.4375,  14,   0.1063, 0.0933, 20, 0.1187, 0.1090),
    ("1/2",    0.5000,  13,   0.1419, 0.1257, 20, 0.1599, 0.1486),
    ("9/16",   0.5625,  12,   0.1820, 0.1620, 18, 0.2030, 0.1890),
    ("5/8",    0.6250,  11,   0.2260, 0.2020, 18, 0.2560, 0.2400),
    ("3/4",    0.7500,  10,   0.3340, 0.3020, 16, 0.3730, 0.3510),
    ("7/8",    0.8750,   9,   0.4620, 0.4190, 14, 0.5090, 0.4800),
    ("1",       1.0000,   8,   0.6060, 0.5510, 12, 0.6630, 0.6250),
    ("1-1/8",   1.1250,   7,   0.7630, 0.6930, 12, 0.8560, 0.8120),
    ("1-1/4",   1.2500,   7,   0.9690, 0.8900, 12, 1.0730, 1.0240),
    ("1-3/8",   1.3750,   6,   1.1550, 1.0540, 12, 1.3150, 1.2600),
    ("1-1/2",   1.5000,   6,   1.4050, 1.2940, 12, 1.5810, 1.5210),
    ("1-3/4",   1.7500,   5,   1.9000, 1.7400, None, None, None),
    ("2",       2.0000,  4.5,  2.5000, 2.3000, None, None, None),
]

def bolt_lookup(return_property, value, lookup_property, thread_type):

    thread_type = thread_type.upper()
    lookup_property = lookup_property.lower()
    return_property = return_property.lower()

    if thread_type not in ("UNC", "UNF"):
        raise ValueError('thread_type must be "UNC" or "UNF"')

    if lookup_property not in ("name", "area", "diameter"):
        raise ValueError(
            'lookup_property must be "name", "area", or "diameter"'
        )

    # Build valid bolt list
    bolts = []

    for bolt in bolt_data:

        name = bolt[0]
        diameter = bolt[1]

        if thread_type == "UNC":
            tpi = bolt[2]
            area = bolt[3]
            minor_area = bolt[4]
        else:
            tpi = bolt[5]
            area = bolt[6]
            minor_area = bolt[7]

        if area is None:
            continue

        bolts.append({
            "name": name,
            "diameter": diameter,
            "tpi": tpi,
            "area": area,
            "minor_area": minor_area,
            "thread": thread_type,
        })

    # NAME LOOKUP
    if lookup_property == "name":
        value = str(value)

        for bolt in bolts:
            if bolt["name"] == value:
                selected_bolt = bolt
                break
        else:
            raise ValueError(
                f'Bolt "{value}" not found for {thread_type}'
            )

    # NUMERICAL LOOKUP — ALWAYS ROUND UP
    else:
        bolts.sort(key=lambda x: x[lookup_property])

        selected_bolt = None

        for bolt in bolts:
            if bolt[lookup_property] >= value:
                selected_bolt = bolt
                break

        if selected_bolt is None:
            raise ValueError(
                f"No {thread_type} bolt is large enough for "
                f"{lookup_property} = {value}"
            )

    # Return all properties
    if return_property == "":
        return selected_bolt

    if return_property not in selected_bolt:
        raise ValueError(
            f"Unknown return property '{return_property}'. "
            f"Available properties: {list(selected_bolt.keys())}"
        )

    return selected_bolt[return_property]

class FlangeSizer:
    def __init__(self, G, P, m, y, b, Sa, Sb, t):
        #define gasket characteristics
        #operating conditions
        self.P = P # internal design pressure [Pa]
        
        #gasket related
        self.m = m # gasket factor []
        self.y = y # gasket unit seating load [Pa]
        self.b = b # effective gasket seating width [m]
        self.G = G # diameter at location of gasket load [m]

        #bolt related
        self.Sa = Sa # Allowable stress for bolt at atmospheric temp [Pa]
        self.Sb = Sb # Allowable stress for bolt at design temp [Pa]
        self.t = t # flange thickness [m]

        # determine gasket operating load
        self.Wm1 = self.operating_load()

        # determine gasket seating load
        self.Wm2 = self.seating_load()
        
        self.Am = self.bolt_area()
        print(self.Wm1*N2LBF, self.Wm2*N2LBF, self.Am/(IN2M**2))
        # choose bolt configuration
        self.bolt_iteration()

        return
    
    def operating_load(self):
        """
        Calculates operating load of a gasket
            Arguments: 
                G (float) diameter at location of gasket load [m]
                P (float) internal design pressure [Pa]
                b (float) effective gasket seating width [m]
                m (float) gasket factor []

            Return: 
                Wm1 (float) Gasket operating load [N]
        """
        H = 0.785 * self.G**2 * self.P # Hydrostatic End Force
        Hp = 2 * self.b * np.pi * self.G * self.m * self.P # Joint Contact Force 
        Wm1 = H + Hp
        return Wm1

    def seating_load(self):
        """
        Calculates seating load of a gasket
            Arguments: 
                G (float) diameter at location of gasket load [m]
                b (float) effective gasket seating width [m]
                y (float) gasket unit seating load [Pa]

            Return: 
                Wm2 (float) Gasket seating load [N]
        """
        Wm2 = np.pi * self.b * self.G * self.y
        return Wm2

    def bolt_area(self):
        """
        Calculates required bolt area for gasket operation
            Arguments: 
                Wm1 (float) Gasket operating load [N]
                Sb (float) Allowable stress for bolt at design temp [Pa]
                Wm2 (float) Gasket seating load [N]
                Sa (float) Allowable stress for bolt at atmospheric temp [Pa]


            Return: 
                Am (float) required bolt area [m^2]
        """
        Am1 = self.Wm1 / self.Sb # Bolt area required for operating stress
        Am2 = self.Wm2 / self.Sa # Bolt area required for seating stress
        Am = max(Am1, Am2)
        if Am1 > Am2:
            print("Bolt Operating Stress is driving.")
        else:
            print("Bolt Seating Stress is driving.")

        return Am

    def bolt_spacing(self, a):
        """
        Calculates maximum bolt spacing for "lethal service"
            Arguments: 
                a (float) nominal bolt diameter [m]
                t (float) flange thickness [m]
                m (float) gasket factor []

            Return: 
                Bs (float) maximum bolt spacing [m]
        """
        Bs = 2 * a + (6 * self.t) / (self.m + 0.5)
        return Bs

    def bolt_iteration(self):
        """Iterates on bolt size to find valid configs of bolt size/count"""
     
        circum = np.pi * (self.G + self.b) # aproximate bolt circle circumference

        print("Possible Bolt Configurations:")

        for boltName in boltNames:
            actualSpacing = 0
            requiredSpacing = 1
            validFlag = True

            # determine max bolt spacing 
            boltArea = bolt_lookup("area", boltName, "name", "UNF") * IN2M**2
            boltDiam = bolt_lookup("diameter", boltName, "name", "UNF") * IN2M
            requiredSpacing = self.bolt_spacing(boltDiam)
            
            minCount = int(np.ceil(self.Am / boltArea))
            count = minCount - 1
            actualSpacing = circum / minCount
            
            while actualSpacing > requiredSpacing: #increase count until spacing is satisfied
                count += 1
                actualSpacing = circum / count

            if validFlag:
                print(f"  {count}x {boltName}")
        
        return boltName, count
  


# For each bolt
#   while spacing is invalid
        # increase count

def main():
    
    G = 3.826 * IN2M
    P = 440 * PSI2PA
    m = 3
    y = 2500 * PSI2PA
    b = 1 * IN2M
    Sa = 30000 * PSI2PA
    Sb = Sa
    t = 1/2 * IN2M

    FlangeSizer(G, P, m, y, b, Sa, Sb, t)
    return 

if __name__ == "__main__":
    main()