bolt_data = [
    # name, diameter, UNC TPI, UNC tensile area, UNC minor area,
    #                  UNF TPI, UNF tensile area, UNF minor area,
    #                  clearance_close, clearance_free
    ("#0",     0.0600,  None, None,   None, 80,   0.0018, 0.0015, 0.0635, 0.0700),
    ("#2",     0.0860,  56,   0.0037, 0.0031, 64, 0.0039, 0.0034, 0.0890, 0.0960),
    ("#4",     0.1120,  40,   0.0060, 0.0050, 48, 0.0066, 0.0057, 0.1160, 0.1285),
    ("#5",     0.1250,  40,   0.0080, 0.0067, 44, 0.0083, 0.0072, 0.1285, 0.1360),
    ("#6",     0.1380,  32,   0.0091, 0.0075, 40, 0.0102, 0.0087, 0.1440, 0.1495),
    ("#8",     0.1640,  32,   0.0140, 0.0120, 36, 0.0147, 0.0129, 0.1695, 0.1770),
    ("#10",    0.1900,  24,   0.0175, 0.0145, 32, 0.0200, 0.0175, 0.1960, 0.2010),
    ("1/4",    0.2500,  20,   0.0318, 0.0269, 28, 0.0364, 0.0326, 0.2570, 0.2660),
    ("5/16",   0.3125,  18,   0.0524, 0.0454, 24, 0.0580, 0.0524, 0.3230, 0.3320),
    ("3/8",    0.3750,  16,   0.0775, 0.0678, 24, 0.0878, 0.0809, 0.3860, 0.3970),
    ("7/16",   0.4375,  14,   0.1063, 0.0933, 20, 0.1187, 0.1090, 0.4531, 0.4687),
    ("1/2",    0.5000,  13,   0.1419, 0.1257, 20, 0.1599, 0.1486, 0.5156, 0.5312),
    ("9/16",   0.5625,  12,   0.1820, 0.1620, 18, 0.2030, 0.1890, 0.5781, 0.5938),
    ("5/8",    0.6250,  11,   0.2260, 0.2020, 18, 0.2560, 0.2400, 0.6406, 0.6562),
    ("3/4",    0.7500,  10,   0.3340, 0.3020, 16, 0.3730, 0.3510, 0.7656, 0.7812),
    ("7/8",    0.8750,   9,   0.4620, 0.4190, 14, 0.5090, 0.4800, 0.8906, 0.9062),
    ("1",       1.0000,   8,   0.6060, 0.5510, 12, 0.6630, 0.6250, 1.0156, 1.0313),

    # Clearance-hole data not provided in your table:
    ("1-1/8",   1.1250,   7,   0.7630, 0.6930, 12, 0.8560, 0.8120, None, None),
    ("1-1/4",   1.2500,   7,   0.9690, 0.8900, 12, 1.0730, 1.0240, None, None),
    ("1-3/8",   1.3750,   6,   1.1550, 1.0540, 12, 1.3150, 1.2600, None, None),
    ("1-1/2",   1.5000,   6,   1.4050, 1.2940, 12, 1.5810, 1.5210, None, None),
    ("1-3/4",   1.7500,   5,   1.9000, 1.7400, None, None, None, None, None),
    ("2",        2.0000,  4.5, 2.5000, 2.3000, None, None, None, None, None),
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
            "clearance_close": bolt[8],
            "clearance_free": bolt[9],
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

    # NUMERICAL LOOKUP — rounds  up
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
