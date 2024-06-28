from numpy import sqrt

def pSBC(p, c0, c1=None, l=False):
    def parse_color(d):
        n = len(d)
        x = {}
        if n > 9:
            d = d.split(",")
            r, g, b, *a = d
            x['r'] = int(r[3:] if r[3] == 'a' else r[4:])
            x['g'] = int(g)
            x['b'] = int(b)
            x['a'] = float(a[0]) if a else -1
        else:
            if n == 8 or n == 6 or n < 4:
                return None
            if n < 6:
                d = "#" + d[1] * 2 + d[2] * 2 + d[3] * 2 + (d[4] * 2 if n > 4 else "")
            d = int(d[1:], 16)
            if n == 9 or n == 5:
                x['r'] = d >> 24 & 255
                x['g'] = d >> 16 & 255
                x['b'] = d >> 8 & 255
                x['a'] = round((d & 255) / 0.255) / 1000
            else:
                x['r'] = d >> 16
                x['g'] = d >> 8 & 255
                x['b'] = d & 255
                x['a'] = -1
        return x

    def round_to_int(val):
        return round(val)

    if not isinstance(p, (int, float)) or p < -1 or p > 1 or not isinstance(c0, str) or (c0[0] not in 'r#') or (c1 and not isinstance(c1, str)):
        return None

    h = len(c0) > 9
    a = isinstance(c1, str)
    h = a and (len(c1) > 9 if c1 != "c" else not h)

    f = parse_color(c0)
    P = p < 0
    t = parse_color(c1) if c1 and c1 != "c" else ({"r": 0, "g": 0, "b": 0, "a": -1} if P else {"r": 255, "g": 255, "b": 255, "a": -1})
    p = -p if P else p
    P = 1 - p

    if not f or not t:
        return None

    if l:
        r = round_to_int(P * f['r'] + p * t['r'])
        g = round_to_int(P * f['g'] + p * t['g'])
        b = round_to_int(P * f['b'] + p * t['b'])
    else:
        r = round_to_int((P * f['r'] ** 2 + p * t['r'] ** 2) ** 0.5)
        g = round_to_int((P * f['g'] ** 2 + p * t['g'] ** 2) ** 0.5)
        b = round_to_int((P * f['b'] ** 2 + p * t['b'] ** 2) ** 0.5)

    a1 = f['a']
    a2 = t['a']
    a = a1 >= 0 or a2 >= 0
    a = a1 if a1 < 0 else a2 if a2 < 0 else a1 * P + a2 * p

    if h:
        return "rgba({}, {}, {}, {:.3f})".format(r, g, b, round(a * 1000) / 1000) if a else "rgb({}, {}, {})".format(r, g, b)
    else:
        return "#{:08x}".format((4294967296 + r * 16777216 + g * 65536 + b * 256 + (round(a * 255) if a else 0)) & 0xFFFFFFFF) if a else "#{:06x}".format((r << 16) + (g << 8) + b)

def hex_to_rgb(hex_color):
    '''Converts a hex color code to RGB values from 0 to 1'''
    # Remove the '#' if it exists
    hex_color = hex_color.lstrip('#')
    
    # Parse hex color code into RGB values
    rgb = list(int(hex_color[i:i+2], 16) for i in (0, 2, 4))
    # Convert to 0-1 scale
    return [i/255 for i in rgb]

def brightness_calculation(hex_color):
    # Convert hex color code to RGB
    rgb = hex_to_rgb(hex_color)
    
    # Calculate brightness
    return sqrt(0.299 * rgb[0]**2 + 0.587 * rgb[1]**2 + 0.114 * rgb[2]**2)