def quaternion_to_euler(qx, qy, qz, qw):
    # taken from http://tech-artists.org/forum/showthread.php?5356-Quaternion-to-Euler-Rotations
    sqw = qw * qw
    sqx = qx * qx
    sqy = qy * qy
    sqz = qz * qz

    normal = math.sqrt(sqw + sqx + sqy + sqz)
    pole_result = (qx * qz) + (qy * qw)

    if (pole_result > (0.5 * normal)): # singularity at north pole
        ry = math.pi/2 #heading/yaw?
        rz = 0 #attitude/roll?
        rx = 2 * math.atan2(qx, qw) #bank/pitch?
        return Euler(rx, ry, rz)
    if (pole_result < (-0.5 * normal)): # singularity at south pole
        ry = -math.pi/2
        rz = 0
        rx = -2 * math.atan2(qx, qw)
        return Euler(rx, ry, rz)

    r11 = 2*(qx*qy + qw*qz)
    r12 = sqw + sqx - sqy - sqz
    r21 = -2*(qx*qz - qw*qy)
    r31 = 2*(qy*qz + qw*qx)
    r32 = sqw - sqx - sqy + sqz

    rx = math.atan2( r31, r32 )
    ry = math.asin ( r21 )
    rz = math.atan2( r11, r12 )

    return (rx, ry, rz)

def get_tuples(line):
    l = line.replace("  X= ","")
    v = l.split(" ")
    v = [float(w) for w in v]
    r = []
    for i in range(0,len(v)/7):
        r.append(v[i*7:i*7+7])
    pos = [[v[0],v[1],v[2]] for v in r]
    orientation = [quaternion_to_euler(v[3],v[4],v[5],v[6]) for v in r]
    return (pos,orientation)

def get_orientation(line):
    l = line.replace("  X= ","")
    v = l.split(" ")
    v = [float(w) for w in v]
    r = []
    for i in range(0,len(v)/7):
        r.append(v[i*7:i*7+7])
    orientation = [quaternion_to_euler(v[3],v[4],v[5],v[6]) for v in r]
    return orientation

def get_position(line):
    l = line.replace("  X= ","")
    v = l.split(" ")
    v = [float(w) for w in v]
    r = []
    for i in range(0,len(v)/7):
        r.append(v[i*7:i*7+7])
    pos = [[v[0],v[1],v[2]] for v in r]
    return pos

def convert_to_csv(lst):
    s = [str(v[0])+","+str(v[1])+","+str(v[2]) for v in lst]
    r = ""
    for a in s:
        r += "," + a
    return r[1:-1]


