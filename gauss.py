import math


def bl2xy(lon: float, lat: float):
    """ 大地2000,经纬度转平面坐标,3度带
    Param:
        lon (float): 经度
        lat (float): 纬度
    Returns:
        (x , y) : x坐标对应经度,y坐标对应纬度
    """

    # 3.1415926535898/180.0
    iPI = 0.0174532925199433
    # 3度带
    zoneWide = 3
    # 长半轴
    a = 6378137
    # 扁率
    f = 1/298.257222101

    projNo = int(lon/zoneWide)
    longitude0 = projNo*3
    longitude0 = longitude0 * iPI
    longitude1 = lon * iPI
    latitude1 = lat * iPI

    e2 = 2 * f - f * f
    ee = e2 * (1.0 - e2)
    NN = a / math.sqrt(1.0 - e2 * math.sin(latitude1) * math.sin(latitude1))
    T = math.tan(latitude1) * math.tan(latitude1)
    C = ee * math.cos(latitude1) * math.cos(latitude1)
    A = (longitude1 - longitude0) * math.cos(latitude1)
    M = a * ((1 - e2 / 4 - 3 * e2 * e2 / 64 - 5 * e2 * e2 * e2 / 256) * latitude1 - (3 * e2 / 8 + 3 * e2 * e2 / 32 + 45 * e2 * e2 * e2 / 1024) *
             math.sin(2 * latitude1) + (15 * e2 * e2 / 256 + 45 * e2 * e2 * e2 / 1024) * math.sin(4 * latitude1) - (35 * e2 * e2 * e2 / 3072) * math.sin(6 * latitude1))

    xval = NN * (A + (1 - T + C) * A * A * A / 6 + (5 - 18 * T +
                                                    T * T + 72 * C - 58 * ee) * A * A * A * A * A / 120)
    yval = M + NN * math.tan(latitude1) * (A * A / 2 + (5 - T + 9 * C + 4 * C * C) * A * A *
                                           A * A / 24 + (61 - 58 * T + T * T + 600 * C - 330 * ee) * A * A * A * A * A * A / 720)

    X0 = 1000000 * projNo + 500000
    Y0 = 0
    xval = xval + X0
    yval = yval + Y0
    return (xval, yval)


def xy2bl(X: float, Y: float):
    """ 大地2000,平面坐标转经纬度,3度带
    Param:
        X (float): 坐标X,对应经度
        Y (float): 坐标Y,对应纬度
    Returns:
        (lon , lat): (经度,纬度)
    """

    # 3.1415926535898/180.0
    iPI = 0.0174532925199433
    # 3度带
    zoneWide = 3
    # 长半轴
    a = 6378137
    # 扁率
    f = 1/298.257222101

    projNo = int(X / 1000000)
    longitude0 = projNo*3
    longitude0 = longitude0 * iPI

    X0 = projNo * 1000000 + 500000
    Y0 = 0
    xval = X - X0
    yval = Y - Y0
    e2 = 2 * f - f * f
    e1 = (1.0 - math.sqrt(1 - e2)) / (1.0 + math.sqrt(1 - e2))
    ee = e2 / (1 - e2)
    M = yval
    u = M / (a * (1 - e2 / 4 - 3 * e2 * e2 / 64 - 5 * e2 * e2 * e2 / 256))
    fai = u + (3 * e1 / 2 - 27 * e1 * e1 * e1 / 32) * math.sin(2 * u) + (21 * e1 * e1 / 16 - 55 * e1 * e1 * e1 * e1 / 32) * \
        math.sin(4 * u) + (151 * e1 * e1 * e1 / 96) * math.sin(6 * u) + \
        (1097 * e1 * e1 * e1 * e1 / 512) * math.sin(8 * u)
    C = ee * math.cos(fai) * math.cos(fai)
    T = math.tan(fai) * math.tan(fai)
    NN = a / math.sqrt(1.0 - e2 * math.sin(fai) * math.sin(fai))
    R = a * (1 - e2) / math.sqrt((1 - e2 * math.sin(fai) * math.sin(fai)) * (1 -
                                                                             e2 * math.sin(fai) * math.sin(fai)) * (1 - e2 * math.sin(fai) * math.sin(fai)))
    D = xval / NN
    longitude1 = longitude0 + (D - (1 + 2 * T + C) * D * D * D / 6 + (5 - 2 * C + 28 *
                                                                      T - 3 * C * C + 8 * ee + 24 * T * T) * D * D * D * D * D / 120) / math.cos(fai)
    latitude1 = fai - (NN * math.tan(fai) / R) * (D * D / 2 - (5 + 3 * T + 10 * C - 4 * C * C - 9 * ee) * D *
                                                  D * D * D / 24 + (61 + 90 * T + 298 * C + 45 * T * T - 256 * ee - 3 * C * C) * D * D * D * D * D * D / 720)
    return (longitude1 / iPI, latitude1 / iPI)


# 示例
print(bl2xy(120.37437907087538, 31.580046209339972))
print(xy2bl(40535536.66, 3495347.439))
