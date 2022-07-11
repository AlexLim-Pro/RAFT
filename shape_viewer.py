import shapefile as shp
import matplotlib.pyplot as plt

sf = shp.Reader("../San_Gaud_data/catchment_Reg12/catchment_Reg12.shp")

plt.figure()
for shape in sf.shapeRecords():
    x = [i[0] for i in shape.shape.points[:]]
    y = [i[1] for i in shape.shape.points[:]]
    plt.plot(x, y, color="#ffffff")

sf2 = shp.Reader("../San_Gaud_data/NHDFlowline_San_Guad/NHDFlowline_San_Guad.shp")

plt.figure()
for shape in sf2.shapeRecords():
    x = [i[0] for i in shape.shape.points[:]]
    y = [i[1] for i in shape.shape.points[:]]
    plt.plot(x, y, color="#ffffff")

plt.show()
